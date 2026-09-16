#!/usr/bin/env python3
"""HydroTuring adapter. Standard library only, to show that the /io contract
needs no scientific Python stack and could be written in any language."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "gwex", "mrso", "snw", "canopy", "channel", "stage"]

MODEL = {"name": "reference_bucket", "version": "1.0.0"}


EVAP_SHAPE = 0.5  # soil moisture at which evaporation reaches its potential rate

SECONDS_PER_DAY = 86400.0

# Default reach geometry, used when the catchment does not hand one over.
DEFAULT_WIDTH_M = 18.0
DEFAULT_SLOPE = 0.0015
DEFAULT_MANNING_N = 0.035

# Steps the contract can name, as a fraction of a day. Forcing and reported
# fluxes are rates in mm per day at every step; the depth moved in one step
# is the rate times this.
TIMESTEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}


def stage_of(runoff_rate_mm_day: float, static: dict) -> float:
    """The level a gauge in the reach would read, in metres.

    The bucket does not route: runoff leaves the stores and the catchment in
    the same step, so it holds no water in transit and its channel store is
    identically zero. A reach that a flow passes straight through still has a
    stage, though, and that is what a gauge in it would read — Manning normal
    depth in the catchment's channel, strictly increasing in the flow.

    The gauge reads the flow it is reporting and nothing else, so the rating
    is a single-valued function of `dis` by construction. That is the honest
    statement for this model: it has no second time constant, because it has
    no store the water waits in, so its two limbs cannot separate by anything
    a survey would resolve. Reading the gauge off the store instead would
    report a constant zero, which is a number that carries no information
    rather than a measurement.
    """
    return _manning(max(runoff_rate_mm_day, 0.0), static)


def _manning(flow_rate_mm_day: float, static: dict) -> float:
    area_km2 = float(static.get("area_km2", 0.0))
    width_m = float(static.get("width_m", DEFAULT_WIDTH_M))
    slope = float(static.get("slope", DEFAULT_SLOPE))
    manning_n = float(static.get("manning_n", DEFAULT_MANNING_N))
    if area_km2 <= 0.0 or width_m <= 0.0 or slope <= 0.0:
        return 0.0
    q_m3s = max(flow_rate_mm_day, 0.0) * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY
    if q_m3s <= 0.0:
        return 0.0
    return (q_m3s * manning_n / (width_m * slope ** 0.5)) ** 0.6


def discharge_m3s(runoff_rate_mm_day: float, static: dict) -> float:
    """The flow through the reach, in m3/s, from the catchment's runoff.

    The model holds no water in transit, so the flow through its reach is the
    flow it generated. Reporting it is what lets a rating be drawn against a
    discharge rather than against a store that is identically zero.
    """
    area_km2 = float(static.get("area_km2", 0.0))
    return max(runoff_rate_mm_day, 0.0) * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY


def simulate(forcing, static, dt_days=1.0):
    """A conceptual bucket that conserves water exactly by construction.

    Interception, degree-day snow, saturation-excess runoff, linear baseflow,
    and soil-moisture-limited evaporation. Every flux is removed from the
    store it actually came from, so the budget closes to floating point.

    `dt_days` is the length of one forcing row. Rates are turned into depths
    with it on the way in and back into rates on the way out, so the same
    catchment integrates the same water whatever step the weather arrives
    at. At a daily step every factor is exactly 1.0 and the arithmetic is
    bit for bit what it was before the step was a parameter.

    When the forcing carries an `abstr` column (mm/day, net of return flow),
    the prescribed withdrawal is taken from the soil store first and any
    remainder from the day's runoff before it leaves, and whatever was
    actually removed is declared as a negative `gwex`. Absent the column the
    model is bit for bit the original bucket.
    """
    soil_cap = static["soil_capacity_mm"]
    canopy_cap = static["canopy_capacity_mm"]
    ddf = static["degree_day_factor_mm_per_C_day"]
    k_base = static["baseflow_coefficient"]
    t_snow = static["snow_threshold_degC"]

    soil = 0.5 * soil_cap
    swe = 0.0
    canopy = 0.0
    rows = []

    for step in forcing:
        pr_rate, tas, pet_rate = step["pr"], step["tas"], step["pet"]
        abstr_rate = step.get("abstr", 0.0)
        pr = pr_rate * dt_days
        pet = pet_rate * dt_days

        # The human term, first call on the store.
        want = max(abstr_rate, 0.0) * dt_days
        removed = min(soil, want)
        soil -= removed

        snowfall = pr if tas < t_snow else 0.0
        rain = 0.0 if tas < t_snow else pr

        swe += snowfall
        melt = min(swe, ddf * max(tas - t_snow, 0.0) * dt_days)
        swe -= melt

        water_in = rain + melt
        intercepted = min(canopy_cap - canopy, water_in)
        canopy += intercepted
        throughfall = water_in - intercepted

        canopy_evap = min(canopy, pet)
        canopy -= canopy_evap
        pet_left = pet - canopy_evap

        soil += throughfall
        surface = max(0.0, soil - soil_cap)
        soil -= surface
        baseflow = k_base * soil * dt_days
        soil -= baseflow
        soil_evap = min(soil, pet_left * min(1.0, soil / (EVAP_SHAPE * soil_cap)))
        soil -= soil_evap

        # The human term, second call: the day's outflow.
        rest = want - removed
        divert = min(surface, rest)
        surface -= divert
        removed += divert
        rest -= divert
        divert = min(baseflow, rest)
        baseflow -= divert
        removed += divert

        runoff = (surface + baseflow) / dt_days
        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (canopy_evap + soil_evap) / dt_days,
            "mrro": runoff,
            "dis": discharge_m3s(runoff, static),
            "gwex": -removed / dt_days,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
            # Runoff leaves the stores and the catchment in the same step: no
            # routing, so the water in transit is identically zero. Reported,
            # not omitted, because it is a statement about the model.
            "channel": 0.0,
            "stage": stage_of(runoff, static),
        })
    return rows

def read_request(path: Path) -> tuple[dict, Path]:
    request = json.loads(path.read_text())
    return request, path.parent


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet", "abstr"):
            if key in row:
                row[key] = float(row[key])
    return rows


def write_result(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    request, io_dir = read_request(request_path)
    forcing = read_forcing(io_dir / request["input"]["forcing"])
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    timestep = request.get("timestep", "PT1D")
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    rows = simulate(forcing, static, TIMESTEP_DAYS[timestep])

    write_result(io_dir / request["output"]["table"], rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows)}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

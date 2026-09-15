#!/usr/bin/env python3
"""HydroTuring adapter for reference_driven_exchange.

The exact bucket with a boundary exchange driven by the external head the case
prescribes: the flux is proportional to the difference between the `gwh`
column and a fixed catchment reference head, so the catchment gains water
whenever the external head stands above its own and loses it whenever the head
falls below. It reads the driver, declares that it does, and follows it.

It is the positive control for mass/exchange-response and exists to
exercise the one distinction the probe makes: because the prescribed head
oscillates inside the rainless windows, this exchange **reverses direction
there** — as often as the head crosses its mean, and no more — and it must
pass. A control whose exchange happened to be one-signed where it is scored
would show nothing.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLUMNS = ["time", "pr", "evspsbl", "mrro", "gwex", "mrso", "snw", "canopy", "channel"]

MODEL = {"name": "reference_driven_exchange", "version": "1.0.0"}

# The catchment's own head, and the conductance of its boundary: a metre of
# head difference moves this many millimetres a day.
REFERENCE_HEAD_M = 10.0
CONDUCTANCE_MM_PER_M_DAY = 2.0


EVAP_SHAPE = 0.5  # soil moisture at which evaporation reaches its potential rate

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

        # The prescribed boundary head, followed as given. Positive into the
        # catchment when the external head stands above the reference; a loss
        # is limited to the water the soil actually holds.
        head = step.get("gwh", REFERENCE_HEAD_M)
        exchange = CONDUCTANCE_MM_PER_M_DAY * (head - REFERENCE_HEAD_M) * dt_days
        exchange = max(exchange, -soil)
        soil += exchange

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

        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (canopy_evap + soil_evap) / dt_days,
            "mrro": (surface + baseflow) / dt_days,
            "gwex": (exchange - removed) / dt_days,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
            # Runoff leaves the stores and the catchment in the same step: no
            # routing, so the water in transit is identically zero. Reported,
            # not omitted, because it is a statement about the model.
            "channel": 0.0,
        })
    return rows

def read_request(path: Path) -> tuple[dict, Path]:
    request = json.loads(path.read_text())
    return request, path.parent


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet", "abstr", "gwh"):
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

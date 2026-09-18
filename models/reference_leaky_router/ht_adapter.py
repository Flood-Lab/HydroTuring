#!/usr/bin/env python3
"""HydroTuring adapter for reference_leaky_router: the exact bucket with one
deliberate fault, sized to the smallest leak `momentum/routing-conservation`
exists to catch.

The kernel sums to 0.999, so a tenth of a percent of every day's runoff enters
the channel store and never leaves. `channel` is reported as generated minus
released, which is what the contract asks, and it grows without bound — but
slowly: the residue is under a millimetre on this probe's case, while the
proportional part of the allowance is hundreds of millimetres during a storm.

Why this model exists next to `reference_stuck_router`. That one loses a tenth
of the flow and is caught three orders of magnitude clear of the allowance, so
it is insensitive to where the allowance sits: raising `min_allowance_mm` by a
factor of twenty still fails it. This one is caught at 4.1 to 4.4 times the
allowance, so it is the case that actually governs the number — put it behind
the gate and a change that moves the allowance out from under the calibration
turns the gate red instead of merely turning a unit test red.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLUMNS = ['time', 'pr', 'evspsbl', 'mrro', 'mrso', 'snw', 'canopy', 'channel']

MODEL = {"name": "reference_leaky_router", "version": "1.0.0"}


EVAP_SHAPE = 0.5  # soil moisture at which evaporation reaches its potential rate

# The share of each day's runoff the kernel releases. The unit test that
# constructs the same fault uses `[0.5, 1/3, 1/6]` normalised and scaled to
# this, so the two are the same arithmetic on the same case.
ROUTED_SHARE = 0.999
_RAW_WEIGHTS = (0.5, 1.0 / 3.0, 1.0 / 6.0)
_RAW_TOTAL = _RAW_WEIGHTS[0] + _RAW_WEIGHTS[1] + _RAW_WEIGHTS[2]
WEIGHTS = tuple(w / _RAW_TOTAL * ROUTED_SHARE for w in _RAW_WEIGHTS)

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
    """The exact bucket, routed through a kernel that keeps a tenth of a percent.

    Identical to `reference_stuck_router` except for the kernel's sum. Every
    flux is removed from the store it actually came from, so the budget closes
    to floating point: what the kernel holds back is in `channel`, and
    generated minus released is exactly the store it reports.
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
    generated = []
    cum_gen = cum_out = 0.0

    for step in forcing:
        pr_rate, tas, pet_rate = step["pr"], step["tas"], step["pet"]
        pr = pr_rate * dt_days
        pet = pet_rate * dt_days

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

        generated.append(surface + baseflow)
        n = len(generated)
        routed = sum(WEIGHTS[k] * generated[n - 1 - k] for k in range(min(len(WEIGHTS), n)))
        cum_gen += surface + baseflow
        cum_out += routed
        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (canopy_evap + soil_evap) / dt_days,
            "mrro": routed / dt_days,
            "channel": cum_gen - cum_out,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
        })
    return rows


def read_request(path: Path) -> tuple[dict, Path]:
    request = json.loads(path.read_text())
    return request, path.parent


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet"):
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

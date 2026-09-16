#!/usr/bin/env python3
"""HydroTuring adapter for reference_unreported_inflow: the bucket that routes
exactly and then, on rainless steps, quietly fills its reach from nowhere.

The physics is `reference_bucket`'s, unchanged: interception, degree-day snow,
saturation-excess runoff, linear baseflow, soil-moisture-limited evaporation,
and every flux removed from the store it came from. The routing is exact — the
runoff generated in a step is the runoff released in it, so `mrro` is the
bucket's and the whole-catchment budget closes to the floating point.

The one thing added is a **bounded, unreported inflow to the reach**. On every
step the sky did not feed, the channel receives a small fixed amount of water
out of a store the model never declares:

    channel += min(0.05 mm - channel, 6e-5 mm)     on a rainless step

The amount is deliberately tiny. At 0.05 mm it is four orders of magnitude
below the runoff of a storm, and — this is the point — it never approaches the
bound `momentum/routing-conservation` enforces, which allows a channel to hold
`max_lag_days` of the largest runoff of the preceding window. On this probe's
record the peak runoff decays through each long dry spell, the allowance
shrinks with it, and the store still stays well inside it.

So the two criteria part company on this model:

* `momentum/routing-conservation` passes it: the store is non-negative and
  never reaches its bound, because the error does not accumulate.
* `momentum/channel-routing-mass` fails it: on step after step where no rain
  fell, the store rose, which is the only thing a reach with nothing entering
  it is not allowed to do.

That is the gap the direction test exists to cover. A bound can only see an
error that grows; a one-sided statement about the sign of the change sees an
error of any size, as long as it is applied where nothing is coming in to mask
it.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLUMNS = ['time', 'pr', 'evspsbl', 'mrro', 'mrso', 'snw', 'canopy', 'channel']

MODEL = {"name": "reference_unreported_inflow", "version": "1.0.0"}

EVAP_SHAPE = 0.5  # soil moisture at which evaporation reaches its potential rate

# Rain below this rate, in mm/day, counts as no rain, so the step is one on
# which the reach is not being fed from the sky.
DRY_PR_MM_PER_DAY = 0.05
# What the reach quietly receives on such a step, and the size of the store it
# comes out of. Both are small enough that the fifteen-day bound stays an order
# of magnitude away, which is what makes this model the one that separates the
# direction test from the bound.
INFLOW_MM_PER_DAY = 6e-5
INFLOW_TOTAL_MM = 0.05

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


def increase(channel: float, pr_rate_mm_day: float, dt_days: float) -> float:
    """The channel store after one step, including the unreported inflow.

    Zero on a step it rained, because the fault is only visible where nothing
    is entering to mask it; otherwise the store moves towards its ceiling at a
    fixed small rate. Split out as a function so the fault is one readable
    line rather than a condition buried in the integration loop.
    """
    if pr_rate_mm_day > DRY_PR_MM_PER_DAY:
        return channel
    return min(INFLOW_TOTAL_MM, channel + INFLOW_MM_PER_DAY * dt_days)


def simulate(forcing, static, dt_days=1.0):
    """The exact bucket, routing exactly, with the reach quietly refilled.

    Identical to `reference_bucket` except for the `channel` column: the water
    budget closes to the floating point and `mrro` is the bucket's, so a
    failure on any criterion is a failure of the reach and not of the
    bookkeeping.
    """
    soil_cap = static["soil_capacity_mm"]
    canopy_cap = static["canopy_capacity_mm"]
    ddf = static["degree_day_factor_mm_per_C_day"]
    k_base = static["baseflow_coefficient"]
    t_snow = static["snow_threshold_degC"]

    soil = 0.5 * soil_cap
    swe = 0.0
    canopy = 0.0
    channel = 0.0
    rows = []

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

        channel = increase(channel, pr_rate, dt_days)

        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (canopy_evap + soil_evap) / dt_days,
            # Routing is exact: what was generated this step is released this
            # step, and the reported channel is not derived from it at all.
            "mrro": (surface + baseflow) / dt_days,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
            "channel": channel,
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
        json.dumps({
            "status": "ok",
            "model": MODEL,
            "n_steps": len(rows),
            "notes": {
                "channel": (
                    f"the reach receives up to {INFLOW_TOTAL_MM} mm from a store "
                    "this model never reports, on the steps where no rain fell"
                ),
            },
        }, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

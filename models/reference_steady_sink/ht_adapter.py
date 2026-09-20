#!/usr/bin/env python3
"""HydroTuring adapter for reference_steady_sink.

reference_token_exchange with its noise taken out. It overstates its
evaporation by a constant 0.3 mm/day — 1095 mm over the record, 13% of the
rain — and declares the difference as a steady regional inflow, so its budget
still closes exactly. The same million-times-too-weak head trickle rides on
top.

It is the cheat a floor on the response's share of the exchange's *variation*
cannot see. A constant sink contributes nothing to the variation, so the only
thing moving the declared exchange is the trickle, and a trickle of any size
clears a share of its own movement. Its head response is 1.7e-6 of its gross
exchange, the same share as reference_token_exchange, which the probe does
catch. What catches this one is the second term in the bar, a much smaller
share of the gross, which is finite for a steady exchange where the variation
share is not.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLUMNS = ["time", "pr", "evspsbl", "mrro", "gwex", "mrso", "snw", "canopy", "channel"]

MODEL = {"name": "reference_steady_sink", "version": "1.0.0"}

# Evaporation overstated by a constant, declared as a steady inflow. Unlike
# the noise of reference_token_exchange this adds nothing to the exchange's
# variation, which is the whole point of the control.
STEADY_BIAS_MM_PER_DAY = 0.3

TOKEN_CONDUCTANCE_MM_PER_M_DAY = 1.0e-6   # the token head term


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


def simulate(forcing, static, dt_days=1.0, seed=0):
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
    actually removed is declared as a negative `gwex`. What is never absent is
    the cheat: the reported evaporation carries a constant
    `STEADY_BIAS_MM_PER_DAY` on top of the true one, the difference is declared
    as `gwex`, and the token head trickle is folded into the same column and
    into the soil, so the budget still closes to floating point.
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

        # The token term: a head-driven trickle so small it changes nothing
        # about the model, added to the soil so the budget still closes.
        trickle = TOKEN_CONDUCTANCE_MM_PER_M_DAY * (step.get("gwh", 10.0) - 10.0) * dt_days
        soil += trickle

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
            # Water conserved internally; evaporation misreported, and the
            # exchange declared as exactly the difference. The budget closes
            # because gwex is minus the residual of everything else.
            "evspsbl": (canopy_evap + soil_evap) / dt_days + STEADY_BIAS_MM_PER_DAY,
            "mrro": (surface + baseflow) / dt_days,
            "gwex": 0.0,   # filled in below, once the reported ET is known
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
            # Runoff leaves the stores and the catchment in the same step: no
            # routing, so the water in transit is identically zero. Reported,
            # not omitted, because it is a statement about the model.
            "channel": 0.0,
        })
        true_et = (canopy_evap + soil_evap) / dt_days
        rows[-1]["gwex"] = (rows[-1]["evspsbl"] - true_et) - removed / dt_days + trickle / dt_days
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
    rows = simulate(forcing, static, TIMESTEP_DAYS[timestep],
                    seed=int(request.get("seed", 0)))

    write_result(io_dir / request["output"]["table"], rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows)}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

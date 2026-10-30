#!/usr/bin/env python3
"""HydroTuring adapter for the Tank Hydrologic Model (nzahasan/tank-model).

The physics is upstream's: `tank_core.computation_helpers.compute_project`
runs `tank_core.tank_basin.tank_discharge` on one Subbasin, at the pinned
commit, unmodified. This file is the units layer and the bookkeeping. README.md
states every parameter and choice; the three that need the code to follow are
these.

Rates. Upstream's coefficients are fractions per row (`del_t` only enters the
conversion to m3/s), and the published values are per day. Each coefficient is
rescaled to the case's step before the model sees it, k -> 1 - (1 - k)^dt.
Outlet heights and initial storages are depths and are not scaled.

The end of a row. Upstream updates `S[t+1] = max(S[t] + P[t+1] - E[t+1] - out[t], 0)`
(Paik et al. 2005, eq. 4-5), so `tank_storage[t]` is the store after row t's
rain and before row t's outflow. The contract's budget differences states
across a row, with the row's fluxes, so what is reported is the store after
row t's outflow: `S_end[t] = S[t] + in[t] - out[t]`, computed from upstream's
own flow arrays. That is the same trajectory, read half a step later.

Evaporation. Upstream returns no evaporation: it subtracts the demand from the
top tank and resets a negative result to zero. What left is therefore
`E - clamp`, the demand less the amount the reset gave back, and that is what
is reported under `evspsbl`. It is not free to absorb a leak: the clamp is
computed from the arrays and checked against the demand, the same reset in the
three lower tanks would be water from nowhere, and both totals are written to
run.json. `pr` on row 0 never enters the model (upstream sets `S[0]` to the
initial storage), so row 0's evaporation is zero; row 0 is always in the spinup.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tank_core import computation_helpers as ch

COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "mrso", "gw"]
MODEL = {"name": "tank", "version": "3.0.0b3-aea105c.1", "revision": "aea105ca32d0957fee76fbcfe3a5b8993938c9a1"}

TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0}
SUBBASIN = "catchment"

# The midpoint of each range in the "Sugawara" column of the table shipped as
# docs/ParameterRange.pdf in the upstream repository (Table 5-3-2). No
# calibration exists for a synthetic catchment. Coefficients per day, depths mm.
PARAMS_PER_DAY = {
    "t0_soc_uo": 0.3,      # S11, 0.1-0.5, upper side outlet of the top tank
    "t0_soc_lo": 0.3,      # S12, 0.1-0.5, lower side outlet of the top tank
    "t0_boc": 0.3,         # B11, 0.1-0.5, bottom outlet of the top tank
    "t1_soc": 0.065,       # S21, 0.03-0.1
    "t1_boc": 0.055,       # B21, 0.01-0.1
    "t2_soc": 0.0075,      # S31, 0.005-0.01
    "t2_boc": 0.0075,      # B31, 0.005-0.01
    "t3_soc": 0.00275,     # S41, 0.0005-0.005
}
HEIGHTS_MM = {
    "t0_soh_uo": 45.0,     # H11, 30-60
    "t0_soh_lo": 20.0,     # H12, 10-30
    "t1_soh": 25.0,        # H21, 0-50
    "t2_soh": 15.0,        # H31, 0-30
}
# IST3 and IST4 are in the same column; it gives no range for the top two
# tanks, so theirs are the midpoint of the table's "This study" column, 0-50.
INITIAL_MM = {"t0_is": 25.0, "t1_is": 25.0, "t2_is": 55.0, "t3_is": 300.0}


def per_step(fraction_per_day: float, dt_days: float) -> float:
    return 1.0 - (1.0 - fraction_per_day) ** dt_days


def parameters(dt_days: float) -> dict:
    scaled = {name: per_step(k, dt_days) for name, k in PARAMS_PER_DAY.items()}
    return {**scaled, **HEIGHTS_MM, **INITIAL_MM}


def simulate(forcing: list[dict], area_km2: float, dt_days: float) -> tuple[list[dict], dict]:
    """Run upstream's model and turn its arrays into the contract's rows."""
    index = pd.Index([row["time"] for row in forcing], name="Time")
    pr_rate = np.array([row["pr"] for row in forcing], dtype=float)
    pet_rate = np.array([row["pet"] for row in forcing], dtype=float)
    pr_mm = pr_rate * dt_days
    pet_mm = pet_rate * dt_days

    basin = {
        "root_node": [SUBBASIN],
        "basin_def": {SUBBASIN: {"type": "Subbasin", "area": area_km2, "parameters": parameters(dt_days)}},
    }
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        discharge, states = ch.compute_project(
            basin,
            pd.DataFrame({SUBBASIN: pr_mm}, index=index),
            pd.DataFrame({SUBBASIN: pet_mm}, index=index),
            del_t=dt_days * 24.0,
        )
    tank = states[SUBBASIN]["tank_storage"]
    side = states[SUBBASIN]["side_outlet_flow"]
    bottom = states[SUBBASIN]["bottom_outlet_flow"]

    # Each tank after row t's flows. Tank 3 has no bottom outlet.
    inflow = np.zeros_like(tank)
    inflow[:, 1:] = bottom
    drained = np.zeros_like(tank)
    drained[:, :3] = bottom
    end = tank + inflow - drained - side

    # What upstream's reset added back on row t+1, tank by tank: the store it
    # computed, S_end[t] + (P - E on the top tank), was negative by this much.
    applied = np.zeros_like(tank)
    applied[1:, 0] = pr_mm[1:] - pet_mm[1:]
    clamp = np.zeros_like(tank)
    clamp[1:] = tank[1:] - (end[:-1] + applied[1:])

    evaporated_mm = pet_mm - clamp[:, 0]
    evaporated_mm[0] = 0.0
    runoff_mm = side.sum(axis=1)

    rows = [
        {
            "time": forcing[i]["time"],
            "pr": forcing[i]["pr"],
            "evspsbl": evaporated_mm[i] / dt_days,
            "mrro": runoff_mm[i] / dt_days,
            "dis": float(discharge[SUBBASIN].iloc[i]),
            "mrso": end[i, 0] + end[i, 1],
            "gw": end[i, 2] + end[i, 3],
        }
        for i in range(len(forcing))
    ]

    # The checks that keep the derived evaporation honest. A clamp larger than
    # the demand would make `evspsbl` negative, and a clamp in a lower tank
    # would be water from nowhere; neither can be hidden from the budget, so
    # both are counted and reported rather than assumed.
    total = end.sum(axis=1)
    residual = pr_mm[1:] - evaporated_mm[1:] - runoff_mm[1:] - np.diff(total)
    checks = {
        "et_limited_steps": int((clamp[1:, 0] > 0).sum()),
        "unmet_demand_mm": float(clamp[1:, 0].sum()),
        "demand_mm": float(pet_mm[1:].sum()),
        "clamp_beyond_demand_mm": float(np.clip(clamp[1:, 0] - pet_mm[1:], 0, None).sum()),
        "lower_tank_clamp_mm": float(clamp[1:, 1:].sum()),
        "max_row_budget_residual_mm": float(np.abs(residual).max()) if len(residual) else 0.0,
        "upstream_warnings": sum(1 for line in captured.getvalue().splitlines() if line.startswith("WARNING")),
        "tank_range_mm": {
            f"tank{i}": [float(end[:, i].min()), float(end[:, i].max())] for i in range(4)
        },
    }
    return rows, checks


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "pet"):
            row[key] = float(row[key])
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent
    forcing = read_forcing(io_dir / request["input"]["forcing"])
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    timestep = request.get("timestep", "PT1D")
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    dt_days = TIMESTEP_DAYS[timestep]
    rows, checks = simulate(forcing, float(static["area_km2"]), dt_days)

    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "n_steps": len(rows),
        "not_reported": {
            "snw": "the model has no snow store; all precipitation enters the top tank as liquid",
            "canopy": "the model has no interception store",
            "channel": "one Subbasin has no reach; upstream's optional Muskingum routing is not configured (README.md)",
            "hfls, hfss, hfg and other energy terms": "the model has no energy balance",
        },
        "notes": {
            "parameters_per_day": PARAMS_PER_DAY,
            "heights_mm": HEIGHTS_MM,
            "initial_mm": INITIAL_MM,
            "mrso": "top tank + second tank, after the row's outflows",
            "gw": "third tank + fourth tank, after the row's outflows",
            "evspsbl": "demand minus what upstream's non-negative reset of the top tank gave back",
            "seed": "unused; the model is deterministic",
        },
        "checks": checks,
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

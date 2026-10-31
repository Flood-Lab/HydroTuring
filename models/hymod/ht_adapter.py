#!/usr/bin/env python3
"""HydroTuring adapter for HYMOD, from hydromodel v0.4.0 (commit 89d7a8e).

The physics is upstream's, unmodified: `hydromodel/models/hymod.py`. This file
is the units layer and the bookkeeping.

`hymod()` returns the streamflow series but only the final stores, and it
never returns evaporation: the second item it returns is the effective
rainfall. So the reported `mrro` is what `hymod()` itself returned, and the
stores come from driving the same module's `excess` and `linres` row by row.
The two are compared on every run and the adapter refuses to report if they
disagree.

Evaporation is not taken as a residual of the budget, which would close the
budget by construction. `excess` is called twice on each row, once with the
row's demand and once with none; the difference between the two end states is
exactly what its evaporation line removed. An error in its runoff accounting
would not be absorbed by that number.

Stores: `mrso` is the soil store (`x_loss`), `gw` the slow reservoir, `channel`
the three quick-flow reservoirs. HYMOD has no snow or interception store, so
`snw` and `canopy` are identically zero; see `zero_stores` below.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from hydromodel.models.hymod import excess, hymod, linres

COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "mrso", "snw", "canopy", "gw", "channel"]
MODEL = {
    "name": "hymod",
    "version": "0.4.0-89d7a8e.2",
    "source": "OuyangWenyu/hydromodel@89d7a8ed1d72ce4fffbbd9897490b089382ecbac",
}
TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0}
SECONDS_PER_DAY = 86400.0

# The SCE-UA calibration of this same five-parameter form on USGS 01666500,
# from the sceua package's documentation. See README.md for why this set.
PUBLISHED = {"cmax": 273.302, "bexp": 0.402, "alpha": 0.706, "ks": 0.031, "kq": 0.675}

# Stores this model does not have. They hold no water, so they are reported as
# identically zero rather than left out, and the budget probes can ask the question.
ZERO_STORES = {
    "snw": "HYMOD has no snow module: all precipitation enters the soil store as liquid, whatever the temperature",
    "canopy": "HYMOD has no interception store",
}

NOT_REPORTED = {
    "gwex": "HYMOD has no exchange with the outside",
    "stage": "HYMOD has no channel geometry or hydraulics, so a stage would be invented by the adapter",
}

# Largest disagreement tolerated between upstream's series and the row-by-row
# drive of the same functions, in mm per step.
AGREEMENT_MM = 1e-9


def parameters(static: dict, dt_days: float) -> dict:
    """The parameter vector for this catchment, in the units of this step.

    `cmax` is the largest point capacity; the store's catchment-mean capacity
    is cmax / (bexp + 1), so the catchment's capacity sets cmax through it. The
    reservoir coefficients are per-day fractions and become fractions of the
    row, 1 - (1 - k)^dt. Identity at a daily step.
    """
    out = dict(PUBLISHED)
    if "soil_capacity_mm" in static:
        out["cmax"] = float(static["soil_capacity_mm"]) * (out["bexp"] + 1.0)
    out["ks"] = 1.0 - (1.0 - out["ks"]) ** dt_days
    out["kq"] = 1.0 - (1.0 - out["kq"]) ** dt_days
    return out


def run_upstream(pr: np.ndarray, pet: np.ndarray, par: dict) -> np.ndarray:
    """Streamflow in mm per step from upstream's own entry point."""
    p_and_e = np.stack([pr, pet], axis=1)[:, np.newaxis, :]
    vector = np.array([[par["cmax"], par["bexp"], par["alpha"], par["ks"], par["kq"]]])
    flow, _effective_rain = hymod(p_and_e, vector, warmup_length=0, normalized_params=False)
    return flow[:, 0, 0]


def drive(pr: np.ndarray, pet: np.ndarray, par: dict) -> dict:
    """The same loop as `hymod()`, keeping what that loop discards.

    The initial stores are upstream's: all empty.
    """
    cmax, bexp, alpha, ks, kq = par["cmax"], par["bexp"], par["alpha"], par["ks"], par["kq"]
    x_loss, x_slow = 0.0, 0.0
    x_quick = [0.0, 0.0, 0.0]
    n = len(pr)
    out = {name: np.zeros(n) for name in ("flow", "evap", "soil", "slow", "quick")}
    for t in range(n):
        p, e = float(pr[t]), float(pet[t])
        er1, er2, after = excess(x_loss, cmax, bexp, p, e)
        # The store before the evaporation line: the same call without demand.
        _, _, before_evap = excess(x_loss, cmax, bexp, p, 0.0)
        x_loss = float(after)
        effective = float(er1) + float(er2)
        x_slow, slow_flow = linres(x_slow, (1.0 - alpha) * effective, ks)
        inflow = alpha * effective
        for i in range(3):
            x_quick[i], inflow = linres(x_quick[i], inflow, kq)
        out["flow"][t] = slow_flow + inflow
        out["evap"][t] = float(before_evap) - x_loss
        out["soil"][t] = x_loss
        out["slow"][t] = x_slow
        out["quick"][t] = sum(x_quick)
    return out


def simulate(forcing: list[dict], static: dict, timestep: str) -> tuple[list[dict], dict]:
    dt_days = TIMESTEP_DAYS[timestep]
    par = parameters(static, dt_days)
    pr_rate = np.array([row["pr"] for row in forcing], dtype=float)
    pet_rate = np.array([row["pet"] for row in forcing], dtype=float)
    pr, pet = pr_rate * dt_days, pet_rate * dt_days

    flow = run_upstream(pr, pet, par)
    state = drive(pr, pet, par)
    gap = float(np.max(np.abs(flow - state["flow"]))) if len(flow) else 0.0
    if gap > AGREEMENT_MM:
        raise RuntimeError(f"row-by-row drive differs from hymod() by {gap} mm per step")

    area_km2 = float(static["area_km2"])
    rows = []
    for t, step in enumerate(forcing):
        mrro = float(flow[t]) / dt_days
        rows.append({
            "time": step["time"],
            "pr": float(pr_rate[t]),
            "evspsbl": float(state["evap"][t]) / dt_days,
            "mrro": mrro,
            "dis": mrro * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY,
            "mrso": float(state["soil"][t]),
            "snw": 0.0,
            "canopy": 0.0,
            "gw": float(state["slow"][t]),
            "channel": float(state["quick"][t]),
        })
    notes = {
        "parameters_at_this_step": par,
        "parameters_published": PUBLISHED,
        "cmax_source": "soil_capacity_mm * (bexp + 1)" if "soil_capacity_mm" in static else "published",
        "agreement_with_hymod_mm": gap,
        "zero_stores": ZERO_STORES,
        "not_reported": NOT_REPORTED,
    }
    return rows, notes


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
    rows, notes = simulate(forcing, static, timestep)
    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows), "notes": notes}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

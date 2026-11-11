#!/usr/bin/env python3
"""HydroTuring adapter for the Xin'anjiang model, from hydromodel v0.4.0.

The physics is upstream's, unmodified: `hydromodel/models/xaj.py` at commit
89d7a8e, with its defaults (`name="xaj"`, `source_type="sources"`,
`source_book="HF"`). This file is the units layer and the bookkeeping.

`xaj()` returns streamflow and evaporation for every row but only the final
stores. So the reported `mrro` and `evspsbl` are what `xaj()` itself returned,
and the stores come from driving the same module's `generation`, `sources` and
`linear_reservoir` row by row. The two are compared on every run and the
adapter refuses to report if they disagree.

Stores, catchment-average, weighted as upstream weights the runoff they feed:

  mrso     (1 - IM) * (WU + WL + WD + S * FR), the three tension layers and the
           free water, which is S deep over the runoff-producing fraction FR.
  gw       the interflow and groundwater linear reservoirs. Upstream carries
           each as its outflow Q alone; the water it stands for is the storage
           of a linear reservoir, Q * C / (1 - C), which balances upstream's
           recursion exactly.
  channel  the channel routing: runoff inside the lag plus the recession store.

XAJ has no snow or interception store, so `snw` and `canopy` are identically
zero; see `zero_stores` below.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from hydromodel.models.xaj import generation, linear_reservoir, sources, xaj

COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "mrso", "snw", "canopy", "gw", "channel"]
MODEL = {
    "name": "xaj",
    "version": "0.4.0-89d7a8e.3",
    "source": "OuyangWenyu/hydromodel@89d7a8ed1d72ce4fffbbd9897490b089382ecbac",
}
TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0}
SECONDS_PER_DAY = 86400.0

# configs/example_xaj_params.yaml at the pinned commit, the file the proposal
# runs: upstream's "typical calibrated values for demonstration purposes".
PUBLISHED = {
    "K": 0.75, "B": 0.25, "IM": 0.06, "UM": 18.0, "LM": 80.0, "DM": 95.0, "C": 0.12,
    "SM": 45.0, "EX": 1.3, "KI": 0.35, "KG": 0.45, "CS": 0.5, "L": 5.5, "CI": 0.65, "CG": 0.985,
}
# Upstream's order for the parameter vector.
ORDER = ["K", "B", "IM", "UM", "LM", "DM", "C", "SM", "EX", "KI", "KG", "CS", "L", "CI", "CG"]

# Stores this model does not have. They hold no water, so they are reported as
# identically zero rather than left out, and the budget probes can ask the question.
ZERO_STORES = {
    "snw": "XAJ has no snow module: all precipitation is treated as liquid, whatever the temperature",
    "canopy": "XAJ has no interception store",
}

NOT_REPORTED = {
    "gwex": "XAJ has no declared exchange with the outside",
    "stage": "XAJ has no channel geometry or hydraulics, so a stage would be invented by the adapter",
}

# Largest disagreement tolerated between upstream's series and the row-by-row
# drive of the same functions, in mm per step.
AGREEMENT_MM = 1e-9

# Upstream starts both linear reservoirs at an outflow of 0.1, in mm per row
# whatever the row is. Read as mm per day it is the same state at every step;
# left per row, the water it stands for, Q * C / (1 - C), goes from 6.6 mm at a
# daily step to 158 mm at an hourly one, and that water drains into the record.
INITIAL_Q_MM_PER_DAY = 0.1


def parameters(static: dict, dt_days: float) -> dict:
    """The parameter vector for this catchment, in the units of this step.

    The catchment's capacity sets UM + LM + DM, keeping their ratios; the free
    water capacity SM keeps its value, as SAC-SMA's lower free-water stores do
    in `models/sacsma_snow17`.

    Upstream's recession constants are per row and its lag is in rows. Per-day
    constants become C^dt. KI and KG become the per-row pair of upstream's own
    conversion in `sources5mm`. L keeps the whole days upstream's int() leaves
    of it, so the lag is the same length of time at every step. All are the
    identity at a daily step.
    """
    out = dict(PUBLISHED)
    if "soil_capacity_mm" in static:
        factor = float(static["soil_capacity_mm"]) / (out["UM"] + out["LM"] + out["DM"])
        for key in ("UM", "LM", "DM"):
            out[key] *= factor
    if dt_days != 1.0:
        ki, kg = out["KI"], out["KG"]
        out["KI"] = (1.0 - (1.0 - (ki + kg)) ** dt_days) / (1.0 + kg / ki)
        out["KG"] = out["KI"] * kg / ki
        for key in ("CS", "CI", "CG"):
            out[key] = out[key] ** dt_days
        out["L"] = round(int(out["L"]) / dt_days)
    return out


def run_upstream(p_and_e: np.ndarray, par: dict, dt_days: float) -> tuple[np.ndarray, np.ndarray]:
    """Streamflow and evaporation in mm per step from upstream's entry point."""
    vector = np.array([[par[name] for name in ORDER]])
    q0 = INITIAL_Q_MM_PER_DAY * dt_days
    flow, evap = xaj(p_and_e, vector, warmup_length=0, normalized_params=False, name="xaj",
                     source_type="sources", source_book="HF", initial_states={"qi": q0, "qg": q0})
    return flow[:, 0, 0], evap[:, 0, 0]


def drive(p_and_e: np.ndarray, par: dict, dt_days: float) -> dict:
    """The same loop as `xaj()`, keeping what that loop discards.

    The initial stores are upstream's for `warmup_length=0`, with the two
    reservoir outflows in the units of the step.
    """
    a = {name: np.array([par[name]]) for name in ORDER}
    im = a["IM"]
    w = (0.5 * a["UM"], 0.5 * a["LM"], 0.5 * a["DM"])
    s = 0.5 * a["SM"]
    fr = np.full(1, 0.1)
    qi = np.full(1, INITIAL_Q_MM_PER_DAY * dt_days)
    qg = np.full(1, INITIAL_Q_MM_PER_DAY * dt_days)
    n = p_and_e.shape[0]
    out = {name: np.zeros(n) for name in ("evap", "direct", "qi", "qg", "tension", "free")}
    for t in range(n):
        (r, rim, e, pe), w = generation(p_and_e[t, :, :], a["K"], a["B"], im, a["UM"], a["LM"], a["DM"], a["C"], *w)
        (rs, ri, rg), (s, fr) = sources(pe, r, a["SM"], a["EX"], a["KI"], a["KG"], s, fr, book="HF")
        qi = linear_reservoir(ri * (1 - im), a["CI"], qi)
        qg = linear_reservoir(rg * (1 - im), a["CG"], qg)
        out["evap"][t] = e[0]
        out["direct"][t] = (rs * (1 - im) + rim)[0]
        out["qi"][t] = qi[0]
        out["qg"][t] = qg[0]
        out["tension"][t] = (w[0] + w[1] + w[2])[0]
        out["free"][t] = (s * fr)[0]
    return out


def route(inflow: np.ndarray, cs: float, lag: int) -> tuple[np.ndarray, np.ndarray]:
    """Upstream's channel routing, and the water it holds at the end of a row.

    For rows from `lag` on, release is cs * previous release + (1 - cs) * the
    inflow of `lag` rows ago, so the store is the last `lag` inflows plus the
    recession store, which balances to round-off. Upstream passes the first
    `lag` rows straight through and then lags those same rows again, which
    releases them twice; the store is zero while it passes them through and the
    step up on row `lag` is that water. It falls inside the spinup.
    """
    n = len(inflow)
    lag = min(lag, n - 1)  # as xaj() limits it
    released, store = np.zeros(n), np.zeros(n)
    recession = cs / (1.0 - cs)
    for t in range(n):
        if t < lag:
            released[t] = inflow[t]
        else:
            released[t] = cs * released[t - 1] + (1.0 - cs) * inflow[t - lag]
            store[t] = float(np.sum(inflow[t - lag + 1 : t + 1])) + released[t] * recession
    return released, store


def simulate(forcing: list[dict], static: dict, timestep: str) -> tuple[list[dict], dict]:
    dt_days = TIMESTEP_DAYS[timestep]
    par = parameters(static, dt_days)
    pr_rate = np.array([row["pr"] for row in forcing], dtype=float)
    pet_rate = np.array([row["pet"] for row in forcing], dtype=float)
    p_and_e = np.stack([pr_rate * dt_days, pet_rate * dt_days], axis=1)[:, np.newaxis, :]

    flow, evap = run_upstream(p_and_e, par, dt_days)
    state = drive(p_and_e, par, dt_days)
    inflow = state["direct"] + state["qi"] + state["qg"]
    released, channel = route(inflow, par["CS"], int(par["L"]))
    gap = max(float(np.max(np.abs(flow - released))), float(np.max(np.abs(evap - state["evap"]))))
    if gap > AGREEMENT_MM:
        raise RuntimeError(f"row-by-row drive differs from xaj() by {gap} mm per step")

    pervious = 1.0 - par["IM"]
    soil = pervious * (state["tension"] + state["free"])
    ci, cg = par["CI"], par["CG"]
    ground = state["qi"] * ci / (1.0 - ci) + state["qg"] * cg / (1.0 - cg)

    area_km2 = float(static["area_km2"])
    rows = []
    for t, step in enumerate(forcing):
        mrro = float(flow[t]) / dt_days
        rows.append({
            "time": step["time"],
            "pr": float(pr_rate[t]),
            "evspsbl": float(evap[t]) / dt_days,
            "mrro": mrro,
            "dis": mrro * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY,
            "mrso": float(soil[t]),
            "snw": 0.0,
            "canopy": 0.0,
            "gw": float(ground[t]),
            "channel": float(channel[t]),
        })
    notes = {
        "parameters_at_this_step": par,
        "parameters_published": PUBLISHED,
        "capacity_source": "soil_capacity_mm sets UM+LM+DM" if "soil_capacity_mm" in static else "published",
        "options": {"name": "xaj", "source_type": "sources", "source_book": "HF"},
        "agreement_with_xaj_mm": gap,
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
        json.dumps(
            {
                "status": "ok",
                "model": MODEL,
                "n_steps": len(rows),
                "zero_stores": ZERO_STORES,
                "not_reported": NOT_REPORTED,
                "notes": notes,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

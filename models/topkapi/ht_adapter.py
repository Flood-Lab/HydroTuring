#!/usr/bin/env python3
"""HydroTuring adapter for TOPKAPI (sahg/PyTOPKAPI 0.4.0).

The physics is upstream's: `pytopkapi.run` at tag v0.4.0, unmodified, reading
the files its own `pretreatment.read_global_parameters` and
`read_cell_parameters` read and the HDF5 forcing it reads, and writing the HDF5
results its own continuity tests difference. This file writes those inputs and
turns those results into the contract's rows. README.md states every parameter
and choice.

The grid is upstream's smallest: the four cells of its continuity test
`pytopkapi/tests/test_continuity/4cells` (Parak, 2006), one hillslope cell
draining into three channel cells in a line. It is a representative hillslope
and reach and not a map of the stated catchment, so the stated area enters only
`dis`, as in models/wflow_sbm and models/lisflood. Forcing is uniform over it.

Rows. Upstream writes the initial state to row 0 of each result and row t's
outcome to row t+1, so result row t+1 is the end of forcing row t.

Budget. Volumes over the grid's area, in mm:
`mrso` is the soil stores; `channel` is the overland and channel stores, water
generated as runoff and not yet out of the grid; `evspsbl` is upstream's cell
evapotranspiration plus its channel evaporation; `mrro` is everything crossing
the outlet, the channel's flow plus the outlet cell's non-channel drainage,
which are the two outflow terms in upstream's own continuity test.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pytopkapi

COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "mrso", "snw", "canopy", "channel"]
MODEL = {"name": "topkapi", "version": "0.4.0-f7b953e.2", "revision": "f7b953eebcb201403f80970846b429bb1081739c"}

SECONDS = {"PT1D": 86400.0, "PT1H": 3600.0, "PT15M": 900.0, "PT5M": 300.0, "PT1M": 60.0}
GROUP = "ht"

# global_parameters.dat of the 4cells test; Dt is the case's. The grid's cell
# side and the three exponents are upstream's values.
X_M = 1000.0
ALPHA_S = 3.0
ALPHA_O = 1.66666667
ALPHA_C = 1.66666667
A_THRES_M2 = 1000000.0
W_MIN_M = 1.0
W_MAX_M = 10.0

# cell_parameters.dat of the 4cells test, verbatim.
# label x y channel Xc dam tan_beta tan_beta_c L Ks theta_r theta_s n_o n_c
# cell_down pVs_t0 Vo_t0 Qc_t0 kc psi_b lambda
CELLS = """\
0 1 4 0 0 0 0.09 0 0.5 0.001 0.04 0.45 0.1 0.03 1 2.4 0 0 1 146.600006 0.322
1 1 3 1 1000 0 0.08 0.08 0.75 0.001 0.05 0.5 0.2 0.04 2 2.4 0 0.001 1 146.600006 0.322
2 1 2 1 1000 0 0.07 0.07 1 0.001 0.06 0.55 0.3 0.05 3 2.4 0 0.001 1 146.600006 0.322
3 1 1 1 1000 0 0.05 0.05 1.25 0.001 0.07 0.6 0.4 0.06 -99 2.4 0 0.001 1 146.600006 0.322
"""
OUTLET = 3
N_CELLS = 4
GRID_AREA_M2 = N_CELLS * X_M ** 2

# numerical_options of upstream's worked example: Todini's quasi-analytical
# solution of each store's ODE, falling back to Runge-Kutta-Fehlberg where it
# fails. README.md says why not the 4cells test's RKF everywhere (solve = 0).
SOLVER = 1


def grid_capacity_mm() -> float:
    table = np.loadtxt(io.StringIO(CELLS))
    return float(((table[:, 11] - table[:, 10]) * table[:, 8] * 1000.0).mean())


def write_inputs(root: Path, forcing: list[dict], dt_s: float, fac_l: float) -> tuple[Path, Path]:
    (root / "global.dat").write_text(
        "X Dt Alpha_s Alpha_o Alpha_c A_thres W_min W_max\n"
        f"{X_M} {dt_s} {ALPHA_S} {ALPHA_O} {ALPHA_C} {A_THRES_M2} {W_MIN_M} {W_MAX_M}\n"
    )
    (root / "cells.dat").write_text(CELLS)

    # Depths per row, the same over every cell. ETo, upstream's open-water
    # demand on the channel, has no separate source here: both are `pet`.
    days = dt_s / 86400.0
    pr_mm = np.array([row["pr"] for row in forcing], dtype=float) * days
    pet_mm = np.array([row["pet"] for row in forcing], dtype=float) * days
    with h5py.File(root / "rain.h5", "w") as fh:
        fh[f"/{GROUP}/rainfall"] = np.repeat(pr_mm[:, None], N_CELLS, axis=1)
    with h5py.File(root / "et.h5", "w") as fh:
        fh[f"/{GROUP}/ETr"] = np.repeat(pet_mm[:, None], N_CELLS, axis=1)
        fh[f"/{GROUP}/ETo"] = np.repeat(pet_mm[:, None], N_CELLS, axis=1)

    results = root / "results.h5"
    ini = root / "run.ini"
    ini.write_text(
        "[input_files]\n"
        f"file_global_param={root}/global.dat\n"
        f"file_cell_param={root}/cells.dat\n"
        f"file_rain={root}/rain.h5\n"
        f"file_ET={root}/et.h5\n"
        "[output_files]\n"
        f"file_out={results}\n"
        "append_output=False\n"
        "[groups]\n"
        f"group_name={GROUP}\n"
        "[external_flow]\n"
        "external_flow=False\n"
        "[numerical_options]\n"
        f"solve_s={SOLVER}\nsolve_o={SOLVER}\nsolve_c={SOLVER}\n"
        "only_channel_output=False\n"
        "[calib_params]\n"
        f"fac_L={fac_l!r}\nfac_Ks=1.0\nfac_n_o=1.0\nfac_n_c=1.0\n"
    )
    return ini, results


def simulate(forcing: list[dict], static: dict, timestep: str) -> tuple[list[dict], dict]:
    dt_s = SECONDS[timestep]
    days = dt_s / 86400.0
    # Upstream's own multiplier on soil depth carries the catchment's capacity:
    # the grid's mean (theta_s - theta_r) * L becomes `soil_capacity_mm`.
    fac_l = 1.0
    if "soil_capacity_mm" in static:
        fac_l = float(static["soil_capacity_mm"]) / grid_capacity_mm()

    with tempfile.TemporaryDirectory() as tmp:
        ini, results = write_inputs(Path(tmp), forcing, dt_s, fac_l)
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            # Serial, because a process pool is neither needed for four cells
            # nor allowed to depend on the host's core count.
            pytopkapi.run(str(ini), quiet=True, parallel_exec=False)
        with h5py.File(results, "r") as fh:
            vs = fh["/Soil/V_s"][...]
            vo = fh["/Overland/V_o"][...]
            vc = fh["/Channel/V_c"][...]
            qc = fh["/Channel/Qc_out"][...]
            q_down = fh["/Q_down"][...]
            et_mm = fh["/ET_out"][...]
            ec_m3 = fh["/Channel/Ec_out"][...]

    to_mm = 1000.0 / GRID_AREA_M2
    soil = vs.sum(axis=1) * to_mm
    transit = (vo.sum(axis=1) + vc.sum(axis=1)) * to_mm
    # `to_mm` turns a volume over the grid into a depth; ET_out is a depth over
    # each cell, so its catchment mean is the mean over the cells.
    evap = et_mm[1:].mean(axis=1) + ec_m3[1:].sum(axis=1) * to_mm
    out = (qc[1:, OUTLET] + q_down[1:, OUTLET]) * dt_s * to_mm

    area_km2 = float(static["area_km2"])
    rows = [
        {
            "time": forcing[i]["time"],
            "pr": forcing[i]["pr"],
            "evspsbl": evap[i] / days,
            "mrro": out[i] / days,
            # The stated area enters here and nowhere else.
            "dis": out[i] / days * area_km2 / 86.4,
            "mrso": soil[i + 1],
            "snw": 0.0,
            "canopy": 0.0,
            "channel": transit[i + 1],
        }
        for i in range(len(forcing))
    ]

    pr_mm = np.array([row["pr"] for row in forcing], dtype=float) * days
    total = soil + transit
    residual = pr_mm - evap - out - np.diff(total)
    checks = {
        "initial_storage_mm": {"soil": float(soil[0]), "transit": float(transit[0])},
        "grid_capacity_mm_before_fac_L": grid_capacity_mm(),
        "fac_L": fac_l,
        "budget_residual_mm": float(residual.sum()),
        "budget_residual_share_of_rain": float(residual.sum() / pr_mm.sum()) if pr_mm.sum() > 0 else None,
        "max_row_budget_residual_mm": float(np.abs(residual).max()),
        "upstream_lines_on_stdout": len(captured.getvalue().splitlines()),
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

    timestep = request.get("timestep", "PT1H")
    if timestep not in SECONDS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    rows, checks = simulate(forcing, static, timestep)

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
        "zero_stores": {
            "snw": "the model has no snow store; all precipitation is rain",
            "canopy": "the model has no interception store",
        },
        "not_reported": {
            "gw": "the model has no store below its soil store",
            "stage": "the channel depth belongs to the 4 km2 grid, not to the stated area that `dis` is scaled to",
            "hfls, hfss, hfg and other energy terms": "the model has no energy balance",
        },
        "notes": {
            "grid": "upstream's 4cells continuity-test fixture; 4 km2, independent of the stated area",
            "mrso": "soil stores over the grid",
            "channel": "overland and channel stores over the grid",
            "evspsbl": "/ET_out (cell mean) plus /Channel/Ec_out over the grid",
            "mrro": "/Channel/Qc_out plus /Q_down at the outlet cell",
            "dis": "mrro * area_km2 / 86.4",
            "ETo": "upstream's open-water demand is the same `pet` as ETr",
            "seed": "unused; the model is deterministic",
        },
        "checks": checks,
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

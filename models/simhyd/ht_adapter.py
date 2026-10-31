#!/usr/bin/env python3
"""HydroTuring adapter for MARRMoT's SIMHYD (m_18_simhyd_7p_3s), standard library only.

The model is MARRMoT v2.1.2's own MATLAB run under Octave by simhyd_run.m, unmodified.
This file is the units layer and the record of what was fixed before evaluation.

Every process in the model is a rate in mm/d against delta_t in days, so the step is
handed over as delta_t and nothing is rescaled. MARRMoT reports fluxes in mm per step;
they are turned back into mm/day here. Stores are absolute, in mm, at the end of the step.

No calibration exists for the synthetic catchment. The parameters are the midpoints of
the ranges m_18_simhyd_7p_3s.m declares, which is what MARRMoT_model.calibrate and
workflow_example_4.m start from (par_ini = mean(parRanges,2)), and the initial stores are
zero, as in workflow_example_3.m and workflow_example_4.m. The two capacities are taken
from static.json when the case provides them. README.md has the sensitivity.
"""

from __future__ import annotations

import argparse
import array
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

MODEL = {"name": "simhyd", "version": "2.1.2-ht.2"}
COLUMNS = ["time", "pr", "evspsbl", "mrro", "mrso", "snw", "canopy", "gw"]
STEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}
RUNNER = os.environ.get("SIMHYD_RUNNER", "/model/simhyd_run.m")

# The order simhyd_run.m writes them in.
OUTPUTS = [
    "Ea", "Q", "Ei", "EXC", "INF", "INT", "REC", "Et", "GWF", "BAS", "SRUN", "Qt",
    "S1", "S2", "S3", "resnorm", "solver", "iter",
]

# obj.parRanges in m_18_simhyd_7p_3s.m at v2.1.2, in theta's order.
PARAMETER_RANGES = {
    "INSC": (0.0, 5.0),     # maximum interception capacity, mm
    "COEFF": (0.0, 600.0),  # maximum infiltration loss, mm/d
    "SQ": (0.0, 15.0),      # infiltration loss exponent, -
    "SMSC": (1.0, 2000.0),  # maximum soil moisture capacity, mm
    "SUB": (0.0, 1.0),      # interflow proportionality constant, -
    "CRAK": (0.0, 1.0),     # recharge proportionality constant, -
    "K": (0.0, 1.0),        # slow flow time scale, 1/d
}
INITIAL_STORES = {"S1": 0.0, "S2": 0.0, "S3": 0.0}


def parameters(static: dict) -> dict:
    """The parameter vector, and where each value came from."""
    theta = {}
    for name, (low, high) in PARAMETER_RANGES.items():
        theta[name] = {"value": 0.5 * (low + high), "source": "midpoint of MARRMoT parRanges"}
    for name, key in (("INSC", "canopy_capacity_mm"), ("SMSC", "soil_capacity_mm")):
        if static.get(key) is not None:
            theta[name] = {"value": float(static[key]), "source": f"static.json {key}"}
    return theta


def run_simhyd(forcing, dt_days, theta, seed):
    """Rows of OUTPUTS, and what the model printed."""
    climate = array.array("d")
    for row in forcing:
        climate.append(float(row["pr"]))
        climate.append(float(row["pet"]))
        # Carried by MARRMoT's climate structure, read by none of this model's equations.
        climate.append(float(row.get("tas", 0.0) or 0.0))
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        with open(tmp_path / "in.bin", "wb") as fh:
            climate.tofile(fh)
        argv = [
            "octave-cli", "--norc", "--no-history", "--no-window-system", RUNNER,
            str(tmp_path / "in.bin"), str(tmp_path / "out.bin"), repr(dt_days), str(int(seed)),
            ",".join(repr(theta[name]["value"]) for name in PARAMETER_RANGES),
            ",".join(repr(value) for value in INITIAL_STORES.values()),
        ]
        done = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if done.returncode != 0:
            raise SystemExit(
                f"octave exited {done.returncode}: {done.stderr[-2000:]} {done.stdout[-500:]}"
            )
        values = array.array("d")
        size = os.path.getsize(tmp_path / "out.bin")
        with open(tmp_path / "out.bin", "rb") as fh:
            values.fromfile(fh, size // values.itemsize)
    width = len(OUTPUTS)
    rows = [values[i:i + width].tolist() for i in range(0, len(values), width)]
    if len(rows) != len(forcing):
        raise SystemExit(f"simhyd_run wrote {len(rows)} rows for {len(forcing)} forcing rows")
    return rows, done.stdout


def water_budget(rows, forcing, dt_days, stdout):
    """What the reported columns do with the water, and what MARRMoT says it did."""
    index = {name: i for i, name in enumerate(OUTPUTS)}
    rain = sum(float(r["pr"]) * dt_days for r in forcing)
    et = sum(r[index["Ea"]] for r in rows)
    q = sum(r[index["Q"]] for r in rows)
    stored = sum(rows[-1][index[s]] for s in INITIAL_STORES) - sum(INITIAL_STORES.values())
    marrmot = [
        line.split()[1] for line in stdout.splitlines() if line.startswith("MARRMOT_WATER_BALANCE_MM")
    ]
    solver = [r[index["solver"]] for r in rows]
    return {
        "reported_columns": {
            "precipitation_mm": rain,
            "evspsbl_mm": et,
            "mrro_mm": q,
            "store_change_mm": stored,
            "residual_mm": rain - et - q - stored,
        },
        "marrmot_water_balance_mm": float(marrmot[0]) if marrmot else None,
        "solver_steps": {
            "newton_raphson": sum(1 for s in solver if s == 1),
            "fsolve": sum(1 for s in solver if s == 2),
            "lsqnonlin": sum(1 for s in solver if s == 3),
            "worst_resnorm": max(r[index["resnorm"]] for r in rows),
        },
        "smallest_store_mm": {s: min(r[index[s]] for r in rows) for s in INITIAL_STORES},
        "largest_store_mm": {s: max(r[index[s]] for r in rows) for s in INITIAL_STORES},
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    io_dir = request_path.parent
    request = json.loads(request_path.read_text())
    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        forcing = list(csv.DictReader(fh))
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    timestep = request.get("timestep", "PT1D")
    if timestep not in STEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    dt_days = STEP_DAYS[timestep]

    theta = parameters(static)
    rows, stdout = run_simhyd(forcing, dt_days, theta, request.get("seed", 0))

    index = {name: i for i, name in enumerate(OUTPUTS)}
    table = io_dir / request["output"]["table"]
    table.parent.mkdir(parents=True, exist_ok=True)
    with open(table, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(COLUMNS)
        for source, row in zip(forcing, rows):
            writer.writerow([
                source["time"],
                source["pr"],
                f"{row[index['Ea']] / dt_days:.15g}",
                f"{row[index['Q']] / dt_days:.15g}",
                f"{row[index['S2']]:.15g}",
                0.0,
                f"{row[index['S1']]:.15g}",
                f"{row[index['S3']]:.15g}",
            ])

    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "n_steps": len(forcing),
        "delta_t_days": dt_days,
        "theta": theta,
        "initial_stores_mm": INITIAL_STORES,
        "mapping": {
            "mrro": "Q (native Qt = SRUN + INT + BAS)",
            "evspsbl": "Ea (native Ei + Et)",
            "canopy": "S1",
            "mrso": "S2",
            "gw": "S3",
        },
        "zero_stores": {
            "snw": "the model has no snow store; all precipitation enters as rain",
        },
        "not_reported": {
            "channel": "the model has no routing; Qt leaves in the step it is generated",
            "dis": "the model produces a depth of runoff, not a routed discharge",
        },
        "water_budget": water_budget(rows, forcing, dt_days, stdout),
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

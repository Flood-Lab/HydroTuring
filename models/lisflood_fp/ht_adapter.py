#!/usr/bin/env python3
"""HydroTuring adapter for LISFLOOD-FP v8.1 (CPU FV1).

Writes the raster, boundary files and parameter file upstream's executable reads,
runs it, and reads its mass file back. It computes no hydraulics. README.md explains
each choice. Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path

MODEL = {"name": "lisflood_fp", "version": "8.1-ht.3"}
EXECUTABLE = Path("/model/lisflood")
WORK = Path("/tmp/lfp")

STEP_SECONDS = {"PT1H": 3600}
# How often the mass file is written. A row's discharge is the mean of the samples
# inside it, because each sample is one solver step's boundary flux and not a mean.
MASS_INTERVAL_S = 60.0
# Upstream's default for the largest solver step.
MAX_TSTEP_S = 10.0


def grid(static: dict) -> tuple:
    """One row of square cells: as wide as the reach, and fitting its length exactly.

    A raster cell is square, so the cell size is the declared width rounded to the
    nearest value that divides the length into whole cells. The wide-channel solver
    works in discharge per unit width, so this does not change the depth or speed of
    the flow: see `unit_discharge`, and the scaling on the way back in `main`.
    """
    length = float(static["reach_length_m"])
    width = float(static["width_m"])
    target = float(os.environ.get("HT_LFP_DX", width))
    n_cells = max(2, int(round(length / target)))
    return n_cells, length / n_cells


def write_inputs(forcing: list, static: dict, dt: float) -> dict:
    n_cells, dx = grid(static)
    slope = float(static["slope"])
    length = float(static["reach_length_m"])
    width = float(static["width_m"])

    # The bed falls to zero at the outfall. Nothing here reads an absolute level.
    bed = [slope * (length - (i + 0.5) * dx) for i in range(n_cells)]
    with open(WORK / "dem.asc", "w") as fh:
        fh.write(f"ncols {n_cells}\nnrows 1\nxllcorner 0\nyllcorner 0\n")
        fh.write(f"cellsize {dx!r}\nNODATA_value -9999\n")
        fh.write(" ".join(repr(z) for z in bed) + "\n")

    # QVAR on the west face, FREE on the east face. The north and south faces default
    # to closed, which are the reach's banks.
    (WORK / "bc.bci").write_text(f"W 0 {dx!r} QVAR inflow\nE 0 {dx!r} FREE\n")

    # q_in is a discharge and QVAR takes one per unit width, so it is divided by the
    # reach's width and by nothing else. Upstream interpolates linearly between the
    # points it is given, and a row's q_in holds across the row, so each row is two
    # points. The ramp between rows sits at the start of the later one, so no row ever
    # reads a later row's value.
    with open(WORK / "bc.bdy", "w") as fh:
        fh.write("HydroTuring prescribed inflow\ninflow\n")
        fh.write(f"{2 * len(forcing)} seconds\n")
        for i, row in enumerate(forcing):
            unit_discharge = float(row["q_in"]) / width
            fh.write(f"{unit_discharge!r} {i * dt + (1.0e-3 if i else 0.0)!r}\n")
            fh.write(f"{unit_discharge!r} {(i + 1) * dt!r}\n")

    sim_time = len(forcing) * dt
    (WORK / "run.par").write_text(
        "fv1\n"
        "DEMfile dem.asc\n"
        "resroot res\n"
        "dirroot results\n"
        f"sim_time {sim_time!r}\n"
        f"initial_tstep {MAX_TSTEP_S!r}\n"
        f"massint {MASS_INTERVAL_S!r}\n"
        f"saveint {2.0 * sim_time!r}\n"
        f"fpfric {float(static['manning_n'])!r}\n"
        "bcifile bc.bci\n"
        "bdyfile bc.bdy\n"
    )
    return {"n_cells": n_cells, "cell_size_m": dx}


def read_mass(path: Path, n_rows: int, dt: float) -> tuple:
    """Per row: the mean of the outflow samples, and the last sampled volume."""
    totals, counts, volumes = [0.0] * n_rows, [0] * n_rows, [0.0] * n_rows
    with open(path) as fh:
        next(fh)
        for line in fh:
            parts = line.split()
            t, volume, qout = float(parts[0]), float(parts[5]), float(parts[8])
            row = min(n_rows - 1, int(t / dt))
            totals[row] += qout
            counts[row] += 1
            volumes[row] = volume
    empty = [i for i, c in enumerate(counts) if not c]
    if empty:
        raise RuntimeError(f"no mass samples in rows {empty[:5]}")
    return [a / c for a, c in zip(totals, counts)], volumes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent

    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        forcing = list(csv.DictReader(fh))
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    shape = str(static["cross_section_shape"]).strip().lower()
    if shape != "rectangular":
        raise ValueError(f"the adapter maps a rectangular section only, not {shape!r}")
    timestep = request["timestep"]
    dt = float(STEP_SECONDS[timestep])
    times = [datetime.fromisoformat(row["time"]) for row in forcing]
    if any(b - a != timedelta(seconds=dt) for a, b in zip(times, times[1:])):
        raise ValueError(f"forcing rows are not {timestep} apart")

    WORK.mkdir(parents=True, exist_ok=True)
    layout = write_inputs(forcing, static, dt)
    proc = subprocess.run(
        [str(EXECUTABLE), "run.par"], cwd=WORK,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    mass_file = WORK / "results" / "res.mass"
    if proc.returncode != 0 or not mass_file.exists():
        sys.stderr.write(proc.stdout[-4000:])
        raise RuntimeError(f"lisflood exited {proc.returncode}")

    qout, volume = read_mass(mass_file, len(forcing), dt)
    # The solver's fluxes and volume are over a row one cell wide. Per unit width they
    # are the reach's, so they are carried to the declared width.
    width = float(static["width_m"])
    scale = width / layout["cell_size_m"]
    area_m2 = float(static["area_km2"]) * 1.0e6

    columns = ["time", "dis", "channel"]
    output = io_dir / request["output"]["table"]
    with open(output, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        for i, row in enumerate(forcing):
            channel = volume[i] * scale / area_m2 * 1000.0
            writer.writerow([row["time"], repr(qout[i] * scale), repr(channel)])

    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "solver": "fv1",
        "grid": layout,
        "declared_width_m": width,
        "mass_interval_s": MASS_INTERVAL_S,
        "not_reported": {
            # Outside the control volume, which starts at the inflow: unknown, not zero.
            "snw": "this configuration takes q_in and never sees weather, so snow is outside its control volume",
            "canopy": "outside the control volume, as snw",
            "mrso": "outside the control volume, as snw",
            "gw": "outside the control volume, as snw",
        },
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

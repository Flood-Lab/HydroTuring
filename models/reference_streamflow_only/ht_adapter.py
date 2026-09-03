#!/usr/bin/env python3
"""HydroTuring adapter. Standard library only, to show that the /io contract
needs no scientific Python stack and could be written in any language."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

MODEL = {"name": "reference_streamflow_only", "version": "1.0.0"}

COLUMNS = ["time", "mrro", "dis"]


def simulate(forcing, static):
    """Stands in for the published rainfall-runoff models of today.

    It predicts discharge and nothing else. It is not dishonest and it is not
    broken. It simply never says enough about its own water budget to be
    checked, so its verdict is FAIL with reason INCOMPLETE rather than
    VIOLATION. That distinction is the whole reason the report carries a
    reason field next to the verdict.
    """
    area_m2 = static["area_km2"] * 1.0e6
    store = 60.0
    rows = []
    for step in forcing:
        store += step["pr"]
        runoff = 0.06 * store
        store -= runoff
        rows.append({
            "time": step["time"],
            "mrro": runoff,
            "dis": runoff * 1.0e-3 * area_m2 / 86400.0,
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent

    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        forcing = list(csv.DictReader(fh))
    for row in forcing:
        for key in ("pr", "tas", "pet"):
            if key in row:
                row[key] = float(row[key])
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    rows = simulate(forcing, static)

    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows)}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

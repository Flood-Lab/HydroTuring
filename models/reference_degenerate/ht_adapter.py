#!/usr/bin/env python3
"""HydroTuring adapter. Standard library only, to show that the /io contract
needs no scientific Python stack and could be written in any language."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLUMNS = ["time", "pr", "evspsbl", "mrro", "mrso", "snw", "canopy"]

MODEL = {"name": "reference_degenerate", "version": "1.0.0"}


def simulate(forcing, static):
    """The trivial solution to a pure closure test.

    Send every drop of precipitation straight back to the atmosphere, produce
    no runoff, and never change any storage. The water balance closes exactly
    for any forcing whatsoever, so no amount of randomisation touches it.

    It is caught because it is not doing hydrology: the runoff ratio is zero,
    and the implied evaporation exceeds what the atmosphere could ever demand.
    """
    rows = []
    for step in forcing:
        rows.append({
            "time": step["time"],
            "pr": step["pr"],
            "evspsbl": step["pr"],
            "mrro": 0.0,
            "mrso": 100.0,
            "snw": 0.0,
            "canopy": 0.0,
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

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

MODEL = {"name": "reference_cheater", "version": "1.0.0"}


def simulate(forcing, static):
    """Closure by construction, which is the attack randomisation cannot see.

    Evaporation and runoff are picked by simple independent rules that have
    nothing to do with the storage. Storage is then set to whatever makes the
    water balance close exactly, on every step, for every seed, forever.

    The closure criterion cannot catch this: the residual is zero by
    definition. What gives it away is that the invented storage is not a
    storage. It drifts without bound and leaves the range a real soil column
    can occupy, which is why the contract demands absolute states and why
    state_bounds runs alongside closure rather than instead of it.
    """
    et_fraction = 0.55
    runoff_fraction = 0.35

    soil = 0.5 * static["soil_capacity_mm"]
    rows = []

    for step in forcing:
        pr, pet = step["pr"], step["pet"]

        evap = et_fraction * pet
        runoff = runoff_fraction * pr
        soil = soil + pr - evap - runoff  # solve for storage, do no hydrology

        rows.append({
            "time": step["time"],
            "pr": pr,
            "evspsbl": evap,
            "mrro": runoff,
            "mrso": soil,
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

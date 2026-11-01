"""Synthetic positive control for mass/interception-response.

A conserving canopy-less model: intercepted water evaporates the same day and
no canopy store is carried or reported, so the probe must score it rather than
treat it as INCOMPLETE. It exists to exercise the representation the README
allows (interception loss routed straight to evaporation).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

MODEL = {"name": "reference_interception_no_canopy", "version": "1.0.0"}
COLUMNS = ["time", "pr", "evspsbl", "mrro", "mrso"]


def simulate(forcing, static):
    cap = float(static["canopy_capacity_mm"])
    soil = 250.0
    rows = []
    for step in forcing:
        pr = max(float(step["pr"]), 0.0)
        pet = max(float(step["pet"]), 0.0)

        # Interception loss is evaporation in the same step: no store.
        intercept = min(pr, cap, pet)
        through = pr - intercept

        soil += through
        soil_evap = min(
            soil, max(pet - intercept, 0.0) * 0.5 * min(1.0, soil / 100.0)
        )
        soil -= soil_evap
        runoff = max(0.0, soil - 200.0)
        soil -= runoff

        rows.append({
            "time": step["time"],
            "pr": float(step["pr"]),
            "evspsbl": intercept + soil_evap,
            "mrro": runoff,
            "mrso": max(soil, 0.0),
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
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    rows = simulate(forcing, static)
    out_path = io_dir / request["output"]["table"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows)}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

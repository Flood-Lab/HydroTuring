#!/usr/bin/env python3
"""Adapter template.

Only `simulate` is yours. Everything else is contract plumbing that does not
change between models. See AGENTS.md for the full contract.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

MODEL = {"name": "_template", "version": "0.1.0"}

# Must match `emits` in model.yaml, plus "time".
COLUMNS = ["time", "pr", "evspsbl", "mrro", "mrso", "snw", "canopy"]


def simulate(forcing: list[dict], static: dict) -> list[dict]:
    """Run the model over the forcing and return one row per forcing row.

    Three rules that are easy to get wrong:

      1. States (mrso, snw, canopy) are ABSOLUTE storages, never tendencies.
         The harness differences them itself.
      2. Echo `pr` exactly as it was given. Rescaling it fails
         forcing_fidelity, even if the rescaling is a correct unit change.
      3. Emit one row per input row, spinup included. Do not trim the spinup.

    Units: fluxes mm/day, states mm, dis m3/s.
    """
    rows = []
    for step in forcing:
        # TODO: replace this with a call into your model.
        raise NotImplementedError("implement simulate() for your model")
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

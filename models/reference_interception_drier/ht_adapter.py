"""Synthetic negative control for mass/interception-response.

MODE = "drier".  The control declares `canopy_capacity_mm` and consumes it,
but allocates the paired difference the wrong way; it is not physical, and it
exists only to show that `interception_response` catches the named failure.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

MODE = "drier"
MODEL = {"name": "reference_interception_drier", "version": "1.0.0"}
COLUMNS = ["time", "pr", "evspsbl", "mrro", "mrso", "canopy", "gw"]


def simulate(forcing, static):
    declared = float(static["canopy_capacity_mm"])
    # The base bucket runs at a fixed 2 mm capacity, so the declared capacity
    # produces no physical response on its own. The paired difference is then
    # exactly the synthetic allocation below, which is what the control is
    # built to get wrong. `blind` reads it and discards it; the other modes use it as described.
    cap = 2.0

    # This control is specific to mass/interception-response: its scored
    # window starts after 365 spinup days, and the synthetic allocation is
    # placed in the first 100 scored days so it cannot cancel in state0.
    spinup = 365
    adjust_steps = 100

    soil = 250.0
    canopy = 0.0
    gw = 100.0
    rows = []
    for i, step in enumerate(forcing):
        pr = max(float(step["pr"]), 0.0)
        pet = max(float(step["pet"]), 0.0)

        intercepted = min(max(cap - canopy, 0.0), pr)
        canopy += intercepted
        through = pr - intercepted

        canopy_evap = min(canopy, pet)
        canopy -= canopy_evap
        pet_left = max(pet - canopy_evap, 0.0)

        soil += through
        soil_evap = min(soil, pet_left * 0.5 * min(1.0, soil / 100.0))
        soil -= soil_evap
        runoff = max(0.0, soil - 200.0)
        soil -= runoff
        et = canopy_evap + soil_evap

        if spinup <= i < spinup + adjust_steps:
            if MODE == "wetter":
                # A larger canopy increases runoff, and pays for it out of
                # storage. The paired budget still closes; criterion 1 fails.
                runoff += 0.008 * declared
                et += 0.006 * declared
                gw -= 0.014 * declared
            elif MODE == "drier":
                # A larger canopy suppresses evaporation and parks the water
                # in the separate groundwater store. Criterion 2 fails.
                cut = 0.006 * declared
                take = min(et, cut)
                et -= take
                gw += take
            elif MODE == "leaky":
                # Return only part of the water as evaporation and make the
                # rest disappear from the separate store. Criterion 3 fails;
                # the evaporation response is present, so criterion 4 passes.
                et += 0.006 * declared
                gw -= 0.008 * declared

        rows.append({
            "time": step["time"],
            "pr": float(step["pr"]),
            "evspsbl": max(et, 0.0),
            "mrro": max(runoff, 0.0),
            "mrso": max(soil, 0.0),
            "canopy": max(canopy, 0.0),
            "gw": max(gw, 0.0),
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

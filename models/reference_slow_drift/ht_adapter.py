#!/usr/bin/env python3
"""Exact bucket physics with a deliberately drifting storage reporter.

The bucket equations match reference_bucket. Keeping this adapter self-contained
lets the subprocess sandbox run it without importing files outside its model.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

BIAS_MM_PER_DAY = 0.012
COLUMNS = ["time", "pr", "evspsbl", "mrro", "mrso", "snw", "canopy", "channel"]
MODEL = {"name": "reference_slow_drift", "version": "1.0.0"}


def simulate(forcing, static):
    soil_cap = static["soil_capacity_mm"]
    canopy_cap = static["canopy_capacity_mm"]
    soil, swe, canopy = 0.5 * soil_cap, 0.0, 0.0
    rows = []
    for index, step in enumerate(forcing):
        pr, tas, pet = step["pr"], step["tas"], step["pet"]
        snowfall = pr if tas < static["snow_threshold_degC"] else 0.0
        rain = pr - snowfall
        swe += snowfall
        melt = min(swe, static["degree_day_factor_mm_per_C_day"]
                   * max(tas - static["snow_threshold_degC"], 0.0))
        swe -= melt
        intercepted = min(canopy_cap - canopy, rain + melt)
        canopy += intercepted
        throughfall = rain + melt - intercepted
        canopy_evap = min(canopy, pet)
        canopy -= canopy_evap
        pet_left = pet - canopy_evap
        soil += throughfall
        surface = max(0.0, soil - soil_cap)
        soil -= surface
        baseflow = static["baseflow_coefficient"] * soil
        soil -= baseflow
        soil_evap = min(soil, pet_left * min(1.0, soil / (0.5 * soil_cap)))
        soil -= soil_evap
        # The error begins at input day one, including spinup, without feedback.
        rows.append({
            "time": step["time"], "pr": pr,
            "evspsbl": canopy_evap + soil_evap,
            "mrro": surface + baseflow - BIAS_MM_PER_DAY,
            "mrso": soil + BIAS_MM_PER_DAY * (index + 1),
            "snw": swe, "canopy": canopy, "channel": 0.0,
        })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    if request.get("timestep", "PT1D") != "PT1D":
        raise SystemExit("reference_slow_drift supports PT1D only")
    io_dir = request_path.parent
    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        forcing = list(csv.DictReader(fh))
    for row in forcing:
        for key in ("pr", "tas", "pet"):
            row[key] = float(row[key])
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    rows = simulate(forcing, static)
    output = io_dir / request["output"]["table"]
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows)})
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

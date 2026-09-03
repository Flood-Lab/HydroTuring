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

MODEL = {"name": "reference_calendar", "version": "1.0.0"}

# A drift term keyed to the calendar, of the kind a model picks up when the
# year is among its features and its training record happens to trend. It
# conserves water perfectly; it simply believes the date means something.
EPOCH_YEAR = 2000
DRIFT_PER_YEAR = 0.02
DRIFT_RANGE = (0.25, 3.0)

EVAP_SHAPE = 0.5  # soil moisture at which evaporation reaches its potential rate


def drift(timestamp: str) -> float:
    """Recession multiplier for a date, from the year alone."""
    try:
        year = int(str(timestamp)[:4])
    except ValueError:
        return 1.0
    factor = 1.0 + DRIFT_PER_YEAR * (year - EPOCH_YEAR)
    return min(max(factor, DRIFT_RANGE[0]), DRIFT_RANGE[1])


def simulate(forcing, static):
    """The bucket model whose recession depends on what year it thinks it is.

    Water is conserved to floating point and every storage stays inside its
    range, so closure and bounds have nothing to report. The failure is only
    visible by running the same weather twice under different dates, which is
    what an invariance probe does and what nothing else in the suite does.
    """
    soil_cap = static["soil_capacity_mm"]
    canopy_cap = static["canopy_capacity_mm"]
    ddf = static["degree_day_factor_mm_per_C_day"]
    k_base = static["baseflow_coefficient"]
    t_snow = static["snow_threshold_degC"]

    soil = 0.5 * soil_cap
    swe = 0.0
    canopy = 0.0
    rows = []

    for step in forcing:
        pr, tas, pet = step["pr"], step["tas"], step["pet"]

        snowfall = pr if tas < t_snow else 0.0
        rain = 0.0 if tas < t_snow else pr

        swe += snowfall
        melt = min(swe, ddf * max(tas - t_snow, 0.0))
        swe -= melt

        water_in = rain + melt
        intercepted = min(canopy_cap - canopy, water_in)
        canopy += intercepted
        throughfall = water_in - intercepted

        canopy_evap = min(canopy, pet)
        canopy -= canopy_evap
        pet_left = pet - canopy_evap

        soil += throughfall
        surface = max(0.0, soil - soil_cap)
        soil -= surface

        # The only line that differs from reference_bucket.
        baseflow = min(soil, k_base * drift(step["time"]) * soil)
        soil -= baseflow

        soil_evap = min(soil, pet_left * min(1.0, soil / (EVAP_SHAPE * soil_cap)))
        soil -= soil_evap

        rows.append({
            "time": step["time"],
            "pr": pr,
            "evspsbl": canopy_evap + soil_evap,
            "mrro": surface + baseflow,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
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

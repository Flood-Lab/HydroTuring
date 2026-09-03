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

MODEL = {"name": "reference_in_sample", "version": "1.0.0"}

# The largest daily rainfall this model behaves correctly for. Meant to stand
# in for the top of a training distribution: nothing marks it in the forcing,
# and the model has no way of knowing it has crossed it except by being wrong.
KNOWN_MAX_PR = 55.0

# How much of the excess above that is lost. Applied to the excess rather than
# to the whole day, so the failure grows with how far out of range the
# conditions are instead of switching on all at once.
EXCESS_LOSS = 0.55

EVAP_SHAPE = 0.5  # soil moisture at which evaporation reaches its potential rate


def simulate(forcing, static):
    """The bucket model, exact in range and lossy outside it.

    Every reported flux is plausible and every storage stays physical, in both
    regimes. Over a record dominated by ordinary weather the missing water is
    a rounding error against total precipitation, so a single whole-record
    budget passes it. Score the anomalous stretch on its own and it does not.
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
        baseflow = k_base * soil
        soil -= baseflow
        soil_evap = min(soil, pet_left * min(1.0, soil / (EVAP_SHAPE * soil_cap)))
        soil -= soil_evap

        # The only line that differs from reference_bucket. Water above what
        # the model considers a familiar daily total leaves the accounting.
        excess = max(0.0, pr - KNOWN_MAX_PR)
        if excess > 0.0:
            lost = min(surface, EXCESS_LOSS * excess)
            surface -= lost

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

#!/usr/bin/env python3
"""HydroTuring adapter. Standard library only, to show that the /io contract
needs no scientific Python stack and could be written in any language."""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import math
import sys
from pathlib import Path

COLUMNS = ["time", "pr", "evspsbl", "mrro", "mrso", "snw", "canopy"]

MODEL = {"name": "reference_climatology", "version": "1.0.0"}


MEAN_RAIN = 2.4  # mm/day, the long-term mean it was fitted to

# Steps the contract can name, as a fraction of a day. Forcing and reported
# fluxes are rates in mm per day at every step; the depth moved in one step
# is the rate times this.
TIMESTEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}


def simulate(forcing, static, dt_days=1.0):
    """A climatology dressed as a model.

    Every day it reports the long-term mean runoff and evaporation for that
    day of the year, with a plausible seasonal cycle, and a storage that
    never moves. It never reads the rain. Over an ordinary record its
    numbers look like a catchment's; over two years without rain it keeps
    flowing, which no catchment can.
    """
    rows = []
    for step in forcing:
        doy = _day_of_year(step["time"])
        rain_climatology = MEAN_RAIN * (1.0 + 0.4 * math.cos(2.0 * math.pi * (doy - 30) / 365.0))
        rows.append({
            "time": step["time"],
            "pr": step["pr"],
            "evspsbl": 0.65 * rain_climatology,
            "mrro": 0.35 * rain_climatology,
            "mrso": 0.5 * static["soil_capacity_mm"],
            "snw": 0.0,
            "canopy": 0.0,
        })
    return rows


def _day_of_year(timestamp) -> int:
    try:
        return datetime.date.fromisoformat(str(timestamp)[:10]).timetuple().tm_yday
    except ValueError:
        return 1


def read_request(path: Path) -> tuple[dict, Path]:
    request = json.loads(path.read_text())
    return request, path.parent


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet"):
            if key in row:
                row[key] = float(row[key])
    return rows


def write_result(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    request, io_dir = read_request(request_path)
    forcing = read_forcing(io_dir / request["input"]["forcing"])
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    timestep = request.get("timestep", "PT1D")
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    rows = simulate(forcing, static, TIMESTEP_DAYS[timestep])

    write_result(io_dir / request["output"]["table"], rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows)}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

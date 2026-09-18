#!/usr/bin/env python3
"""Steady reach whose gauge silently uses the wrong bed slope."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

MODEL = {"name": "reference_wrong_slope", "version": "1.0.0"}
SECONDS_PER_DAY = 86400.0
WRONG_SLOPE = 0.01


def _capacity(depth: float, width: float, roughness: float, slope: float) -> float:
    area = width * depth
    radius = area / (width + 2.0 * depth)
    return area * radius ** (2.0 / 3.0) * math.sqrt(slope) / roughness


def normal_depth(discharge: float, width: float, roughness: float, slope: float) -> float:
    if discharge <= 0.0:
        return 0.0
    lo, hi = 0.0, 1.0
    while _capacity(hi, width, roughness, slope) < discharge:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _capacity(mid, width, roughness, slope) < discharge:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def simulate(forcing: list[dict], static: dict) -> list[dict]:
    area_km2 = float(static["area_km2"])
    width = float(static["width_m"])
    bed = float(static["bed_elevation_m"])
    roughness = float(static["manning_n"])
    rows = []
    for step in forcing:
        effective = max(float(step["pr"]) - float(step["pet"]), 0.0)
        discharge = effective * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY
        depth = normal_depth(discharge, width, roughness, WRONG_SLOPE)
        rows.append({"time": step["time"], "dis": discharge, "stage": bed + depth})
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
    output = io_dir / request["output"]["table"]
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["time", "dis", "stage"])
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok", "model": MODEL, "n_steps": len(rows),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

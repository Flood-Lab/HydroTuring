#!/usr/bin/env python3
"""Negative routing control: all effective rain exits with zero lag."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


MODEL = {"name": "reference_instant_router", "version": "1.0.0"}
COLUMNS = ["time", "mrro"]
RUNOFF_COEFFICIENT = 0.65
GEOMETRY_KEYS = (
    "area_km2",
    "main_channel_length_km",
    "centroid_channel_length_km",
)
TIMESTEP_DAYS = {"PT1D": 1.0}


def validate_geometry(static: dict) -> tuple[float, float, float]:
    """Validate the same public geometry interface used by the positive model."""
    missing = [name for name in GEOMETRY_KEYS if name not in static]
    if missing:
        raise ValueError(f"missing static catchment attributes: {missing}")
    try:
        area_km2, length_km, centroid_length_km = (
            float(static[name]) for name in GEOMETRY_KEYS
        )
    except (TypeError, ValueError):
        raise ValueError("routing geometry must be numeric") from None
    if not all(
        math.isfinite(value)
        for value in (area_km2, length_km, centroid_length_km)
    ):
        raise ValueError("routing geometry must be finite")
    if area_km2 <= 0.0 or length_km <= 0.0 or centroid_length_km <= 0.0:
        raise ValueError("area and channel lengths must be positive")
    if centroid_length_km > length_km:
        raise ValueError("centroid channel length cannot exceed main-channel length")
    return area_km2, length_km, centroid_length_km


def simulate(
    forcing: list[dict], static: dict, dt_days: float
) -> tuple[list[dict], tuple[float, float, float]]:
    geometry = validate_geometry(static)
    if not math.isfinite(dt_days) or dt_days <= 0.0:
        raise ValueError("timestep must be finite and positive")
    rows = []
    for step in forcing:
        rain_rate = float(step["pr"])
        if not math.isfinite(rain_rate) or rain_rate < 0.0:
            raise ValueError("precipitation must be finite and non-negative")
        rows.append({"time": step["time"], "mrro": RUNOFF_COEFFICIENT * rain_rate})
    return rows, geometry


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent
    timestep = str(request.get("timestep", "PT1D"))
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")

    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        forcing = list(csv.DictReader(fh))
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    rows, geometry = simulate(forcing, static, TIMESTEP_DAYS[timestep])

    output = io_dir / request["output"]["table"]
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    (io_dir / request["output"]["run"]).write_text(
        json.dumps(
            {
                "status": "ok",
                "model": MODEL,
                "n_steps": len(rows),
                "routing": {
                    "method": "none",
                    "lag_hours": 0.0,
                    "validated_geometry": dict(zip(GEOMETRY_KEYS, geometry)),
                },
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

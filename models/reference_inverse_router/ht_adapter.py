#!/usr/bin/env python3
"""Negative control with one deliberate reversal in routing lag by area."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


MODEL = {"name": "reference_inverse_router", "version": "1.0.0"}
COLUMNS = ["time", "mrro"]
RUNOFF_COEFFICIENT = 0.65
SNYDER_SI_CONVERSION = 0.75
SNYDER_CT = 4.0
SNYDER_STANDARD_DURATION_RATIO = 5.5
DESIGN_STORM_DURATION_HOURS = 24.0
GEOMETRY_KEYS = (
    "area_km2",
    "main_channel_length_km",
    "centroid_channel_length_km",
)
TARGET_LAGS_DAYS = (
    (30.0, 1.0),
    (300.0, 2.0),
    (3000.0, 1.0),
    (10000.0, 3.0),
)
TIMESTEP_DAYS = {"PT1D": 1.0}


def validate_geometry(static: dict) -> tuple[float, float, float]:
    """Read and validate every public geometry field used by the probe."""
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


def corrected_snyder_lag_hours(length_km: float, centroid_length_km: float) -> float:
    raw_lag = (
        SNYDER_SI_CONVERSION
        * SNYDER_CT
        * (length_km * centroid_length_km) ** 0.3
    )
    standard_duration = raw_lag / SNYDER_STANDARD_DURATION_RATIO
    return raw_lag - (standard_duration - DESIGN_STORM_DURATION_HOURS) / 4.0


def target_lag_days(area_km2: float) -> float:
    """Select the deliberate 1, 2, 1, 3 day sequence by ascending area."""
    for target_area, lag_days in TARGET_LAGS_DAYS:
        if math.isclose(area_km2, target_area, rel_tol=0.0, abs_tol=1e-9):
            return lag_days
    expected = [area for area, _lag in TARGET_LAGS_DAYS]
    raise ValueError(f"unsupported synthetic area {area_km2:g}; expected one of {expected}")


def discrete_delay(values: list[float], lag_days: float, dt_days: float) -> list[float]:
    """Apply a causal, conservative whole-step delay."""
    if not math.isfinite(lag_days) or lag_days < 0.0:
        raise ValueError("routing lag must be finite and non-negative")
    lag_steps_float = lag_days / dt_days
    lag_steps = int(round(lag_steps_float))
    if not math.isclose(lag_steps_float, lag_steps, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("target routing lag must be an integer number of steps")
    routed = [0.0] * len(values)
    for target in range(lag_steps, len(values)):
        routed[target] = values[target - lag_steps]
    return routed


def simulate(
    forcing: list[dict], static: dict, dt_days: float
) -> tuple[list[dict], dict]:
    area_km2, length_km, centroid_length_km = validate_geometry(static)
    physical_lag_hours = corrected_snyder_lag_hours(length_km, centroid_length_km)
    applied_lag_days = target_lag_days(area_km2)

    effective_depth = []
    for step in forcing:
        rain_rate = float(step["pr"])
        if not math.isfinite(rain_rate) or rain_rate < 0.0:
            raise ValueError("precipitation must be finite and non-negative")
        effective_depth.append(RUNOFF_COEFFICIENT * rain_rate * dt_days)
    routed_depth = discrete_delay(effective_depth, applied_lag_days, dt_days)
    rows = [
        {"time": step["time"], "mrro": routed_depth[i] / dt_days}
        for i, step in enumerate(forcing)
    ]
    routing = {
        "method": "deliberately reversed discrete Snyder-Hack lag sequence",
        "physical_lag_hours": physical_lag_hours,
        "applied_lag_hours": 24.0 * applied_lag_days,
        "runoff_coefficient": RUNOFF_COEFFICIENT,
    }
    return rows, routing


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
    rows, routing = simulate(forcing, static, TIMESTEP_DAYS[timestep])

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
                "routing": routing,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

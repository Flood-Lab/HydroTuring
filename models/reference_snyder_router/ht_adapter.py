#!/usr/bin/env python3
"""Snyder-Hack positive reference for momentum/routing-lag-consistency.

A fixed share of gross rainfall becomes effective rainfall. It is routed by a
causal triangular unit hydrograph whose peak follows the public geometry's
duration-corrected Snyder lag from the effective-rainfall centroid. Integrating
the triangle over model steps gives a non-negative discrete kernel that sums
to one.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


MODEL = {"name": "reference_snyder_router", "version": "1.0.0"}
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


def snyder_lag_hours(
    static: dict, *, ct: float = SNYDER_CT
) -> tuple[float, float, float]:
    """Return corrected lag, raw lag and standard rain duration, in hours."""
    _area_km2, length_km, centroid_length_km = validate_geometry(static)
    if not math.isfinite(ct) or ct <= 0.0:
        raise ValueError("Snyder Ct must be finite and positive")
    raw_lag = (
        SNYDER_SI_CONVERSION
        * ct
        * (length_km * centroid_length_km) ** 0.3
    )
    standard_duration = raw_lag / SNYDER_STANDARD_DURATION_RATIO
    corrected_lag = raw_lag - (
        standard_duration - DESIGN_STORM_DURATION_HOURS
    ) / 4.0
    if not math.isfinite(corrected_lag) or corrected_lag <= 0.0:
        raise ValueError("Snyder geometry produced a non-positive corrected lag")
    return corrected_lag, raw_lag, standard_duration


def triangular_unit_hydrograph(mode_steps: float) -> list[float]:
    """Integrate a unit-area symmetric triangle over regular model steps.

    The continuous support is ``[0, 2 * mode_steps]`` and its mode is exactly
    ``mode_steps``. Bin integration preserves the sub-step mode without using
    a two-cell fractional-delay approximation.
    """
    if not math.isfinite(mode_steps) or mode_steps <= 0.0:
        raise ValueError("unit-hydrograph mode must be finite and positive")

    def cdf(position: float) -> float:
        if position <= 0.0:
            return 0.0
        if position < mode_steps:
            return 0.5 * (position / mode_steps) ** 2
        if position < 2.0 * mode_steps:
            return 1.0 - 0.5 * ((2.0 * mode_steps - position) / mode_steps) ** 2
        return 1.0

    n_bins = max(1, int(math.ceil(2.0 * mode_steps)))
    weights = [cdf(float(i + 1)) - cdf(float(i)) for i in range(n_bins)]
    total = math.fsum(weights)
    if total <= 0.0 or any(weight < -1e-15 for weight in weights):
        raise ValueError("invalid triangular unit-hydrograph weights")
    return [max(0.0, weight) / total for weight in weights]


def route(values: list[float], kernel: list[float]) -> list[float]:
    """Apply a causal unit-hydrograph kernel on the supplied regular grid."""
    routed = [0.0] * len(values)
    for source, value in enumerate(values):
        for offset, weight in enumerate(kernel):
            target = source + offset
            if target >= len(routed):
                break
            routed[target] += value * weight
    return routed


def simulate(
    forcing: list[dict],
    static: dict,
    dt_days: float,
    *,
    ct: float = SNYDER_CT,
) -> tuple[list[dict], dict]:
    corrected_lag, raw_lag, standard_duration = snyder_lag_hours(static, ct=ct)
    effective_depth = []
    for step in forcing:
        rain_rate = float(step["pr"])
        if not math.isfinite(rain_rate) or rain_rate < 0.0:
            raise ValueError("precipitation must be finite and non-negative")
        effective_depth.append(RUNOFF_COEFFICIENT * rain_rate * dt_days)

    # Snyder lag starts at the centroid of the excess-rainfall block.  The
    # kernel, however, is indexed from the beginning of that block.  For this
    # one-day rectangular design storm those origins differ by 12 hours.
    rainfall_centroid_offset = 0.5 * DESIGN_STORM_DURATION_HOURS
    kernel_mode_hours = corrected_lag + rainfall_centroid_offset
    kernel = triangular_unit_hydrograph(kernel_mode_hours / (24.0 * dt_days))
    routed_depth = route(effective_depth, kernel)
    rows = [
        {"time": step["time"], "mrro": routed_depth[i] / dt_days}
        for i, step in enumerate(forcing)
    ]
    routing = {
        "method": "causal conservative triangular Snyder unit hydrograph",
        "raw_snyder_lag_hours": raw_lag,
        "standard_excess_rain_duration_hours": standard_duration,
        "design_storm_duration_hours": DESIGN_STORM_DURATION_HOURS,
        "snyder_ct": ct,
        "lag_hours": corrected_lag,
        "rainfall_centroid_offset_hours": rainfall_centroid_offset,
        "kernel_mode_from_storm_start_hours": kernel_mode_hours,
        "kernel_weights": kernel,
        "kernel_sum": math.fsum(kernel),
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

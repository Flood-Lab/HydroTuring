"""Steady uniform-flow consistency for a rectangular open channel.

At steady, uniform flow the Saint-Venant momentum equation reduces to one
local balance: the friction slope equals the bed slope.  A model that reports
both discharge and stage therefore has to make those two diagnostics belong
to the same reach geometry and roughness.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL,
    PASS,
    CriterionResult,
    criterion,
    make_window,
    segments,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


def _failure(message: str, **diagnostics) -> CriterionResult:
    return CriterionResult(
        "uniform_flow_friction",
        FAIL,
        message,
        diagnostics=diagnostics,
    )


def _cv(values: np.ndarray) -> float:
    mean = float(np.mean(values))
    if not np.isfinite(mean) or abs(mean) <= np.finfo(float).eps:
        return float("inf")
    return float(np.std(values) / abs(mean))


def _relative_quarter_shift(values: np.ndarray) -> float:
    """First-to-last-quarter movement, normalized by the block mean."""
    quarter = max(1, len(values) // 4)
    scale = abs(float(np.mean(values)))
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        return float("inf")
    return float(
        abs(float(np.mean(values[-quarter:])) - float(np.mean(values[:quarter])))
        / scale
    )


@criterion("uniform_flow_friction")
def uniform_flow_friction(
    run: RunResult,
    probe: ProbeSpec,
    params: dict,
) -> CriterionResult:
    """Require Manning friction slope to reproduce the declared bed slope.

    For a rectangular section,

        d = stage - z_b
        A = w d
        R_h = A / (w + 2 d)
        S_f = (n Q / (A R_h^(2/3)))^2.

    The verdict is the mean absolute value of ``S_f / S_0 - 1`` over the
    declared final steady block.  A separate response gate refuses to score a
    block whose discharge or depth is still moving, rather than letting a
    transient hydrograph masquerade as uniform flow.
    """
    window = make_window(run, probe)
    required = ("dis", "stage")
    missing = [name for name in required if name not in window.table]
    if missing:
        return _failure(f"missing required output(s): {', '.join(missing)}")

    static = run.case.static
    static_names = (
        "width_m",
        "bed_elevation_m",
        "slope",
        "manning_n",
        "cross_section_shape",
    )
    absent = [name for name in static_names if name not in static]
    if absent:
        raise ValueError(f"case is missing static value(s): {', '.join(absent)}")

    shape = str(static["cross_section_shape"]).strip().lower()
    if shape != "rectangular":
        raise ValueError(
            "uniform_flow_friction requires cross_section_shape='rectangular', "
            f"not {static['cross_section_shape']!r}"
        )

    try:
        width = float(static["width_m"])
        bed = float(static["bed_elevation_m"])
        slope = float(static["slope"])
        roughness = float(static["manning_n"])
    except (TypeError, ValueError):
        raise ValueError("reach geometry and roughness must be numeric") from None

    geometry = np.asarray([width, bed, slope, roughness], dtype=float)
    if not np.isfinite(geometry).all():
        raise ValueError("reach geometry and roughness must be finite")
    if width <= 0.0 or slope <= 0.0 or roughness <= 0.0:
        raise ValueError("width_m, slope and manning_n must be positive")

    steady_days = float(params.get("steady_days", 90.0))
    if not np.isfinite(steady_days) or steady_days <= 0.0:
        raise ValueError("steady_days must be positive")
    steady_steps = max(1, int(round(steady_days / window.dt_days)))
    max_cv = float(params.get("max_cv", 0.01))
    max_shift = float(params.get("max_relative_trend", 0.01))
    tolerance = float(params.get("tolerance", 0.05))
    if (
        not np.isfinite(max_cv)
        or max_cv < 0.0
        or not np.isfinite(max_shift)
        or max_shift < 0.0
        or not np.isfinite(tolerance)
        or tolerance < 0.0
    ):
        raise ValueError("criterion tolerances must be finite and non-negative")

    expected = [str(value) for value in params.get(
        "plateaus", ["low", "medium", "high"]
    )]
    blocks = segments(window, "_plateau")
    by_label = {label: (start, stop) for label, start, stop in blocks}
    missing_plateaus = [label for label in expected if label not in by_label]
    if missing_plateaus:
        raise ValueError(
            "uniform-flow generator is missing plateau(s): "
            + ", ".join(missing_plateaus)
        )

    failures: list[str] = []
    plateau_diagnostics: dict[str, dict[str, float | int]] = {}
    mean_errors: list[float] = []
    p95_errors: list[float] = []
    max_errors: list[float] = []
    for label in expected:
        start, stop = by_label[label]
        if stop - start < steady_steps:
            raise ValueError(
                f"plateau {label!r} needs {steady_steps} rows, but the generator "
                f"supplied {stop - start}"
            )
        block = window.table.iloc[stop - steady_steps:stop]
        discharge = np.asarray(block["dis"], dtype=float)
        stage = np.asarray(block["stage"], dtype=float)
        depth = stage - bed
        if not np.isfinite(discharge).all() or not np.isfinite(depth).all():
            return _failure(
                f"{label} plateau: discharge and stage must be finite",
                plateau=label,
            )
        if np.any(discharge <= 0.0):
            return _failure(
                f"{label} plateau is wet, but reported discharge is non-positive",
                plateau=label,
                min_discharge_m3s=float(np.min(discharge)),
            )
        if np.any(depth <= 0.0):
            return _failure(
                f"{label} plateau is wet, but stage is at or below the bed",
                plateau=label,
                min_depth_m=float(np.min(depth)),
            )

        q_cv = _cv(discharge)
        depth_cv = _cv(depth)
        q_shift = _relative_quarter_shift(discharge)
        depth_shift = _relative_quarter_shift(depth)
        if max(q_cv, depth_cv) > max_cv or max(q_shift, depth_shift) > max_shift:
            failures.append(
                f"{label} is not steady (CV Q/depth {q_cv:.2%}/{depth_cv:.2%}; "
                f"quarter shift {q_shift:.2%}/{depth_shift:.2%})"
            )

        area = width * depth
        hydraulic_radius = area / (width + 2.0 * depth)
        friction_slope = (
            roughness * discharge / (area * hydraulic_radius ** (2.0 / 3.0))
        ) ** 2
        velocity = discharge / area
        froude = velocity / np.sqrt(9.80665 * depth)
        normalized = np.abs(friction_slope - slope) / slope
        mean_error = float(np.mean(normalized))
        p95_error = float(np.percentile(normalized, 95))
        max_error = float(np.max(normalized))
        if mean_error > tolerance:
            failures.append(
                f"{label} residual {mean_error:.2%} exceeds {tolerance:.2%}"
            )
        mean_errors.append(mean_error)
        p95_errors.append(p95_error)
        max_errors.append(max_error)
        plateau_diagnostics[label] = {
            "mean_absolute_normalized_residual": mean_error,
            "p95_absolute_normalized_residual": p95_error,
            "max_absolute_normalized_residual": max_error,
            "discharge_cv": q_cv,
            "depth_cv": depth_cv,
            "discharge_quarter_shift": q_shift,
            "depth_quarter_shift": depth_shift,
            "min_depth_m": float(np.min(depth)),
            "max_depth_m": float(np.max(depth)),
            "mean_friction_slope": float(np.mean(friction_slope)),
            "mean_froude_number": float(np.mean(froude)),
            "max_froude_number": float(np.max(froude)),
            "steady_steps": steady_steps,
        }

    worst_mean = max(mean_errors)
    summary = ", ".join(
        f"{label} {plateau_diagnostics[label]['mean_absolute_normalized_residual']:.2%}"
        for label in expected
    )
    status = FAIL if failures else PASS
    message = (
        "; ".join(failures)
        if failures
        else f"three steady plateaus satisfy |S_f - S_0| / S_0: {summary} "
             f"(limit {tolerance:.2%})"
    )
    return CriterionResult(
        "uniform_flow_friction",
        status,
        message,
        value=worst_mean,
        threshold=tolerance,
        diagnostics={
            "worst_mean_absolute_normalized_residual": worst_mean,
            "worst_p95_absolute_normalized_residual": max(p95_errors),
            "worst_max_absolute_normalized_residual": max(max_errors),
            "bed_slope": slope,
            "steady_steps": steady_steps,
            "plateaus": plateau_diagnostics,
        },
    )

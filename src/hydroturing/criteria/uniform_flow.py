"""Steady uniform-flow consistency for a rectangular open channel.

At steady, uniform flow the Saint-Venant momentum equation reduces to one
local balance: the friction slope equals the bed slope.  A model that reports
both discharge and stage therefore has to make those two diagnostics belong
to the same reach geometry and roughness.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
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
    static_names = ("width_m", "bed_elevation_m", "slope", "manning_n")
    absent = [name for name in static_names if name not in static]
    if absent:
        return _failure(f"case is missing static value(s): {', '.join(absent)}")

    try:
        width = float(static["width_m"])
        bed = float(static["bed_elevation_m"])
        slope = float(static["slope"])
        roughness = float(static["manning_n"])
    except (TypeError, ValueError):
        return _failure("reach geometry and roughness must be numeric")

    geometry = np.asarray([width, bed, slope, roughness], dtype=float)
    if not np.isfinite(geometry).all():
        return _failure("reach geometry and roughness must be finite")
    if width <= 0.0 or slope <= 0.0 or roughness <= 0.0:
        return _failure("width_m, slope and manning_n must be positive")

    steady_days = float(params.get("steady_days", 365.0))
    if not np.isfinite(steady_days) or steady_days <= 0.0:
        return _failure("steady_days must be positive")
    steady_steps = max(1, int(round(steady_days / window.dt_days)))
    if len(window.table) < steady_steps:
        return _failure(
            f"steady block needs {steady_steps} rows, but the scored window has "
            f"{len(window.table)}"
        )

    block = window.table.iloc[-steady_steps:]
    discharge = np.asarray(block["dis"], dtype=float)
    stage = np.asarray(block["stage"], dtype=float)
    depth = stage - bed
    if not np.isfinite(discharge).all() or not np.isfinite(depth).all():
        return _failure("discharge and stage must be finite on the steady block")
    if np.any(discharge <= 0.0):
        return _failure(
            "the generated reach is wet, but reported discharge is non-positive",
            min_discharge_m3s=float(np.min(discharge)),
        )
    if np.any(depth <= 0.0):
        return _failure(
            "the generated reach is wet, but stage is at or below the bed",
            min_depth_m=float(np.min(depth)),
        )

    max_cv = float(params.get("max_cv", 0.01))
    q_cv = _cv(discharge)
    depth_cv = _cv(depth)
    if q_cv > max_cv or depth_cv > max_cv:
        return _failure(
            f"final block is not steady: discharge CV {q_cv:.3%}, depth CV "
            f"{depth_cv:.3%} (limit {max_cv:.3%})",
            discharge_cv=q_cv,
            depth_cv=depth_cv,
        )

    area = width * depth
    hydraulic_radius = area / (width + 2.0 * depth)
    friction_slope = (
        roughness * discharge / (area * hydraulic_radius ** (2.0 / 3.0))
    ) ** 2
    normalized = np.abs(friction_slope - slope) / slope
    mean_error = float(np.mean(normalized))
    p95_error = float(np.percentile(normalized, 95))
    max_error = float(np.max(normalized))
    tolerance = float(params.get("tolerance", 0.05))
    if not np.isfinite(tolerance) or tolerance < 0.0:
        return _failure("tolerance must be finite and non-negative")

    status = PASS if mean_error <= tolerance else FAIL
    message = (
        f"mean |S_f - S_0| / S_0 is {mean_error:.2%} "
        f"(limit {tolerance:.2%}); p95 {p95_error:.2%}, max {max_error:.2%}"
    )
    return CriterionResult(
        "uniform_flow_friction",
        status,
        message,
        value=mean_error,
        threshold=tolerance,
        diagnostics={
            "mean_absolute_normalized_residual": mean_error,
            "p95_absolute_normalized_residual": p95_error,
            "max_absolute_normalized_residual": max_error,
            "discharge_cv": q_cv,
            "depth_cv": depth_cv,
            "min_depth_m": float(np.min(depth)),
            "max_depth_m": float(np.max(depth)),
            "mean_friction_slope": float(np.mean(friction_slope)),
            "bed_slope": slope,
            "steady_steps": steady_steps,
        },
    )

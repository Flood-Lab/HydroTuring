"""Transient flood-wave celerity against a Manning/kinematic expectation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


class _ResponseFailure(ValueError):
    """A model answer that cannot demonstrate a measurable routing response."""

    def __init__(self, variant: str, message: str, diagnostics: dict[str, Any] | None = None):
        super().__init__(message)
        self.variant = variant
        self.diagnostics = {"variant": variant, **(diagnostics or {})}


@dataclass(frozen=True)
class _Measurement:
    state: str
    c_obs: float
    c_kin: float
    residual: float
    dt_hours: float
    short_centroid_h: float
    long_centroid_h: float
    short_variance_h2: float
    long_variance_h2: float
    prescribed_base_m3s: float
    short_base_m3s: float
    long_base_m3s: float
    short_base_variation: float
    long_base_variation: float
    diffusivity_obs_m2_s: float | None


def _positive(static: dict[str, Any], key: str, variant: str) -> float:
    try:
        value = float(static[key])
    except (KeyError, TypeError, ValueError):
        raise ValueError(f"variant '{variant}' needs numeric static '{key}'") from None
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"variant '{variant}' needs positive finite static '{key}'")
    return value


def _normal_depth(q: float, width: float, slope: float, manning_n: float) -> float:
    """Invert rectangular Manning discharge by monotone bisection."""
    if q <= 0:
        raise ValueError("wave celerity needs positive base discharge")

    def discharge(h: float) -> float:
        area = width * h
        perimeter = width + 2.0 * h
        radius = area / perimeter
        return (1.0 / manning_n) * area * radius ** (2.0 / 3.0) * slope ** 0.5

    lo, hi = 1.0e-6, 1.0
    while discharge(hi) < q:
        hi *= 2.0
        if hi > 1.0e4:
            raise ValueError("could not bracket Manning normal depth")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if discharge(mid) < q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _kinematic_celerity(q: float, static: dict[str, Any], variant: str) -> float:
    shape = str(static.get("cross_section_shape", "")).strip().lower()
    if shape != "rectangular":
        raise ValueError(
            f"variant '{variant}' needs cross_section_shape='rectangular'"
        )
    width = _positive(static, "width_m", variant)
    slope = _positive(static, "slope", variant)
    manning_n = _positive(static, "manning_n", variant)
    h = _normal_depth(q, width, slope, manning_n)
    area = width * h
    perimeter = width + 2.0 * h
    # Q(A) = n^-1 A (A/P)^(2/3) S^(1/2).
    # Hence dQ/dA = Q [5/(3A) - 4/(3 b P)] for a rectangular section.
    celerity = q * (5.0 / (3.0 * area) - 4.0 / (3.0 * width * perimeter))
    if not np.isfinite(celerity) or celerity <= 0:
        raise ValueError(f"variant '{variant}' produced invalid kinematic celerity")
    return float(celerity)


def _centroid(
    run: RunResult,
    probe: ProbeSpec,
    params: dict,
    variant: str,
    event_column: str,
) -> tuple[float, float, float, float, float]:
    window = make_window(run, probe)
    if "dis" not in window.table:
        raise _ResponseFailure(variant, f"variant '{variant}' does not report dis")
    discharge = pd.to_numeric(window.table["dis"], errors="coerce").to_numpy(float)
    if len(discharge) != len(window.forcing):
        raise _ResponseFailure(
            variant,
            f"variant '{variant}' returned {len(discharge)} discharge rows for "
            f"{len(window.forcing)} forcing rows",
            {"discharge_rows": len(discharge), "forcing_rows": len(window.forcing)},
        )
    nonfinite = int((~np.isfinite(discharge)).sum())
    if nonfinite:
        raise _ResponseFailure(
            variant,
            f"variant '{variant}' has {nonfinite} non-finite discharge values",
            {"nonfinite_count": nonfinite},
        )

    marker = event_column
    if marker not in window.forcing:
        raise ValueError(f"wave celerity needs hidden forcing column '{marker}'")
    event = pd.to_numeric(window.forcing[marker], errors="coerce").to_numpy(float)
    indices = np.flatnonzero(event > 0)
    if not len(indices):
        raise ValueError("wave-celerity scored window has no pulse")
    first, last = int(indices[0]), int(indices[-1])
    if np.any(np.diff(indices) != 1):
        raise ValueError("wave-celerity pulse must be contiguous")

    baseline_hours = float(params.get("baseline_hours", 12.0))
    baseline_steps = max(2, int(round(baseline_hours / (window.dt_days * 24.0))))
    if first < baseline_steps:
        raise ValueError("wave-celerity pulse has too little pre-event baseline")
    base_slice = discharge[first - baseline_steps:first]
    base = float(np.mean(base_slice))

    # q_in is the externally prescribed hydraulic operating point.  The
    # theory target must not be chosen by the model under test: first verify
    # that the routed steady discharge actually represents that prescribed
    # state, then evaluate dQ/dA at q_in rather than at a model-selected Q.
    if "q_in" not in window.forcing:
        raise ValueError("wave celerity needs prescribed forcing column 'q_in'")
    inflow = pd.to_numeric(window.forcing["q_in"], errors="coerce").to_numpy(float)
    prescribed_slice = inflow[first - baseline_steps:first]
    if np.any(~np.isfinite(prescribed_slice)):
        raise ValueError("wave-celerity prescribed q_in baseline is non-finite")
    prescribed_base = float(np.mean(prescribed_slice))
    if not np.isfinite(prescribed_base) or prescribed_base <= 0.0:
        raise ValueError("wave-celerity prescribed q_in baseline must be positive")
    input_span = float(np.max(np.abs(prescribed_slice - prescribed_base)))
    if input_span > 1.0e-10 * max(1.0, abs(prescribed_base)):
        raise ValueError("wave-celerity prescribed q_in must be steady before the pulse")

    base_tolerance = float(params.get("base_flow_relative_tolerance", 0.01))
    if not np.isfinite(base_tolerance) or base_tolerance < 0.0:
        raise ValueError("base_flow_relative_tolerance must be finite and non-negative")
    base_error = abs(base - prescribed_base) / prescribed_base
    base_variation = float(np.max(np.abs(base_slice - base))) / prescribed_base
    if base_error > base_tolerance + 1.0e-12:
        raise _ResponseFailure(
            variant,
            f"variant '{variant}' settled base discharge {base:.6g} m3/s "
            f"does not match prescribed q_in {prescribed_base:.6g} m3/s "
            f"within {100*base_tolerance:.1f}%",
            {
                "baseline_discharge_m3s": base,
                "prescribed_base_inflow_m3s": prescribed_base,
                "base_relative_error": base_error,
                "base_flow_relative_tolerance": base_tolerance,
            },
        )
    if base_variation > base_tolerance + 1.0e-12:
        raise _ResponseFailure(
            variant,
            f"variant '{variant}' is not settled before the pulse "
            f"(baseline variation {100*base_variation:.2f}% of q_in; "
            f"limit {100*base_tolerance:.1f}%)",
            {
                "baseline_discharge_m3s": base,
                "prescribed_base_inflow_m3s": prescribed_base,
                "base_relative_variation": base_variation,
                "base_flow_relative_tolerance": base_tolerance,
            },
        )

    response_hours = float(params.get("response_hours", 72.0))
    if not np.isfinite(response_hours) or response_hours <= 0:
        raise ValueError("wave-celerity response_hours must be finite and positive")
    response_steps = max(1, int(np.ceil(response_hours / (window.dt_days * 24.0))))
    response_stop = min(len(discharge), first + response_steps)
    response = np.zeros_like(discharge)
    response[first:response_stop] = np.clip(
        discharge[first:response_stop] - base, 0.0, None
    )
    total = float(response.sum())
    min_fraction = float(params.get("min_response_fraction", 1.0e-4))
    pulse_depth = float(np.sum(event) * window.dt_days)
    if pulse_depth <= 0:
        raise ValueError("wave-celerity pulse has zero integrated depth")
    if total <= 0:
        raise _ResponseFailure(
            variant,
            f"variant '{variant}' has no measurable positive response",
            {"baseline_discharge_m3s": base, "response_peak_m3s": 0.0},
        )
    # Use response magnitude relative to the base flow as a simple non-degeneracy guard.
    response_peak = float(response.max())
    response_ratio = response_peak / max(abs(base), 1.0e-12)
    if response_ratio < min_fraction:
        raise _ResponseFailure(
            variant,
            f"variant '{variant}' response is too small to time "
            f"({response_ratio:.3g} of base; minimum {min_fraction:g})",
            {
                "baseline_discharge_m3s": base,
                "response_peak_m3s": response_peak,
                "response_to_base_ratio": response_ratio,
                "min_response_fraction": min_fraction,
            },
        )

    # The first-moment identity is valid only for a contained response.  A
    # slowly decaying tail truncated by response_hours would bias the centroid,
    # so refuse to emit a celerity when the positive excess at the window edge
    # remains material relative to the event peak.
    max_tail_fraction = float(params.get("max_tail_fraction", 0.01))
    if not np.isfinite(max_tail_fraction) or max_tail_fraction < 0.0:
        raise ValueError("max_tail_fraction must be finite and non-negative")
    tail_excess = max(float(discharge[response_stop - 1] - base), 0.0)
    tail_fraction = tail_excess / max(response_peak, 1.0e-12)
    if tail_fraction > max_tail_fraction + 1.0e-12:
        raise _ResponseFailure(
            variant,
            f"variant '{variant}' response is not contained in the response window "
            f"(tail {tail_fraction:.3g} of peak; maximum {max_tail_fraction:g})",
            {
                "baseline_discharge_m3s": base,
                "response_peak_m3s": response_peak,
                "response_tail_m3s": tail_excess,
                "response_tail_fraction": tail_fraction,
                "max_tail_fraction": max_tail_fraction,
                "response_hours": response_hours,
            },
        )

    centres_h = (np.arange(len(response), dtype=float) + 0.5) * window.dt_days * 24.0
    centroid_h = float(np.dot(centres_h, response) / total)
    variance_h2 = float(np.dot((centres_h - centroid_h) ** 2, response) / total)
    return centroid_h, variance_h2, base, prescribed_base, base_variation


def _pair(runs: dict[str, RunResult], probe: ProbeSpec, params: dict, state: str) -> _Measurement:
    short_name, long_name = "short", "long"
    if short_name not in runs or long_name not in runs:
        raise ValueError("wave celerity needs short and long variants")
    short, long = runs[short_name], runs[long_name]

    # Within a hydraulic state, forcing and every static property except length
    # must be byte-identical.
    if not short.case.forcing.reset_index(drop=True).equals(
        long.case.forcing.reset_index(drop=True)
    ):
        raise ValueError(f"{state} short/long variants must have identical forcing")
    keys = set(short.case.static) | set(long.case.static)
    for key in keys - {"reach_length_m"}:
        if short.case.static.get(key) != long.case.static.get(key):
            raise ValueError(f"{state} short/long variants differ in static '{key}'")

    length_s = _positive(short.case.static, "reach_length_m", short_name)
    length_l = _positive(long.case.static, "reach_length_m", long_name)
    if length_l <= length_s:
        raise ValueError(f"{state} long reach must exceed short reach")

    event_columns = params.get("event_columns", {})
    if not isinstance(event_columns, dict) or state not in event_columns:
        raise ValueError(f"wave celerity needs an event column for state '{state}'")
    event_column = str(event_columns[state])
    t_short, var_short_h2, q_short, q_in_short, base_var_short = _centroid(
        short, probe, params, short_name, event_column
    )
    t_long, var_long_h2, q_long, q_in_long, base_var_long = _centroid(
        long, probe, params, long_name, event_column
    )
    if not np.isclose(q_in_short, q_in_long, rtol=0.0, atol=1.0e-10):
        raise ValueError(f"{state} short/long prescribed q_in baselines disagree")
    base_tolerance = float(params.get("base_flow_relative_tolerance", 0.01))
    if not np.isclose(q_short, q_long, rtol=base_tolerance, atol=1.0e-9):
        raise _ResponseFailure(
            "pair",
            f"{state} short/long settled base discharges disagree",
            {
                "state": state,
                "short_base_discharge_m3s": q_short,
                "long_base_discharge_m3s": q_long,
                "prescribed_base_inflow_m3s": q_in_short,
            },
        )
    dt_hours = t_long - t_short
    timing_floor = float(params.get("timing_floor_hours", 0.05))
    if not np.isfinite(dt_hours) or dt_hours <= timing_floor:
        raise _ResponseFailure(
            "pair",
            f"{state}: propagation is not downstream/resolved "
            f"(short/long centroid difference {dt_hours:.6g} h; "
            f"minimum {timing_floor:g} h)",
            {
                "state": state,
                "delta_t_h": float(dt_hours) if np.isfinite(dt_hours) else None,
                "short_centroid_h": t_short,
                "long_centroid_h": t_long,
                "timing_floor_hours": timing_floor,
            },
        )

    # Both model outputs are read at the outlet, so the paired propagation
    # distance is the full difference in reach lengths.
    dx = length_l - length_s
    c_obs = dx / (dt_hours * 3600.0)
    c_kin = _kinematic_celerity(q_in_short, short.case.static, short_name)
    residual = (c_obs - c_kin) / c_kin

    # For the constant-parameter diffusion-wave (Hayami) kernel,
    # Var[T_L] = 2 D L / c^3. Under the same paired-convolution assumptions
    # used by the centroid estimator, subtracting variances cancels the common
    # upstream response just as subtracting centroids does. Keep this as a
    # diagnostic only: #148's verdict is about celerity, not diffusivity.
    delta_variance_s2 = (var_long_h2 - var_short_h2) * 3600.0**2
    diffusivity_obs = None
    if delta_variance_s2 > 0.0:
        diffusivity_obs = float(c_obs**3 * delta_variance_s2 / (2.0 * dx))

    return _Measurement(
        state,
        float(c_obs),
        c_kin,
        float(residual),
        float(dt_hours),
        t_short,
        t_long,
        var_short_h2,
        var_long_h2,
        q_in_short,
        q_short,
        q_long,
        base_var_short,
        base_var_long,
        diffusivity_obs,
    )


@criterion("wave_celerity_bounds", paired=True)
def wave_celerity_bounds(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    """Celerity is downstream, near dQ/dA, and rises with hydraulic state."""
    states = tuple(params.get("states", ["low", "medium", "high"]))
    if len(states) < 2:
        raise ValueError("wave_celerity_bounds needs at least two hydraulic states")
    expected = {"short", "long"}
    if set(runs) != expected:
        raise ValueError(
            f"wave-celerity variants mismatch: expected {sorted(expected)}, got {sorted(runs)}"
        )

    tolerance = float(params.get("relative_tolerance", 0.05))
    separation = float(params.get("ordering_margin_fraction", 0.0))
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("relative_tolerance must be finite and non-negative")
    if not np.isfinite(separation) or separation < 0:
        raise ValueError("ordering_margin_fraction must be finite and non-negative")

    try:
        measured = [_pair(runs, probe, params, state) for state in states]
    except _ResponseFailure as exc:
        return CriterionResult(
            name="wave_celerity_bounds",
            status=FAIL,
            value=1.0,
            threshold=0.0,
            message=str(exc),
            diagnostics=dict(exc.diagnostics),
        )

    failures: list[str] = []
    worst = 0.0
    diagnostics: dict[str, Any] = {}
    for item in measured:
        if item.c_obs <= 0:
            failures.append(f"{item.state}: propagation is not downstream/resolved")
            violation = 1.0
        else:
            violation = max(0.0, abs(item.residual) - tolerance)
            if abs(item.residual) > tolerance + 1.0e-12:
                failures.append(
                    f"{item.state}: c_obs={item.c_obs:.3f} m/s differs from "
                    f"c_kin={item.c_kin:.3f} m/s by {100*item.residual:+.1f}%"
                )
        worst = max(worst, violation)
        diagnostics[item.state] = {
            "c_obs_m_s": item.c_obs,
            "c_kin_m_s": item.c_kin,
            "relative_residual": item.residual,
            "short_centroid_h": item.short_centroid_h,
            "long_centroid_h": item.long_centroid_h,
            "delta_t_h": item.dt_hours,
            "short_response_variance_h2": item.short_variance_h2,
            "long_response_variance_h2": item.long_variance_h2,
            "delta_response_variance_h2": (
                item.long_variance_h2 - item.short_variance_h2
            ),
            "prescribed_base_inflow_m3s": item.prescribed_base_m3s,
            "short_base_discharge_m3s": item.short_base_m3s,
            "long_base_discharge_m3s": item.long_base_m3s,
            "short_base_relative_variation": item.short_base_variation,
            "long_base_relative_variation": item.long_base_variation,
            "paired_diffusivity_m2_s": item.diffusivity_obs_m2_s,
        }

    c = [m.c_obs for m in measured]
    ordering_failures = []
    for left, right, left_name, right_name in zip(c, c[1:], states, states[1:]):
        required = max(0.0, separation * max(abs(left), abs(right)))
        if right - left <= required:
            ordering_failures.append(
                f"{left_name}->{right_name}: {left:.3f}->{right:.3f} m/s "
                f"does not exceed the {required:.3f} m/s resolution margin"
            )
    if ordering_failures:
        failures.extend(ordering_failures)
        worst = max(worst, 1.0)

    message = (
        "celerity is downstream, within the Manning/kinematic allowance, and "
        + " < ".join(f"{s}={m.c_obs:.3f} m/s" for s, m in zip(states, measured))
        if not failures else "; ".join(failures)
    )
    return CriterionResult(
        name="wave_celerity_bounds",
        status=PASS if not failures else FAIL,
        value=float(worst),
        threshold=0.0,
        message=message,
        diagnostics={
            "relative_tolerance": tolerance,
            "base_flow_relative_tolerance": float(
                params.get("base_flow_relative_tolerance", 0.01)
            ),
            "ordering_margin_fraction": separation,
            "states": diagnostics,
        },
    )

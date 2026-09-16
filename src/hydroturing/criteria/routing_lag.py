"""Rainfall-runoff timing against a synthetic geomorphic expectation.

The routing-lag probe serves one isolated design storm to catchments of
different size.  The weather is byte-identical between variants; only the
catchment geometry changes.  These paired criteria recover the storm response
from reported runoff, measure its peak lag from the added-rainfall centroid,
and ask two separate questions:

* ``lag_time_bounds``: is each lag within the broad range predicted by a
  Snyder synthetic-unit-hydrograph scaling?
* ``scaling_monotonicity``: is the lag non-decreasing, within daily sampling
  uncertainty, and observably longer across the full area ladder?

The separation is deliberate.  An instantaneous router and a router whose
delay shrinks with basin size are different counterexamples and the acceptance
gate should be able to say which assertion caught each one.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

_STATIC_KEYS = (
    "area_km2",
    "main_channel_length_km",
    "centroid_channel_length_km",
)
_SNYDER_SI_CONVERSION = 0.75
_SNYDER_CT = 4.0
_SNYDER_STANDARD_DURATION_RATIO = 5.5
_SNYDER_MIN_AREA_KM2 = 10.0 * 2.589988110336
_SNYDER_MAX_AREA_KM2 = 10_000.0 * 2.589988110336


@dataclass(frozen=True)
class _LagMeasurement:
    variant: str
    area_km2: float
    expected_hours: float
    observed_hours: float
    response_fraction: float
    diagnostics: dict[str, Any]


class _ResponseFailure(ValueError):
    """A finite model answer that cannot demonstrate a storm response.

    Generator and paired-case contract errors remain ordinary ``ValueError``
    instances and become harness errors.  This subclass is caught by the
    criteria and becomes a scientific FAIL: a model that emits zero runoff,
    non-finite runoff, or a still-rising response has answered the case, but
    has not supplied a physically measurable routing lag.
    """

    def __init__(self, message: str, diagnostics: dict[str, Any] | None = None):
        super().__init__(message)
        self.diagnostics = diagnostics or {}


def _number(static: dict[str, Any], key: str, variant: str) -> float:
    if key not in static:
        raise ValueError(
            f"variant '{variant}' is missing static field '{key}'; "
            f"routing lag needs {list(_STATIC_KEYS)}"
        )
    try:
        value = float(static[key])
    except (TypeError, ValueError):
        raise ValueError(
            f"variant '{variant}' has non-numeric static field '{key}'"
        ) from None
    if not np.isfinite(value):
        raise ValueError(
            f"variant '{variant}' has non-finite static field '{key}'"
        )
    return value


def _snyder_lag_hours(
    static: dict[str, Any], variant: str, event_duration_hours: float
) -> tuple[float, dict[str, float]]:
    """Return a reproducible Snyder lag derived from drainage area.

    Snyder's SI relation is ``t_lag = 0.75 Ct (L Lc)^0.3`` hours for lengths
    in kilometres. ``Ct=4`` sits inside the documented 0.4--8 regional range;
    the probe's factor-of-two band deliberately covers the common 1.8--2.2
    values as well.  Snyder's standard excess-rain duration is ``t_lag/5.5``;
    the published non-standard-duration correction is applied for the
    probe's one-day storm.

    The generator derives both lengths from area with Hack's published
    length-area relation.  All three geometry values are supplied to the
    model, which avoids giving basin area two contradictory meanings across
    this probe and the suite's area-unit invariance experiment.
    """
    area = _number(static, "area_km2", variant)
    length_km = _number(static, "main_channel_length_km", variant)
    centroid_length_km = _number(
        static, "centroid_channel_length_km", variant
    )
    if area <= 0 or length_km <= 0 or centroid_length_km <= 0:
        raise ValueError(
            f"variant '{variant}' needs positive area and channel lengths"
        )
    if not _SNYDER_MIN_AREA_KM2 <= area <= _SNYDER_MAX_AREA_KM2:
        raise ValueError(
            f"variant '{variant}' has area_km2={area:g}, outside Snyder's "
            f"published 10--10,000 mi2 range "
            f"({_SNYDER_MIN_AREA_KM2:.2f}--{_SNYDER_MAX_AREA_KM2:.2f} km2)"
        )
    if centroid_length_km > length_km:
        raise ValueError(
            f"variant '{variant}' has centroid channel length longer than its "
            "main channel"
        )

    standard_lag = (
        _SNYDER_SI_CONVERSION
        * _SNYDER_CT
        * (length_km * centroid_length_km) ** 0.3
    )
    standard_duration = standard_lag / _SNYDER_STANDARD_DURATION_RATIO
    if not np.isfinite(event_duration_hours) or event_duration_hours <= 0:
        raise ValueError("the design-storm duration must be finite and positive")
    expected = standard_lag - (standard_duration - event_duration_hours) / 4.0
    if not np.isfinite(expected) or expected <= 0:
        raise ValueError(
            f"variant '{variant}' produces an invalid Snyder lag ({expected!r} hours)"
        )
    return float(expected), {
        "area_km2": area,
        "main_channel_length_km": length_km,
        "centroid_channel_length_km": centroid_length_km,
        "snyder_ct": _SNYDER_CT,
        "snyder_standard_lag_hours": standard_lag,
        "snyder_standard_rain_duration_hours": standard_duration,
        "design_storm_duration_hours": event_duration_hours,
    }


def _positive_param(params: dict, key: str, default: float, *, allow_zero: bool = False) -> float:
    value = float(params.get(key, default))
    valid = np.isfinite(value) and (value >= 0 if allow_zero else value > 0)
    if not valid:
        relation = "non-negative" if allow_zero else "positive"
        raise ValueError(f"routing lag parameter '{key}' must be finite and {relation}")
    return value


def _same_forcing(runs: dict[str, RunResult], ordered: list[str]) -> None:
    """The geometry transform may not also change the weather or its window."""
    reference = runs[ordered[0]].case
    for variant in ordered[1:]:
        case = runs[variant].case
        if case.timestep != reference.timestep:
            raise ValueError(
                "routing-lag variants use different timesteps "
                f"({reference.timestep} and {case.timestep})"
            )
        if case.spinup_steps != reference.spinup_steps:
            raise ValueError(
                "routing-lag variants carry different numbers of spinup steps"
            )
        if not case.forcing.reset_index(drop=True).equals(
            reference.forcing.reset_index(drop=True)
        ):
            raise ValueError(
                f"variant '{variant}' does not carry exactly the same forcing as "
                f"variant '{ordered[0]}'; only static geometry may change"
            )


def _ordered_variants(runs: dict[str, RunResult], probe: ProbeSpec) -> list[str]:
    if len(runs) < 2:
        raise ValueError("routing lag needs at least two catchment variants")
    ordered = list(probe.variants) if probe.variants else list(runs)
    missing = [name for name in ordered if name not in runs]
    extra = [name for name in runs if name not in ordered]
    if missing or extra:
        raise ValueError(
            "routing-lag run variants do not match the probe declaration "
            f"(missing={missing}, extra={extra})"
        )
    return ordered


def _event_definition(
    run: RunResult,
    probe: ProbeSpec,
    params: dict,
    variant: str,
) -> tuple[np.ndarray, np.ndarray, int, int, float]:
    """Return added rain, precipitation, event bounds and its centroid in days."""
    window = make_window(run, probe)
    event_column = str(params.get("event_column", "_event_pr"))
    if event_column not in window.forcing.columns:
        raise ValueError(
            f"routing lag needs hidden forcing column '{event_column}' to identify "
            "the added design storm"
        )
    if "pr" not in window.forcing.columns:
        raise ValueError("routing lag needs 'pr' in the forcing")

    try:
        event = pd.to_numeric(window.forcing[event_column], errors="raise").to_numpy(float)
        precipitation = pd.to_numeric(window.forcing["pr"], errors="raise").to_numpy(float)
    except (TypeError, ValueError):
        raise ValueError(
            f"variant '{variant}' has non-numeric pr or {event_column} forcing"
        ) from None
    if not np.isfinite(event).all() or not np.isfinite(precipitation).all():
        raise ValueError(
            f"variant '{variant}' has non-finite pr or {event_column} forcing"
        )
    if np.any(event < 0):
        raise ValueError(f"variant '{variant}' has negative added-event precipitation")
    forcing_atol = _positive_param(
        params, "forcing_atol", 1e-9, allow_zero=True
    )
    if np.any(precipitation < -forcing_atol):
        raise ValueError(f"variant '{variant}' has negative precipitation")
    if np.any(event - precipitation > forcing_atol):
        raise ValueError(
            f"variant '{variant}' marks more added rain than total precipitation"
        )

    indices = np.flatnonzero(event > forcing_atol)
    if len(indices) == 0:
        raise ValueError("the scored forcing contains no added design storm")
    if np.any(np.diff(indices) != 1):
        raise ValueError("the added design storm must be one contiguous event")
    first, last = int(indices[0]), int(indices[-1])
    background_rain = precipitation - event
    if np.any(background_rain > forcing_atol):
        raise ValueError(
            "the scored routing-lag record must be dry outside its annotated "
            "design storm"
        )
    dt_days = window.dt_days
    min_pre = _positive_param(
        params, "min_pre_event_days", 10.0, allow_zero=True
    )
    min_post = _positive_param(
        params, "min_post_event_days", 30.0, allow_zero=True
    )
    pre_days = first * dt_days
    post_days = (len(event) - 1 - last) * dt_days
    if pre_days + 1e-12 < min_pre or post_days + 1e-12 < min_post:
        raise ValueError(
            "the scored window truncates the routing experiment: it keeps "
            f"{pre_days:g} pre-event and {post_days:g} post-event days, but "
            f"needs at least {min_pre:g} and {min_post:g}"
        )

    weights = event * dt_days
    event_depth = float(weights.sum())
    if event_depth <= 0:
        raise ValueError("the added design storm has zero integrated depth")
    centres = (np.arange(len(event), dtype=float) + 0.5) * dt_days
    centroid_days = float(np.dot(centres, weights) / event_depth)
    return event, precipitation, first, last, centroid_days


def _measure_one(
    run: RunResult,
    probe: ProbeSpec,
    params: dict,
    variant: str,
) -> _LagMeasurement:
    window = make_window(run, probe)
    event, _precipitation, first, last, rain_centroid_days = _event_definition(
        run, probe, params, variant
    )
    runoff_name = str(params.get("runoff", "mrro"))
    if runoff_name not in window.table.columns:
        raise _ResponseFailure(
            f"variant '{variant}' does not report '{runoff_name}'",
            {"variant": variant, "runoff": runoff_name},
        )
    try:
        runoff = pd.to_numeric(window.table[runoff_name], errors="coerce").to_numpy(float)
    except (TypeError, ValueError):  # pragma: no cover - to_numeric handles this
        runoff = np.full(len(window.table), np.nan)
    if len(runoff) != len(event):
        raise _ResponseFailure(
            f"variant '{variant}' returned {len(runoff)} scored runoff rows for "
            f"{len(event)} forcing rows",
            {"variant": variant, "runoff_rows": len(runoff), "forcing_rows": len(event)},
        )
    bad = int((~np.isfinite(runoff)).sum())
    if bad:
        raise _ResponseFailure(
            f"variant '{variant}' has {bad} non-finite runoff values",
            {"variant": variant, "nonfinite_count": bad},
        )

    baseline_days = _positive_param(params, "baseline_days", 5.0)
    baseline_steps = max(1, int(np.ceil(baseline_days / window.dt_days)))
    baseline_start = max(0, first - baseline_steps)
    baseline_values = runoff[baseline_start:first]
    if len(baseline_values) == 0:
        raise ValueError(
            "the design storm has no pre-event rows from which to estimate baseflow"
        )
    baseline = float(np.median(baseline_values))
    response = runoff - baseline
    search = response[first:]
    peak = float(np.max(search))

    event_depth = float(np.sum(event) * window.dt_days)
    positive_volume = float(np.clip(search, 0.0, None).sum() * window.dt_days)
    response_fraction = positive_volume / event_depth
    min_response = _positive_param(
        params, "min_response_fraction", 0.01, allow_zero=True
    )
    if peak <= 0 or response_fraction < min_response:
        raise _ResponseFailure(
            f"variant '{variant}' has no measurable storm response "
            f"({response_fraction:.3g} of the added rain; minimum {min_response:g})",
            {
                "variant": variant,
                "response_fraction": response_fraction,
                "min_response_fraction": min_response,
                "peak_response": peak,
                "baseline_runoff": baseline,
            },
        )

    peak_rtol = _positive_param(
        params, "peak_tie_rtol", 1e-12, allow_zero=True
    )
    peak_atol = _positive_param(
        params, "peak_tie_atol", 1e-12, allow_zero=True
    )
    tied_local = np.flatnonzero(np.isclose(search, peak, rtol=peak_rtol, atol=peak_atol))
    tied = tied_local + first
    if int(tied[-1]) == len(response) - 1:
        raise _ResponseFailure(
            f"variant '{variant}' runoff is still at its peak on the final row; "
            "the routing lag is outside the scored response window",
            {"variant": variant, "peak_on_final_row": True},
        )
    runoff_peak_days = float(
        np.mean((tied.astype(float) + 0.5) * window.dt_days)
    )
    observed_hours = (runoff_peak_days - rain_centroid_days) * 24.0
    event_duration_hours = (last - first + 1) * window.dt_days * 24.0
    expected_hours, geometry = _snyder_lag_hours(
        run.case.static, variant, event_duration_hours
    )
    diagnostics: dict[str, Any] = {
        **geometry,
        "expected_lag_hours": expected_hours,
        "expected_lag_days": expected_hours / 24.0,
        "observed_lag_hours": observed_hours,
        "observed_lag_days": observed_hours / 24.0,
        "rain_centroid_day": rain_centroid_days,
        "runoff_peak_day": runoff_peak_days,
        "tied_peak_steps": int(len(tied)),
        "baseline_runoff": baseline,
        "response_fraction": response_fraction,
        "added_rain_mm": event_depth,
    }
    return _LagMeasurement(
        variant=variant,
        area_km2=geometry["area_km2"],
        expected_hours=expected_hours,
        observed_hours=observed_hours,
        response_fraction=response_fraction,
        diagnostics=diagnostics,
    )


def _measure_all(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> list[_LagMeasurement]:
    ordered = _ordered_variants(runs, probe)
    _same_forcing(runs, ordered)
    measurements = [_measure_one(runs[name], probe, params, name) for name in ordered]
    areas = [m.area_km2 for m in measurements]
    if len(set(areas)) != len(areas):
        raise ValueError("routing-lag variants must have distinct catchment areas")
    return measurements


def _response_fail(name: str, exc: _ResponseFailure, threshold: float | None = None):
    return CriterionResult(
        name=name,
        status=FAIL,
        threshold=threshold,
        message=str(exc),
        diagnostics=dict(exc.diagnostics),
    )


@criterion("lag_time_bounds", paired=True)
def lag_time_bounds(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    """Every measured peak lag must lie inside its Snyder plausibility band."""
    lower_ratio = _positive_param(params, "lower_ratio", 0.5, allow_zero=True)
    upper_ratio = _positive_param(params, "upper_ratio", 2.0)
    if lower_ratio > upper_ratio:
        raise ValueError("lag_time_bounds lower_ratio exceeds upper_ratio")
    tolerance_hours = 24.0 * _positive_param(
        params, "discretization_tolerance_days", 0.5, allow_zero=True
    )
    try:
        measurements = _measure_all(runs, probe, params)
    except _ResponseFailure as exc:
        return _response_fail("lag_time_bounds", exc)

    # The rain centroid is a float quotient, so a lag of a whole number of steps
    # can carry about 1e-13 h of rounding; that must not decide a bound.
    timing_atol_hours = 1e-6
    failures: list[str] = []
    worst_deviation = 0.0
    diagnostics: dict[str, Any] = {}
    for measured in measurements:
        low = max(0.0, lower_ratio * measured.expected_hours - tolerance_hours)
        high = upper_ratio * measured.expected_hours + tolerance_hours
        observed = measured.observed_hours
        if observed < low - timing_atol_hours:
            deviation = (low - observed) / measured.expected_hours
        elif observed > high + timing_atol_hours:
            deviation = (observed - high) / measured.expected_hours
        else:
            deviation = 0.0
        worst_deviation = max(worst_deviation, float(deviation))
        item = dict(measured.diagnostics)
        item.update(
            {
                "lower_bound_hours": low,
                "upper_bound_hours": high,
                "observed_to_expected_ratio": observed / measured.expected_hours,
                "normalized_bound_violation": float(deviation),
            }
        )
        diagnostics[measured.variant] = item
        if deviation > 0:
            failures.append(
                f"{measured.variant} lag {observed / 24.0:.2f} d is outside "
                f"[{low / 24.0:.2f}, {high / 24.0:.2f}] d"
            )

    ok = not failures
    ranges = ", ".join(
        f"{m.variant} {m.observed_hours / 24.0:.2f}/"
        f"{m.expected_hours / 24.0:.2f} d"
        for m in measurements
    )
    return CriterionResult(
        name="lag_time_bounds",
        status=PASS if ok else FAIL,
        value=worst_deviation,
        threshold=0.0,
        message=(
            f"observed/Snyder lags are plausible ({ranges})"
            if ok
            else "; ".join(failures)
        ),
        diagnostics={
            "lower_ratio": lower_ratio,
            "upper_ratio": upper_ratio,
            "discretization_tolerance_hours": tolerance_hours,
            "variants": diagnostics,
        },
    )


@criterion("scaling_monotonicity", paired=True)
def scaling_monotonicity(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    """Lag may quantise at one day, but may not reverse with basin area.

    Requiring a positive jump at every adjacent scale would reject two true
    sub-day lags that land on the same daily row.  We therefore combine a
    small tolerance on each adjacent comparison with a resolvable end-to-end
    span across the full area ladder.
    """
    reversal_tolerance_days = _positive_param(
        params, "reversal_tolerance_days", 0.5, allow_zero=True
    )
    min_span_days = _positive_param(
        params, "min_span_days", 2.0, allow_zero=True
    )
    try:
        measurements = sorted(
            _measure_all(runs, probe, params), key=lambda measured: measured.area_km2
        )
    except _ResponseFailure as exc:
        return _response_fail("scaling_monotonicity", exc, threshold=0.0)

    for left, right in zip(measurements, measurements[1:]):
        if right.expected_hours <= left.expected_hours:
            raise ValueError(
                "Snyder expected lag does not increase with the generated catchment "
                f"areas ({left.variant} to {right.variant}); fix the static geometry"
            )

    increments: dict[str, float] = {}
    reversal_excesses: list[float] = []
    failures: list[str] = []
    for left, right in zip(measurements, measurements[1:]):
        increase = (right.observed_hours - left.observed_hours) / 24.0
        key = f"{left.variant}->{right.variant}"
        increments[key] = float(increase)
        reversal_excess = max(0.0, -increase - reversal_tolerance_days)
        reversal_excesses.append(reversal_excess)
        if reversal_excess > 1e-12:
            failures.append(
                f"lag decreases by {-increase:.2f} d from {left.variant} to "
                f"{right.variant} (allowed sampling reversal "
                f"{reversal_tolerance_days:g} d)"
            )

    minimum = min(increments.values())
    span_days = (
        measurements[-1].observed_hours - measurements[0].observed_hours
    ) / 24.0
    span_shortfall_days = max(0.0, min_span_days - span_days)
    if span_shortfall_days > 1e-12:
        failures.append(
            f"lag span from {measurements[0].variant} to "
            f"{measurements[-1].variant} is {span_days:.2f} d "
            f"(minimum {min_span_days:g} d)"
        )
    ok = not failures
    violation_days = max(
        span_shortfall_days,
        max(reversal_excesses, default=0.0),
    )
    return CriterionResult(
        name="scaling_monotonicity",
        status=PASS if ok else FAIL,
        value=violation_days,
        threshold=0.0,
        message=(
            "rainfall-runoff lag is non-decreasing with catchment area "
            f"(full span {span_days:.2f} d)"
            if ok
            else "; ".join(failures)
        ),
        diagnostics={
            "reversal_tolerance_days": reversal_tolerance_days,
            "min_span_days": min_span_days,
            "span_days": span_days,
            "span_shortfall_days": span_shortfall_days,
            "maximum_reversal_excess_days": max(
                reversal_excesses, default=0.0
            ),
            "minimum_adjacent_increment_days": minimum,
            "increments_days": increments,
            "ordered_variants": [m.variant for m in measurements],
            "variants": {m.variant: dict(m.diagnostics) for m in measurements},
        },
    )

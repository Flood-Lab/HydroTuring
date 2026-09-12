"""Create extreme rainfall beyond fitted 100-year depths to stress-test water closure.

The aim is to test whether water-balance closure generalizes to extreme
rainfall a model may not have encountered during training. Closure learned
as a statistical regularity may fail under unfamiliar forcing, even when
it appears reliable within the training distribution. Extremeness is
defined against the synthetic baseline climate; the model's actual
training distribution is not assumed known.

Every seed supplies independent weather and 100-year calibration streams.
After a 365-day spinup, select the median-wet year from twenty 365-day
scored blocks: the annual depth closest to their median, with chronological
tie-breaking. This confines the stress to one representative year while
retaining the full continuous record for event and whole-window budgets.
GEV distributions fitted by L-moments supply 1/3/7-day construction targets.
Whole rainfall events are overlapped to exceed at least one of these
duration-specific 100-year depths. Each chronological group stops at its
first exceedance of any target; the remaining year-end group is retained
even below target. Groups are then placed chronologically with one dry day
between them, starting on May 1 in the selected scored year.
These thresholds control construction only:
event_water_closure scores all complete events.

All calibration code lives here. No additional dependency or data download
is needed. Underscore columns and DataFrame attrs stay on the host; models
receive only time, pr, tas and pet, plus the unchanged catchment attributes.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

PERIOD_YEARS = 20
SPINUP_DAYS = 365
N_STEPS = int(PERIOD_YEARS * 365) + SPINUP_DAYS
TARGET_RETURN_PERIOD = 100
INTER_EVENT_DRY_DAYS = 1
STORM_START_MONTH_DAY = (5, 1)
CALIBRATION_YEARS = 100
REPORT_DURATIONS = (1, 3, 7)
RETURN_PERIODS = (2, 5, 10, 20, 50, 100, 200, 500)
EULER_GAMMA = 0.5772156649015329

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 40.0,
}


def _rng(seed: int, stream: int) -> np.random.Generator:
    """Fixed child 0 = model weather; child 1 = frequency calibration."""
    children = np.random.SeedSequence(seed).spawn(2)
    return np.random.Generator(np.random.PCG64(children[stream]))


def _rainfall(rng: np.random.Generator, day: np.ndarray) -> np.ndarray:
    """The catchment-closure rainfall law, also used for the calibration."""
    doy = day % 365
    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(len(day)) < p_wet
    return np.round(np.where(wet, rng.gamma(shape=0.7, scale=13.0, size=len(day)), 0.0), 6)


def generate_baseline(seed: int) -> tuple[pd.DataFrame, dict]:
    """Unmodified catchment-closure weather, drawn from this seed's child 0.

    The climate equations and static attributes are unchanged. Snow remains
    possible; neither temperature nor PET is altered to activate a model's
    defect. A separate child stream supplies the DDF calibration rainfall.
    """
    rng = _rng(seed, 0)
    day = np.arange(N_STEPS)
    doy = day % 365
    pr = _rainfall(rng, day)

    seasonal = 9.0 - 12.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(N_STEPS)
    innovation = rng.normal(0.0, 2.6, N_STEPS)
    for t in range(1, N_STEPS):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    tas = seasonal + noise

    daylength = 1.0 + 0.35 * np.cos(2 * np.pi * (doy - 172) / 365)
    pet = np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
        }
    )
    return forcing, dict(STATIC)


def rainfall_events(pr: np.ndarray) -> list[tuple[int, int]]:
    """Maximal positive-rainfall runs as (start, stop), with stop exclusive.

    Runs touching record edges are returned too; callers decide whether the
    preceding and following dry days are observed and the event is complete.
    """
    rain = np.asarray(pr, dtype=float)
    if rain.ndim != 1 or not np.isfinite(rain).all() or (rain < 0).any():
        raise ValueError("rainfall must be a finite, nonnegative one-dimensional array")
    edges = np.diff(np.r_[False, rain > 0, False].astype(int))
    return list(zip(np.flatnonzero(edges == 1).tolist(), np.flatnonzero(edges == -1).tolist()))


def event_duration_maxima(pr: np.ndarray, durations=REPORT_DURATIONS) -> dict[int, float]:
    """Maximum D-day depths attributable to an isolated event, in mm.

    Outside the supplied hyetograph rainfall is zero, so an event shorter
    than D contributes its full depth. Longer events use sliding D-day
    sums. Adjacent events never help a group reach its construction target.
    Depths use the six-decimal precision of the generated daily forcing.
    """
    rain = np.asarray(pr, dtype=float)
    rainfall_events(rain)
    if not len(rain):
        raise ValueError("event depth needs a nonempty hyetograph")
    return {
        d: float(np.round(rain.sum() if len(rain) < d else
                          np.convolve(rain, np.ones(d), mode="valid").max(), 6))
        for d in _duration_tuple(durations)
    }


def overlap_event_groups(
    pr: np.ndarray, start: int, stop: int | None = None, *, thresholds: dict[int, float],
) -> tuple[np.ndarray, list[dict]]:
    """Greedily overlap complete events until any duration target is exceeded.

    Add whole source events, with starts aligned and shorter events padded
    with zeros. Stop at the first strict exceedance of any supplied D-day
    depth, including a first source that already exceeds a target. Restart
    with the next unused source. Combine all remaining sources at year end
    even if below target. Keep window-crossing and record-edge events intact.
    Place each group at its first start and clear its original source span;
    every eligible event is used exactly once and the time axis is unchanged.
    """
    rain = np.asarray(pr, dtype=float)
    events = rainfall_events(rain)
    stop = len(rain) if stop is None else stop
    if (any(isinstance(v, bool) or not isinstance(v, (int, np.integer)) for v in (start, stop))
            or not 0 <= start <= stop <= len(rain)):
        raise ValueError("overlap bounds must satisfy 0 <= start <= stop <= record length")
    durations = _duration_tuple(thresholds)
    if any(not np.isfinite(thresholds[d]) or thresholds[d] <= 0 for d in durations):
        raise ValueError("construction thresholds must be finite and positive")
    eligible = [
        (a, b) for a, b in events
        if start <= a and b <= stop and a > 0 and b < len(rain)
    ]
    result = rain.copy()
    groups = []
    sources = []
    combined = np.array([], dtype=float)
    for index, (a, b) in enumerate(eligible):
        sources.append((a, b))
        combined = np.pad(combined, (0, max(0, b - a - len(combined))))
        combined[:b - a] += rain[a:b]
        combined = np.round(combined, 6)
        depths = event_duration_maxima(combined, durations)
        exceeded = [int(d) for d in durations if depths[d] > thresholds[d]]
        if exceeded or index == len(eligible) - 1:
            first = sources[0][0]
            end = first + len(combined)
            result[first:b] = 0.0
            result[first:end] = combined
            groups.append({
                "event_id": len(groups) + 1, "start": first, "stop": end,
                "source_events": tuple(sources), "source_event_count": len(sources),
                "depths_mm": {str(d): depths[d] for d in durations},
                "exceeded_durations_days": exceeded, "target_exceeded": bool(exceeded),
                "stop_reason": "threshold_exceeded" if exceeded else "end_of_year",
            })
            sources = []
            combined = np.array([], dtype=float)
    return result, groups


def pack_event_groups(
    pr: np.ndarray, groups: list[dict], gap_days: int = INTER_EVENT_DRY_DAYS,
    *, anchor: int | None = None, window: tuple[int, int] | None = None,
) -> tuple[np.ndarray, list[dict]]:
    """Place constructed groups consecutively with an exact dry-day gap.

    Use the supplied anchor, or the first group's existing start by default.
    Preserve each hyetograph, order and target status, including the unmet
    remainder. Clear old locations before writing the new ones. The packed
    sequence must fit inside the destination window (stop exclusive) and
    have observed dry neighbors in the record, without overwriting or
    touching an ungrouped wet event. No rainfall is truncated to make it fit.
    Source-event indices remain original; start/stop become the placed dates.
    """
    rain = np.asarray(pr, dtype=float)
    rainfall_events(rain)
    if isinstance(gap_days, bool) or not isinstance(gap_days, (int, np.integer)) or gap_days < 1:
        raise ValueError("gap_days must be a positive integer")
    result = rain.copy()
    if not groups:
        return result, []
    hyetographs = []
    previous_stop = -1
    for group in groups:
        a, b = group["start"], group["stop"]
        if (any(isinstance(v, bool) or not isinstance(v, (int, np.integer)) for v in (a, b))
                or not 0 < a < b < len(rain) or a <= previous_stop):
            raise ValueError("groups must be ordered, separate and inside the record")
        if rain[a - 1] != 0 or rain[b] != 0 or not (rain[a:b] > 0).all():
            raise ValueError("each group must be a complete wet event with dry neighbors")
        hyetographs.append(rain[a:b].copy())
        result[a:b] = 0.0
        previous_stop = b
    first = groups[0]["start"] if anchor is None else anchor
    if window is None:
        lower, upper = 0, len(rain)
    else:
        if not isinstance(window, (tuple, list)) or len(window) != 2:
            raise ValueError("window must contain start and exclusive stop bounds")
        lower, upper = window
    if (any(isinstance(v, bool) or not isinstance(v, (int, np.integer))
            for v in (first, lower, upper))
            or not 0 <= lower <= first < upper <= len(rain) or first == 0):
        raise ValueError("anchor and window must identify an interior destination in the record")
    packed_stop = first + sum(len(h) for h in hyetographs) + gap_days * (len(groups) - 1)
    if packed_stop > upper or packed_stop >= len(rain):
        raise ValueError("packed events must fit inside the destination window with a following dry day")
    if (result[first - 1:packed_stop + 1] > 0).any():
        raise ValueError("packing would overwrite or touch ungrouped rainfall")
    packed = []
    cursor = first
    for group, hyetograph in zip(groups, hyetographs):
        end = cursor + len(hyetograph)
        result[cursor:end] = hyetograph
        packed.append({
            **group, "original_start": group["start"], "original_stop": group["stop"],
            "start": cursor, "stop": end,
        })
        cursor = end + gap_days
    return result, packed


def select_median_wet_year(pr: np.ndarray) -> dict:
    """Choose the scored 365-day block nearest the median annual rainfall.

    Exclude spinup and use all twenty complete scored blocks. Daily depths
    are expressed as integer micro-mm at the forcing's six-decimal precision
    before summation. Twice the median is the sum of the two central annual
    totals, so distances and ties are compared exactly without floating-point
    ambiguity. Among equally close years choose the earliest chronological
    block, also used to rank equal annual totals. The selected fraction is
    descriptive, not a guarantee of a whole-record residual below 5%.
    """
    rain = np.asarray(pr, dtype=float)
    rainfall_events(rain)
    if len(rain) != N_STEPS:
        raise ValueError(f"year selection expects the full {N_STEPS}-row case")
    blocks = rain[SPINUP_DAYS:].reshape(PERIOD_YEARS, 365)
    totals = [sum(int(round(float(value) * 1_000_000)) for value in block)
              for block in blocks]
    ordered = sorted(totals)
    twice_median = ordered[(PERIOD_YEARS - 1) // 2] + ordered[PERIOD_YEARS // 2]
    distances = [abs(2 * total - twice_median) for total in totals]
    index = min(range(PERIOD_YEARS), key=lambda i: (distances[i], i))
    ranked = sorted(range(PERIOD_YEARS), key=lambda i: (totals[i], i))
    total_precip = sum(totals)
    start = SPINUP_DAYS + index * 365
    return {
        "year_number": index + 1, "start_row": start, "stop_row": start + 365,
        "baseline_precip_mm": totals[index] / 1_000_000,
        "annual_precip_mm": [total / 1_000_000 for total in totals],
        "median_precip_mm": twice_median / 2_000_000,
        "distance_from_median_mm": distances[index] / 2_000_000,
        "annual_rank_ascending": ranked.index(index) + 1,
        "scored_precip_mm": total_precip / 1_000_000,
        "selected_fraction_of_scored_precip": totals[index] / total_precip if total_precip else None,
        "selection": "365-day scored block closest to median annual precipitation",
        "tie_break": "earliest chronological block, using integer micro-mm distances",
    }


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    """Generate one continuous case with host-only overlap annotations."""
    forcing, static = generate_baseline(seed)
    baseline_rain = forcing["pr"].to_numpy().copy()
    selection = select_median_wet_year(baseline_rain)
    start, stop = selection["start_row"], selection["stop_row"]
    selection.update({
        "start_time": str(forcing["time"].iloc[start]),
        "end_time": str(forcing["time"].iloc[stop - 1]),
    })
    calibration = calibrate_ddf(seed)
    thresholds = {d: calibration["return_levels_mm"][str(d)][str(TARGET_RETURN_PERIOD)]
                  for d in REPORT_DURATIONS}
    rain, groups = overlap_event_groups(baseline_rain, start, stop, thresholds=thresholds)
    dates = pd.to_datetime(forcing["time"])
    month, day = STORM_START_MONTH_DAY
    candidates = np.flatnonzero(
        ((dates.dt.month == month) & (dates.dt.day == day)).to_numpy()[start:stop]
    )
    if len(candidates) != 1:
        raise ValueError("the selected scored block must contain exactly one May 1 anchor")
    anchor = start + int(candidates[0])
    rain, groups = pack_event_groups(rain, groups, anchor=anchor, window=(start, stop))
    forcing["pr"] = np.round(rain, 6)
    row = np.arange(N_STEPS)
    forcing["_regime"] = np.where(
        (row >= start) & (row < stop), "anomaly", "ordinary"
    )
    event_id = np.zeros(N_STEPS, dtype=int)
    event_start = np.zeros(N_STEPS, dtype=bool)
    event_end = np.zeros(N_STEPS, dtype=bool)
    for event in groups:
        a, b = event["start"], event["stop"]
        event_id[a:b] = event["event_id"]
        # Construction diagnostics only; closure detects every complete wet
        # event from supplied precipitation, including unmodified events.
        event_start[a] = True
        event_end[b - 1] = True
    forcing["_event_id"] = event_id
    forcing["_event_start"] = event_start
    forcing["_event_end"] = event_end
    forcing.attrs["rainfall_diagnostics"] = {
        "seed": int(seed), "group_count": len(groups),
        "target_exceeded_group_count": sum(g["target_exceeded"] for g in groups),
        "selection": selection, "calibration": calibration,
        "construction": {
            "target_return_period_years": TARGET_RETURN_PERIOD,
            "thresholds_mm": {str(d): thresholds[d] for d in REPORT_DURATIONS},
            "window_convention": "maximum D-day sum of the isolated event, zero outside it",
            "inter_event_dry_days": INTER_EVENT_DRY_DAYS,
            "anchor_row": anchor, "anchor_time": str(forcing["time"].iloc[anchor]),
            "placement": "chronological groups from May 1 in the selected scored year, separated by one dry day",
        },
        "groups": [{**g, "source_events": [list(pair) for pair in g["source_events"]]}
                   for g in groups],
        "events": describe_constructed_events(forcing, groups, calibration),
        "modified_year": describe_modified_year(forcing, calibration, start, stop),
    }
    return forcing, static


def event_summary(forcing: pd.DataFrame) -> pd.DataFrame:
    """One descriptive row per merged wet event; no frequency-based selection."""
    columns = [
        "event_id", "start", "end", "duration_days", "precip_mm",
        "mean_mm_per_day", "peak_mm_per_day",
    ]
    records = []
    for event_id, event in forcing.loc[forcing["_event_id"] > 0].groupby("_event_id"):
        records.append({
            "event_id": int(event_id),
            "start": event["time"].iloc[0],
            "end": event["time"].iloc[-1],
            "duration_days": len(event),
            "precip_mm": float(event["pr"].sum()),
            "mean_mm_per_day": float(event["pr"].mean()),
            "peak_mm_per_day": float(event["pr"].max()),
        })
    return pd.DataFrame.from_records(records, columns=columns)


def sample_l_moments(sample) -> tuple[float, float, float]:
    """Unbiased sample L1, L2 and L-skewness from ranked observations.

    Probability-weighted moments b_r use binomial(i, r)/binomial(n-1, r)
    for zero-based rank i. See Hosking's lmom SAMLMU/PELGEV conventions:
    https://cran.r-universe.dev/lmom/doc/manual.html
    """
    x = np.asarray(sample, dtype=float)
    if x.ndim != 1 or len(x) < 3 or not np.isfinite(x).all():
        raise ValueError("L-moments need at least three finite observations in one dimension")
    x = np.sort(x)
    n = len(x)
    rank = np.arange(n, dtype=float)
    b0 = float(x.mean())
    b1 = float(np.mean(x * rank / (n - 1)))
    b2 = float(np.mean(x * rank * (rank - 1) / ((n - 1) * (n - 2))))
    l2 = 2 * b1 - b0
    if l2 <= 0 or not np.isfinite([b0, b1, b2, l2]).all():
        raise ValueError("L-moments require a finite, positive L2 (nonconstant sample)")
    tau3 = (6 * b2 - 6 * b1 + b0) / l2
    if not -1 < tau3 < 1:
        raise ValueError("L-skewness must lie strictly between -1 and 1")
    return b0, l2, float(tau3)


def fit_gev_lmom(sample) -> dict:
    """GEV fitted to L-moments; Hosking shape k = -conventional xi.

    Solve tau3 = 2*(1-3**(-k))/(1-2**(-k))-3 by deterministic bisection.
    The root is unique for k > -1. No likelihood fit, random initialisation,
    automatic distribution switch or approximation-only shape estimate.
    """
    l1, l2, tau3 = sample_l_moments(sample)
    log2, log3 = math.log(2), math.log(3)
    gumbel_tau = 2 * log3 / log2 - 3

    def skew(k):
        if abs(k) < 1e-10:
            return gumbel_tau
        return 2 * math.expm1(-k * log3) / math.expm1(-k * log2) - 3

    if abs(tau3 - gumbel_tau) < 1e-10:
        k = 0.0
    else:
        lo, hi = -1.0, 1.0
        while skew(hi) > tau3:
            hi *= 2
        for _ in range(100):
            mid = (lo + hi) / 2
            if skew(mid) > tau3:
                lo = mid
            else:
                hi = mid
        k = (lo + hi) / 2
    if abs(k) < 1e-8:
        k = 0.0
        scale = l2 / log2
        location = l1 - EULER_GAMMA * scale
    else:
        log_gamma = math.lgamma(1 + k)
        scale = l2 * k / (-math.expm1(-k * log2)) * math.exp(-log_gamma)
        location = l1 + scale * math.expm1(log_gamma) / k
    if scale <= 0 or not np.isfinite([location, scale, k]).all():
        raise ValueError("GEV L-moment fit produced invalid parameters")
    return {
        "location": location, "scale": scale, "shape_k": k,
        "l1": l1, "l2": l2, "tau3": tau3, "n": len(sample),
    }


def _gev_parameters(fit: dict) -> tuple[float, float, float]:
    location, scale, k = (float(fit[v]) for v in ("location", "scale", "shape_k"))
    if scale <= 0 or not np.isfinite([location, scale, k]).all():
        raise ValueError("GEV parameters must be finite with positive scale")
    return location, scale, k


def gev_return_level(fit: dict, return_period: float) -> float:
    """Rainfall depth with annual exceedance probability 1 / return_period."""
    location, scale, k = _gev_parameters(fit)
    t = float(return_period)
    if not math.isfinite(t) or t <= 1:
        raise ValueError("return_period must be finite and greater than one year")
    log_a = math.log(-math.log1p(-1 / t))
    return location - scale * (log_a if k == 0 else math.expm1(k * log_a) / k)


def gev_return_period(fit: dict, depth: float) -> float:
    """Inverse annual exceedance probability, with explicit support limits.

    Infinity means beyond a fitted upper endpoint or representable tail; it
    is not a defensible numeric return-period estimate. Reports convert it
    to null and an explanatory status, rather than cap it at an invented T.
    """
    location, scale, k = _gev_parameters(fit)
    if not math.isfinite(depth):
        raise ValueError("depth must be finite")
    z = (float(depth) - location) / scale
    if k and 1 - k * z <= 0:
        return math.inf if k > 0 else 1.0
    log_u = -z if k == 0 else math.log1p(-k * z) / k
    if log_u > 700:
        return 1.0
    u = math.exp(log_u)
    sf = -math.expm1(-u)
    return 1 / sf if sf > 0 else math.inf


def _duration_tuple(durations) -> tuple[int, ...]:
    values = tuple(durations)
    if not values or any(
        isinstance(d, bool) or not isinstance(d, (int, np.integer)) or d < 1 for d in values
    ) or len(set(values)) != len(values):
        raise ValueError("durations must be distinct positive integers")
    return tuple(int(d) for d in values)


def annual_maxima(
    pr, years: int, durations=REPORT_DURATIONS, days_per_year: int = 365,
) -> dict[int, np.ndarray]:
    """Rolling-depth maxima in complete years, assigned by window end.

    Input includes max(durations)-1 prefix days, then years*days_per_year
    days. Thus every first-year ending day also has a full rolling window.
    Dry days and windows crossing a year boundary are retained.
    """
    durations = _duration_tuple(durations)
    if any(isinstance(v, bool) or not isinstance(v, (int, np.integer)) or v < 1
           for v in (years, days_per_year)):
        raise ValueError("years and days_per_year must be positive integers")
    rain = np.asarray(pr, dtype=float)
    rainfall_events(rain)  # shared finite, nonnegative, 1-D validation
    prefix = max(durations) - 1
    if len(rain) != prefix + years * days_per_year:
        raise ValueError("calibration rain must include the full years and rolling-window prefix")
    maxima = {}
    for d in durations:
        rolling = np.convolve(rain, np.ones(d), mode="valid")
        # Valid convolution row 0 ends on input row d-1.
        annual = rolling[prefix - d + 1:].reshape(years, days_per_year)
        maxima[d] = annual.max(axis=1)
    return maxima


def calibrate_ddf(seed: int, years: int = CALIBRATION_YEARS, durations=REPORT_DURATIONS) -> dict:
    """100 synthetic years on child 1, separate from model weather on child 0."""
    durations = _duration_tuple(durations)
    if isinstance(years, bool) or not isinstance(years, (int, np.integer)) or years < 3:
        raise ValueError("calibration needs at least three complete years")
    prefix = max(durations) - 1
    day = np.arange(-prefix, years * 365)
    rain = _rainfall(_rng(seed, 1), day)
    maxima = annual_maxima(rain, years, durations)
    fits = {str(d): fit_gev_lmom(maxima[d]) for d in durations}
    levels = {
        str(d): {str(t): gev_return_level(fits[str(d)], t) for t in RETURN_PERIODS}
        for d in durations
    }
    ordered = sorted(durations)
    crossings = [
        t for t in RETURN_PERIODS
        if any(levels[str(a)][str(t)] > levels[str(b)][str(t)]
               for a, b in zip(ordered[:-1], ordered[1:]))
    ]
    return {
        "years": int(years), "days_per_year": 365, "prefix_days": prefix,
        "method": "annual maxima; GEV; unbiased sample L-moments; deterministic bisection",
        "shape_convention": "Hosking k = -xi (same sign as scipy genextreme c)",
        "random_stream": 1, "year_assignment": "365-day blocks, by rolling-window end",
        "durations_days": list(durations), "fits": fits,
        "annual_maxima_mm": {str(d): maxima[d].tolist() for d in durations},
        "return_levels_mm": levels, "crossing_return_levels_years": crossings,
        "notes": (
            "Return periods describe the specified synthetic precipitation climate. "
            "They are estimates from finite samples, not exact recurrence times or "
            "external model training limits. Values above the calibration record "
            "length are extrapolations. Confidence intervals are not computed. "
            "Independent duration fits are reported without artificial smoothing."
        ),
    }


def _frequency_estimate(depth: float, duration: int, calibration: dict) -> dict:
    """JSON-safe marginal annual frequency diagnostics for one depth."""
    fit = calibration["fits"][str(duration)]
    t = gev_return_period(fit, depth)
    above_support = fit["shape_k"] > 0 and depth >= fit["location"] + fit["scale"] / fit["shape_k"]
    status = "estimated" if math.isfinite(t) else (
        "above_fitted_upper_endpoint" if above_support else "tail_not_representable"
    )
    return {
        "max_depth_mm": depth, "return_period_years": t if math.isfinite(t) else None,
        "return_period_status": status, "extrapolates_record": t > calibration["years"],
        "exceeds_calibration_max": depth > max(calibration["annual_maxima_mm"][str(duration)]),
    }


def describe_constructed_events(forcing: pd.DataFrame, groups: list[dict], calibration: dict) -> list[dict]:
    """Per-group duration depths and construction thresholds for human review.

    These are isolated-event maxima, not annual maxima or a joint return
    period across durations. Below-target year-end groups remain visible.
    """
    records = []
    for group in groups:
        start, stop = group["start"], group["stop"]
        event = forcing["pr"].iloc[start:stop].to_numpy(dtype=float)
        depths = event_duration_maxima(event, REPORT_DURATIONS)
        for d, depth in depths.items():
            threshold = calibration["return_levels_mm"][str(d)][str(TARGET_RETURN_PERIOD)]
            records.append({
                "event_id": group["event_id"], "source_event_count": group["source_event_count"],
                "start_time": str(forcing["time"].iloc[start]),
                "end_time": str(forcing["time"].iloc[stop - 1]),
                "event_duration_days": stop - start, "precip_mm": float(np.round(event.sum(), 6)),
                "duration_days": d, "threshold_mm": threshold,
                "exceeds_threshold": depth > threshold, "target_exceeded": group["target_exceeded"],
                "stop_reason": group["stop_reason"], **_frequency_estimate(depth, d, calibration),
            })
    return records


def describe_modified_year(
    forcing: pd.DataFrame, calibration: dict, start: int, stop: int,
) -> list[dict]:
    """Maximum 1/3/7-day rainfall ending in the selected 365-day scored block.

    As in calibration, windows are assigned by their end and may include
    preceding days, dry gaps or multiple wet events. Maxima at different
    durations need not be from the same event. No extreme label is assigned.
    """
    rain = forcing["pr"].to_numpy(dtype=float)
    rainfall_events(rain)
    if len(rain) != N_STEPS:
        raise ValueError(f"modified-year description expects the full {N_STEPS}-row case")
    if (any(isinstance(v, bool) or not isinstance(v, (int, np.integer)) for v in (start, stop))
            or not SPINUP_DAYS <= start < stop <= N_STEPS
            or stop - start != 365 or (start - SPINUP_DAYS) % 365):
        raise ValueError("description bounds must identify one complete post-spinup 365-day year")
    records = []
    for d in calibration["durations_days"]:
        rolling = np.convolve(rain, np.ones(d), mode="valid")
        first = start - d + 1
        offset = int(np.argmax(rolling[first:stop - d + 1]))
        window_end = start + offset
        window_start = window_end - d + 1
        depth = float(rolling[first + offset])
        records.append({
            "duration_days": d, "start_time": str(forcing["time"].iloc[window_start]),
            "end_time": str(forcing["time"].iloc[window_end]),
            **_frequency_estimate(depth, d, calibration),
        })
    return records


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    seeds = parser.add_mutually_exclusive_group()
    seeds.add_argument("--seed", type=int, default=20260912)
    seeds.add_argument("--gate-seeds", action="store_true", help="describe all five fixed acceptance seeds")
    parser.add_argument("--output-dir", type=Path, help="write report.json, events.csv, modified_year.csv and ddf.csv")
    args = parser.parse_args()
    if args.gate_seeds:
        from hydroturing.seeds import gate_seeds
        selected_seeds = gate_seeds("mass/extreme-event-closure", 5)
    else:
        selected_seeds = [args.seed]
    reports, summaries, ddf_rows, event_rows = [], [], [], []
    for seed in selected_seeds:
        frame, _ = generate(seed)
        report = frame.attrs["rainfall_diagnostics"]
        reports.append(report)
        selected = report["selection"]
        event_rows.extend({"seed": seed, **row} for row in report["events"])
        summaries.extend({
            "seed": seed, "year_number": selected["year_number"],
            "year_start_time": selected["start_time"], "year_end_time": selected["end_time"],
            "baseline_annual_precip_mm": selected["baseline_precip_mm"],
            "median_annual_precip_mm": selected["median_precip_mm"],
            "scored_precip_mm": selected["scored_precip_mm"],
            "selected_fraction_of_scored_precip": selected["selected_fraction_of_scored_precip"],
            **row,
        } for row in report["modified_year"])
        for d, levels in report["calibration"]["return_levels_mm"].items():
            ddf_rows.extend({
                "seed": seed, "duration_days": int(d), "return_period_years": int(t),
                "depth_mm": depth, "extrapolates_record": int(t) > CALIBRATION_YEARS,
                "duration_fit_crossing": int(t) in report["calibration"]["crossing_return_levels_years"],
            } for t, depth in levels.items())
        print(f"seed={seed}; selected median-wet scored year={selected['year_number']} "
              f"({selected['start_time']} to {selected['end_time']}); "
              f"baseline precipitation={selected['baseline_precip_mm']:.3f} mm; "
              f"{report['target_exceeded_group_count']}/{report['group_count']} groups "
              f"exceed at least one {TARGET_RETURN_PERIOD}-year duration threshold")
        fraction = selected["selected_fraction_of_scored_precip"]
        if fraction is not None:
            print(f"  median annual precipitation={selected['median_precip_mm']:.3f} mm; "
                  f"selected year's share of scored precipitation={fraction:.3%}")
        print(f"  source events per group: {[g['source_event_count'] for g in report['groups']]}")
        if report["groups"]:
            print(f"  packed storms: {frame['time'].iloc[report['groups'][0]['start']]} to "
                  f"{frame['time'].iloc[report['groups'][-1]['stop'] - 1]}; "
                  f"{INTER_EVENT_DRY_DAYS} dry day between groups")
        if report["calibration"]["crossing_return_levels_years"]:
            print(f"  independent duration fits cross at return periods: "
                  f"{report['calibration']['crossing_return_levels_years']} years")
    summary = pd.DataFrame(summaries)
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"Fitted {TARGET_RETURN_PERIOD}-year depths control construction only; closure scores all complete events.")
    print("Return periods refer to individual durations under the synthetic climate, not a joint frequency.")
    print("A value above 100 years extrapolates the calibration record; no confidence intervals are computed.")
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "report.json").write_text(
            json.dumps(reports, indent=2, allow_nan=False) + "\n", encoding="utf-8",
        )
        summary.to_csv(args.output_dir / "modified_year.csv", index=False)
        pd.DataFrame(event_rows).to_csv(args.output_dir / "events.csv", index=False)
        pd.DataFrame(ddf_rows).to_csv(args.output_dir / "ddf.csv", index=False)
        print(f"Wrote diagnostics to {args.output_dir}")

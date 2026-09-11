"""Stress criteria: push a model somewhere the training record never went.

Each of these states a limit that any correct physical model satisfies
exactly, and that a model which learned a hydrograph rather than the
hydrology has no reason to. A response cannot precede its cause. Without
rain, runoff can only fall, and only as much water can drain as was held.
Under constant weather the catchment settles, and what runs off plus what
evaporates is what fell. More rain cannot mean less runoff, and it cannot
mean more runoff than the rain that was added.

None of these is a closure statement, which is why they can be scored on a
model that reports runoff alone.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.criteria.response import pick
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

STATE_VARS = ("mrso", "snw", "canopy")


def _capacity(run: RunResult, keys: list[str]) -> float:
    """The most water the catchment could have been holding, from its
    attributes plus any snow the model itself reported at the start."""
    static = run.case.static
    missing = [k for k in keys if k not in static]
    if missing:
        raise ValueError(f"the catchment attributes lack {missing}, needed for the storage bound")
    bound = float(sum(float(static[k]) for k in keys))
    return bound


@criterion("causality", paired=True)
def causality(runs: dict[str, RunResult], probe: ProbeSpec, params: dict) -> CriterionResult:
    """Nothing may change before the cause, and something must change after.

    The perturbed variant carries the control's weather with one storm
    added on one day. Up to that day the two records are identical, so
    every reported variable must be identical too, to floating point. A
    model that differs before the storm has used the storm: it saw the
    future, through a bidirectional pass, a centred filter, an attention
    window, or statistics taken over the record it was handed. After the
    storm the runoff has to respond, so that a model cannot pass by
    ignoring the rain altogether.
    """
    driver = str(params.get("driver", "pr"))
    rtol = float(params.get("rtol", 1e-6))
    min_share = float(params.get("min_share", 0.02))
    response_var = str(params.get("response", "mrro"))
    variables = list(params.get("variables", ["mrro", "evspsbl", *STATE_VARS]))

    control = make_window(pick(runs, params, "control", "control"), probe)
    pulsed = make_window(pick(runs, params, "perturbed", "pulse"), probe)
    if len(control.table) != len(pulsed.table):
        raise ValueError("the variants must have the same number of scored steps")
    if driver not in control.forcing.columns:
        raise ValueError(f"causality needs '{driver}' in the forcing")

    diff = pulsed.forcing[driver].to_numpy(dtype=float) - control.forcing[driver].to_numpy(dtype=float)
    changed = np.nonzero(np.abs(diff) > 1e-9)[0]
    if len(changed) == 0:
        raise ValueError("the variants carry the same driver; there is no cause to look for")
    cause = int(changed[0])
    if cause == 0:
        raise ValueError("the cause is the first scored step; put it inside the window so there is a before")
    added = float(control.volume(diff).sum())
    if added <= 0:
        raise ValueError(f"the perturbed variant removes {driver}; a cause must add it")
    when = str(control.forcing["time"].iloc[cause])

    present = [v for v in variables if v in control.table.columns and v in pulsed.table.columns]
    if not present:
        raise ValueError(f"causality found none of {variables} in the model result")

    deviations: dict[str, float] = {}
    for var in present:
        a = control.table[var].to_numpy(dtype=float)
        b = pulsed.table[var].to_numpy(dtype=float)
        scale = max(float(np.abs(a).mean()), 1e-12)
        deviations[var] = float(np.abs(b[:cause] - a[:cause]).max() / scale)
    worst_var, worst = max(deviations.items(), key=lambda kv: kv[1])

    failures = []
    if worst > rtol:
        failures.append(
            f"{worst_var} differs from the control by {worst:.3e} (relative) before "
            f"the storm of {when}; a model may not respond before the rain it responds to"
        )
    share = None
    if response_var in present:
        a = control.table[response_var].to_numpy(dtype=float)[cause:]
        b = pulsed.table[response_var].to_numpy(dtype=float)[cause:]
        share = float(control.volume(b - a).sum() / added)
        if share < min_share:
            failures.append(
                f"the {added:.0f} mm storm of {when} changed {response_var} by "
                f"{share:+.3f} of itself afterwards; at least {min_share:g} is expected"
            )

    ok = not failures
    return CriterionResult(
        name="causality",
        status=PASS if ok else FAIL,
        value=worst,
        threshold=rtol,
        message=(
            f"nothing moves before the storm of {when} (worst {worst:.1e} relative) and "
            f"{response_var} answers it afterwards ({share:+.3f} of the added rain)"
            if ok
            else "; ".join(failures)
        ),
        diagnostics={
            "cause_step": cause,
            "cause_time": when,
            "added_mm": added,
            "pre_cause_deviations": deviations,
            "response_share": share,
        },
    )


def _blocks(series: np.ndarray, size: int) -> np.ndarray:
    n = len(series) // size
    if n == 0:
        return np.array([series.mean()])
    return series[: n * size].reshape(n, size).mean(axis=1)


@criterion("dry_down")
def dry_down(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Without rain, runoff can only fall, and only stored water can drain.

    Two things a rainless record fixes exactly. The total that runs off
    cannot exceed what the catchment held when the rain stopped, which is
    bounded by its storage capacities. And nothing arrives to raise the
    flow, so runoff, and every storage, can only decrease. The decrease is
    judged on block means after a short lag, so that routing delay and a
    model's own numerical noise are not mistaken for water from nowhere.
    """
    driver = str(params.get("driver", "pr"))
    keys = list(params.get("capacity", ["soil_capacity_mm", "canopy_capacity_mm"]))
    lag = int(params.get("lag_days", 7))
    block_days = int(params.get("aggregate_days", 7))
    tolerance = float(params.get("tolerance", 0.05))

    w = make_window(run, probe)
    if driver not in w.forcing.columns:
        raise ValueError(f"dry_down needs '{driver}' in the forcing")
    if float(w.forcing[driver].sum()) > 1e-9:
        raise ValueError("dry_down needs a rainless scored record; the generator let rain in")
    if "mrro" not in w.table.columns:
        raise ValueError("dry_down needs 'mrro' in the model result")

    runoff = w.table["mrro"].to_numpy(dtype=float)
    total = float(w.volume(runoff).sum())
    # What the catchment held when the rain stopped: the capacities of its
    # bounded stores, plus whatever it reports in the stores that have no
    # capacity to name, snow, groundwater and water in transit. A slow
    # groundwater reservoir can legitimately hold and drain more than the
    # soil column's capacity; a model that reports it is not penalised for
    # having it, and one that does not report it is held to the capacities.
    bound = _capacity(run, keys)
    for var in ("snw", "gw", "channel"):
        if var in w.state0.index:
            bound += max(0.0, float(w.state0[var]))
    days = len(runoff) * w.dt_days

    failures = []
    if runoff.min() < -1e-9:
        failures.append(f"runoff goes negative ({runoff.min():.3g} mm/day)")
    if total > bound:
        failures.append(
            f"{total:.0f} mm ran off in {days:g} rainless days, more than the "
            f"{bound:.0f} mm the catchment could have held"
        )

    lag_rows = int(round(lag / w.dt_days))
    size = max(1, int(round(block_days / w.dt_days)))
    rises: dict[str, float] = {}
    series_by_var = {"mrro": runoff}
    for var in STATE_VARS:
        if var in w.table.columns:
            series_by_var[var] = w.table[var].to_numpy(dtype=float)
    for var, series in series_by_var.items():
        blocks = _blocks(series[lag_rows:], size)
        reference = max(abs(float(blocks[0])), 1e-12)
        increases = np.diff(blocks)
        worst_rise = float(increases.max() / reference) if len(increases) else 0.0
        rises[var] = worst_rise
        if worst_rise > tolerance:
            at = int(np.argmax(increases)) + 1
            failures.append(
                f"{var} rose by {worst_rise:.1%} of its post-lag level between "
                f"{block_days}-day blocks {at} and {at + 1} with no rain to raise it"
            )

    ok = not failures
    return CriterionResult(
        name="dry_down",
        status=PASS if ok else FAIL,
        value=total / bound if bound > 0 else None,
        threshold=1.0,
        message=(
            f"drains without rain: {total:.0f} mm over {days:g} dry days against a "
            f"{bound:.0f} mm bound, and nothing rises"
            if ok
            else "; ".join(failures)
        ),
        diagnostics={
            "runoff_total_mm": total,
            "storage_bound_mm": bound,
            "dry_days": days,
            "worst_block_rise": rises,
        },
    )


@criterion("steady_state")
def steady_state(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Under constant weather the catchment settles, and the budget balances.

    Hold rain, temperature and demand constant long enough and every store
    reaches the level at which what comes in equals what goes out. From
    then on nothing reported may vary, runoff cannot exceed the rain, and
    where evaporation is reported the two together must equal the rain,
    because storage is no longer changing. A recurrent state that drifts,
    oscillates or keys on the calendar shows here and nowhere else.
    """
    driver = str(params.get("driver", "pr"))
    last_days = int(params.get("last_days", 365))
    tolerance = float(params.get("tolerance", 0.01))
    floor = float(params.get("floor", 0.05))
    state_floor = float(params.get("state_floor", 1.0))

    w = make_window(run, probe)
    for column in (driver, "tas", "pet"):
        if column in w.forcing.columns and float(w.forcing[column].std()) > 1e-9:
            raise ValueError(f"steady_state needs constant forcing; '{column}' varies")
    rain = float(w.forcing[driver].iloc[0])
    rows = max(1, int(round(last_days / w.dt_days)))
    if rows >= len(w.table):
        raise ValueError("the scored record is not longer than the stretch judged as settled")

    ranges: dict[str, float] = {}
    means: dict[str, float] = {}
    failures = []
    for var in ("mrro", "evspsbl", *STATE_VARS):
        if var not in w.table.columns:
            continue
        tail = w.table[var].to_numpy(dtype=float)[-rows:]
        scale = max(float(np.abs(tail).mean()), state_floor if var in STATE_VARS else floor)
        rel = float((tail.max() - tail.min()) / scale)
        ranges[var], means[var] = rel, float(tail.mean())
        if rel > tolerance:
            failures.append(
                f"{var} still varies by {rel:.1%} of its level in the last {last_days} "
                "days of constant weather"
            )
    if "mrro" not in means:
        raise ValueError("steady_state needs 'mrro' in the model result")
    if means["mrro"] > rain * (1.0 + tolerance):
        failures.append(
            f"runs off {means['mrro']:.3f} mm/day under {rain:.3f} mm/day of rain at steady state"
        )
    residual = None
    if "evspsbl" in means:
        declared = 0.0
        if "gwex" in w.table.columns:
            declared = float(w.table["gwex"].to_numpy(dtype=float)[-rows:].mean())
        residual = rain + declared - means["mrro"] - means["evspsbl"]
        if abs(residual) > tolerance * rain:
            failures.append(
                f"runoff and evaporation leave {residual:+.3f} mm/day of the "
                f"{rain:.3f} mm/day rain unaccounted for at steady state"
            )

    worst_var, worst = max(ranges.items(), key=lambda kv: kv[1])
    ok = not failures
    return CriterionResult(
        name="steady_state",
        status=PASS if ok else FAIL,
        value=worst,
        threshold=tolerance,
        message=(
            f"settles: runoff {means['mrro']:.3f} mm/day under {rain:.3f} mm/day of rain, "
            f"nothing varies by more than {worst:.2%} ({worst_var}) in the last {last_days} days"
            if ok
            else "; ".join(failures)
        ),
        diagnostics={
            "rain_mm_per_day": rain,
            "last_days": last_days,
            "means": means,
            "relative_ranges": ranges,
            "budget_residual": residual,
        },
    )


@criterion("monotone_response", paired=True)
def monotone_response(runs: dict[str, RunResult], probe: ProbeSpec, params: dict) -> CriterionResult:
    """More rain cannot mean less runoff, nor more runoff than was added.

    The variants are the same record with one storm scaled up a ladder.
    Along it the integrated runoff may not fall, and each rung may not add
    more runoff than it added rain. Across the whole ladder the extra rain
    must show up: evaporation is bounded by demand and storage by capacity,
    so of a storm far larger than either, most has to run off. And on every
    rung, what runs off cannot exceed what fell plus what the catchment
    could have held. A learned response that flattens or turns down beyond
    the training range fails the first two; one that amplifies fails the
    third.
    """
    driver = str(params.get("driver", "pr"))
    var = str(params.get("variable", "mrro"))
    min_share = float(params.get("min_share", 0.1))
    max_share = float(params.get("max_share", 1.05))
    slack = float(params.get("tolerance", 0.01))
    keys = list(params.get("capacity", ["soil_capacity_mm", "canopy_capacity_mm"]))

    if len(runs) < 2:
        raise ValueError("monotone_response needs at least two rungs")
    ladder = []
    for name, run in runs.items():
        w = make_window(run, probe)
        if driver not in w.forcing.columns or var not in w.table.columns:
            raise ValueError(f"monotone_response needs '{driver}' in the forcing and '{var}' in the result")
        ladder.append((name, float(w.volume(w.forcing[driver]).sum()), float(w.volume(w.table[var]).sum()), run))
    ladder.sort(key=lambda r: r[1])

    failures = []
    rungs = []
    for (lo, p_lo, q_lo, _), (hi, p_hi, q_hi, run_hi) in zip(ladder, ladder[1:]):
        added = p_hi - p_lo
        if added <= 0:
            raise ValueError(f"variants '{lo}' and '{hi}' carry the same {driver}; the ladder has a flat rung")
        share = (q_hi - q_lo) / added
        rungs.append({"from": lo, "to": hi, "added_mm": added, "share": share})
        if share < -slack:
            failures.append(
                f"{var} fell by {q_lo - q_hi:.0f} mm when {added:.0f} mm more rain fell ('{lo}' to '{hi}')"
            )
        elif share > max_share:
            failures.append(
                f"{var} rose by {share:.2f} of the {added:.0f} mm added from '{lo}' to '{hi}', "
                "more than the rain that was added"
            )
    for name, p_total, q_total, run in ladder:
        bound = p_total + _capacity(run, keys)
        w0 = make_window(run, probe)
        if "snw" in w0.state0.index:
            bound += max(0.0, float(w0.state0["snw"]))
        if q_total > bound * (1.0 + slack):
            failures.append(
                f"'{name}' runs off {q_total:.0f} mm, more than the {p_total:.0f} mm that fell "
                f"plus the {bound - p_total:.0f} mm the catchment could have held"
            )

    (_, p_bottom, q_bottom, _), (top, p_top, q_top, _) = ladder[0], ladder[-1]
    overall = (q_top - q_bottom) / (p_top - p_bottom)
    if overall < min_share:
        failures.append(
            f"scaling the storm to '{top}' added {p_top - p_bottom:.0f} mm of rain and only "
            f"{q_top - q_bottom:.0f} mm of {var} ({overall:.2f}); a catchment cannot absorb "
            f"that much"
        )

    ok = not failures
    return CriterionResult(
        name="monotone_response",
        status=PASS if ok else FAIL,
        value=overall,
        threshold=min_share,
        message=(
            f"{var} rises with the storm on every rung and returns {overall:.2f} of the "
            f"{p_top - p_bottom:.0f} mm added at the top"
            if ok
            else "; ".join(failures)
        ),
        diagnostics={
            "rungs": rungs,
            "totals": {name: {"rain_mm": p, var: q} for name, p, q, _ in ladder},
            "overall_share": overall,
        },
    )


@criterion("response_nonnegativity", paired=True)
def response_nonnegativity(runs: dict[str, RunResult], probe: ProbeSpec, params: dict) -> CriterionResult:
    """More water can never lower the flow, at any moment.

    The extreme-rain probe checks that added rain adds runoff in total; this
    checks it step by step. From the storm onward the perturbed run's runoff
    may not fall below the control's by more than a small share of what was
    added, on any step. A network's impulse response can dip negative while
    its integral is fine, and every other criterion integrates.
    """
    driver = str(params.get("driver", "pr"))
    var = str(params.get("variable", "mrro"))
    slack = float(params.get("tolerance", 0.001))  # share of the added rain, per day

    control = make_window(pick(runs, params, "control", "control"), probe)
    pulsed = make_window(pick(runs, params, "perturbed", "pulse"), probe)
    if len(control.table) != len(pulsed.table):
        raise ValueError("the variants must have the same number of scored steps")
    diff = pulsed.forcing[driver].to_numpy(dtype=float) - control.forcing[driver].to_numpy(dtype=float)
    changed = np.nonzero(np.abs(diff) > 1e-9)[0]
    if len(changed) == 0:
        raise ValueError("the variants carry the same driver; there is no perturbation to look for")
    cause = int(changed[0])
    added = float(control.volume(diff).sum())
    if added <= 0:
        raise ValueError(f"the perturbed variant removes {driver}; it must add it")

    a = control.table[var].to_numpy(dtype=float)[cause:]
    b = pulsed.table[var].to_numpy(dtype=float)[cause:]
    dip = a - b  # positive where the perturbed run is lower
    worst = float(dip.max())
    limit = slack * added
    ok = worst <= limit
    i = int(dip.argmax()) + cause
    return CriterionResult(
        name="response_nonnegativity",
        status=PASS if ok else FAIL,
        value=worst / added if added > 0 else None,
        threshold=slack,
        message=(
            f"{var} never falls below the control after the added {added:.0f} mm "
            f"(largest dip {worst / added:.1e} of it)"
            if ok
            else f"{var} falls {worst:.3f} mm/day below the control on "
            f"{control.forcing['time'].iloc[i]} after {added:.0f} mm of rain was added; "
            "more water cannot lower the flow"
        ),
        diagnostics={"cause_step": cause, "added_mm": added, "worst_dip_mm_per_day": worst,
                     "worst_step": i},
    )


@criterion("antecedent_monotonicity", paired=True)
def antecedent_monotonicity(runs: dict[str, RunResult], probe: ProbeSpec, params: dict) -> CriterionResult:
    """The same storm on a wetter catchment yields more runoff, and no more
    than the extra water that was there.

    The `wet` variant adds rain before the storm and nothing during or after
    it. The storm is the first rain after the last step on which the
    variants differ, however many rainless days a probe leaves in between,
    and the window opens on the storm and runs `window_days` from it. It
    cannot open any earlier: the wet catchment is still draining the
    antecedent rain in those quiet days, and that recession is not runoff
    from a storm that has not yet fallen. Over the window the wet run's
    runoff must exceed the dry run's by at least a share of the storm, the
    unbroken run of rain the window opens on up to its first dry step, and
    by no more than the antecedent rain that was added. A
    memoryless model, one that runs off a fixed share of each day's rain,
    answers both storms identically and fails the first; a model that
    manufactures water fails the second.
    """
    driver = str(params.get("driver", "pr"))
    var = str(params.get("variable", "mrro"))
    min_share = float(params.get("min_share", 0.02))
    window_days = float(params.get("window_days", 30.0))

    dry = make_window(pick(runs, params, "control", "dry"), probe)
    wet = make_window(pick(runs, params, "perturbed", "wet"), probe)
    if len(dry.table) != len(wet.table):
        raise ValueError("the variants must have the same number of scored steps")
    rain = dry.forcing[driver].to_numpy(dtype=float)
    diff = wet.forcing[driver].to_numpy(dtype=float) - rain
    added_steps = np.nonzero(np.abs(diff) > 1e-9)[0]
    if len(added_steps) == 0:
        raise ValueError("the variants carry the same driver; there is no antecedent rain")
    antecedent = float(dry.volume(diff).sum())
    if antecedent <= 0:
        raise ValueError("the wet variant must add rain, not remove it")
    # The variants carry the same rain from here on, so the first of it is
    # the storm they share. It need not come on the next step.
    after = int(added_steps[-1]) + 1
    falls = np.nonzero(rain[after:] > 1e-9)[0]
    if len(falls) == 0:
        raise ValueError("no storm follows the antecedent rain in the scored window")
    storm = after + int(falls[0])
    n = max(1, int(round(window_days / dry.dt_days)))
    stop = min(len(dry.table), storm + n)
    lull = np.nonzero(rain[storm:stop] <= 1e-9)[0]
    storm_end = storm + int(lull[0]) if len(lull) else stop
    storm_mm = float(dry.volume(rain[storm:storm_end]).sum())

    q_dry = float(dry.volume(dry.table[var].to_numpy(dtype=float)[storm:stop]).sum())
    q_wet = float(wet.volume(wet.table[var].to_numpy(dtype=float)[storm:stop]).sum())
    extra = q_wet - q_dry
    failures = []
    if extra < min_share * storm_mm:
        failures.append(
            f"the storm of {storm_mm:.0f} mm ran off {extra:+.1f} mm more on the wetter catchment; "
            f"at least {min_share:g} of the storm is expected"
        )
    if extra > antecedent + 1e-6:
        failures.append(
            f"the wetter catchment ran off {extra:.1f} mm more but only {antecedent:.0f} mm more "
            "had fallen on it"
        )
    ok = not failures
    return CriterionResult(
        name="antecedent_monotonicity",
        status=PASS if ok else FAIL,
        value=extra / storm_mm,
        threshold=min_share,
        message=(
            f"the wetter catchment runs off {extra:.1f} mm more from the {storm_mm:.0f} mm storm "
            f"({extra / storm_mm:.2f} of it), within the {antecedent:.0f} mm it had been given"
            if ok else "; ".join(failures)
        ),
        diagnostics={"antecedent_mm": antecedent, "storm_step": storm,
                     "storm_time": str(dry.forcing["time"].iloc[storm]), "storm_mm": storm_mm,
                     "runoff_dry_mm": q_dry, "runoff_wet_mm": q_wet},
    )


@criterion("phase_invariance", paired=True)
def phase_invariance(runs: dict[str, RunResult], probe: ProbeSpec, params: dict) -> CriterionResult:
    """Rain that would have been snow is still the same water.

    The `warm` variant lifts every sub-freezing day above the snow threshold
    and changes nothing else, demand included. Timing moves; mass does not:
    integrated over the record, runoff (and evaporation, where reported) may
    differ between the variants by at most a share of the rain. A snow store
    that loses water it never reports shows here as a difference that only
    exists when it snows.
    """
    threshold = float(params.get("threshold", 0.05))
    volumes = list(params.get("volumes", ["mrro", "evspsbl"]))
    control = make_window(pick(runs, params, "control", "control"), probe)
    warm = make_window(pick(runs, params, "perturbed", "warm"), probe)
    if len(control.table) != len(warm.table):
        raise ValueError("the variants must have the same number of scored steps")
    if not np.allclose(control.forcing["pr"], warm.forcing["pr"]):
        raise ValueError("the variants must carry the same precipitation")
    rain = float(control.volume(control.forcing["pr"]).sum())
    if rain <= 0:
        raise ValueError("no rain fell")
    shares = {}
    for var in volumes:
        if var in control.table.columns and var in warm.table.columns:
            a = float(control.volume(control.table[var]).sum())
            b = float(warm.volume(warm.table[var]).sum())
            shares[var] = (b - a) / rain
    if not shares:
        raise ValueError(f"phase_invariance found none of {volumes} in the model result")
    worst_var, worst = max(shares.items(), key=lambda kv: abs(kv[1]))
    ok = abs(worst) <= threshold
    return CriterionResult(
        name="phase_invariance",
        status=PASS if ok else FAIL,
        value=abs(worst),
        threshold=threshold,
        message=(
            f"turning snow into rain moves the integrated volumes by at most {abs(worst):.1%} of "
            f"the rain ({worst_var})"
            if ok
            else f"{worst_var} differs by {worst:+.1%} of the rain when snow falls as rain "
            f"(limit {threshold:.0%}); the same water fell either way"
        ),
        diagnostics={"shares": shares, "rain_mm": rain},
    )

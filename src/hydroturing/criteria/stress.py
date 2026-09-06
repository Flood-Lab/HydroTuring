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
        residual = rain - means["mrro"] - means["evspsbl"]
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

"""Physical bounds on states and fluxes.

This is the criterion that catches closure by construction, which is the
attack randomised forcing cannot touch. A model that solves for storage as
whatever balances the budget will pass any closure check, forever, on every
seed. What it cannot do is keep that storage physical: the invented storage
drifts, and over a long random sequence it leaves the range a real soil
column or snowpack can occupy.

That is why the model contract requires absolute storage states rather than
tendencies. Reporting dS/dt would hide exactly the quantity this checks.
"""

from __future__ import annotations

import math

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

EPS = 1e-6  # mm, absorbs floating point noise without hiding a real excursion

# A capacity in static.json describes the catchment to a model that reads it. A
# model that learns its own, as a differentiable HBV does, cannot be held to a
# number it never used, so for it that upper bound is N/A. Dropping the bound
# outright would let invented storage through just by not declaring the input,
# so a store past the capacity must still not run away: the means of its first
# and last whole years, up to a third of the record, may differ by at most this
# share of the range it covered before those last years. dHBV2 and the MCP models
# reach 0.162 on their gate seeds. reference_cheater has a lowest of 0.399 and
# 0.556 over 205 seeds on the closure and ungauged probes, and is still caught on
# the spin-up probe's gate seeds.
RUNAWAY_SHARE = 0.25
# Under two years a seasonal or storm-driven rise looks like a drift, so a record
# that short cannot show the store is stationary, and the bound stays. Otherwise
# invented storage would escape on the default flood window, which a model gets
# by asking for nothing. A model that wants the bound dropped asks for the record:
# `window_days: full`.
RUNAWAY_MIN_DAYS = 730.0
YEAR_DAYS = 365.25
RATIO_TOLERANCE = 1e-6


def _resolve(bound, static: dict) -> float:
    """Bounds may be a number, null (unbounded), or a key into static.json."""
    if bound is None:
        return math.inf
    if isinstance(bound, (int, float)):
        return float(bound)
    if isinstance(bound, str):
        if bound not in static:
            raise ValueError(f"state bound '{bound}' is not in the case attributes")
        return float(static[bound])
    raise ValueError(f"cannot interpret state bound {bound!r}")


def _unread_capacity(spec, run: RunResult) -> str | None:
    """The static key behind an upper bound, if the model does not consume it.

    A run with no manifest, as the tests build by hand, is held to the bound.
    """
    if run.model is None or not isinstance(spec, (list, tuple)):
        return None
    key = spec[1]
    if isinstance(key, str) and key not in run.model.needs_static + run.model.uses_static:
        return key
    return None


def _runaway_share(values: np.ndarray, dt_days: float) -> float:
    """How far a store's last years sit from its first, against what it did before them.

    Whole years, so a season cannot read as a drift, and the range before the last
    block, so a store that accelerates cannot hide behind the height it reached.
    """
    if not np.isfinite(values).all():
        return math.inf
    years = max(1, int(len(values) * dt_days / 3 / YEAR_DAYS))
    block = round(years * YEAR_DAYS / dt_days)
    earlier = values[:-block]
    drift = abs(float(values[-block:].mean() - values[:block].mean()))
    span = float(earlier.max() - earlier.min())
    if drift == 0.0:
        return 0.0
    return drift / span if span > 0.0 else math.inf


@criterion("state_bounds")
def state_bounds(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Every reported storage must stay inside its physical range."""
    # When a paired probe requests all variants, score the labelled
    # evaluation cycle rather than only the control's full record.
    w = make_window(run, probe, params.get("phase"))
    static = run.case.static
    violations = []
    not_applied = []
    worst_excess, worst_var = 0.0, None

    for var, spec in params.items():
        if var not in w.table.columns:
            continue
        lo = _resolve(spec[0], static) if isinstance(spec, (list, tuple)) else 0.0
        hi = _resolve(spec[1], static) if isinstance(spec, (list, tuple)) else math.inf
        values = np.asarray(w.table[var], dtype=float)

        runaway = ""
        unread = _unread_capacity(spec, run)
        if unread and values.max() > hi + EPS:
            if len(values) * w.dt_days < RUNAWAY_MIN_DAYS:
                runaway = (
                    f", and {len(values) * w.dt_days:g} days cannot show it does not run "
                    f"away though the model does not consume {unread}"
                )
            elif (share := _runaway_share(values, w.dt_days)) > RUNAWAY_SHARE:
                runaway = (
                    f", and runs away by {share:.2f} of its earlier range "
                    f"(limit {RUNAWAY_SHARE:g}) though it does not consume {unread}"
                )
            else:
                not_applied.append(
                    f"{var} upper bound N/A, the model does not consume {unread} "
                    f"(peak {values.max():.3g} against {hi:g})"
                )
                hi = math.inf

        below = float(np.maximum(lo - values.min(), 0.0))
        above = float(np.maximum(values.max() - hi, 0.0)) if math.isfinite(hi) else 0.0
        excess = max(below, above)

        if excess > EPS:
            n_bad = int(((values < lo - EPS) | (values > hi + EPS)).sum())
            violations.append(
                f"{var} leaves [{lo:g}, {hi:g}] on {n_bad} steps "
                f"(range {values.min():.3g} to {values.max():.3g}){runaway}"
            )
            if excess > worst_excess:
                worst_excess, worst_var = excess, var

    ok = not violations
    return CriterionResult(
        name="state_bounds",
        status=PASS if ok else FAIL,
        value=worst_excess,
        threshold=EPS,
        message=(
            "; ".join(["all storages stay physical", *not_applied])
            if ok else "; ".join(violations)
        ),
        diagnostics={
            "worst_variable": worst_var,
            "n_violations": len(violations),
            "bounds_not_applied": not_applied,
        },
    )


@criterion("et_plausible")
def et_plausible(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Cumulative evapotranspiration must be non-negative and bounded by potential ET.

    Cheap, and it catches the degenerate escape of setting ET equal to
    precipitation so that the water balance closes with no runoff at all. A row
    may go negative: `evspsbl` is net of dew and frost where the model nets them.
    """
    max_ratio = float(params.get("max_ratio_to_pet", 1.0))
    pet_var = params.get("pet_variable", "pet")
    w = make_window(run, probe)

    if "evspsbl" not in w.table.columns:
        raise ValueError("et_plausible needs 'evspsbl' in the model result")
    et = np.asarray(w.table["evspsbl"], dtype=float)

    if et.sum() < -EPS:
        return CriterionResult(
            name="et_plausible", status=FAIL, value=float(et.sum()), threshold=0.0,
            message="cumulative evapotranspiration is negative",
        )

    if pet_var not in w.forcing.columns:
        return CriterionResult(
            name="et_plausible", status=PASS,
            message="cumulative ET non-negative; no potential ET in the forcing to bound it against",
        )

    pet = np.asarray(w.forcing[pet_var], dtype=float)
    cum_ratio = float(et.sum() / max(pet.sum(), 1e-12))
    # A model capped at PET in single precision can sum to a few parts in 1e8 over it.
    ok = cum_ratio <= max_ratio + RATIO_TOLERANCE

    return CriterionResult(
        name="et_plausible",
        status=PASS if ok else FAIL,
        value=cum_ratio,
        threshold=max_ratio,
        message=(
            f"cumulative ET is {cum_ratio:.3f} of potential ET "
            f"(limit {max_ratio:g})"
        ),
        diagnostics={
            "n_steps_over_pet": int((et > pet + EPS).sum()),
            "cumulative_et_mm": float(et.sum()),
            "cumulative_pet_mm": float(pet.sum()),
        },
    )

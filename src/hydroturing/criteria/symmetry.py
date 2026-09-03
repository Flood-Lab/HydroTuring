"""Invariance: the answer must not depend on how the question was written.

Physics has symmetries. A catchment does not care what year the calendar says
it is, and doubling its area doubles the discharge in cubic metres per second
while leaving every depth in millimetres exactly where it was. A model that
has learned the physics inherits those symmetries for free. A model that has
learned the training set does not, and breaks them in ways that are obvious
once you look and invisible if you never run the same case twice.

This is cheap to check and unusually hard to game, because there is no
tolerance to tune towards: the two runs either agree or they do not.

Two kinds of expectation, declared per variable:

  unchanged   the transform must not move it at all
  scaled      it must move by exactly the declared factor

Comparison is relative to the mean magnitude of the control run, so a variable
that is legitimately near zero for most of the record is not judged against
its own noise.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL,
    PASS,
    CriterionResult,
    criterion,
    make_window,
)
from hydroturing.criteria.response import pick
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


def _deviation(control: np.ndarray, other: np.ndarray, factor: float) -> float:
    """Worst relative departure of `other` from `factor` times `control`."""
    expected = factor * control
    scale = max(float(np.abs(expected).mean()), 1e-12)
    return float(np.abs(other - expected).max() / scale)


@criterion("invariance", paired=True)
def invariance(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    """The transformed run must reproduce the control, up to the declared factor."""
    rtol = float(params.get("rtol", 1e-6))
    unchanged = list(params.get("unchanged", []))
    scaled = dict(params.get("scaled", {}))

    if not unchanged and not scaled:
        raise ValueError(
            "invariance asserts nothing: declare `unchanged` variables, "
            "`scaled` ones, or both"
        )

    control = make_window(pick(runs, params, "control", "control"), probe)
    transformed = make_window(pick(runs, params, "transformed", "transformed"), probe)

    if len(control.table) != len(transformed.table):
        raise ValueError(
            f"the two variants produced windows of different length "
            f"({len(control.table)} and {len(transformed.table)}); an "
            "invariance transform must preserve the number of steps"
        )

    expectations = [(var, 1.0) for var in unchanged]
    expectations += [(var, float(f)) for var, f in scaled.items()]

    worst_var, worst = None, 0.0
    deviations: dict[str, float] = {}
    for var, factor in expectations:
        if var not in control.table.columns or var not in transformed.table.columns:
            raise ValueError(f"invariance needs '{var}' in the model result")
        deviation = _deviation(
            np.asarray(control.table[var], dtype=float),
            np.asarray(transformed.table[var], dtype=float),
            factor,
        )
        deviations[var] = deviation
        if deviation > worst:
            worst_var, worst = var, deviation

    ok = worst <= rtol
    return CriterionResult(
        name="invariance",
        status=PASS if ok else FAIL,
        value=worst,
        threshold=rtol,
        message=(
            f"the transform leaves the answer unchanged (worst departure "
            f"{worst:.3e})"
            if ok
            else (
                f"'{worst_var}' departs from its expected value by {worst:.3e} "
                f"(relative, limit {rtol:g})"
            )
        ),
        diagnostics={
            "worst_variable": worst_var,
            "deviations": deviations,
        },
    )


@criterion("resolution_invariance", paired=True)
def resolution_invariance(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    """Integrated volumes must not depend on the step the weather was given at.

    The same weather twice: once at a fine step, once aggregated to a coarse
    one so that every hour or day carries the same water. Water is
    conserved under that aggregation, so over the same stretch of time the
    runoff volume, the evapotranspiration volume and the storage a model
    ends with must agree between the two runs, to within a share of the
    precipitation that fell.

    What is not asserted matters as much. Peaks, stages, velocity heads and
    momentum fluxes are nonlinear in the flow and legitimately change with
    the step; by Jensen's inequality aggregation can only attenuate them,
    never sharpen them, and that ordering is a separate test that needs a
    routed channel. This criterion is about mass alone.

    The failure it exists to catch: a model whose arithmetic assumes a
    particular step. Given minutes it drains its stores as if each minute
    were an hour, and its volumes come out different at every resolution it
    is run at.
    """
    threshold = float(params.get("threshold", 0.05))
    volumes = list(params.get("volumes", ["mrro", "evspsbl"]))
    states = list(params.get("states", probe.requires_states))
    driver = str(params.get("driver", "pr"))

    fine = make_window(pick(runs, params, "fine", "fine"), probe)
    coarse = make_window(pick(runs, params, "coarse", "coarse"), probe)

    fine_days = len(fine.table) * fine.dt_days
    coarse_days = len(coarse.table) * coarse.dt_days
    if abs(fine_days - coarse_days) > max(fine.dt_days, coarse.dt_days):
        raise ValueError(
            f"the two variants cover different stretches of time ({fine_days:g} "
            f"and {coarse_days:g} days); the coarse variant must be the fine one "
            "aggregated, not a different record"
        )
    if driver not in fine.forcing.columns or driver not in coarse.forcing.columns:
        raise ValueError(f"resolution_invariance needs '{driver}' in the forcing")

    total_fine = float(fine.volume(fine.forcing[driver]).sum())
    total_coarse = float(coarse.volume(coarse.forcing[driver]).sum())
    if total_fine <= 0:
        return CriterionResult(
            name="resolution_invariance", status=FAIL,
            message=f"no {driver} fell in the scored window; case is degenerate",
        )
    if abs(total_fine - total_coarse) > 1e-4 * total_fine:
        raise ValueError(
            f"the variants do not carry the same water ({total_fine:.4f} and "
            f"{total_coarse:.4f} mm of {driver}); aggregation must preserve totals"
        )

    deviations: dict[str, float] = {}
    totals: dict[str, dict[str, float]] = {"fine": {}, "coarse": {}}
    for var in volumes:
        if var not in fine.table.columns or var not in coarse.table.columns:
            raise ValueError(f"resolution_invariance needs '{var}' in the model result")
        v_fine = float(fine.volume(fine.table[var]).sum())
        v_coarse = float(coarse.volume(coarse.table[var]).sum())
        totals["fine"][var], totals["coarse"][var] = v_fine, v_coarse
        deviations[var] = abs(v_fine - v_coarse) / total_fine
    for var in states:
        if var not in fine.table.columns or var not in coarse.table.columns:
            continue
        s_fine = float(fine.table[var].iloc[-1])
        s_coarse = float(coarse.table[var].iloc[-1])
        totals["fine"][f"{var}_end"], totals["coarse"][f"{var}_end"] = s_fine, s_coarse
        deviations[f"{var}_end"] = abs(s_fine - s_coarse) / total_fine

    worst_var, worst = max(deviations.items(), key=lambda kv: kv[1])
    ok = worst <= threshold
    fine_step = pick(runs, params, "fine", "fine").case.timestep
    coarse_step = pick(runs, params, "coarse", "coarse").case.timestep
    return CriterionResult(
        name="resolution_invariance",
        status=PASS if ok else FAIL,
        value=worst,
        threshold=threshold,
        message=(
            f"integrated volumes agree between {fine_step} and {coarse_step} "
            f"(worst {worst_var} {worst:.3%} of {driver})"
            if ok
            else (
                f"{worst_var} differs by {worst:.2%} of {driver} between "
                f"{fine_step} and {coarse_step} (limit {threshold:.0%})"
            )
        ),
        diagnostics={
            "fine_step": fine_step,
            "coarse_step": coarse_step,
            "driver_total_mm": total_fine,
            "deviations": deviations,
            "totals": totals,
        },
    )

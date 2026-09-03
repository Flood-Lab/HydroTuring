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

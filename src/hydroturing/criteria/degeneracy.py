"""Anti-degeneracy: the model has to actually be doing hydrology.

A pure closure test has a trivial solution. Set runoff to zero, set
evapotranspiration equal to precipitation, hold storage constant, and the
water balance closes exactly on every seed forever. Such a model conserves
water perfectly and knows nothing.

Randomised forcing does not help here either, because the degenerate answer
is correct for every possible forcing. The defence has to assert that the
model produces a non-trivial partition and responds to what it is given.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@criterion("non_degenerate")
def non_degenerate(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Partition, variability, and response to forcing must all be non-trivial."""
    ratio_lo, ratio_hi = params.get("runoff_ratio", [0.02, 0.98])
    min_cv = float(params.get("min_flux_cv", 0.1))
    # `null` turns the response check off. It assumes runoff is driven by
    # recent rainfall, which is not true of a snow-dominated catchment, where
    # runoff follows melt timing instead. A probe spanning that range has to be
    # able to say the test does not apply rather than fail honest models with
    # it.
    min_response = params.get("min_response", 0.05)
    min_response = None if min_response is None else float(min_response)
    cv_vars = params.get("cv_variables", ["mrro", "evspsbl"])

    w = make_window(run, probe)
    failures = []
    diagnostics: dict[str, float] = {}

    pr = np.asarray(w.forcing["pr"], dtype=float)
    total_pr = float(pr.sum())

    # The runoff ratio is a climatological statistic. It says something about
    # the partition only once storage has cycled, which takes a year; over a
    # shorter scored window it measures the storage change instead. A melt
    # flood returns several times the rain that fell in that month and a
    # dry-down returns almost none, and neither is degenerate. So the bounds
    # are applied to a window of at least a year and reported, not judged,
    # on a flood-event window.
    scored_days = len(w.table) * w.dt_days
    if "mrro" in w.table.columns and total_pr > 0:
        runoff_ratio = float(np.asarray(w.table["mrro"], dtype=float).sum() / total_pr)
        diagnostics["runoff_ratio"] = runoff_ratio
        if scored_days < 365:
            diagnostics["runoff_ratio_check"] = (
                f"not applied: scored window is {scored_days:g} days, under a year"
            )
        elif not (ratio_lo <= runoff_ratio <= ratio_hi):
            failures.append(
                f"runoff ratio {runoff_ratio:.4f} outside [{ratio_lo:g}, {ratio_hi:g}]"
            )

    for var in cv_vars:
        if var not in w.table.columns:
            continue
        values = np.asarray(w.table[var], dtype=float)
        mean = float(values.mean())
        cv = float(values.std() / mean) if abs(mean) > 1e-12 else 0.0
        diagnostics[f"cv_{var}"] = cv
        if cv < min_cv:
            failures.append(f"{var} is nearly constant (cv {cv:.4f} < {min_cv:g})")

    # A model that ignores its forcing entirely is degenerate even if its
    # fluxes happen to vary. Correlate runoff against precipitation smoothed
    # over a short window, so that legitimate baseflow lag is not punished.
    if min_response is not None and "mrro" in w.table.columns and len(pr) > 30:
        import pandas as pd

        smooth_pr = pd.Series(pr).rolling(7, min_periods=1).mean().to_numpy()
        runoff = np.asarray(w.table["mrro"], dtype=float)
        if runoff.std() > 1e-12 and smooth_pr.std() > 1e-12:
            response = float(np.corrcoef(smooth_pr, runoff)[0, 1])
        else:
            response = 0.0
        diagnostics["runoff_precip_correlation"] = response
        if response < min_response:
            failures.append(
                f"runoff barely responds to precipitation (r {response:.3f} < {min_response:g})"
            )

    ok = not failures
    return CriterionResult(
        name="non_degenerate",
        status=PASS if ok else FAIL,
        value=diagnostics.get("runoff_ratio"),
        message="partition and variability are non-trivial" if ok else "; ".join(failures),
        diagnostics=diagnostics,
    )

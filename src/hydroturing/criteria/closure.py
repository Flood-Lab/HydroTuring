"""Conservation closure, and the fidelity check that protects its denominator.

The 5 percent engineering rule is the accept threshold: if the budget closes
to better than 5 percent of the driving flux, it counts as closed. Two notes
that belong in the methods section of any paper using this.

First, 5 percent is loose next to what a physics model achieves, which is
nearer 1e-6. The rule states engineering acceptability, not numerical rigour,
and the report always carries the measured residual alongside the verdict.

Second, a percentage needs a denominator that does not pass through zero.
Precipitation and channel inflow are strictly non-negative, so they work
directly. Net radiation changes sign every night, so the energy form
accumulates |Rn| and additionally applies an absolute floor.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL, PASS, CriterionResult, criterion, make_window, reported_states,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

DENOMINATORS = {
    "sum_pr": ("pr", False),
    "sum_abs_rn": ("rn", True),
    "sum_inflow": ("q_in", False),
}


@criterion("closure")
def closure(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Cumulative budget residual as a share of the driving flux."""
    threshold = float(params.get("threshold", 0.05))
    floor = params.get("floor")
    denom_key = params.get("denominator", "sum_pr")
    if denom_key not in DENOMINATORS:
        raise ValueError(f"unknown denominator '{denom_key}'")
    forcing_var, take_abs = DENOMINATORS[denom_key]

    w = make_window(run, probe)
    if forcing_var not in w.forcing.columns:
        raise ValueError(
            f"closure denominator '{denom_key}' needs forcing column "
            f"'{forcing_var}', which this probe's generator does not produce"
        )

    # The driver is taken from the forcing, never from what the model echoed
    # back. A model that quietly rescales its input is caught separately by
    # forcing_fidelity, not by silently changing the denominator here.
    drive = w.volume(w.forcing[forcing_var])
    if take_abs:
        drive = np.abs(drive)

    sinks = params.get("sinks", ["evspsbl", "mrro"])
    outflow = np.zeros(len(w.table))
    for var in sinks:
        if var not in w.table.columns:
            raise ValueError(f"closure needs '{var}' in the model result")
        outflow += w.volume(w.table[var])

    # A model with an explicit exchange with the outside, a regional
    # groundwater term, an inter-basin transfer, may declare it as `gwex`
    # (positive into the catchment). Declared, it is a source in the budget
    # and the budget can close; hidden, it is the residual. The denominator
    # stays the rain, so a declared source does not dilute the residual.
    sources = params.get("sources", ["gwex"])
    declared = np.zeros(len(w.table))
    for var in sources:
        if var in w.table.columns:
            declared += w.volume(w.table[var])
    drive = drive + declared

    states = reported_states(w, probe)
    storage = w.storage(states)
    storage_change = float(storage[-1]) - w.storage_initial(states)

    step_residual = drive - outflow - np.diff(storage, prepend=w.storage_initial(states))
    cumulative = float(drive.sum() - outflow.sum() - storage_change)

    total_drive = float((drive - declared).sum())
    if total_drive <= 0:
        return CriterionResult(
            name="closure", status=FAIL,
            message=f"denominator {denom_key} accumulated to zero; case is degenerate",
        )

    relative = abs(cumulative) / total_drive
    ok = relative <= threshold
    if floor is not None and not ok:
        # The absolute floor rescues a case where the denominator is small but
        # the residual is physically negligible, which is the nighttime
        # problem in the energy budget.
        mean_abs = float(np.abs(step_residual).mean())
        if mean_abs <= float(floor):
            ok = True

    # Closing to machine precision on every step is not physics, it is
    # arithmetic: a model that solves for one budget term as the residual
    # produces exactly this. Flag it rather than reward it.
    suspicious = bool(relative < 1e-10 and np.abs(step_residual).max() < 1e-9)

    return CriterionResult(
        name="closure",
        status=PASS if ok else FAIL,
        value=relative,
        threshold=threshold,
        message=(
            f"cumulative residual {relative:.4%} of {denom_key} "
            f"(limit {threshold:.1%})"
        ),
        diagnostics={
            "cumulative_residual": cumulative,
            "denominator_total": total_drive,
            "declared_sources_total": float(declared.sum()),
            "storage_change": storage_change,
            "max_step_residual": float(np.abs(step_residual).max()),
            "mean_step_residual": float(np.abs(step_residual).mean()),
            "suspicious_exact": suspicious,
        },
    )


@criterion("paired_closure", paired=True)
def paired_closure(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    """Require closure in each named variant, without changing single-run scoring.

    The ordinary ``closure`` criterion still judges the control run alone.
    Paired probes opt into this wrapper when both budgets must be checked;
    the existing residual, denominator and floor rules apply independently
    to every selected run.
    """
    variants = params.get("variants", probe.variants)
    if not isinstance(variants, (list, tuple)) or len(variants) < 2:
        raise ValueError("paired_closure needs at least two variant names")
    if any(not isinstance(name, str) or not name for name in variants):
        raise ValueError("paired_closure variant names must be nonempty strings")
    if len(set(variants)) != len(variants):
        raise ValueError("paired_closure variant names must be unique")
    missing = [name for name in variants if name not in runs]
    if missing:
        raise ValueError(f"paired_closure is missing requested variants: {missing}")

    closure_params = {key: value for key, value in params.items() if key != "variants"}
    results = {name: closure(runs[name], probe, closure_params) for name in variants}
    failed = [name for name, result in results.items() if not result.passed]
    values = [result.value for result in results.values() if result.value is not None]
    return CriterionResult(
        name="paired_closure",
        status=FAIL if failed else PASS,
        value=max(values, default=None),
        threshold=float(params.get("threshold", 0.05)),
        message=(
            "; ".join(f"{name}: {results[name].message}" for name in failed)
            if failed
            else f"water budgets close in all variants: {', '.join(variants)}"
        ),
        diagnostics={
            "variants": {
                name: {
                    "status": result.status,
                    "value": result.value,
                    "threshold": result.threshold,
                    "message": result.message,
                    "diagnostics": result.diagnostics,
                }
                for name, result in results.items()
            },
            "failed_variants": failed,
            "suspicious_exact": any(
                result.diagnostics.get("suspicious_exact", False)
                for result in results.values()
            ),
        },
    )


@criterion("forcing_fidelity")
def forcing_fidelity(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """The model must report back the forcing it was actually given.

    Without this, a model could rescale its input precipitation and close a
    budget against its own private version of the driver.
    """
    rtol = float(params.get("rtol", 1e-6))
    variables = params.get("variables", ["pr"])
    w = make_window(run, probe)

    worst_var, worst = None, 0.0
    for var in variables:
        if var not in w.table.columns or var not in w.forcing.columns:
            continue
        given = np.asarray(w.forcing[var], dtype=float)
        echoed = np.asarray(w.table[var], dtype=float)
        scale = max(float(np.abs(given).mean()), 1e-12)
        deviation = float(np.abs(echoed - given).max() / scale)
        if deviation > worst:
            worst_var, worst = var, deviation

    ok = worst <= rtol
    return CriterionResult(
        name="forcing_fidelity",
        status=PASS if ok else FAIL,
        value=worst,
        threshold=rtol,
        message=(
            "reported forcing matches the input"
            if ok
            else f"'{worst_var}' deviates from the given forcing by {worst:.3e} (relative)"
        ),
        diagnostics={"worst_variable": worst_var},
    )

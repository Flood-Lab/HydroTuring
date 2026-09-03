"""Does conservation survive outside the range the model has seen?

A model can learn closure as a statistical regularity of its training
distribution rather than as a structural property, and nothing in a
single-window budget check will notice: the ordinary years dominate the
cumulative residual and the one record storm disappears into it.

`regime_transfer` scores the budget separately over labelled stretches of the
record and compares them. The generator marks each step with a `_regime`
column — `_` because the harness strips those columns before staging, so the
model is never handed a label saying which part it is being tested on. It has
to be told apart by its physics, from the forcing alone.

Two conditions, both necessary. The out-of-range stretch must close on its own
terms, and it must not close conspicuously worse than the in-range stretch. The
first alone lets a model with a large uniform error hide; the second alone is
meaningless for a model whose reference residual is already floating-point
noise, which is why the degradation test carries a floor.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL,
    PASS,
    CriterionResult,
    criterion,
    make_window,
    segments,
    storage_at,
)
from hydroturing.criteria.closure import DENOMINATORS
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


def _block_residual(w, states, forcing_var, take_abs, sinks, start, stop):
    """Cumulative budget residual and driver total over one block of steps."""
    drive = w.volume(w.forcing[forcing_var].to_numpy()[start:stop])
    if take_abs:
        drive = np.abs(drive)

    outflow = np.zeros(stop - start)
    for var in sinks:
        if var not in w.table.columns:
            raise ValueError(f"regime_transfer needs '{var}' in the model result")
        outflow += w.volume(w.table[var].to_numpy()[start:stop])

    storage_end = float(w.storage(states)[stop - 1]) if states else 0.0
    storage_start = storage_at(w, states, start)

    residual = float(drive.sum() - outflow.sum() - (storage_end - storage_start))
    return residual, float(drive.sum())


@criterion("regime_transfer")
def regime_transfer(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """The budget must close out of range as well as it closes in range."""
    reference_label = str(params.get("reference", "reference"))
    extrapolation_label = str(params.get("extrapolation", "extrapolation"))
    threshold = float(params.get("threshold", 0.05))
    max_degradation = float(params.get("max_degradation", 3.0))
    # Without this floor the degradation test is unusable: an exact model has a
    # reference residual near 1e-15, and three times nothing is nothing.
    degradation_floor = float(params.get("degradation_floor", 0.005))
    sinks = params.get("sinks", ["evspsbl", "mrro"])

    denom_key = params.get("denominator", "sum_pr")
    if denom_key not in DENOMINATORS:
        raise ValueError(f"unknown denominator '{denom_key}'")
    forcing_var, take_abs = DENOMINATORS[denom_key]

    w = make_window(run, probe)
    if forcing_var not in w.forcing.columns:
        raise ValueError(
            f"regime_transfer denominator '{denom_key}' needs forcing column "
            f"'{forcing_var}', which this probe's generator does not produce"
        )

    states = probe.requires_states
    blocks = segments(w)
    found = {label for label, _, _ in blocks}
    for needed in (reference_label, extrapolation_label):
        if needed not in found:
            raise ValueError(
                f"regime_transfer expects a '{needed}' regime in the scored "
                f"window; the generator labelled {sorted(found)}"
            )

    totals = {reference_label: [0.0, 0.0], extrapolation_label: [0.0, 0.0]}
    for label, start, stop in blocks:
        if label not in totals:
            continue
        residual, drive = _block_residual(
            w, states, forcing_var, take_abs, sinks, start, stop
        )
        # Residuals are accumulated in absolute value across blocks so that a
        # gain in one stretch cannot cancel a loss in another and report a
        # model as conservative when it is merely inconsistent.
        totals[label][0] += abs(residual)
        totals[label][1] += drive

    ref_residual, ref_drive = totals[reference_label]
    ext_residual, ext_drive = totals[extrapolation_label]

    if ext_drive <= 0:
        return CriterionResult(
            name="regime_transfer", status=FAIL,
            message=(
                f"the '{extrapolation_label}' regime accumulated no "
                f"{denom_key}; there is nothing out of range to test against"
            ),
        )

    ext_relative = ext_residual / ext_drive
    ref_relative = ref_residual / ref_drive if ref_drive > 0 else 0.0
    allowed = max(max_degradation * ref_relative, degradation_floor)

    failures = []
    if ext_relative > threshold:
        failures.append(
            f"residual out of range is {ext_relative:.4%} of {denom_key} "
            f"(limit {threshold:.1%})"
        )
    if ext_relative > allowed:
        ratio = ext_relative / ref_relative if ref_relative > 0 else float("inf")
        failures.append(
            f"closure degrades {ratio:.1f}x out of range "
            f"({ref_relative:.4%} in range, {ext_relative:.4%} out, "
            f"limit {allowed:.4%})"
        )

    ok = not failures
    return CriterionResult(
        name="regime_transfer",
        status=PASS if ok else FAIL,
        value=ext_relative,
        threshold=min(threshold, allowed),
        message=(
            f"closure holds out of range ({ref_relative:.4%} in range, "
            f"{ext_relative:.4%} out)"
            if ok
            else "; ".join(failures)
        ),
        diagnostics={
            "reference_relative": ref_relative,
            "extrapolation_relative": ext_relative,
            "degradation_allowed": allowed,
            "reference_denominator": ref_drive,
            "extrapolation_denominator": ext_drive,
            "n_blocks": len(blocks),
        },
    )

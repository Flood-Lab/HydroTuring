"""Storage-trend checks for repeated-forcing water probes."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from hydroturing.criteria.base import (
    FAIL,
    PASS,
    CriterionResult,
    criterion,
    make_window,
    reported_states,
    storage_at,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@criterion("total_storage_drift")
def total_storage_drift(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Check net change in all reported stores over the final repeated block."""
    unknown = set(params) - {"block_days", "relative_to_precipitation", "absolute_tolerance_mm"}
    if unknown:
        raise ValueError(f"total_storage_drift: unknown parameters {sorted(unknown)}")

    block_days = float(params.get("block_days", 1825))
    relative = float(params.get("relative_to_precipitation", 2e-4))
    absolute = float(params.get("absolute_tolerance_mm", 1e-6))
    if not math.isfinite(block_days) or block_days <= 0:
        raise ValueError("total_storage_drift block_days must be finite and positive")
    if not math.isfinite(relative) or relative < 0:
        raise ValueError("total_storage_drift relative_to_precipitation must be finite and nonnegative")
    if not math.isfinite(absolute) or absolute < 0:
        raise ValueError("total_storage_drift absolute_tolerance_mm must be finite and nonnegative")

    w = make_window(run, probe)
    if "pr" not in w.forcing.columns:
        raise ValueError("total_storage_drift needs supplied forcing column 'pr'")

    block_steps_float = block_days / w.dt_days
    block_steps = int(round(block_steps_float))
    if not math.isclose(block_steps_float, block_steps, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("total_storage_drift block_days must align with the run timestep")
    if block_steps <= 0:
        raise ValueError("total_storage_drift block must contain at least one step")
    if len(w.table) < block_steps:
        return CriterionResult(
            name="total_storage_drift",
            status=FAIL,
            message=(
                f"scored window has {len(w.table)} steps, shorter than "
                f"the {block_steps}-step drift block"
            ),
        )

    states = reported_states(w, probe)
    variables = list(states)
    invalid = []
    for var in variables:
        values = pd.to_numeric(w.table[var], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all():
            invalid.append(var)
    start = len(w.table) - block_steps
    if start == 0:
        initial = pd.to_numeric(w.state0[list(states)], errors="coerce").to_numpy(dtype=float)
        invalid.extend(f"initial {v}" for v, x in zip(states, initial) if not np.isfinite(x))
    precipitation = pd.to_numeric(w.forcing["pr"].iloc[start:], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(precipitation).all() or (precipitation < 0).any():
        invalid.append("forcing pr")
    if invalid:
        return CriterionResult(
            name="total_storage_drift",
            status=FAIL,
            message=f"non-finite storage-drift data: {', '.join(invalid)}",
            diagnostics={"non_finite_variables": invalid},
        )

    w.table = w.table.copy()
    w.table[variables] = w.table[variables].astype(float)
    storage_start = storage_at(w, states, start)
    storage_end = float(w.storage(states)[-1])
    delta = storage_end - storage_start
    block_precip = float(w.volume(precipitation).sum())
    allowance = relative * block_precip + absolute
    ok = abs(delta) <= allowance

    return CriterionResult(
        name="total_storage_drift",
        status=PASS if ok else FAIL,
        value=abs(delta),
        threshold=allowance,
        message=(
            f"total reported storage changes {delta:.6g} mm over final "
            f"{block_days:g}-day block (allowance {allowance:.6g} mm, "
            f"{relative:.1e} of block precipitation + {absolute:g} mm)"
        ),
        diagnostics={
            "states": list(states),
            "block_start": start,
            "block_stop": len(w.table),
            "block_days": block_days,
            "block_steps": block_steps,
            "block_precipitation_mm": block_precip,
            "storage_start_mm": storage_start,
            "storage_end_mm": storage_end,
            "storage_delta_mm": delta,
            "relative_tolerance": relative,
            "absolute_tolerance_mm": absolute,
        },
    )

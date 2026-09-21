"""Whether repeated seasonal forcing reaches one repeatable state cycle."""

from __future__ import annotations

from itertools import combinations
from typing import Any

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec, STATE_VARS


def _evaluation(run: RunResult, label: str) -> tuple[np.ndarray, Any]:
    phase = run.case.forcing.get("_phase")
    if phase is None:
        raise ValueError("spinup_cycle_invariance needs host-side '_phase' annotations")
    indices = np.flatnonzero(phase.astype(str).to_numpy() == "evaluation")
    if not len(indices):
        raise ValueError(f"variant '{label}' has no evaluation cycle")
    if not np.array_equal(indices, np.arange(indices[0], indices[0] + len(indices))):
        raise ValueError(f"variant '{label}' evaluation cycle is not contiguous")
    return indices, run.case.forcing.iloc[indices].reset_index(drop=True)


def _same_forcing(left, right, left_name: str, right_name: str) -> None:
    """Require byte-identical visible drivers for each selected pair."""
    left_visible = {
        name for name in left.columns if not name.startswith("_") and name != "time"
    }
    right_visible = {
        name for name in right.columns if not name.startswith("_") and name != "time"
    }
    if left_visible != right_visible:
        raise ValueError(
            "evaluation forcing has different visible columns: "
            f"{left_name}={sorted(left_visible)}, {right_name}={sorted(right_visible)}"
        )
    for name in sorted(left_visible):
        if not np.array_equal(left[name].to_numpy(), right[name].to_numpy()):
            raise ValueError(
                f"evaluation forcing differs for '{name}' between {left_name} and {right_name}"
            )


@criterion("spinup_cycle_invariance", paired=True)
def spinup_cycle_invariance(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict[str, Any]
) -> CriterionResult:
    """The same periodic year must not depend on extra prior cycles.

    The normal form compares all names in ``params.variants`` pairwise. The
    two-name ``short``/``long`` form remains accepted for older out-of-tree
    probes.
    """
    configured = params.get("variants")
    if configured:
        names = [str(name) for name in configured]
    else:
        names = [
            str(params.get("short", probe.control or "short")),
            str(params.get("long", "long")),
        ]
    if len(names) < 2 or len(set(names)) != len(names):
        raise ValueError("spinup_cycle_invariance needs at least two distinct variants")
    missing = [name for name in names if name not in runs]
    if missing:
        raise ValueError(f"spinup_cycle_invariance missing variants: {missing}")

    threshold = float(params.get("threshold", 0.05))
    flux_floor = float(params.get("flux_floor_mm_per_day", 0.05))
    state_floor = float(params.get("state_floor_mm", 1.0))
    required = list(params.get("variables", ["evspsbl", "mrro", *probe.requires_states]))
    optional = list(params.get("optional", ["gwex", "gw", "channel"]))
    if threshold < 0.0 or flux_floor <= 0.0 or state_floor <= 0.0:
        raise ValueError("spinup_cycle_invariance needs non-negative threshold and positive floors")

    evaluations = {name: _evaluation(runs[name], name) for name in names}
    for left_name, right_name in combinations(names, 2):
        left_run, right_run = runs[left_name], runs[right_name]
        left_rows, left_forcing = evaluations[left_name]
        right_rows, right_forcing = evaluations[right_name]
        if len(left_rows) != len(right_rows):
            raise ValueError(
                f"{left_name}- and {right_name}-spinup evaluations have different lengths "
                f"({len(left_rows)} and {len(right_rows)})"
            )
        if left_rows[0] == right_rows[0]:
            raise ValueError(f"{left_name} and {right_name} selected the same evaluation cycle")
        if left_run.case.dt_days != right_run.case.dt_days:
            raise ValueError("spinup variants must use the same timestep")
        _same_forcing(left_forcing, right_forcing, left_name, right_name)

    variables = list(required)
    for var in optional:
        present = [var in runs[name].table.columns for name in names]
        if any(present) and not all(present):
            raise ValueError(f"optional output '{var}' is present in only some variants")
        if all(present):
            variables.append(var)

    state_variables = [
        var for var in STATE_VARS
        if all(var in runs[name].table.columns for name in names)
    ]
    if state_variables:
        variables.append("total_reported_storage")

    deviations: dict[str, float] = {}
    failures: list[str] = []
    for left_name, right_name in combinations(names, 2):
        left_run, right_run = runs[left_name], runs[right_name]
        left_rows, right_rows = evaluations[left_name][0], evaluations[right_name][0]
        for var in variables:
            if var == "total_reported_storage":
                a = left_run.table[state_variables].to_numpy(dtype=float)[left_rows].sum(axis=1)
                b = right_run.table[state_variables].to_numpy(dtype=float)[right_rows].sum(axis=1)
            else:
                if var not in left_run.table.columns or var not in right_run.table.columns:
                    raise ValueError(f"spinup_cycle_invariance needs '{var}' in every result")
                a = left_run.table[var].to_numpy(dtype=float)[left_rows]
                b = right_run.table[var].to_numpy(dtype=float)[right_rows]
            key = f"{left_name}<->{right_name}:{var}"
            if not np.isfinite(a).all() or not np.isfinite(b).all():
                failures.append(f"'{var}' contains non-finite values in {left_name}/{right_name}")
                continue
            floor = state_floor if var in STATE_VARS or var == "total_reported_storage" else flux_floor
            scale = max(float(np.abs(a).mean()), floor)
            deviations[key] = float(np.abs(b - a).max() / scale)

    worst_var, worst = max(deviations.items(), key=lambda item: item[1], default=(None, 0.0))
    if worst > threshold:
        failures.append(
            f"'{worst_var}' differs by {worst:.2%} between the paired evaluations "
            f"(limit {threshold:.2%}). Outputs still differ after the prescribed spin-up; "
            "insufficient spin-up is one possible cause, so this result alone does not "
            "establish a physical violation."
        )

    return CriterionResult(
        name="spinup_cycle_invariance",
        status=FAIL if failures else PASS,
        value=worst,
        threshold=threshold,
        message=(
            "; ".join(failures)
            if failures
            else (
                "the repeated forcing reaches the same evaluation-year response for "
                f"all {len(names)} selected spin-up cycles "
                f"(worst {worst_var} departure {worst:.2%})"
            )
        ),
        diagnostics={
            "evaluation_start_rows": {
                name: int(evaluations[name][0][0]) for name in names
            },
            "variants": names,
            "evaluation_rows": int(len(evaluations[names[0]][0])),
            "reported_storage_variables": state_variables,
            "deviations": deviations,
        },
    )

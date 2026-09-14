"""Whether repeated seasonal forcing reaches one repeatable state cycle.

The input weather repeats exactly every configured cycle.  The same adapter is
run twice with identical visible inputs; host-only annotations select the year
seen after N or N+K repetitions.  The criterion compares those selected years,
so the adapter cannot key off a variant name or a different seed.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion
from hydroturing.criteria.response import pick
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


def _same_forcing(short, long) -> None:
    # The comparison must cover both schemas.  Checking only short's columns
    # would let an adapter or generator add a long-only driver unnoticed.
    short_visible = {
        name for name in short.columns if not name.startswith("_") and name != "time"
    }
    long_visible = {
        name for name in long.columns if not name.startswith("_") and name != "time"
    }
    if short_visible != long_visible:
        raise ValueError(
            "evaluation forcing has different visible columns: "
            f"short={sorted(short_visible)}, long={sorted(long_visible)}"
        )
    for name in sorted(short_visible):
        a = short[name].to_numpy()
        b = long[name].to_numpy()
        if not np.array_equal(a, b):
            raise ValueError(f"evaluation forcing differs for '{name}'")


@criterion("spinup_cycle_invariance", paired=True)
def spinup_cycle_invariance(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict[str, Any]
) -> CriterionResult:
    """The same periodic year must not depend on four extra prior cycles.

    This is a state-space test, not a claim that one arbitrary year is enough
    for all catchments.  The generator supplies N and N+K copies of the same
    365-row forcing cycle, then compares their following identical cycle.  A
    physical reference that has reached its periodic attractor produces the
    same seasonal hydrograph and reported stores in both selected years.
    """
    short_name = str(params.get("short", probe.control or "short"))
    long_name = str(params.get("long", "long"))
    threshold = float(params.get("threshold", 0.05))
    flux_floor = float(params.get("flux_floor_mm_per_day", 0.05))
    state_floor = float(params.get("state_floor_mm", 1.0))
    required = list(params.get("variables", ["evspsbl", "mrro", *probe.requires_states]))
    optional = list(params.get("optional", ["gwex", "gw", "channel"]))
    if threshold < 0.0 or flux_floor <= 0.0 or state_floor <= 0.0:
        raise ValueError(
            "spinup_cycle_invariance needs non-negative threshold and positive floors"
        )

    short_run = pick(runs, params, "short", short_name)
    long_run = pick(runs, params, "long", long_name)
    short_rows, short_forcing = _evaluation(short_run, short_name)
    long_rows, long_forcing = _evaluation(long_run, long_name)
    if len(short_rows) != len(long_rows):
        raise ValueError(
            "short- and long-spinup evaluations have different lengths "
            f"({len(short_rows)} and {len(long_rows)})"
        )
    if short_rows[0] == long_rows[0]:
        raise ValueError("short and long variants selected the same evaluation cycle")
    if short_run.case.dt_days != long_run.case.dt_days:
        raise ValueError("spinup variants must use the same timestep")
    _same_forcing(short_forcing, long_forcing)

    variables = list(required)
    for var in optional:
        in_short = var in short_run.table.columns
        in_long = var in long_run.table.columns
        if in_short != in_long:
            raise ValueError(f"optional output '{var}' is present in only one variant")
        if in_short:
            variables.append(var)

    # Score the aggregate inventory as well as individual stores.  Several
    # small stores can each sit below the 1 mm floor while their combined
    # drift is physically significant.
    state_variables = [
        var for var in STATE_VARS
        if var in short_run.table.columns and var in long_run.table.columns
    ]
    if state_variables:
        variables.append("total_reported_storage")

    deviations: dict[str, float] = {}
    failures: list[str] = []
    for var in variables:
        if var == "total_reported_storage":
            a = short_run.table[state_variables].to_numpy(dtype=float)[short_rows].sum(axis=1)
            b = long_run.table[state_variables].to_numpy(dtype=float)[long_rows].sum(axis=1)
        else:
            if var not in short_run.table.columns or var not in long_run.table.columns:
                raise ValueError(f"spinup_cycle_invariance needs '{var}' in both results")
            a = short_run.table[var].to_numpy(dtype=float)[short_rows]
            b = long_run.table[var].to_numpy(dtype=float)[long_rows]
        if not np.isfinite(a).all() or not np.isfinite(b).all():
            failures.append(f"'{var}' contains non-finite values")
            continue
        floor = (
            state_floor
            if var in STATE_VARS or var == "total_reported_storage"
            else flux_floor
        )
        scale = max(float(np.abs(a).mean()), floor)
        deviations[var] = float(np.abs(b - a).max() / scale)

    worst_var, worst = max(
        deviations.items(), key=lambda item: item[1], default=(None, 0.0)
    )
    if worst > threshold:
        failures.append(
            f"'{worst_var}' differs by {worst:.2%} between the N- and N+K-cycle evaluations "
            f"(limit {threshold:.2%})"
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
                "the repeated forcing reaches the same evaluation-year response after "
                f"N and N+K cycles "
                f"(worst {worst_var} departure {worst:.2%})"
            )
        ),
        diagnostics={
            "short_evaluation_start_row": int(short_rows[0]),
            "long_evaluation_start_row": int(long_rows[0]),
            "evaluation_rows": int(len(short_rows)),
            "reported_storage_variables": state_variables,
            "deviations": deviations,
        },
    )

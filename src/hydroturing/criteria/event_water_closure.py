"""Event water closure for each complete precipitation event, independently.

An event is a maximal run of positive supplied precipitation. At the daily
step, one dry day separates events; at a finer step, one dry step does so.
Only wet steps belong to the budget. Water still in snow, soil, groundwater
or routing at the end is accounted for by the storage change, so a recession
tail need not be appended. Dry-period closure is a separate assertion.

This criterion measures event closure, not rainfall extremity. It does not
use generator labels or assume anything about a model's training climate.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from hydroturing.criteria.base import (
    FAIL, PASS, CriterionResult, criterion, make_window, reported_states, storage_at,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@criterion("event_water_closure")
def event_water_closure(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Require |P + GWex - ET - Q - delta S| / P <= threshold for every event.

    Diagnostic start/stop indices are relative to the post-spinup window,
    with stop exclusive. Skipped events crossing spinup can have a negative
    start. Dates identify the first and last wet steps, inclusively.
    """
    unknown = set(params) - {"threshold"}
    if unknown:
        raise ValueError(f"event_water_closure: unknown parameters {sorted(unknown)}")
    threshold = float(params.get("threshold", 0.05))
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("event_water_closure threshold must be finite and nonnegative")

    forcing = run.case.forcing
    if "pr" not in forcing:
        raise ValueError("event_water_closure needs supplied forcing column 'pr'")
    rain = forcing["pr"].to_numpy(dtype=float)
    if not np.isfinite(rain).all() or (rain < 0).any():
        raise ValueError("event_water_closure needs finite, nonnegative precipitation")

    w = make_window(run, probe)
    spinup = run.case.spinup_steps
    # Detect on the full record, including spinup, to avoid manufacturing a
    # new event where an existing wet spell crosses the scoring boundary.
    edges = np.diff(np.r_[False, rain > 0, False].astype(int))
    bounds, skipped = [], []
    for a, b in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)):
        a, b = int(a), int(b)
        if b <= spinup:
            continue
        reason = None
        if a < spinup:
            reason = "crosses_spinup_boundary"
        elif a == 0:
            reason = "no_preceding_dry_step"
        elif b == len(rain):
            reason = "no_following_dry_step"
        if reason:
            skipped.append({"start": a - spinup, "stop": b - spinup, "reason": reason})
        else:
            bounds.append((a - spinup, b - spinup))

    diagnostics = {
        "event_count": len(bounds), "skipped_events": skipped,
        "denominator": "sum_pr", "dt_days": w.dt_days,
    }
    if not bounds:
        return CriterionResult(
            name="event_water_closure", status=FAIL, threshold=threshold,
            message="no complete precipitation events after spinup; nothing to score",
            diagnostics={**diagnostics, "events": []},
        )

    states = reported_states(w, probe)
    variables = ["evspsbl", "mrro", *states]
    if "gwex" in w.table:
        variables.append("gwex")
    missing = [v for v in variables if v not in w.table]
    if missing:
        raise ValueError(f"event_water_closure needs model result columns {missing}")

    # Protocol validation does not check optional gwex/gw/channel. In
    # particular, pandas storage sums would otherwise silently skip NaN.
    invalid = []
    for var in variables:
        values = pd.to_numeric(w.table[var], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all():
            invalid.append(var)
    if bounds[0][0] == 0:
        initial = pd.to_numeric(w.state0[list(states)], errors="coerce").to_numpy(dtype=float)
        invalid.extend(f"initial {v}" for v, x in zip(states, initial) if not np.isfinite(x))
    if invalid:
        return CriterionResult(
            name="event_water_closure", status=FAIL, threshold=threshold,
            message=f"non-finite model water-budget data: {', '.join(invalid)}",
            diagnostics={**diagnostics, "non_finite_variables": invalid, "events": []},
        )

    # Keep numeric coercion local to this criterion, including optional stores.
    w.table = w.table.copy()
    w.table[variables] = w.table[variables].astype(float)

    with np.errstate(over="ignore", invalid="ignore"):
        precipitation = w.volume(w.forcing["pr"])
        evap = w.volume(w.table["evspsbl"])
        runoff = w.volume(w.table["mrro"])
        exchange = w.volume(w.table["gwex"]) if "gwex" in w.table else np.zeros(len(w.table))
        storage = w.storage(states)

        events = []
        for event_id, (a, b) in enumerate(bounds, start=1):
            p = float(precipitation[a:b].sum())
            e = float(evap[a:b].sum())
            q = float(runoff[a:b].sum())
            g = float(exchange[a:b].sum())
            s0 = storage_at(w, states, a)
            s1 = float(storage[b - 1])
            change = s1 - s0
            residual = p + g - e - q - change
            # Finite individual values can still overflow during integration.
            if p <= 0 or not np.isfinite([p, e, q, g, s0, s1, change, residual]).all():
                return CriterionResult(
                    name="event_water_closure", status=FAIL, threshold=threshold,
                    message=f"event {event_id} has a non-finite or degenerate water budget",
                    diagnostics={**diagnostics, "invalid_event": {"start": a, "stop": b}},
                )
            relative = abs(residual) / p
            if not np.isfinite(relative):
                return CriterionResult(
                    name="event_water_closure", status=FAIL, threshold=threshold,
                    message=f"event {event_id} has a non-finite relative residual",
                    diagnostics={**diagnostics, "invalid_event": {"start": a, "stop": b}},
                )
            events.append({
                "event_id": event_id, "start": a, "stop": b,
                "start_time": str(w.forcing["time"].iloc[a]),
                "end_time": str(w.forcing["time"].iloc[b - 1]),
                "duration_days": (b - a) * w.dt_days,
                "precip_mm": p, "gwex_mm": g, "evap_mm": e, "runoff_mm": q,
                "storage_start_mm": s0, "storage_end_mm": s1,
                "storage_change_mm": change, "residual_mm": residual,
                "relative_residual": relative, "passed": relative <= threshold,
            })

    worst = max(events, key=lambda event: event["relative_residual"])
    n_failed = sum(not event["passed"] for event in events)
    return CriterionResult(
        name="event_water_closure", status=FAIL if n_failed else PASS,
        value=worst["relative_residual"], threshold=threshold,
        message=(
            f"worst event residual {worst['relative_residual']:.4%} of event precipitation "
            f"({n_failed}/{len(events)} events fail; limit {threshold:.1%})"
        ),
        diagnostics={
            **diagnostics, "states": list(states), "n_failed_events": n_failed,
            "events": events, "worst_event": worst,
        },
    )

"""Criteria for a groundwater control volume and river exchange."""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@criterion("groundwater_balance")
def groundwater_balance(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Check recharge + net river-aquifer exchange (positive into the
    aquifer) equals the change in groundwater storage. The exchange is a
    named flux (default gw_sw_exchange), not gwex: this moves water between
    two in-catchment stores rather than crossing the catchment boundary.

    A model with its own boundary term to the outside of this control
    volume (a GHB or WEL package, a regional groundwater exchange) may
    declare it under `sources` (default `["gwex"]`, matching closure.py);
    declared, it is added to the budget rather than left to show up as an
    unexplained residual.
    """
    w = make_window(run, probe)
    recharge_name = params.get("recharge", "gw_recharge")
    exchange_name = params.get("exchange", "gw_sw_exchange")
    storage_name = params.get("storage", "gw")
    for name in (recharge_name, exchange_name, storage_name):
        if name not in (w.forcing.columns if name == recharge_name else w.table.columns):
            raise ValueError(f"groundwater_balance needs '{name}'")

    storage = w.table[storage_name].to_numpy(dtype=float)
    recharge = w.volume(w.forcing[recharge_name])
    exchange = w.volume(w.table[exchange_name])
    sources = params.get("sources", ["gwex"])
    declared = np.zeros(len(w.table))
    for var in sources:
        if var in w.table.columns:
            declared += w.volume(w.table[var])
    storage_change = np.diff(storage, prepend=float(w.state0[storage_name]))
    residual = recharge + exchange + declared - storage_change
    total_recharge = float(np.sum(recharge))
    max_abs = float(np.max(np.abs(residual)))
    cumulative_residual = float(np.sum(residual))
    rel_tol = float(params.get("rel_tol", 0.01))
    abs_tol = float(params.get("abs_tol_mm", 1.0e-5))
    # A rounded storage state contributes a difference-of-roundings step
    # error that scales with the state's own magnitude, not with recharge or
    # exchange, and that scale is the model's choice of datum, not a probe
    # constant. Scale a second floor to the state so a model reporting `gw`
    # at a different offset is not penalized for its own rounding.
    storage_floor = 1.0e-5 * float(np.max(np.abs(storage)))
    step_floor = max(abs_tol, storage_floor)
    denominator = np.maximum(
        np.maximum(np.abs(recharge), np.abs(exchange)), 1.0e-12
    )
    allowed = np.maximum(rel_tol * denominator, step_floor)
    step_ok = bool(np.all(np.abs(residual) <= allowed))
    # Rounding noise in the storage difference cancels out over the record;
    # a systematic leak does not. The per-step floor alone cannot tell them
    # apart, so also gate on the cumulative residual against a tolerance
    # that does not inherit the per-step floor's storage-scaled slack.
    cumulative_allowed = max(rel_tol * total_recharge, 2.0 * step_floor)
    cumulative_ok = abs(cumulative_residual) <= cumulative_allowed
    ok = step_ok and cumulative_ok
    relative_total = float(np.abs(cumulative_residual) / max(total_recharge, 1.0e-12))
    return CriterionResult(
        name="groundwater_balance",
        status=PASS if ok else FAIL,
        value=relative_total,
        threshold=rel_tol,
        message=(
            f"groundwater balance closes; max step residual {max_abs:.6g} mm"
            if ok
            else (
                f"groundwater balance residual is {max_abs:.6g} mm per step"
                if not step_ok
                else f"groundwater balance drifts by {cumulative_residual:.6g} mm cumulative"
            )
        ),
        diagnostics={
            "total_recharge_mm": total_recharge,
            "cumulative_residual_mm": cumulative_residual,
            "max_abs_step_residual_mm": max_abs,
            "step_floor_mm": step_floor,
            "cumulative_allowed_mm": cumulative_allowed,
        },
    )


@criterion("exchange_directions")
def exchange_directions(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Require both groundwater-to-river and river-to-groundwater flow.

    Per the repository's convention for gw_sw_exchange (positive into the
    aquifer), gw_to_sw (aquifer losing to the river) is <= 0 and sw_to_gw
    (river losing to the aquifer) is >= 0.
    """
    w = make_window(run, probe)
    for name in ("gw_to_sw", "sw_to_gw"):
        if name not in w.table.columns:
            raise ValueError(f"exchange_directions needs '{name}'")
    forward = float(-np.minimum(w.volume(w.table["gw_to_sw"]), 0.0).sum())
    reverse = float(np.maximum(w.volume(w.table["sw_to_gw"]), 0.0).sum())
    minimum = float(params.get("minimum_gross_mm", 0.1))
    ok = forward >= minimum and reverse >= minimum
    return CriterionResult(
        name="exchange_directions",
        status=PASS if ok else FAIL,
        value=min(forward, reverse),
        threshold=minimum,
        message=(
            f"both exchange directions occur ({forward:.3f} and {reverse:.3f} mm)"
            if ok
            else f"exchange is one-sided ({forward:.3f} and {reverse:.3f} mm)"
        ),
        diagnostics={"gw_to_sw_mm": forward, "sw_to_gw_mm": reverse},
    )

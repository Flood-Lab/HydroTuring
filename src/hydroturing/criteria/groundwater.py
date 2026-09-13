"""Criteria for a groundwater control volume and river exchange."""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@criterion("groundwater_balance")
def groundwater_balance(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Check recharge + net exchange (gwex, positive into the aquifer per
    AGENTS.md) equals the change in groundwater storage."""
    w = make_window(run, probe)
    recharge_name = params.get("recharge", "gw_recharge")
    exchange_name = params.get("exchange", "gwex")
    storage_name = params.get("storage", "gw")
    for name in (recharge_name, exchange_name, storage_name):
        if name not in (w.forcing.columns if name == recharge_name else w.table.columns):
            raise ValueError(f"groundwater_balance needs '{name}'")

    recharge = w.volume(w.forcing[recharge_name])
    exchange = w.volume(w.table[exchange_name])
    storage_change = np.diff(
        w.table[storage_name].to_numpy(dtype=float),
        prepend=float(w.state0[storage_name]),
    )
    residual = recharge + exchange - storage_change
    total_recharge = float(np.sum(recharge))
    max_abs = float(np.max(np.abs(residual)))
    rel_tol = float(params.get("rel_tol", 0.01))
    abs_tol = float(params.get("abs_tol_mm", 1.0e-5))
    denominator = np.maximum(
        np.maximum(np.abs(recharge), np.abs(exchange)), 1.0e-12
    )
    allowed = np.maximum(rel_tol * denominator, abs_tol)
    ok = bool(np.all(np.abs(residual) <= allowed))
    relative_total = float(np.abs(np.sum(residual)) / max(total_recharge, 1.0e-12))
    return CriterionResult(
        name="groundwater_balance",
        status=PASS if ok else FAIL,
        value=relative_total,
        threshold=rel_tol,
        message=(
            f"groundwater balance closes; max step residual {max_abs:.6g} mm"
            if ok
            else f"groundwater balance residual is {max_abs:.6g} mm per step"
        ),
        diagnostics={
            "total_recharge_mm": total_recharge,
            "cumulative_residual_mm": float(np.sum(residual)),
            "max_abs_step_residual_mm": max_abs,
        },
    )


@criterion("exchange_directions")
def exchange_directions(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Require both groundwater-to-river and river-to-groundwater flow.

    Per the repository's gwex convention (positive into the aquifer),
    gw_to_sw (aquifer losing to the river) is <= 0 and sw_to_gw (river
    losing to the aquifer) is >= 0.
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
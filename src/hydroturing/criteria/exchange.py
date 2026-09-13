"""Consistency of signed groundwater-to-surface-water components."""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@criterion("exchange_components")
def exchange_components(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Require directional components to sum to the declared net exchange.

    The repository's convention is positive ``gwex`` into the catchment
    (AGENTS.md). ``sw_to_gw`` (river losing to the aquifer) is the positive
    component and ``gw_to_sw`` (aquifer losing to the river) is the negative
    component, so ``gw_to_sw + sw_to_gw == gwex``.
    """
    w = make_window(run, probe)
    required = ("gw_to_sw", "sw_to_gw", "gwex")
    missing = [name for name in required if name not in w.table.columns]
    if missing:
        raise ValueError(f"exchange_components needs {missing} in the model result")

    gw_to_sw = w.volume(w.table["gw_to_sw"])
    sw_to_gw = w.volume(w.table["sw_to_gw"])
    gwex = w.volume(w.table["gwex"])
    residual = gw_to_sw + sw_to_gw - gwex
    denominator = np.maximum(np.maximum(np.abs(gw_to_sw), np.abs(sw_to_gw)), 1.0e-9)
    relative = float(np.max(np.abs(residual) / denominator))
    sign_violation = bool((gw_to_sw > 1.0e-9).any() or (sw_to_gw < -1.0e-9).any())
    rel_tol = float(params.get("rel_tol", 1.0e-6))
    abs_tol = float(params.get("abs_tol", 1.0e-6))
    max_abs = float(np.max(np.abs(residual)))
    ok = not sign_violation and bool(
        np.all(np.abs(residual) <= np.maximum(rel_tol * denominator, abs_tol))
    )
    return CriterionResult(
        name="exchange_components",
        status=PASS if ok else FAIL,
        value=relative,
        threshold=rel_tol,
        message=(
            "directional exchange components sum to gwex"
            if ok
            else f"directional exchange is inconsistent with gwex (max residual {max_abs:.6g} mm)"
        ),
        diagnostics={
            "max_abs_residual_mm": max_abs,
            "max_relative_residual": relative,
            "minimum_gw_to_sw": float(np.min(gw_to_sw)),
            "maximum_sw_to_gw": float(np.max(sw_to_gw)),
        },
    )
"""Criteria for a groundwater control volume and river exchange."""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

# A directional component may sit this far on the wrong side of zero before
# its sign counts as reversed: float noise, not a share of the flow.
SIGN_TOL_MM = 1.0e-9
# Keeps a ratio finite when a configuration allows no residual at all.
_TINY = float(np.finfo(float).tiny)


@criterion("groundwater_balance")
def groundwater_balance(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Check recharge + net river-aquifer exchange (positive into the
    aquifer) equals the change in groundwater storage. The exchange is a
    named flux (default gw_sw_exchange), not gwex: this moves water between
    two in-catchment stores rather than crossing the catchment boundary.

    A model with its own boundary term acting on the aquifer alone (a GHB
    or WEL package, a regional groundwater exchange) may declare it under
    `sources` (default `[]`); declared, it is added to this control
    volume's budget rather than left to show up as an unexplained
    residual. `gwex` is not a safe default here: it is a whole-catchment
    boundary term that may be taken from any reported store (soil, channel,
    the aquifer itself), so a probe.yaml that wants a source counted here
    must name a flux scoped to the aquifer, such as `gw_boundary`, not
    `gwex` itself.
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
    sources = params.get("sources", [])
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
    # at a different offset is not penalized for its own rounding. `gw` is
    # an absolute storage, not a recharge or exchange volume, so this floor
    # is capped rather than left to scale without bound: an uncapped floor
    # lets a model raise its own tolerance by reporting `gw` at a larger
    # offset, hiding a fixed-size leak or an unexplained storage jump behind
    # a floor sized to a datum the probe never asked for.
    storage_floor = min(1.0e-5 * float(np.max(np.abs(storage))), 0.05)
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
    # Two checks with two different allowances, so no single residual can be
    # set against one threshold. Report each residual as a share of its own
    # allowance and keep the worse: the run fails exactly when this exceeds 1.
    step_ratio = float(np.max(np.abs(residual) / np.maximum(allowed, _TINY)))
    cumulative_ratio = abs(cumulative_residual) / max(cumulative_allowed, _TINY)
    return CriterionResult(
        name="groundwater_balance",
        status=PASS if ok else FAIL,
        value=max(step_ratio, cumulative_ratio),
        threshold=1.0,
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
            "step_ratio": step_ratio,
            "cumulative_ratio": cumulative_ratio,
            "cumulative_residual_share_of_recharge": relative_total,
        },
    )


@criterion("exchange_components")
def exchange_components(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Require directional components to sum to the declared net exchange.

    The net flux is named by ``net`` (default ``gw_sw_exchange``), not
    ``gwex``: river-aquifer exchange moves water between two stores inside
    the control volume, while ``gwex`` is a source or sink crossing the
    catchment boundary (AGENTS.md, closure.py). ``sw_to_gw`` (river losing to
    the aquifer) is the positive component and ``gw_to_sw`` (aquifer losing
    to the river) is the negative component, so
    ``gw_to_sw + sw_to_gw == net``.
    """
    w = make_window(run, probe)
    net_name = params.get("net", "gw_sw_exchange")
    required = ("gw_to_sw", "sw_to_gw", net_name)
    missing = [name for name in required if name not in w.table.columns]
    if missing:
        raise ValueError(f"exchange_components needs {missing} in the model result")

    gw_to_sw = w.volume(w.table["gw_to_sw"])
    sw_to_gw = w.volume(w.table["sw_to_gw"])
    net = w.volume(w.table[net_name])
    residual = gw_to_sw + sw_to_gw - net
    denominator = np.maximum(np.maximum(np.abs(gw_to_sw), np.abs(sw_to_gw)), 1.0e-9)
    relative = float(np.max(np.abs(residual) / denominator))
    # How far a component strays onto the wrong side of zero.
    sign_excess = max(float(np.max(gw_to_sw)), -float(np.min(sw_to_gw)), 0.0)
    sign_violation = sign_excess > SIGN_TOL_MM
    rel_tol = float(params.get("rel_tol", 1.0e-6))
    abs_tol = float(params.get("abs_tol", 1.0e-6))
    max_abs = float(np.max(np.abs(residual)))
    allowed = np.maximum(rel_tol * denominator, abs_tol)
    residual_ok = bool(np.all(np.abs(residual) <= allowed))
    ok = not sign_violation and residual_ok
    # Each arm as a share of its own tolerance, the worse one reported, so
    # the run fails exactly when this exceeds 1.
    residual_ratio = float(np.max(np.abs(residual) / np.maximum(allowed, _TINY)))
    sign_ratio = sign_excess / SIGN_TOL_MM
    return CriterionResult(
        name="exchange_components",
        status=PASS if ok else FAIL,
        value=max(residual_ratio, sign_ratio),
        threshold=1.0,
        message=(
            f"directional exchange components sum to {net_name}"
            if ok
            else f"directional exchange is inconsistent with {net_name} (max residual {max_abs:.6g} mm)"
        ),
        diagnostics={
            "max_abs_residual_mm": max_abs,
            "max_relative_residual": relative,
            "sign_excess_mm": sign_excess,
            "minimum_gw_to_sw": float(np.min(gw_to_sw)),
            "maximum_sw_to_gw": float(np.max(sw_to_gw)),
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

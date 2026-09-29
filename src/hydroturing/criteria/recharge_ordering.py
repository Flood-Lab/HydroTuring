"""Order preservation for a nonnegative direct aquifer-recharge pulse."""
from __future__ import annotations
import numpy as np
from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

@criterion("aquifer_recharge_ordering", paired=True)
def aquifer_recharge_ordering(runs: dict[str, RunResult], probe: ProbeSpec, params: dict) -> CriterionResult:
    """Check storage ordering, exchange ordering and cumulative partitioning.

    E_S is a maximum because one negative excursion violates comparison; E_Q
    is a sum because any reversed increment is water sent in the wrong
    direction; E_B is a maximum so local losses cannot cancel later. The
    storage floor is capped to protect absolute-storage rounding without
    letting a large datum buy unlimited tolerance.
    """
    control_name = str(params.get("control", "control")); perturbed_name = str(params.get("perturbed", "recharge_added"))
    if control_name not in runs or perturbed_name not in runs: raise ValueError("aquifer_recharge_ordering needs control and perturbed variants")
    control = make_window(runs[control_name], probe); perturbed = make_window(runs[perturbed_name], probe)
    recharge = str(params.get("recharge", "gw_recharge")); exchange = str(params.get("exchange", "gw_sw_exchange")); storage = str(params.get("storage", "gw"))
    for window, label in ((control, "control"), (perturbed, "perturbed")):
        if recharge not in window.forcing.columns: raise ValueError(f"{label} forcing needs '{recharge}'")
        for name in (exchange, storage):
            if name not in window.table.columns: raise ValueError(f"{label} result needs '{name}'")
    if len(control.table) != len(perturbed.table): raise ValueError("paired recharge runs have different scored lengths")
    added = perturbed.volume(perturbed.forcing[recharge]) - control.volume(control.forcing[recharge])
    tol = float(params.get("input_tol_mm", 1.0e-10))
    if np.any(added < -tol): raise ValueError("aquifer_recharge_ordering requires nonnegative added recharge at every step")
    cumulative_added = np.cumsum(added); added_total = float(cumulative_added[-1])
    if added_total <= tol: raise ValueError("aquifer_recharge_ordering needs positive added recharge")
    delta_storage = perturbed.table[storage].to_numpy(float) - control.table[storage].to_numpy(float)
    delta_storage -= float(perturbed.state0[storage]) - float(control.state0[storage])
    d_exchange = control.volume(control.table[exchange]) - perturbed.volume(perturbed.table[exchange]); cumulative_exchange = np.cumsum(d_exchange)
    source_delta = np.zeros(len(control.table))
    for source in params.get("sources", []):
        if source in control.table.columns and source in perturbed.table.columns:
            source_delta += perturbed.volume(perturbed.table[source]) - control.volume(control.table[source])
    storage_error = float(np.max(np.maximum(-delta_storage, 0.0))); exchange_error = float(np.sum(np.maximum(-d_exchange, 0.0)))
    budget_error = float(np.max(np.abs(delta_storage + cumulative_exchange + np.cumsum(source_delta) - cumulative_added)))
    relative = float(params.get("rel_tol", 0.01)) * added_total; absolute = float(params.get("abs_tol_mm", 1.0e-3))
    storage_max = max(float(np.max(np.abs(control.table[storage]))), float(np.max(np.abs(perturbed.table[storage])))); storage_floor = min(1.0e-5 * storage_max, 0.05)
    epsilon = max(relative, absolute, storage_floor)
    ratios = {"storage": storage_error / epsilon, "exchange": exchange_error / epsilon, "budget": budget_error / epsilon}; worst_name, worst = max(ratios.items(), key=lambda item: item[1]); ok = worst <= 1.0
    return CriterionResult(name="aquifer_recharge_ordering", status=PASS if ok else FAIL, value=worst, threshold=1.0, message=(f"added recharge preserves storage and exchange ordering (I={added_total:.6g} mm)" if ok else f"recharge ordering fails {worst_name}: {max(storage_error, exchange_error, budget_error):.6g} mm exceeds {epsilon:.6g} mm"), diagnostics={"added_recharge_mm": added_total, "epsilon_mm": epsilon, "E_S_mm": storage_error, "E_Q_mm": exchange_error, "E_B_mm": budget_error, "ratios": ratios})

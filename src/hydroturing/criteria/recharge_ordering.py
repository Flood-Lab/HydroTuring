"""Ordering for a direct aquifer-recharge pulse.

`aquifer_recharge_ordering` compares a control run with one that differs only
by a nonnegative pulse of `gw_recharge` (mass/aquifer-recharge-ordering, #135).
It is the groundwater comparison principle: for a passive aquifer whose
outflow does not fall as its storage rises, adding recharge can never leave
less water in it, and can never make the net exchange less favourable to the
river. Let I_n be the cumulative added recharge, D_n the cumulative
incremental exchange toward the river (control minus perturbed
`gw_sw_exchange`, which is positive into the aquifer), B_n the cumulative
change in any declared aquifer boundary flux (positive into the aquifer, the
sign `groundwater_balance` gives it) and dS_n the paired storage difference.
Three one-sided errors are scored:

    E_S = max_n [-dS_n]_+                    storage never falls below control
    E_Q = sum_i [-d_i]_+ dt                  the exchange increment never reverses
    E_B = max_n |dS_n + D_n - B_n - I_n|     the added water is all accounted for

E_S is a maximum because one excursion below the control breaks the order;
E_Q is a sum because every reversed increment is water sent the wrong way;
E_B is a running maximum so that a loss cannot be cancelled by a later gain.
If each run closes its own budget, which `groundwater_balance` checks on both
variants, E_B is zero to within the two tolerances: it restates that check as
a pair. The independent content of the criterion is E_S and E_Q.

All three share one tolerance scaled to the intervention rather than to the
background storage, eps = max(rel_tol * I_*, abs_tol_mm, min(1e-5 * S_max,
0.05 mm)). The last term absorbs rounding in a large absolute `gw` without
letting a datum buy unlimited slack, the same capped floor `groundwater_balance`
uses.

The order is non-strict: a zero exchange increment passes, so a model that
banks the added recharge and never releases it is not caught here. Asking
that the water drain back out needs the aquifer's response time, and for a
distributed aquifer that depends on a transmissivity and a geometry the case
does not supply.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL,
    PASS,
    CriterionResult,
    criterion,
    make_window,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@criterion("aquifer_recharge_ordering", paired=True)
def aquifer_recharge_ordering(runs: dict[str, RunResult], probe: ProbeSpec, params: dict) -> CriterionResult:
    """Added recharge never lowers storage, never reverses the exchange
    increment, and is partitioned between them without loss or creation."""
    name = "aquifer_recharge_ordering"
    control_name = str(params.get("control", "control"))
    perturbed_name = str(params.get("perturbed", "recharge_added"))
    for variant in (control_name, perturbed_name):
        if variant not in runs:
            raise ValueError(f"{name} needs variant {variant!r}")
    control = make_window(runs[control_name], probe)
    perturbed = make_window(runs[perturbed_name], probe)
    recharge = str(params.get("recharge", "gw_recharge"))
    exchange = str(params.get("exchange", "gw_sw_exchange"))
    storage = str(params.get("storage", "gw"))
    for window, label in ((control, control_name), (perturbed, perturbed_name)):
        if recharge not in window.forcing.columns:
            raise ValueError(f"{name}: {label} forcing needs '{recharge}'")
        for column in (exchange, storage):
            if column not in window.table.columns:
                raise ValueError(f"{name}: {label} result needs '{column}'")
    if len(control.table) != len(perturbed.table):
        raise ValueError(f"{name}: paired runs have different scored lengths")

    added = perturbed.volume(perturbed.forcing[recharge]) - control.volume(control.forcing[recharge])
    input_tol = float(params.get("input_tol_mm", 1.0e-10))
    # The order is only implied by physics when water is added at every step;
    # a variant that removes some is a different experiment.
    if np.any(added < -input_tol):
        raise ValueError(f"{name} requires nonnegative added recharge at every step")
    added_total = float(np.sum(added))
    if added_total <= input_tol:
        raise ValueError(f"{name} needs positive added recharge")

    s_control = control.table[storage].to_numpy(dtype=float)
    s_perturbed = perturbed.table[storage].to_numpy(dtype=float)
    s0 = float(perturbed.state0[storage]) - float(control.state0[storage])
    delta_storage = (s_perturbed - s_control) - s0
    d_exchange = control.volume(control.table[exchange]) - perturbed.volume(perturbed.table[exchange])
    d_boundary = np.zeros(len(control.table))
    for source in params.get("sources", []):
        if source in control.table.columns and source in perturbed.table.columns:
            d_boundary += perturbed.volume(perturbed.table[source]) - control.volume(control.table[source])
    # Checked here rather than left to the comparisons below: the worst arm is
    # picked with max(), which passes over a NaN ratio.
    if not (np.isfinite(delta_storage).all() and np.isfinite(d_exchange).all() and np.isfinite(d_boundary).all()):
        return CriterionResult(name=name, status=FAIL, threshold=1.0, message="non-finite storage, exchange or boundary flux")

    storage_max = float(np.max(np.abs(np.concatenate([s_control, s_perturbed]))))
    epsilon = max(
        float(params.get("rel_tol", 0.01)) * added_total,
        float(params.get("abs_tol_mm", 1.0e-3)),
        min(1.0e-5 * storage_max, 0.05),
    )
    budget = delta_storage + np.cumsum(d_exchange) - np.cumsum(d_boundary) - np.cumsum(added)
    errors = {
        "storage": float(np.max(np.maximum(-delta_storage, 0.0))),
        "exchange": float(np.sum(np.maximum(-d_exchange, 0.0))),
        "budget": float(np.max(np.abs(budget))),
    }
    ratios = {arm: error / epsilon for arm, error in errors.items()}
    worst_arm = max(ratios, key=ratios.get)
    ok = ratios[worst_arm] <= 1.0
    return CriterionResult(
        name=name,
        status=PASS if ok else FAIL,
        value=ratios[worst_arm],
        threshold=1.0,
        message=(
            f"added recharge preserves storage and exchange ordering (I={added_total:.6g} mm)"
            if ok
            else f"recharge ordering fails {worst_arm}: {errors[worst_arm]:.6g} mm exceeds {epsilon:.6g} mm"
        ),
        diagnostics={
            "added_recharge_mm": added_total,
            "epsilon_mm": epsilon,
            "E_S_mm": errors["storage"],
            "E_Q_mm": errors["exchange"],
            "E_B_mm": errors["budget"],
            "ratios": ratios,
        },
    )

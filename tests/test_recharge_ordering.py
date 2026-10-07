"""mass/aquifer-recharge-ordering: each arm caught on its own.

The gate shows that the criterion separates its references, but a must_fail
baseline names a criterion, not an arm of it, and no reference reverses the
exchange. These tests build runs on the probe's own generator and check that
each arm of `aquifer_recharge_ordering` is the one that catches the fault it is
there for; a reversed exchange increment is caught by E_Q alone. Deleting any
one arm, or flipping the boundary's sign, fails a test here.
"""

from __future__ import annotations

import runpy
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.recharge_ordering import aquifer_recharge_ordering
from hydroturing.harness import build_case, evaluate_criteria
from hydroturing.protocol import RunResult
from hydroturing.seeds import gate_seeds

REPO_ROOT = Path(__file__).resolve().parents[1]
PROBE = registry.find_probe("mass/aquifer-recharge-ordering")
SEED = gate_seeds(PROBE.id, 1)[0]


def _reference(mode: str):
    # The three reference_recharge_* adapters are one file that differs only in MODE.
    simulate = runpy.run_path(str(REPO_ROOT / "models/reference_recharge_order_exact/ht_adapter.py"))["simulate"]
    return lambda forcing, static: simulate(forcing, static, mode)


def _tables(mode: str = "exact") -> dict[str, tuple]:
    simulate = _reference(mode)
    out = {}
    for variant in PROBE.variants:
        case = build_case(PROBE, SEED, variant)
        table = pd.DataFrame(simulate(case.forcing.to_dict("records"), case.static)).drop(columns="time")
        out[variant] = (case, table)
    return out


def _runs(tables) -> dict[str, RunResult]:
    runs = {}
    for variant, (case, table) in tables.items():
        q = table["gw_sw_exchange"]
        table = table.assign(gw_to_sw=np.minimum(q, 0.0), sw_to_gw=np.maximum(q, 0.0))
        runs[variant] = RunResult(case=case, table=table, meta={}, wall_seconds=0.0)
    return runs


def _params(name: str) -> dict:
    return dict(next(c.params for c in PROBE.criteria if c.name == name))


def _ordering(runs):
    return aquifer_recharge_ordering(runs, PROBE, _params("aquifer_recharge_ordering"))


def _others(runs) -> dict[str, bool]:
    results = evaluate_criteria(runs, PROBE, control="control")
    return {r.name: r.passed for r in results if r.name != "aquifer_recharge_ordering"}


def test_exact_reference_passes_every_criterion():
    runs = _runs(_tables())
    ordering = _ordering(runs)
    assert ordering.passed
    assert ordering.diagnostics["E_B_mm"] < 1e-9
    assert all(_others(runs).values())


def test_storage_arm_catches_an_overshooting_update():
    runs = _runs(_tables("overshoot"))
    result = _ordering(runs)
    assert not result.passed
    ratios = result.diagnostics["ratios"]
    # E_S, about 20 times the tolerance. E_Q follows at about 1.5 as the river
    # refills the deficit, which is a reversed increment too.
    assert ratios["storage"] > 10.0
    assert ratios["storage"] == max(ratios.values())
    assert ratios["budget"] <= 1.0


def test_exchange_arm_catches_a_reversed_increment_the_budget_cannot_see():
    tables = _tables()
    control = tables["control"][1]
    case, perturbed = tables["recharge_added"]
    # Send the added water's exchange response the wrong way, and move gw by
    # the same amount so this run still closes its own budget.
    q = control["gw_sw_exchange"] - (perturbed["gw_sw_exchange"] - control["gw_sw_exchange"])
    gw = perturbed["gw"] + np.cumsum(q - perturbed["gw_sw_exchange"])
    tables["recharge_added"] = (case, perturbed.assign(gw_sw_exchange=q, gw=gw))
    runs = _runs(tables)
    result = _ordering(runs)
    assert not result.passed
    ratios = result.diagnostics["ratios"]
    assert ratios["exchange"] > 1.0
    assert ratios["storage"] <= 1.0 and ratios["budget"] <= 1.0
    assert all(_others(runs).values())


def test_budget_arm_catches_a_recharge_blind_model():
    runs = _runs(_tables("blind"))
    result = _ordering(runs)
    assert not result.passed
    ratios = result.diagnostics["ratios"]
    assert ratios["budget"] > 1.0
    assert ratios["storage"] <= 1.0 and ratios["exchange"] <= 1.0


def test_the_budget_arm_is_a_running_maximum():
    """A loss that a later gain cancels is still a loss: 1 mm missing from
    the perturbed storage for five days, then back."""
    tables = _tables()
    case, table = tables["recharge_added"]
    gw = table["gw"].copy()
    gw.iloc[case.spinup_steps + 5 : case.spinup_steps + 10] -= 1.0
    tables["recharge_added"] = (case, table.assign(gw=gw))
    result = _ordering(_runs(tables))
    assert not result.passed
    ratios = result.diagnostics["ratios"]
    assert ratios["budget"] > 1.0
    assert ratios["storage"] <= 1.0 and ratios["exchange"] <= 1.0


def test_a_store_that_keeps_the_pulse_is_outside_this_probe():
    """The order is non-strict, so an exchange that ignores the aquifer's head
    passes: it banks the pulse and never releases it. Pinned so that a
    recovery check, when one is added, shows up here as a change of scope."""
    tables = {}
    for variant in PROBE.variants:
        case = build_case(PROBE, SEED, variant)
        s = case.static
        b = s["aquifer_specific_yield"] * 1000.0
        k = s["river_conductance_m2_per_day"] / (s["area_km2"] * 1.0e6) * 1000.0
        q = k * (case.forcing["sw_stage_m"].to_numpy() - s["aquifer_initial_head_m"])
        gw0 = b * (s["aquifer_initial_head_m"] - s["aquifer_bottom_m"])
        gw = gw0 + np.cumsum(case.forcing["gw_recharge"].to_numpy() + q)
        tables[variant] = (case, pd.DataFrame({"gw_sw_exchange": q, "gw": gw}))
    runs = _runs(tables)
    result = _ordering(runs)
    assert result.passed
    assert result.diagnostics["E_Q_mm"] == 0.0
    assert all(_others(runs).values())


def _with_boundary(sign: float):
    """The exact cell plus a head-dependent boundary (a GHB) draining it."""
    tables = {}
    for variant in PROBE.variants:
        case = build_case(PROBE, SEED, variant)
        s = case.static
        b = s["aquifer_specific_yield"] * 1000.0
        k = s["river_conductance_m2_per_day"] / (s["area_km2"] * 1.0e6) * 1000.0
        g = 0.2 * k
        head0 = head = s["aquifer_initial_head_m"]
        rows = []
        for stage, r in zip(case.forcing["sw_stage_m"], case.forcing["gw_recharge"]):
            new = (b * head + r + k * stage + g * head0) / (b + k + g)
            q, boundary = k * (stage - new), g * (head0 - new)
            head = new
            rows.append({"gw_sw_exchange": q, "gw_boundary": sign * boundary, "gw": b * (head - s["aquifer_bottom_m"])})
        tables[variant] = (case, pd.DataFrame(rows))
    return _runs(tables)


def test_a_declared_boundary_is_credited_into_the_aquifer():
    runs = _with_boundary(+1.0)
    assert _ordering(runs).passed
    assert all(_others(runs).values())


def test_a_boundary_reported_with_the_wrong_sign_fails_the_budget():
    result = _ordering(_with_boundary(-1.0))
    assert not result.passed
    assert result.diagnostics["ratios"]["budget"] > 1.0


def test_removing_recharge_is_refused_rather_than_judged():
    tables = _tables()
    case, table = tables["recharge_added"]
    forcing = case.forcing.copy()
    forcing.loc[5, "gw_recharge"] = -1.0
    tables["recharge_added"] = (replace(case, forcing=forcing), table)
    with pytest.raises(ValueError, match="nonnegative"):
        _ordering(_runs(tables))


@pytest.mark.parametrize("column", ["gw", "gw_sw_exchange"])
def test_a_non_finite_output_fails(column):
    tables = _tables()
    case, table = tables["recharge_added"]
    table = table.copy()
    table.loc[10, column] = np.nan
    tables["recharge_added"] = (case, table)
    assert not _ordering(_runs(tables)).passed

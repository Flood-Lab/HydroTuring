"""Independent analytic values, interval units, and isolated broken controls."""

from __future__ import annotations

from dataclasses import replace
import importlib.util

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import CRITERIA
from hydroturing.harness import build_case, run_probe
from hydroturing.protocol import RunResult
from hydroturing.spec import REPO_ROOT
from test_aquifer_equilibrium import STATIC


def _adapter():
    path = REPO_ROOT / "models/reference_confined_aquifer/ht_adapter.py"
    assert path.exists(), "the confined-aquifer reference implementation is missing"
    spec = importlib.util.spec_from_file_location("confined_reference", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _forcing(stages, recharge=None):
    return [dict(time=str(i), sw_stage_m=h, gw_recharge=r)
            for i, (h, r) in enumerate(zip(stages, recharge or [0.0] * len(stages)))]


@pytest.mark.parametrize("method,expected", [
    ("exact", [3000.0, 3000.0, 3000.190325163928, 3000.362538493844]),
    ("backward_euler", [3000.0, 3000.0, 3000.181818181818, 3000.3471074380165]),
])
def test_reference_matches_hand_calculated_daily_solution(method, expected):
    rows = _adapter().simulate(_forcing([95.0, 95.0, 96.0, 96.0]), STATIC, method=method)
    assert [r["gw"] for r in rows] == pytest.approx(expected, abs=1.0e-9)


@pytest.mark.parametrize("method", ["exact", "backward_euler"])
@pytest.mark.parametrize("dt_days", [1.0, 1.0 / 24.0, 1.0 / 1440.0])
def test_interval_mean_exchange_closes_storage_with_recharge(method, dt_days):
    rows = _adapter().simulate(
        _forcing([95.0, 96.0, 94.0, 95.0], [0.0, 0.1, 0.0, 0.2]), STATIC,
        dt_days=dt_days, method=method,
    )
    previous = 3000.0
    for row, recharge in zip(rows, [0.0, 0.1, 0.0, 0.2]):
        assert (recharge + row["gw_sw_exchange"]) * dt_days == pytest.approx(row["gw"] - previous, abs=1.0e-10)
        assert row["gw_to_sw"] <= 0.0 <= row["sw_to_gw"]
        assert row["gw_to_sw"] + row["sw_to_gw"] == row["gw_sw_exchange"]
        previous = row["gw"]


def test_exact_hourly_solution_matches_same_duration_daily_solution():
    simulate = _adapter().simulate
    daily = simulate(_forcing([96.0] * 2), STATIC)
    hourly = simulate(_forcing([96.0] * 48), STATIC, dt_days=1.0 / 24.0)
    assert daily[-1]["gw"] == pytest.approx(hourly[-1]["gw"], abs=1.0e-10)
    assert sum(r["gw_sw_exchange"] for r in daily) == pytest.approx(
        sum(r["gw_sw_exchange"] / 24.0 for r in hourly), abs=1.0e-10,
    )


@pytest.mark.parametrize("kwargs", [{}, {"method": "backward_euler"}])
def test_honest_reference_never_leaves_two_mm_head_capacity(kwargs):
    simulate = _adapter().simulate
    for stage, lower, upper in [(96.0, 3000.0, 3002.0), (94.0, 2998.0, 3000.0)]:
        rows = simulate(_forcing([stage] * 90), STATIC, **kwargs)
        assert all(lower <= r["gw"] <= upper for r in rows)


@pytest.mark.parametrize("kwargs", [{"head_gain": 1.4}, {"storage_gain": 1.4}])
def test_single_factor_mutants_conserve_water_but_exceed_prescribed_capacity(kwargs):
    rows = _adapter().simulate(_forcing([96.0] * 90), STATIC, **kwargs)
    assert rows[-1]["gw"] > 3002.1
    assert sum(r["gw_sw_exchange"] for r in rows) == pytest.approx(rows[-1]["gw"] - 3000.0, abs=1.0e-9)


@pytest.mark.parametrize("kwargs", [{"head_gain": 1.001}, {"storage_gain": 1.001}])
def test_calibrated_policy_catches_conservative_per_mille_faults_that_five_percent_misses(kwargs):
    probe = registry.find_probe("mass/aquifer-equilibrium-bounds")
    simulate = _adapter().simulate
    runs = {}
    for variant in probe.variants:
        case = build_case(probe, 7, variant)
        rows = simulate(case.forcing.to_dict("records"), case.static, **kwargs)
        runs[variant] = RunResult(case, pd.DataFrame(rows), {}, 0.0)
    envelope = next(item for item in probe.criteria if item.name == "aquifer_equilibrium_bounds")
    result = CRITERIA[envelope.name](runs, probe, envelope.params)
    assert not result.passed
    assert result.value > 1.0
    original = dict(envelope.params, rel_tol=0.05, abs_tol_mm=0.001)
    assert CRITERIA[envelope.name](runs, probe, original).passed
    for item in probe.criteria:
        if item.name != envelope.name:
            assert all(CRITERIA[item.name](run, probe, item.params).passed for run in runs.values())


@pytest.mark.parametrize("kwargs", [{}, {"head_gain": 1.4}, {"storage_gain": 1.4}])
def test_recharge_ordering_does_not_establish_the_equilibrium_envelope(kwargs):
    simulate = _adapter().simulate
    probe = registry.find_probe("mass/aquifer-recharge-ordering")
    case0 = build_case(registry.find_probe("mass/aquifer-equilibrium-bounds"), 7, "raised")
    case1 = replace(case0, forcing=case0.forcing.copy(deep=True))
    case1.forcing.loc[10:19, "gw_recharge"] = 0.1
    runs = {}
    for variant, case in [("control", case0), ("recharge_added", case1)]:
        rows = simulate(case.forcing.to_dict("records"), case.static, **kwargs)
        runs[variant] = RunResult(case, pd.DataFrame(rows), {}, 0.0)
    result = CRITERIA["aquifer_recharge_ordering"](runs, probe, {})
    assert result.passed


@pytest.mark.parametrize("model,expected", [
    ("reference_confined_aquifer", []),
    ("reference_confined_aquifer_implicit", []),
    ("reference_amplified_river_head", ["aquifer_equilibrium_bounds"]),
    ("reference_inflated_aquifer_storage", ["aquifer_equilibrium_bounds"]),
])
def test_reference_contract_and_gate_failures_are_isolated(model, expected):
    assert (REPO_ROOT / "models" / model / "model.yaml").exists(), "the reference manifest is missing"
    probe = registry.find_probe("mass/aquifer-equilibrium-bounds")
    result = run_probe(registry.find_model(model), probe, [7])
    assert result.failing == expected
    assert result.verdict == ("FAIL" if expected else "PASS")


@pytest.mark.parametrize("kwargs", [{"dt_days": 0.0}, {"method": "unsupported"}])
def test_reference_rejects_invalid_integration_configuration(kwargs):
    with pytest.raises(ValueError):
        _adapter().simulate(_forcing([96.0]), STATIC, **kwargs)

"""A missing hydraulic response may remove a pass, never a measured failure."""

from __future__ import annotations

import pytest

from hydroturing import registry
from hydroturing.harness import run_probe
from hydroturing.runner import get_runner
from hydroturing.scoring import FAIL, INCOMPATIBLE, NOT_SCORED, PASS
from test_aquifer_equilibrium_validation import _audit


def _evaluate(monkeypatch, tmp_path, factors, *, bad_components=False):
    """Scale actual reference output, keeping its storage/exchange ledger."""
    from hydroturing import harness

    model = registry.find_model("reference_confined_aquifer")
    probe = registry.find_probe("mass/aquifer-equilibrium-bounds")
    real_runner = get_runner(model)

    class ResponseRunner:
        def run(self, model, probe, case, io_dir):
            run = real_runner.run(model, probe, case, io_dir)
            stage = float(case.forcing["sw_stage_m"].iloc[case.spinup_steps])
            h0 = case.static["aquifer_initial_head_m"]
            variant = "raised" if stage > h0 else "lowered" if stage < h0 else "control"
            if variant != "control":
                factor = factors[case.seed].get(variant, 1.0)
                initial = float(run.table["gw"].iloc[case.spinup_steps - 1])
                run.table["gw"] = initial + factor * (run.table["gw"] - initial)
                for column in ("gw_sw_exchange", "gw_to_sw", "sw_to_gw"):
                    run.table[column] *= factor
                if bad_components:
                    run.table.loc[len(run.table) - 1, "sw_to_gw"] += 0.1
            return run

    monkeypatch.setattr(harness, "get_runner", lambda _: ResponseRunner())
    return run_probe(model, probe, list(factors), workdir=tmp_path)


@pytest.mark.parametrize("factor", [0.0, 0.005])
def test_frozen_and_weak_responses_are_incompatible_through_the_harness(monkeypatch, tmp_path, factor):
    result = _evaluate(monkeypatch, tmp_path, {7: {"raised": factor, "lowered": factor}})
    assert result.verdict == NOT_SCORED
    assert result.reason == INCOMPATIBLE
    assert "raised" in result.incompatible[0] and "lowered" in result.incompatible[0]


def test_one_missing_arm_prevents_a_harness_pass(monkeypatch, tmp_path):
    result = _evaluate(monkeypatch, tmp_path, {7: {"raised": 0.0}})
    assert result.verdict == NOT_SCORED
    assert result.reason == INCOMPATIBLE


def test_a_bad_exercised_arm_still_fails_when_the_other_arm_is_frozen(monkeypatch, tmp_path):
    result = _evaluate(monkeypatch, tmp_path, {7: {"raised": 0.0, "lowered": 1.4}})
    assert result.verdict == FAIL
    assert result.failing == ["aquifer_equilibrium_bounds"]


def test_unexercised_response_does_not_hide_an_exchange_component_failure(monkeypatch, tmp_path):
    result = _evaluate(monkeypatch, tmp_path, {7: {"raised": 0.0, "lowered": 0.0}}, bad_components=True)
    assert result.verdict == FAIL
    assert "exchange_components" in result.failing


def test_a_good_seed_cannot_buy_a_pass_when_another_seed_is_unexercised(monkeypatch, tmp_path):
    result = _evaluate(monkeypatch, tmp_path, {7: {"raised": 0.0}, 42: {}})
    assert result.verdict == NOT_SCORED
    assert result.reason == INCOMPATIBLE
    assert "only 1 of 2 seeds" in result.incompatible[0]


def test_a_failing_seed_cannot_be_erased_by_an_unexercised_seed(monkeypatch, tmp_path):
    result = _evaluate(monkeypatch, tmp_path, {7: {"raised": 0.0}, 42: {"raised": 1.4}})
    assert result.verdict == FAIL
    assert result.failing == ["aquifer_equilibrium_bounds"]


def test_slow_aquifer_study_distinguishes_insufficient_response_from_bad_physics():
    result = _audit().response_threshold_study()
    records = result["records"]
    # Across 90 days, tau=3000 gives 2.9554%, tau=9000 gives 0.9950%.
    for method in ("exact", "backward_euler"):
        for sigma in (0.0002, 0.002, 0.01):
            fast = next(r for r in records if r["method"] == method
                        and r["storativity"] == sigma and r["timescale_days"] == 3000)
            slow = next(r for r in records if r["method"] == method
                        and r["storativity"] == sigma and r["timescale_days"] == 9000)
            assert fast["verdict"] == PASS
            assert slow["verdict"] == NOT_SCORED
            assert slow["ordinary_criteria_passed"]
            assert slow["physical_checks_passed"]
            assert all(0.0099 < f < 0.01 for f in slow["response_fractions"].values())
            if method == "exact":
                assert fast["analytic_response_fraction"] == pytest.approx(0.02955446645)
                assert slow["analytic_response_fraction"] == pytest.approx(0.00995016625)
    assert all(r["verdict"] != FAIL for r in records)

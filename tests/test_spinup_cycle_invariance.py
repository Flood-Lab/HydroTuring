"""Tests for periodic spin-up convergence and its deliberately restless control."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.harness import (
    build_case,
    evaluate_criteria,
    resolve_window_days,
    run_probe,
    select_window,
    window_case,
)
from hydroturing.protocol import FORCING_FILE, stage
from hydroturing.scoring import FAIL, PASS, VIOLATION
from hydroturing.seeds import gate_seeds

VALIDATION_SEED = 20260912
CYCLE_DAYS = 365
PERIOD_DAYS = 3652
SPINUP_DAYS = 365


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/spinup-cycle-invariance")


def _evaluation_rows(case):
    return np.flatnonzero(case.forcing["_phase"].eq("evaluation").to_numpy())


@pytest.mark.parametrize("seed", [7, VALIDATION_SEED])
def test_weather_repeats_exactly_and_host_selects_two_spinup_lengths(probe, seed):
    short = build_case(probe, seed, "short")
    long = build_case(probe, seed, "long")
    assert short.n_steps == long.n_steps == PERIOD_DAYS + SPINUP_DAYS
    assert short.spinup_steps == long.spinup_steps == SPINUP_DAYS
    assert short.static == long.static
    pd.testing.assert_frame_equal(short.forcing.drop(columns="_phase"), long.forcing.drop(columns="_phase"))

    # Every calendar year maps back to the same 365-value cycle.  Leap years
    # repeat February 28 rather than shifting all later forcing by one day.
    scored = short.forcing.iloc[SPINUP_DAYS:].copy()
    years = pd.to_datetime(scored["time"]).dt.year
    for column in ("pr", "tas", "pet"):
        yearly = [scored.loc[years == year, column].to_numpy() for year in sorted(years.unique())]
        assert all(len(values) in (365, 366) for values in yearly)
        for values in yearly:
            if len(values) == 366:
                values = np.delete(values, 59)
            np.testing.assert_array_equal(values, yearly[0][:365])
    assert short.forcing["pr"].sum() > 0.0
    assert short.forcing["tas"].min() > short.static["snow_threshold_degC"]

    for case, expected_year in ((short, 2007), (long, 2011)):
        rows = _evaluation_rows(case)
        assert len(rows) == CYCLE_DAYS
        assert pd.to_datetime(case.forcing["time"].iloc[rows[0]]).year == expected_year
        np.testing.assert_array_equal(rows, np.arange(rows[0], rows[0] + CYCLE_DAYS))


def test_adapters_receive_the_same_visible_case_metadata(probe, tmp_path):
    model = registry.find_model("reference_bucket")
    requests = []
    visible = []
    for variant in probe.variants:
        case = build_case(probe, 7, variant)
        io_dir = tmp_path / variant
        request_path = stage(io_dir, case, probe, model)
        requests.append(json.loads(request_path.read_text()))
        visible.append(pd.read_csv(io_dir / FORCING_FILE))
    assert requests[0] == requests[1]
    pd.testing.assert_frame_equal(visible[0], visible[1])
    assert "_phase" not in visible[0].columns


def test_full_record_preconditions_score_both_variants(probe):
    from hydroturing.criteria.base import make_window
    from hydroturing.protocol import RunResult

    for variant in probe.variants:
        case = build_case(probe, VALIDATION_SEED, variant)
        assert case.spinup_steps == SPINUP_DAYS
        table = pd.DataFrame({name: np.zeros(case.n_steps) for name in ("mrso", "snw", "canopy")})
        window = make_window(RunResult(case, table, {}, 0.0), probe)
        assert len(window.table) == PERIOD_DAYS
        assert window.state0.equals(table.iloc[SPINUP_DAYS - 1])


def test_all_variants_keeps_a_long_only_failure_with_no_value(probe, monkeypatch):
    """A failed variant must not be hidden by a passing numeric result."""
    from hydroturing.criteria.base import CriterionResult
    from hydroturing.protocol import RunResult
    from hydroturing import criteria as criteria_mod

    cases = {variant: build_case(probe, VALIDATION_SEED, variant) for variant in probe.variants}
    runs = {
        variant: RunResult(case, pd.DataFrame(), {}, 0.0)
        for variant, case in cases.items()
    }

    def fake_criterion(subject, _probe, _params):
        if isinstance(subject, dict):
            return CriterionResult("spinup_cycle_invariance", "pass", "ok", value=0.0)
        variant = subject.case.probe_id.rsplit("@", 1)[-1]
        if variant == "long":
            return CriterionResult("closure", "fail", "long output invalid", value=None)
        return CriterionResult("closure", "pass", "short output valid", value=0.01)

    monkeypatch.setattr(criteria_mod, "get", lambda _name: fake_criterion)
    results = evaluate_criteria(runs, probe, control="short")
    closure = next(result for result in results if result.name == "closure")
    assert not closure.passed
    assert closure.value is None
    assert closure.diagnostics["worst_variant"] == "long"


@pytest.mark.parametrize(
    "model_name", ["reference_bucket", "flex_lumped", "flex_topo", "sacsma_snow17"]
)
def test_physical_references_reach_the_same_cycle_on_an_independent_seed(probe, model_name):
    assert VALIDATION_SEED not in gate_seeds(probe.id, probe.n_seeds)
    outcome = run_probe(registry.find_model(model_name), probe, [VALIDATION_SEED])
    assert outcome.verdict == PASS, outcome.error or outcome.failing
    result = next(item for item in outcome.criteria if item.name == "spinup_cycle_invariance")
    assert result.passed
    assert result.value < result.threshold


def test_hidden_internal_clock_fails_while_its_budget_and_states_remain_valid(probe):
    outcome = run_probe(registry.find_model("reference_restless"), probe, [VALIDATION_SEED])
    assert outcome.verdict == FAIL and outcome.reason == VIOLATION, outcome.error
    assert outcome.failing == ["spinup_cycle_invariance"]
    results = {item.name: item for item in outcome.criteria}
    assert results["closure"].passed
    assert results["state_bounds"].passed
    assert results["non_degenerate"].passed
    assert results["spinup_cycle_invariance"].value > results["spinup_cycle_invariance"].threshold


def test_exact_closure_flag_survives_per_variant_aggregation(probe):
    outcome = run_probe(registry.find_model("reference_bucket"), probe, [VALIDATION_SEED])
    assert f"suspicious_exact:{probe.id}" in outcome.flags


@pytest.mark.parametrize(
    ("model_name", "criterion"),
    [
        ("reference_leaky", "closure"),
        ("reference_cheater", "state_bounds"),
        ("reference_degenerate", "non_degenerate"),
    ],
)
def test_generic_guards_still_catch_their_own_failures(probe, model_name, criterion):
    outcome = run_probe(registry.find_model(model_name), probe, [VALIDATION_SEED])
    assert outcome.verdict == FAIL and outcome.reason == VIOLATION, outcome.error
    assert criterion in outcome.failing


def test_minimum_window_keeps_both_host_selected_evaluation_cycles(probe):
    short = build_case(probe, VALIDATION_SEED, "short")
    long = build_case(probe, VALIDATION_SEED, "long")
    days = resolve_window_days(registry.find_model("reference_bucket"), probe, override=30)
    bounds = select_window(short, probe, days)
    assert bounds.days == PERIOD_DAYS
    for case in (short, long):
        kept = window_case(case, bounds)
        expected = _evaluation_rows(case)
        np.testing.assert_array_equal(_evaluation_rows(kept), expected)

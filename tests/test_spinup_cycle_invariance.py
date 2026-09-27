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
def test_weather_repeats_exactly_and_host_selects_three_spinup_lengths(probe, seed):
    cases = {variant: build_case(probe, seed, variant) for variant in probe.variants}
    assert all(case.n_steps == PERIOD_DAYS + SPINUP_DAYS for case in cases.values())
    assert all(case.spinup_steps == SPINUP_DAYS for case in cases.values())
    assert len({json.dumps(case.static, sort_keys=True) for case in cases.values()}) == 1
    # Every calendar year maps back to the same 365-value cycle.  Leap years
    # repeat February 28 rather than shifting all later forcing by one day.
    scored = cases["short"].forcing.iloc[SPINUP_DAYS:].copy()
    years = pd.to_datetime(scored["time"]).dt.year
    for column in ("pr", "tas", "pet"):
        yearly = [scored.loc[years == year, column].to_numpy() for year in sorted(years.unique())]
        assert all(len(values) in (365, 366) for values in yearly)
        for values in yearly:
            if len(values) == 366:
                values = np.delete(values, 59)
            np.testing.assert_array_equal(values, yearly[0][:365])
    assert cases["short"].forcing["pr"].sum() > 0.0
    assert cases["short"].forcing["tas"].min() > cases["short"].static["snow_threshold_degC"]

    expected_years = {"short": 2007, "plus3": 2007, "long": 2007}
    for variant, case in cases.items():
        rows = _evaluation_rows(case)
        assert len(rows) == CYCLE_DAYS
        assert pd.to_datetime(case.forcing["time"].iloc[rows[0]]).year == expected_years[variant]
        np.testing.assert_array_equal(rows, np.arange(rows[0], rows[0] + CYCLE_DAYS))


def test_adapters_receive_common_metadata_and_aligned_evaluation_forcing(probe, tmp_path):
    model = registry.find_model("reference_bucket")
    requests = []
    visible = []
    for variant in probe.variants:
        case = build_case(probe, 7, variant)
        io_dir = tmp_path / variant
        request_path = stage(io_dir, case, probe, model)
        requests.append(json.loads(request_path.read_text()))
        visible.append(pd.read_csv(io_dir / FORCING_FILE))
    assert all(request == requests[0] for request in requests)
    # Dates differ intentionally so each history reaches the same evaluation
    # year without comparing leap and non-leap calendar phases. The existing
    # probe contract keeps the row count fixed while preserving that alignment.
    assert len({frame["time"].iloc[0] for frame in visible}) == 3
    assert len({len(frame) for frame in visible}) == 1
    scored = []
    for frame in visible:
        years = pd.to_datetime(frame["time"]).dt.year
        assert set(years[years == 2007]) == {2007}
        scored.append(frame.loc[years == 2007, ["pr", "tas", "pet"]].reset_index(drop=True))
    for frame in scored[1:]:
        pd.testing.assert_frame_equal(scored[0], frame)
    assert "_phase" not in visible[0].columns


def test_full_record_preconditions_score_all_variants(probe):
    from hydroturing.criteria.base import make_window
    from hydroturing.protocol import RunResult

    for variant in probe.variants:
        case = build_case(probe, VALIDATION_SEED, variant)
        assert case.spinup_steps == SPINUP_DAYS
        table = pd.DataFrame({name: np.zeros(case.n_steps) for name in ("mrso", "snw", "canopy")})
        window = make_window(RunResult(case, table, {}, 0.0), probe)
        assert len(window.table) == PERIOD_DAYS
        assert window.state0.equals(table.iloc[SPINUP_DAYS - 1])


def test_phase_scoped_window_uses_the_row_before_the_labelled_block(probe):
    from hydroturing.criteria.base import make_window
    from hydroturing.protocol import RunResult

    case = build_case(probe, VALIDATION_SEED, "short")
    table = pd.DataFrame({"mrso": np.arange(case.n_steps, dtype=float)})
    window = make_window(RunResult(case, table, {}, 0.0), probe, phase="evaluation")
    rows = _evaluation_rows(case)
    assert len(window.table) == CYCLE_DAYS
    assert window.state0["mrso"] == rows[0] - 1


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
    ("value", "expected"),
    [("false", False), ("0", False), ("off", False), ("true", True), ("1", True)],
)
def test_harness_boolean_options_honour_quoted_yaml_values(value, expected):
    """Quoted YAML booleans must keep their written meaning."""
    from hydroturing.harness import _as_bool

    assert _as_bool(value, default=False) is expected


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
    assert (
        "The model did not reproduce the same evaluation cycle after different amounts "
        "of identical prior history. This signature does not identify the mechanism: it "
        "may reflect a hidden state or a physical store that has not yet settled."
        in results["spinup_cycle_invariance"].message
    )


def test_storage_floors_can_be_calibrated_per_variable(probe):
    from hydroturing.criteria.spinup import spinup_cycle_invariance
    from hydroturing.protocol import RunResult

    runs = {}
    for variant in probe.variants:
        case = build_case(probe, VALIDATION_SEED, variant)
        canopy = 0.05 if variant != "plus3" else 0.09
        runs[variant] = RunResult(
            case,
            pd.DataFrame({"canopy": np.full(case.n_steps, canopy)}),
            {},
            0.0,
        )

    result = spinup_cycle_invariance(
        runs,
        probe,
        {
            "variants": ["short", "plus3", "long"],
            "variables": ["canopy"],
            "optional": [],
            "state_floor_mm": 1.0,
            "state_floors_mm": {"canopy": 0.01},
        },
    )
    assert not result.passed
    assert result.diagnostics["state_floors_mm"]["canopy"] == 0.01
    assert result.diagnostics["deviations"]["short<->plus3:canopy"] > 0.05


def test_storage_floor_names_reject_typographical_errors(probe):
    from hydroturing.criteria.spinup import spinup_cycle_invariance
    from hydroturing.protocol import RunResult

    runs = {}
    for variant in probe.variants:
        case = build_case(probe, VALIDATION_SEED, variant)
        runs[variant] = RunResult(
            case,
            pd.DataFrame({"canopy": np.ones(case.n_steps)}),
            {},
            0.0,
        )

    with pytest.raises(ValueError, match="unknown storage variables"):
        spinup_cycle_invariance(
            runs,
            probe,
            {
                "variants": ["short", "plus3", "long"],
                "variables": ["canopy"],
                "optional": [],
                "state_floors_mm": {"canopyy": 0.01},
            },
        )


@pytest.mark.parametrize("period", [2, 3, 4, 6, 8])
def test_pairwise_offsets_catch_hidden_clocks_that_can_alias_one_offset(probe, period):
    """N+3 and N+4 leave no short integer clock aligned in every comparison."""
    from hydroturing.criteria.spinup import spinup_cycle_invariance
    from hydroturing.protocol import RunResult

    cycle_numbers = {"short": 5, "plus3": 8, "long": 9}
    runs = {}
    for variant, cycle in cycle_numbers.items():
        case = build_case(probe, VALIDATION_SEED, variant)
        level = 10.0 * (1.0 + 0.10 * (cycle % period))
        table = pd.DataFrame(
            {
                "mrro": np.full(case.n_steps, level),
                "evspsbl": np.full(case.n_steps, 2.0),
            }
        )
        runs[variant] = RunResult(case, table, {}, 0.0)

    result = spinup_cycle_invariance(
        runs,
        probe,
        {"variants": ["short", "plus3", "long"], "variables": ["mrro"], "optional": []},
    )
    assert not result.passed
    assert any("short<->plus3:mrro" in key for key in result.diagnostics["deviations"])
    assert any("plus3<->long:mrro" in key for key in result.diagnostics["deviations"])


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


def test_minimum_window_keeps_all_host_selected_evaluation_cycles(probe):
    cases = [build_case(probe, VALIDATION_SEED, variant) for variant in probe.variants]
    days = resolve_window_days(registry.find_model("reference_bucket"), probe, override=30)
    bounds = select_window(cases[0], probe, days)
    assert bounds.days == PERIOD_DAYS
    for case in cases:
        kept = window_case(case, bounds)
        expected = _evaluation_rows(case)
        np.testing.assert_array_equal(_evaluation_rows(kept), expected)

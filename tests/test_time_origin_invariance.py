"""Isolate absolute calendar origin from weather, season and model seed."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.closure import closure
from hydroturing.harness import (
    build_case,
    load_generator,
    run_probe,
    select_window,
    window_case,
)
from hydroturing.protocol import stage
from hydroturing.runner import get_runner
from hydroturing.scoring import FAIL, PASS, VIOLATION
from hydroturing.seeds import gate_seeds


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/time-origin-invariance")


@pytest.mark.parametrize("seed", [0, 7, 20260908])
def test_only_the_absolute_calendar_origin_changes(probe, seed):
    control = build_case(probe, seed, "control")
    shifted = build_case(probe, seed, "shifted")
    repeated = build_case(probe, seed, "control")
    assert control.n_steps == shifted.n_steps == 4015
    assert control.spinup_steps == shifted.spinup_steps == 365
    assert control.static == shifted.static == repeated.static
    assert control.forcing.to_csv(index=False) == repeated.forcing.to_csv(index=False)
    assert (
        control.forcing.drop(columns="time").to_csv(index=False)
        == shifted.forcing.drop(columns="time").to_csv(index=False)
    )

    a = pd.DatetimeIndex(control.forcing["time"])
    b = pd.DatetimeIndex(shifted.forcing["time"])
    assert (a.year - b.year == 28).all()
    assert (a != b).all()
    for component in ("month", "day", "dayofyear", "dayofweek", "is_leap_year"):
        np.testing.assert_array_equal(getattr(a, component), getattr(b, component))
    assert ((a.month == 2) & (a.day == 29)).sum() == 3
    assert (a[1:] - a[:-1] == pd.Timedelta(days=1)).all()
    assert (b[1:] - b[:-1] == pd.Timedelta(days=1)).all()
    assert control.forcing["tas"].min() >= 2.0
    assert control.forcing["pr"].iloc[365:].sum() > 0


def test_seed_changes_weather_and_static_is_not_shared(probe):
    generator = load_generator(probe)
    first, static = generator.generate(7)
    second, _ = generator.generate(8)
    assert not first["pr"].equals(second["pr"])
    static["soil_capacity_mm"] = -1
    assert generator.generate(7)[1]["soil_capacity_mm"] == 320.0
    with pytest.raises(ValueError, match="unknown variant"):
        generator.generate(7, "not-a-variant")


def test_adapter_receives_identical_seed_and_metadata(probe, tmp_path):
    model = registry.find_model("reference_bucket")
    requests = []
    for variant in probe.variants:
        io_dir = tmp_path / variant
        request_path = stage(io_dir, build_case(probe, 7, variant), probe, model)
        requests.append(json.loads(request_path.read_text()))
    assert requests[0] == requests[1]
    assert (tmp_path / "control/input/static.json").read_bytes() == (
        tmp_path / "shifted/input/static.json"
    ).read_bytes()


@pytest.mark.parametrize(
    "model_name", ["reference_bucket", "flex_lumped", "flex_topo", "sacsma_snow17"]
)
def test_physical_models_pass_a_separate_validation_seed(probe, model_name):
    seed = 20260908
    assert seed not in gate_seeds(probe.id, probe.n_seeds)
    outcome = run_probe(registry.find_model(model_name), probe, [seed])
    assert outcome.verdict == PASS, outcome.error or outcome.failing


@pytest.mark.parametrize("seed", [11, 99, 20260908])
def test_calendar_dependence_is_detected_despite_closing_both_budgets(probe, seed, tmp_path):
    outcome = run_probe(registry.find_model("reference_calendar"), probe, [seed])
    assert outcome.verdict == FAIL and outcome.reason == VIOLATION, outcome.error
    assert outcome.failing == ["invariance"]
    results = {criterion.name: criterion for criterion in outcome.criteria}
    assert results["invariance"].value > 1e-3  # Large separation from the 1e-9 limit.
    assert results["closure"].passed
    assert results["non_degenerate"].passed
    # The scored closure is the control budget. Check the shifted reference
    # separately to demonstrate why closure alone misses calendar dependence.
    model = registry.find_model("reference_calendar")
    shifted = get_runner(model).run(
        model, probe, build_case(probe, seed, "shifted"), tmp_path / "shifted"
    )
    params = next(c.params for c in probe.criteria if c.name == "closure")
    assert closure(shifted, probe, dict(params)).passed


@pytest.mark.parametrize(
    ("model_name", "expected"),
    [("reference_degenerate", "non_degenerate"), ("reference_leaky", "closure")],
)
def test_invariant_but_invalid_outputs_fail_the_guards(probe, model_name, expected):
    outcome = run_probe(registry.find_model(model_name), probe, [20260908])
    assert outcome.verdict == FAIL and outcome.reason == VIOLATION, outcome.error
    assert expected in outcome.failing
    assert next(c for c in outcome.criteria if c.name == "invariance").passed


def test_repeated_control_run_is_identical(probe, tmp_path):
    model = registry.find_model("reference_calendar")
    case = build_case(probe, 20260908, "control")
    runner = get_runner(model)
    a = runner.run(model, probe, case, tmp_path / "first")
    b = runner.run(model, probe, case, tmp_path / "repeat")
    pd.testing.assert_frame_equal(a.table, b.table, check_exact=True)


def test_event_window_preserves_pairing_and_annual_guard(probe):
    control = build_case(probe, 20260908, "control")
    shifted = build_case(probe, 20260908, "shifted")
    bounds = select_window(control, probe, 365)
    a = window_case(control, bounds)
    b = window_case(shifted, bounds)
    pd.testing.assert_frame_equal(
        a.forcing.drop(columns="time"), b.forcing.drop(columns="time"), check_exact=True
    )
    assert a.spinup_steps == b.spinup_steps
    for name, expected in (("reference_bucket", PASS), ("reference_calendar", FAIL)):
        outcome = run_probe(registry.find_model(name), probe, [20260908], window=30)
        assert outcome.window_days == 365
        assert outcome.verdict == expected, outcome.error or outcome.failing
        if name == "reference_calendar":
            assert outcome.failing == ["invariance"]
        else:
            guard = next(c for c in outcome.criteria if c.name == "non_degenerate")
            assert "runoff_ratio_check" not in guard.diagnostics

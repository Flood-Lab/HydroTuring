"""Hand-calculated event budgets, independently of the synthetic generator."""

from __future__ import annotations

import json
from dataclasses import asdict
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from hydroturing.criteria import get, is_paired
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.protocol import Case, RunResult
from hydroturing.spec import TIMESTEP_DAYS


@pytest.fixture
def probe():
    return SimpleNamespace(requires_states=("mrso", "snw", "canopy"))


def build(pr, *, spinup=1, timestep="PT1D", **outputs):
    """Supply explicit fluxes and absolute endpoint stores; infer no budget terms."""
    rain = np.asarray(pr, dtype=float)
    time = pd.date_range(
        "2001-01-01", periods=len(rain),
        freq=pd.Timedelta(days=TIMESTEP_DAYS[timestep]),
    )
    forcing = pd.DataFrame({"time": time, "pr": rain})
    table = pd.DataFrame({
        "time": time, "pr": rain.copy(), "evspsbl": 0.0, "mrro": 0.0,
        "mrso": 0.0, "snw": 0.0, "canopy": 0.0,
    })
    for variable, values in outputs.items():
        table[variable] = np.asarray(values, dtype=float)
    case = Case("test/event-water-closure", 1, forcing, {}, spinup, timestep)
    return RunResult(case, table, {}, 0.0)


def score(run, probe, **params):
    return get("event_water_closure")(run, probe, params)


def test_hand_calculated_event_at_the_first_scored_row(probe):
    # 30 mm rain - 3 mm ET - 7 mm runoff = 20 mm added to storage.
    run = build(
        [0, 10, 20, 0], evspsbl=[0, 1, 2, 0], mrro=[0, 3, 4, 0],
        mrso=[10, 16, 30, 30],
    )
    result = score(run, probe)
    assert result.status == PASS
    assert result.name == "event_water_closure"
    assert not is_paired(result.name)
    assert result.value == pytest.approx(0.0)
    assert result.threshold == pytest.approx(0.05)
    assert result.diagnostics["event_count"] == 1
    assert result.diagnostics["n_failed_events"] == 0
    event = result.diagnostics["events"][0]
    assert (event["start"], event["stop"]) == (0, 2)
    for name, expected in {
        "duration_days": 2, "precip_mm": 30, "gwex_mm": 0,
        "evap_mm": 3, "runoff_mm": 7, "storage_start_mm": 10,
        "storage_end_mm": 30, "storage_change_mm": 20,
        "residual_mm": 0, "relative_residual": 0,
    }.items():
        assert event[name] == pytest.approx(expected), name
    assert event["passed"]
    assert {"event_id", "start_time", "end_time"} <= event.keys()
    assert pd.Timestamp(event["start_time"]) == run.case.forcing.loc[1, "time"]
    assert isinstance(result.diagnostics["skipped_events"], list)
    json.dumps(asdict(result), allow_nan=False)


def test_opposite_event_errors_pass_cumulative_closure_but_fail_events(probe):
    run = build([0, 10, 0, 10, 0], mrro=[0, 8, 0, 12, 0])
    assert get("closure")(run, probe, {}).status == PASS
    result = score(run, probe)
    assert result.status == FAIL
    assert result.value == pytest.approx(0.2)
    assert result.diagnostics["n_failed_events"] == 2
    assert [e["residual_mm"] for e in result.diagnostics["events"]] == [2, -2]


def test_one_large_event_cannot_dilute_a_small_failing_event(probe):
    run = build([0, 100, 0, 1, 0], mrro=[0, 100, 0, 0.9, 0])
    assert get("closure")(run, probe, {}).status == PASS
    result = score(run, probe)
    assert result.status == FAIL
    assert result.value == pytest.approx(0.1)
    assert result.diagnostics["n_failed_events"] == 1
    assert result.diagnostics["worst_event"]["precip_mm"] == pytest.approx(1)


def test_within_event_signed_errors_are_integrated_as_one_budget(probe):
    # This measures each event's net loss, not the sum of absolute step errors.
    run = build([0, 10, 10, 0], mrro=[0, 8, 12, 0])
    result = score(run, probe)
    assert result.status == PASS
    assert result.diagnostics["event_count"] == 1
    assert result.value == pytest.approx(0.0)


def test_events_crossing_spinup_or_record_edges_are_not_partly_scored(probe):
    # Only the 10 mm event is complete within the scored window. The first
    # event lacks a preceding row, the 12 mm event crosses spinup, and the
    # final 9 mm event has no following dry row.
    run = build(
        [8, 0, 6, 6, 0, 10, 0, 9], spinup=3,
        mrro=[0, 0, 0, 0, 0, 10, 0, 0],
    )
    result = score(run, probe)
    assert result.status == PASS
    assert result.diagnostics["event_count"] == 1
    event = result.diagnostics["events"][0]
    assert (event["start"], event["stop"]) == (2, 3)
    assert event["precip_mm"] == pytest.approx(10)
    assert len(result.diagnostics["skipped_events"]) >= 2


def test_all_reported_stores_are_counted_without_double_counting(probe):
    # Soil +2, snow -1, canopy +1, groundwater +2, routing +1 = +5 mm.
    run = build(
        [0, 10, 0], evspsbl=[0, 2, 0], mrro=[0, 3, 0],
        mrso=[10, 12, 12], snw=[2, 1, 1], canopy=[1, 2, 2],
        gw=[5, 7, 7], channel=[4, 5, 5],
    )
    result = score(run, probe)
    assert result.status == PASS
    event = result.diagnostics["events"][0]
    assert event["storage_start_mm"] == pytest.approx(22)
    assert event["storage_end_mm"] == pytest.approx(27)
    assert event["storage_change_mm"] == pytest.approx(5)


def test_later_event_uses_its_immediate_preceding_storage(probe):
    # The dry gap drains 4 mm, leaving 21 mm before the second event.
    run = build(
        [0, 8, 0, 6, 0], mrro=[0, 3, 4, 3, 0], mrso=[20, 25, 21, 24, 24],
    )
    result = score(run, probe)
    assert result.status == PASS
    first, second = result.diagnostics["events"]
    assert (first["storage_start_mm"], second["storage_start_mm"]) == (20, 21)
    assert second["storage_change_mm"] == pytest.approx(3)


@pytest.mark.parametrize("exchange,runoff", [(2, 12), (-2, 8)])
def test_groundwater_exchange_sign_is_positive_into_the_catchment(probe, exchange, runoff):
    run = build([0, 10, 0], gwex=[0, exchange, 0], mrro=[0, runoff, 0])
    result = score(run, probe)
    assert result.status == PASS
    assert result.diagnostics["events"][0]["gwex_mm"] == pytest.approx(exchange)


def test_declared_source_does_not_enlarge_the_rain_denominator(probe):
    # 1 mm rain + 9 mm exchange - 9.9 mm runoff leaves 10% of rain unaccounted.
    run = build([0, 1, 0], gwex=[0, 9, 0], mrro=[0, 9.9, 0])
    result = score(run, probe)
    assert result.status == FAIL
    assert result.value == pytest.approx(0.1)


def test_events_and_budget_use_given_rain_instead_of_the_models_echo(probe):
    run = build([0, 10, 0], mrro=[0, 10, 0])
    run.table["pr"] = 0.0
    result = score(run, probe)
    assert result.status == PASS
    assert result.diagnostics["events"][0]["precip_mm"] == pytest.approx(10)
    assert get("forcing_fidelity")(run, probe, {}).status == FAIL


@pytest.mark.parametrize("timestep,rate,depth,duration", [
    ("PT1D", 24, 48, 2),
    ("PT1H", 24, 2, 2 / 24),
    ("PT15M", 96, 2, 2 / 96),
])
def test_flux_rates_use_the_case_step_and_one_dry_step_separates_events(
    probe, timestep, rate, depth, duration,
):
    run = build(
        [0, rate, rate, 0, rate, 0], timestep=timestep,
        mrro=[0, rate, rate, 0, rate, 0],
    )
    result = score(run, probe)
    assert result.status == PASS
    assert result.diagnostics["event_count"] == 2
    first, second = result.diagnostics["events"]
    assert first["precip_mm"] == pytest.approx(depth)
    assert first["runoff_mm"] == pytest.approx(depth)
    assert first["duration_days"] == pytest.approx(duration)
    assert second["precip_mm"] == pytest.approx(depth / 2)


def test_sublimation_component_and_discharge_are_not_extra_water_losses(probe):
    run = build(
        [0, 10, 0], evspsbl=[0, 3, 0], mrro=[0, 7, 0],
        sbl=[0, 2, 0], dis=[0, 100, 0],
    )
    assert score(run, probe).status == PASS


def test_dry_day_fluxes_are_outside_the_event_budget(probe):
    run = build(
        [0, 10, 0], evspsbl=[0, 0, 100], mrro=[0, 10, 100],
    )
    assert get("closure")(run, probe, {}).status == FAIL
    assert score(run, probe).status == PASS


@pytest.mark.parametrize("runoff,expected", [(19.01, PASS), (19, PASS), (18.99, FAIL)])
def test_five_percent_boundary_is_inclusive(probe, runoff, expected):
    result = score(build([0, 20, 0], mrro=[0, runoff, 0]), probe)
    assert result.status == expected
    assert result.value == pytest.approx(abs(20 - runoff) / 20)
    assert result.threshold == pytest.approx(0.05)


def test_threshold_override_changes_the_same_budget_decision(probe):
    run = build([0, 20, 0], mrro=[0, 18, 0])
    assert score(run, probe).status == FAIL
    relaxed = score(run, probe, threshold=0.1)
    assert relaxed.status == PASS
    assert relaxed.value == relaxed.threshold == pytest.approx(0.1)


@pytest.mark.parametrize("runoff,expected", [(1e-14, PASS), (0, FAIL)])
def test_tiny_positive_rain_keeps_its_actual_denominator(probe, runoff, expected):
    result = score(build([0, 1e-14, 0], mrro=[0, runoff, 0]), probe)
    assert result.status == expected
    assert result.diagnostics["events"][0]["precip_mm"] == 1e-14
    assert result.value == (0 if expected == PASS else 1)


@pytest.mark.parametrize("rain", [
    [0, 0, 0], [1, 1, 1], [0, 0, 1], [1, 0, 0],
], ids=["no-rain", "whole-record-wet", "unfinished-final-event", "missing-prior-state"])
def test_no_complete_event_is_a_finite_scored_failure(probe, rain):
    result = score(build(rain, spinup=0), probe)
    assert result.status == FAIL
    assert result.value is None
    assert result.diagnostics["event_count"] == 0
    json.dumps(asdict(result), allow_nan=False)


@pytest.mark.parametrize("variable,invalid", [
    ("evspsbl", np.nan), ("mrro", np.inf), ("mrso", -np.inf),
    ("snw", np.nan), ("canopy", np.inf), ("gw", np.nan),
    ("channel", np.inf), ("gwex", -np.inf),
])
def test_nonfinite_required_or_optional_budget_data_is_a_scored_failure(
    probe, variable, invalid,
):
    run = build([0, 10, 0], mrro=[0, 10, 0])
    if variable not in run.table:
        run.table[variable] = 0.0
    run.table.loc[1, variable] = invalid
    result = score(run, probe)
    assert result.status == FAIL
    assert result.value is None
    json.dumps(asdict(result), allow_nan=False)


@pytest.mark.parametrize("variable", ["mrso", "gw", "channel"])
def test_nonfinite_preceding_storage_cannot_be_silently_ignored(probe, variable):
    run = build([0, 10, 0], mrro=[0, 10, 0])
    if variable not in run.table:
        run.table[variable] = 0.0
    run.table.loc[0, variable] = np.nan
    result = score(run, probe)
    assert result.status == FAIL
    assert result.value is None
    json.dumps(asdict(result), allow_nan=False)


@pytest.mark.parametrize("rain", [-1, np.nan, np.inf, -np.inf])
def test_invalid_forcing_rain_is_an_input_error(probe, rain):
    with pytest.raises(ValueError):
        score(build([0, rain, 0]), probe)


@pytest.mark.parametrize("threshold", [-0.01, np.nan, np.inf, -np.inf])
def test_invalid_threshold_is_a_configuration_error(probe, threshold):
    with pytest.raises(ValueError):
        score(build([0, 10, 0], mrro=[0, 10, 0]), probe, threshold=threshold)


@pytest.mark.parametrize("variable", ["evspsbl", "mrro", "mrso", "snw", "canopy"])
def test_missing_required_budget_columns_raise(probe, variable):
    run = build([0, 10, 0], mrro=[0, 10, 0])
    run.table = run.table.drop(columns=variable)
    with pytest.raises(ValueError):
        score(run, probe)

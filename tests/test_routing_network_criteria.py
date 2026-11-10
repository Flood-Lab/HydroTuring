"""Tests for reach-indexed routing criteria."""

import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.criteria.routing import routing_input_fidelity
from hydroturing.protocol import Case, RunResult


@pytest.fixture
def input_fidelity_run():
    times = pd.date_range("2000-01-01", periods=3, freq="D")
    forcing = pd.DataFrame(
        {
            "time": times,
            "q_in_A": [0.0, 2.0, 4.0],
            "q_in_B": [0.0, 3.0, 5.0],
        }
    )
    case = Case(
        probe_id="mass/routing-network-closure",
        seed=1,
        forcing=forcing,
        static={"routing_network": {"reaches": ["A", "B"]}},
        spinup_steps=1,
    )

    rows = []
    for i, time in enumerate(times):
        for reach, forcing_name in (("A", "q_in_A"), ("B", "q_in_B")):
            value = forcing.loc[i, forcing_name]
            rows.append(
                {
                    "time": time,
                    "reach_id": reach,
                    "reach_q_in": value,
                    "q_out": value,
                    "channel_storage": 0.0,
                }
            )

    return RunResult(
        case=case,
        table=pd.DataFrame({"time": times}),
        meta={},
        wall_seconds=0.0,
        routing=pd.DataFrame(rows),
    )


@pytest.mark.parametrize("changed_reach", [None, "A", "B"])
def test_input_fidelity_checks_both_headwaters(input_fidelity_run, changed_reach):
    run = input_fidelity_run
    if changed_reach is not None:
        mask = run.routing["reach_id"] == changed_reach
        run.routing.loc[mask, "reach_q_in"] *= 1.1

    probe = registry.find_probe("mass/routing-network-closure")
    result = routing_input_fidelity(
        run,
        probe,
        {"headwater_forcing": {"A": "q_in_A", "B": "q_in_B"}},
    )

    assert result.status == (PASS if changed_reach is None else FAIL)


@pytest.fixture
def closure_run():
    times = pd.date_range("2000-01-01", periods=3, freq="D")
    case = Case(
        probe_id="mass/routing-network-closure",
        seed=1,
        forcing=pd.DataFrame({"time": times}),
        static={
            "routing_network": {
                "reaches": ["A", "B", "C"],
                "contributing_area_km2": {"A": 1.0, "B": 1.0, "C": 2.0},
                "junctions": [
                    {"id": "J", "upstream": ["A", "B"], "downstream": "C"}
                ],
            }
        },
        spinup_steps=1,
    )

    rows = []
    for time in times:
        for reach, value in (("A", 2.0), ("B", 3.0), ("C", 5.0)):
            rows.append(
                {
                    "time": time,
                    "reach_id": reach,
                    "reach_q_in": value,
                    "q_out": value,
                    "channel_storage": 0.0,
                }
            )

    return RunResult(
        case=case,
        table=pd.DataFrame({"time": times}),
        meta={},
        wall_seconds=0.0,
        routing=pd.DataFrame(rows),
    )


@pytest.mark.parametrize("with_errors", [False, True])
def test_reach_closure_does_not_cancel_errors_over_time(closure_run, with_errors):
    from hydroturing.criteria.routing import routing_reach_closure

    run = closure_run
    if with_errors:
        a_rows = run.routing.index[run.routing["reach_id"] == "A"]
        run.routing.loc[a_rows[1], "q_out"] = 1.0
        run.routing.loc[a_rows[2], "q_out"] = 3.0

    probe = registry.find_probe("mass/routing-network-closure")
    result = routing_reach_closure(run, probe, {})

    assert result.status == (FAIL if with_errors else PASS)


@pytest.mark.parametrize("with_errors", [False, True])
def test_junction_closure_does_not_cancel_errors_over_time(closure_run, with_errors):
    from hydroturing.criteria.routing import routing_junction_closure

    run = closure_run
    if with_errors:
        c_rows = run.routing.index[run.routing["reach_id"] == "C"]
        run.routing.loc[c_rows[1], "reach_q_in"] = 4.0
        run.routing.loc[c_rows[2], "reach_q_in"] = 6.0

    probe = registry.find_probe("mass/routing-network-closure")
    result = routing_junction_closure(run, probe, {})

    assert result.status == (FAIL if with_errors else PASS)


@pytest.mark.parametrize(
    ("storage_m3", "expected"),
    [(49.0, PASS), (51.0, FAIL)],
)
def test_recession_bound_checks_absolute_allowance(closure_run, storage_m3, expected):
    from hydroturing.criteria.routing import routing_recession_bound

    run = closure_run
    run.routing["reach_q_in"] = 0.0
    a_rows = run.routing.index[run.routing["reach_id"] == "A"]
    run.routing.loc[a_rows[2], "channel_storage"] = storage_m3

    probe = registry.find_probe("mass/routing-network-closure")
    result = routing_recession_bound(
        run,
        probe,
        {"min_allowance_mm": 0.05},
    )

    assert result.status == expected


@pytest.mark.parametrize(
    "negative_column",
    [None, "reach_q_in", "q_out", "channel_storage"],
)
def test_routing_state_bounds_rejects_negative_values(closure_run, negative_column):
    from hydroturing.criteria.routing import routing_state_bounds

    run = closure_run
    if negative_column is not None:
        a_rows = run.routing.index[run.routing["reach_id"] == "A"]
        run.routing.loc[a_rows[1], negative_column] = -1.0

    probe = registry.find_probe("mass/routing-network-closure")
    result = routing_state_bounds(run, probe, {})

    assert result.status == (PASS if negative_column is None else FAIL)


@pytest.mark.parametrize("criterion_name", ["reach", "junction"])
def test_zero_reported_flow_does_not_pass(closure_run, criterion_name):
    from hydroturing.criteria.routing import (
        routing_junction_closure,
        routing_reach_closure,
    )

    run = closure_run
    run.routing["reach_q_in"] = 0.0
    run.routing["q_out"] = 0.0

    probe = registry.find_probe("mass/routing-network-closure")
    criterion = (
        routing_reach_closure
        if criterion_name == "reach"
        else routing_junction_closure
    )
    result = criterion(run, probe, {})

    assert result.status == FAIL
    assert "degenerate" in result.message
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


@pytest.mark.parametrize(
    ("lateral_fraction", "expected"),
    [(0.01, PASS), (0.10, FAIL)],
)
def test_junction_closure_counts_lateral_inflow_as_residual(
    closure_run, lateral_fraction, expected
):
    from hydroturing.criteria.routing import routing_junction_closure

    run = closure_run
    c_rows = run.routing["reach_id"] == "C"
    run.routing.loc[c_rows, "reach_q_in"] *= 1.0 + lateral_fraction
    run.routing.loc[c_rows, "q_out"] *= 1.0 + lateral_fraction

    probe = registry.find_probe("mass/routing-network-closure")
    result = routing_junction_closure(run, probe, {"threshold": 0.05})

    assert result.status == expected


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


@pytest.mark.parametrize("extra_flow", [0.0, 0.001])
def test_recession_bound_cannot_be_evaded_by_extra_tail_flow(extra_flow):
    from hydroturing.criteria import get

    times = pd.date_range("2000-01-01", periods=62, freq="D")
    forcing = pd.DataFrame({
        "time": times,
        "q_in_A": [1.0 if i == 1 else 0.0 for i in range(len(times))],
        "q_in_B": [1.0 if i == 1 else 0.0 for i in range(len(times))],
    })

    case = Case(
        probe_id="mass/routing-network-closure",
        seed=1,
        forcing=forcing,
        static={
            "routing_network": {
                "reaches": ["A", "B", "C"],
                "contributing_area_km2": {
                    "A": 1.0,
                    "B": 1.0,
                    "C": 2.0,
                },
                "junctions": [
                    {"id": "J", "upstream": ["A", "B"], "downstream": "C"}
                ],
            }
        },
        spinup_steps=1,
    )

    rows = []
    for i, time in enumerate(times):
        for reach in ("A", "B", "C"):
            flow = (2.0 if reach == "C" else 1.0) if i == 1 else 0.0
            storage = 0.0
            if i == len(times) - 1 and reach == "C":
                flow += extra_flow
                storage = 200.0

            rows.append({
                "time": time,
                "reach_id": reach,
                "reach_q_in": flow,
                "q_out": flow,
                "channel_storage": storage,
            })

    run = RunResult(
        case=case,
        table=pd.DataFrame({"time": times}),
        meta={},
        wall_seconds=0.0,
        routing=pd.DataFrame(rows),
    )

    probe = registry.find_probe("mass/routing-network-closure")
    results = {
        criterion.name: get(criterion.name)(run, probe, criterion.params)
        for criterion in probe.criteria
    }

    assert results["routing_input_fidelity"].status == PASS
    assert results["routing_reach_closure"].status == PASS
    assert results["routing_junction_closure"].status == PASS
    assert results["routing_state_bounds"].status == PASS
    assert results["routing_recession_bound"].status == FAIL


def test_recession_bound_allows_conserved_delayed_junction_inflow():
    from hydroturing.criteria import get

    times = pd.date_range("2000-01-01", periods=62, freq="D")
    forcing = pd.DataFrame({
        "time": times,
        "q_in_A": [1.0 if i == 1 else 0.0 for i in range(len(times))],
        "q_in_B": [1.0 if i == 1 else 0.0 for i in range(len(times))],
    })

    case = Case(
        probe_id="mass/routing-network-closure",
        seed=1,
        forcing=forcing,
        static={
            "routing_network": {
                "reaches": ["A", "B", "C"],
                "contributing_area_km2": {
                    "A": 1.0,
                    "B": 1.0,
                    "C": 2.0,
                },
                "junctions": [
                    {"id": "J", "upstream": ["A", "B"], "downstream": "C"}
                ],
            }
        },
        spinup_steps=1,
    )

    seconds_per_day = 86400.0
    release = 50.0 / seconds_per_day
    rows = []

    for i, time in enumerate(times):
        if i == 0:
            values = {
                "A": (0.0, 0.0, 0.0),
                "B": (0.0, 0.0, 0.0),
                "C": (0.0, 0.0, 0.0),
            }
        elif i == 1:
            values = {
                "A": (1.0, 1.0 - release, 50.0),
                "B": (1.0, 1.0 - release, 50.0),
                "C": (2.0 - 2 * release, 2.0 - 4 * release, 100.0),
            }
        elif i == len(times) - 1:
            values = {
                "A": (0.0, release, 0.0),
                "B": (0.0, release, 0.0),
                "C": (2 * release, 0.0, 200.0),
            }
        else:
            values = {
                "A": (0.0, 0.0, 50.0),
                "B": (0.0, 0.0, 50.0),
                "C": (0.0, 0.0, 100.0),
            }

        for reach, (inflow, outflow, storage) in values.items():
            rows.append({
                "time": time,
                "reach_id": reach,
                "reach_q_in": inflow,
                "q_out": outflow,
                "channel_storage": storage,
            })

    run = RunResult(
        case=case,
        table=pd.DataFrame({"time": times}),
        meta={},
        wall_seconds=0.0,
        routing=pd.DataFrame(rows),
    )

    probe = registry.find_probe("mass/routing-network-closure")
    for criterion in probe.criteria:
        result = get(criterion.name)(run, probe, criterion.params)
        assert result.status == PASS, (
            f"{criterion.name}: {result.status} — {result.message}"
        )
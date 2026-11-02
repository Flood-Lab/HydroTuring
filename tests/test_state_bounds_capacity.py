"""A capacity bounds a store only for a model that reads it, unless the store runs away."""

from __future__ import annotations

from dataclasses import replace
import runpy

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.base import PASS, get
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult
from hydroturing.seeds import gate_seeds
from hydroturing.spec import REPO_ROOT


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/catchment-closure")


@pytest.fixture(scope="module")
def case(probe):
    return build_case(probe, gate_seeds(probe.id, 1)[0])


def simulated(name, case):
    simulate = runpy.run_path(str(REPO_ROOT / f"models/{name}/ht_adapter.py"))["simulate"]
    return pd.DataFrame(simulate(case.forcing.to_dict("records"), case.static))


def state_bounds(run, probe):
    params = next(c.params for c in probe.criteria if c.name == "state_bounds")
    return get("state_bounds")(run, probe, params)


def unread():
    """The exact bucket's manifest with the capacities taken out of what it consumes."""
    return replace(registry.find_model("reference_bucket"), needs_static=(), uses_static=())


def lifted(case):
    """A stationary soil store sitting wholly above the catchment's capacity."""
    table = simulated("reference_bucket", case)
    table["mrso"] += case.static["soil_capacity_mm"]
    return table


def test_a_stationary_store_past_an_unread_capacity_is_not_held_to_it(probe, case):
    result = state_bounds(RunResult(case, lifted(case), {}, 0, model=unread()), probe)
    assert result.status == PASS
    assert "mrso upper bound N/A" in result.message
    assert "soil_capacity_mm" in result.message


def test_a_model_that_reads_the_capacity_is_held_to_it(probe, case):
    model = registry.find_model("reference_bucket")
    assert "soil_capacity_mm" in model.needs_static
    assert not state_bounds(RunResult(case, lifted(case), {}, 0, model=model), probe).passed


def test_a_run_with_no_manifest_is_held_to_the_bound(probe, case):
    assert not state_bounds(RunResult(case, lifted(case), {}, 0), probe).passed


def test_storage_solved_from_the_budget_cannot_escape_by_not_declaring(probe, case):
    model = registry.find_model("reference_cheater")
    assert "soil_capacity_mm" not in model.needs_static + model.uses_static
    run = RunResult(case, simulated("reference_cheater", case), {}, 0, model=model)
    result = state_bounds(run, probe)
    assert not result.passed
    assert "runs away" in result.message


def test_the_lower_bound_still_applies(probe, case):
    table = simulated("reference_bucket", case)
    table["mrso"] -= 2.0 * case.static["soil_capacity_mm"]
    assert not state_bounds(RunResult(case, table, {}, 0, model=unread()), probe).passed


def test_a_store_inside_the_capacity_says_nothing_about_it(probe, case):
    run = RunResult(case, simulated("reference_bucket", case), {}, 0, model=unread())
    result = state_bounds(run, probe)
    assert result.passed
    assert result.message == "all storages stay physical"


def test_a_record_too_short_to_show_the_store_is_stationary_keeps_the_bound(probe, case):
    """Or invented storage escapes on the default flood window."""
    rows = case.spinup_steps + 100
    short = replace(case, forcing=case.forcing.iloc[:rows].reset_index(drop=True))
    table = simulated("reference_bucket", short)
    table["mrso"] = case.static["soil_capacity_mm"] + np.ones(rows)
    result = state_bounds(RunResult(short, table, {}, 0, model=unread()), probe)
    assert not result.passed
    assert "cannot show it does not run away" in result.message


def over_capacity(case, rows, store):
    """The exact bucket on `rows` rows, with its soil store replaced after the spinup."""
    short = replace(case, forcing=case.forcing.iloc[:rows].reset_index(drop=True))
    table = simulated("reference_bucket", short)
    scored = rows - case.spinup_steps
    table["mrso"] = case.static["soil_capacity_mm"] + np.r_[np.ones(case.spinup_steps), store(scored)]
    return RunResult(short, table, {}, 0, model=unread())


def test_a_store_that_accelerates_cannot_hide_behind_the_height_it_reached(probe, case):
    run = over_capacity(case, len(case.forcing), lambda n: 400.0 * 1.005 ** np.arange(n))
    result = state_bounds(run, probe)
    assert not result.passed
    assert "runs away" in result.message


def test_a_store_that_turns_infinite_runs_away(probe, case):
    def store(n):
        values = np.ones(n)
        values[n // 2:] = np.inf
        return values

    result = state_bounds(over_capacity(case, len(case.forcing), store), probe)
    assert not result.passed
    assert "runs away" in result.message


@pytest.mark.parametrize("days", [730, 1095, 1200])
def test_a_season_is_not_a_drift(probe, case, days):
    def store(n):
        return 100.0 + 100.0 * np.sin(2.0 * np.pi * np.arange(n) / 365.0)

    result = state_bounds(over_capacity(case, case.spinup_steps + days, store), probe)
    assert result.passed
    assert "mrso upper bound N/A" in result.message

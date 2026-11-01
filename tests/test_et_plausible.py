"""et_plausible reads the record's total, within what single precision can hold."""

from __future__ import annotations

import runpy

import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.base import FAIL, PASS, get
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


@pytest.fixture
def table(case):
    simulate = runpy.run_path(str(REPO_ROOT / "models/reference_bucket/ht_adapter.py"))["simulate"]
    return pd.DataFrame(simulate(case.forcing.to_dict("records"), case.static))


def scored(table, case, probe):
    params = next(c.params for c in probe.criteria if c.name == "et_plausible")
    return get("et_plausible")(RunResult(case, table, {}, 0), probe, params)


def test_a_row_of_condensation_is_allowed(table, case, probe):
    table.loc[len(table) - 10, "evspsbl"] = -0.5
    assert scored(table, case, probe).status == PASS


def test_a_record_that_takes_water_from_the_air_overall_fails(table, case, probe):
    table["evspsbl"] = -0.1
    result = scored(table, case, probe)
    assert result.status == FAIL
    assert "negative" in result.message


@pytest.mark.parametrize("excess, status", [(5e-8, PASS), (1e-4, FAIL)])
def test_the_ceiling_allows_single_precision_and_nothing_more(table, case, probe, excess, status):
    table["evspsbl"] = case.forcing["pet"].to_numpy(dtype=float) * (1.0 + excess)
    assert scored(table, case, probe).status == status

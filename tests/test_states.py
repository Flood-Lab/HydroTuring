"""A model may carry stores the probe did not name, and the budget must
count them. A groundwater zone or water in transit reported under `gw` or
`channel` is part of the storage change; leaving it out opens the budget by
exactly the water it holds."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult
from hydroturing.runner import get_runner
from hydroturing.criteria.base import PASS


@pytest.fixture(scope="module")
def bucket_run():
    probe = registry.find_probe("mass/catchment-closure")
    model = registry.find_model("reference_bucket")
    case = build_case(probe, 20260903)
    run = get_runner(model).run(model, probe, case, Path(tempfile.mkdtemp()))
    return probe, run


def _closure(probe, run):
    params = dict(next(c.params for c in probe.criteria if c.name == "closure"))
    return get("closure")(run, probe, params)


def test_a_store_reported_under_another_name_still_closes(bucket_run):
    probe, run = bucket_run
    table = run.table.copy()
    # Split the soil column: half stays soil water, half is reported as the
    # groundwater the probe never asked for. The water is all still there.
    table["gw"] = 0.5 * table["mrso"]
    table["mrso"] = 0.5 * table["mrso"]
    split = RunResult(case=run.case, table=table, meta=run.meta, wall_seconds=0.0)
    assert _closure(probe, split).status == PASS


def test_a_store_left_out_opens_the_budget(bucket_run):
    probe, run = bucket_run
    table = run.table.copy()
    table["mrso"] = 0.5 * table["mrso"]  # and the other half goes unreported
    missing = RunResult(case=run.case, table=table, meta=run.meta, wall_seconds=0.0)
    result = _closure(probe, missing)
    exact = _closure(probe, run)
    assert result.value > exact.value

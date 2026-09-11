"""Non-finite model outputs are scored violations, not harness failures."""

from __future__ import annotations

from dataclasses import asdict
import json

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.symmetry import invariance
from hydroturing.harness import run_probe
from hydroturing.protocol import Case, RunResult, read_result
from hydroturing.runner import get_runner
from hydroturing.scoring import FAIL, VIOLATION


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/time-origin-invariance")


def _runs(control, shifted, *, spinup_steps=0):
    runs = {}
    for variant, values in (("control", control), ("shifted", shifted)):
        table = pd.DataFrame(values)
        forcing = pd.DataFrame({"time": pd.date_range("2000-01-01", periods=len(table))})
        case = Case("test/invariance", 0, forcing, {}, spinup_steps)
        runs[variant] = RunResult(case, table, {}, 0.0)
    return runs


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
@pytest.mark.parametrize("variant", ["control", "shifted"])
@pytest.mark.parametrize("variable", ["mrro", "gw"])
def test_nonfinite_output_names_variable_and_actual_variant(probe, value, variant, variable):
    runs = _runs({"mrro": [0.0], "gw": [0.0]}, {"mrro": [0.0], "gw": [0.0]})
    runs[variant].table.loc[0, variable] = value
    result = invariance(
        runs, probe,
        {"transformed": "shifted", "unchanged": ["mrro"], "optional": ["gw"]},
    )
    assert not result.passed
    assert result.name == "invariance" and result.value is None
    assert result.diagnostics == {
        "variable": variable, "variant": variant, "nonfinite_count": 1,
    }
    assert variable in result.message and variant in result.message
    # An undefined deviation must not put NaN or Infinity into score JSON.
    json.dumps(asdict(result), allow_nan=False)


@pytest.mark.parametrize("variable,value", [("mrro", np.inf), ("channel", np.nan)])
def test_nonfinite_output_preserves_full_scorecard(probe, monkeypatch, variable, value):
    model = registry.find_model("reference_bucket")
    runner = get_runner(model)

    class InvalidOutputRunner:
        def run(self, model, probe, case, io_dir):
            run = runner.run(model, probe, case, io_dir)
            if case.probe_id.endswith("@shifted"):
                frame = run.table.copy()
                if variable not in frame:
                    frame[variable] = 0.0
                frame.loc[case.spinup_steps, variable] = value
                frame.to_csv(io_dir / "output/result.csv", index=False)
                # Exercise the actual output reader before harness scoring.
                return read_result(io_dir, case, probe, run.wall_seconds)
            if variable not in run.table:
                run.table[variable] = 0.0
            return run

    monkeypatch.setattr("hydroturing.harness.get_runner", lambda model: InvalidOutputRunner())
    outcome = run_probe(model, probe, [20260908])
    assert outcome.verdict == FAIL and outcome.reason == VIOLATION, outcome.error
    assert outcome.error is None
    assert outcome.failing == ["invariance"]
    results = {criterion.name: criterion for criterion in outcome.criteria}
    assert set(results) == {"invariance", "closure", "non_degenerate"}
    assert results["closure"].passed and results["non_degenerate"].passed
    assert results["invariance"].diagnostics["variable"] == variable
    assert results["invariance"].diagnostics["variant"] == "shifted"
    assert results["invariance"].per_seed[0]["status"] == "fail"
    json.dumps(asdict(outcome), allow_nan=False)


@pytest.mark.parametrize("factor", [np.nan, np.inf, -np.inf])
def test_nonfinite_scale_is_a_configuration_error(probe, factor):
    runs = _runs({"dis": [0.0]}, {"dis": [0.0]})
    with pytest.raises(ValueError, match="finite scaled factor"):
        invariance(runs, probe, {"transformed": "shifted", "scaled": {"dis": factor}})


@pytest.mark.parametrize(
    "control, shifted, factor",
    [
        ([1e308, 1e308], [5e307, 5e307], 1.0),
        ([1e308], [1.0], 2.0),
        ([-1e308], [1e308], 1.0),
    ],
)
def test_finite_inputs_with_overflow_cannot_silently_pass(probe, control, shifted, factor):
    runs = _runs({"dis": control}, {"dis": shifted})
    with pytest.raises(ValueError, match="non-finite intermediate or deviation"):
        invariance(runs, probe, {"transformed": "shifted", "scaled": {"dis": factor}})

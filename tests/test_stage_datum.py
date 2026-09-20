"""Stage diagnostics are invariant to the chosen vertical datum."""

from __future__ import annotations

import importlib.util
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.harness import build_case
from hydroturing.protocol import Case, RunResult
from hydroturing.seeds import gate_seeds


def _run(bed_elevation_m: float) -> RunResult:
    time = pd.date_range("2001-01-01", periods=60, freq="D")
    depth = 1.0 + 0.3 * np.sin(np.linspace(0.0, 4.0 * np.pi, len(time)))
    forcing = pd.DataFrame({"time": time, "pr": np.full(len(time), 2.0)})
    case = Case(
        probe_id="momentum/stage-discharge-monotonic",
        seed=1,
        forcing=forcing,
        static={"bed_elevation_m": bed_elevation_m},
        spinup_steps=0,
        timestep="PT1D",
    )
    table = pd.DataFrame({"time": time, "stage": bed_elevation_m + depth})
    return RunResult(case, table, {}, 0.0)


def test_stage_variability_uses_depth_not_absolute_elevation():
    params = {
        "cv_variables": ["stage"],
        "min_flux_cv": 0.1,
        "min_response": None,
    }
    local = get("non_degenerate")(_run(0.0), None, params)
    elevated = get("non_degenerate")(_run(250.0), None, params)
    assert elevated.status == local.status
    assert elevated.diagnostics["cv_stage"] == pytest.approx(
        local.diagnostics["cv_stage"], abs=1e-12
    )


PROBE_ID = "momentum/stage-discharge-monotonic"


@pytest.fixture(scope="module")
def stage_probe():
    return registry.find_probe(PROBE_ID)


def test_stage_probe_requires_the_datum_its_criteria_form_a_depth_from(stage_probe):
    """`non_degenerate` takes CV(stage), which no datum cancels out of.

    Without the declaration a model carrying its own vertical datum reads as a
    flat gauge and fails `non_degenerate` for a difference of convention, which
    is the `reference_flat_stage` verdict and not the one it has earned.
    """
    assert "bed_elevation_m" in stage_probe.requires_static
    case = build_case(stage_probe, gate_seeds(stage_probe.id, 1)[0])
    assert "bed_elevation_m" in case.static


def _rating_reference(static: dict, case) -> pd.DataFrame:
    path = registry.find_model("reference_rating").path / "ht_adapter.py"
    spec = importlib.util.spec_from_file_location("rating_reference", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    forcing = case.forcing.to_dict("records")
    return pd.DataFrame(module.simulate(forcing, static))


def test_every_criterion_on_this_probe_survives_a_shifted_datum(stage_probe):
    """Move the case's datum 150 m and nothing about the verdict may move.

    The rating criteria work on differences of stage, so a datum cancels there
    by construction; `non_degenerate` forms a ratio of the level and only
    cancels because it subtracts the declared bed. This is the end-to-end
    statement of that: same model, same weather, a datum 150 m higher.
    """
    case = build_case(stage_probe, gate_seeds(stage_probe.id, 1)[0])
    shifted = dict(case.static)
    shifted["bed_elevation_m"] = float(case.static["bed_elevation_m"]) + 150.0

    at_bed = _rating_reference(dict(case.static), case)
    elevated = _rating_reference(shifted, case)
    assert np.allclose(
        elevated["stage"] - shifted["bed_elevation_m"],
        at_bed["stage"] - float(case.static["bed_elevation_m"]),
    )

    for criterion in stage_probe.criteria:
        scored = [
            get(criterion.name)(
                RunResult(replace(case, static=static), table, {}, 0.0),
                stage_probe,
                dict(criterion.params),
            )
            for static, table in ((dict(case.static), at_bed), (shifted, elevated))
        ]
        local, raised = scored
        assert raised.status == local.status, criterion.name
        assert raised.value == pytest.approx(local.value), criterion.name

"""Stage diagnostics are invariant to the chosen vertical datum."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing.criteria import get
from hydroturing.protocol import Case, RunResult


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

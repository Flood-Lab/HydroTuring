"""Focused counterexamples for the steady uniform-flow criterion."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.protocol import Case, RunResult


STATIC = {
    "width_m": 100.0,
    "bed_elevation_m": 75.0,
    "slope": 0.0015,
    "manning_n": 0.035,
}


def _capacity(depth: float, slope: float, roughness: float = 0.035) -> float:
    area = STATIC["width_m"] * depth
    radius = area / (STATIC["width_m"] + 2.0 * depth)
    return area * radius ** (2.0 / 3.0) * math.sqrt(slope) / roughness


def _depth(discharge: float, slope: float, roughness: float = 0.035) -> float:
    lo, hi = 0.0, 1.0
    while _capacity(hi, slope, roughness) < discharge:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _capacity(mid, slope, roughness) < discharge:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _run(discharge, stage, static=None) -> RunResult:
    discharge = np.asarray(discharge, dtype=float)
    stage = np.asarray(stage, dtype=float)
    n = len(discharge)
    time = pd.date_range("2001-01-01", periods=n, freq="D")
    case = Case(
        probe_id="momentum/uniform-flow-friction-consistency",
        seed=1,
        forcing=pd.DataFrame({"time": time}),
        static=dict(STATIC if static is None else static),
        spinup_steps=0,
        timestep="PT1D",
    )
    table = pd.DataFrame({"time": time, "dis": discharge, "stage": stage})
    return RunResult(case, table, {}, 0.0)


def _params(**overrides):
    params = {"tolerance": 0.05, "max_cv": 0.01, "steady_days": 365}
    params.update(overrides)
    return params


def test_exact_rectangular_normal_depth_passes():
    q = np.full(400, 12.0)
    stage = np.full(400, STATIC["bed_elevation_m"] + _depth(12.0, STATIC["slope"]))
    result = get("uniform_flow_friction")(_run(q, stage), None, _params())
    assert result.status == PASS
    assert result.value < 1e-10
    assert result.diagnostics["steady_steps"] == 365


def test_a_rating_built_with_the_wrong_roughness_fails():
    q = np.full(400, 12.0)
    wrong_depth = _depth(12.0, STATIC["slope"], roughness=0.08)
    stage = np.full(400, STATIC["bed_elevation_m"] + wrong_depth)
    result = get("uniform_flow_friction")(_run(q, stage), None, _params())
    assert result.status == FAIL
    assert result.value > 0.5


def test_a_rating_built_with_the_wrong_slope_fails():
    q = np.full(400, 12.0)
    wrong_depth = _depth(12.0, slope=0.01)
    stage = np.full(400, STATIC["bed_elevation_m"] + wrong_depth)
    result = get("uniform_flow_friction")(_run(q, stage), None, _params())
    assert result.status == FAIL
    assert result.value > 1.0


def test_opposite_signed_residuals_cannot_cancel():
    q = np.full(400, 12.0)
    depths = np.concatenate([
        np.full(200, _depth(12.0, slope=0.5 * STATIC["slope"])),
        np.full(200, _depth(12.0, slope=1.5 * STATIC["slope"])),
    ])
    stage = STATIC["bed_elevation_m"] + depths
    result = get("uniform_flow_friction")(
        _run(q, stage), None, _params(max_cv=10.0)
    )
    assert result.status == FAIL
    assert result.value > 0.45


def test_a_transient_final_block_is_refused_before_the_balance_is_scored():
    q = np.linspace(5.0, 15.0, 400)
    stage = STATIC["bed_elevation_m"] + np.array([
        _depth(value, STATIC["slope"]) for value in q
    ])
    result = get("uniform_flow_friction")(_run(q, stage), None, _params())
    assert result.status == FAIL
    assert "not steady" in result.message


def test_dry_or_nonfinite_answers_fail_instead_of_being_skipped():
    q = np.full(400, 12.0)
    wet_stage = STATIC["bed_elevation_m"] + _depth(12.0, STATIC["slope"])

    dry = np.full(400, wet_stage)
    dry[-1] = STATIC["bed_elevation_m"]
    dry_result = get("uniform_flow_friction")(_run(q, dry), None, _params())
    assert dry_result.status == FAIL
    assert "at or below the bed" in dry_result.message

    nonfinite = np.full(400, wet_stage)
    nonfinite[-1] = np.nan
    nonfinite_result = get("uniform_flow_friction")(
        _run(q, nonfinite), None, _params()
    )
    assert nonfinite_result.status == FAIL
    assert "finite" in nonfinite_result.message

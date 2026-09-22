"""Focused counterexamples for the steady uniform-flow criterion."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import depth_series, get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.harness import build_case, run_probe
from hydroturing.protocol import Case, RunResult
from hydroturing.runner import get_runner
from hydroturing.scoring import INCOMPATIBLE, NOT_SCORED, PASS as PROBE_PASS
from hydroturing.seeds import gate_seeds


STATIC = {
    "width_m": 100.0,
    "cross_section_shape": "rectangular",
    "bed_elevation_m": 75.0,
    "slope": 0.0015,
    "manning_n": 0.035,
}
LABELS = ("low", "medium", "high")
ROWS_PER_PLATEAU = 120


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


def _steady_series(
    discharges=(6.0, 12.0, 24.0),
    *,
    slopes=None,
    roughnesses=None,
    rows_per_plateau=ROWS_PER_PLATEAU,
):
    slopes = slopes or (STATIC["slope"],) * 3
    roughnesses = roughnesses or (STATIC["manning_n"],) * 3
    q = np.concatenate([
        np.full(rows_per_plateau, value) for value in discharges
    ])
    stage = np.concatenate([
        np.full(
            rows_per_plateau,
            STATIC["bed_elevation_m"] + _depth(discharge, slope, roughness),
        )
        for discharge, slope, roughness in zip(discharges, slopes, roughnesses)
    ])
    return q, stage


def _run(discharge, stage, static=None, labels=None) -> RunResult:
    discharge = np.asarray(discharge, dtype=float)
    stage = np.asarray(stage, dtype=float)
    n = len(discharge)
    time = pd.date_range("2001-01-01", periods=n, freq="D")
    if labels is None:
        if n % 3:
            raise ValueError("test series must split into three plateaus")
        labels = np.repeat(LABELS, n // 3)
    forcing = pd.DataFrame({"time": time, "_plateau": labels})
    case = Case(
        probe_id="momentum/uniform-flow-friction-consistency",
        seed=1,
        forcing=forcing,
        static=dict(STATIC if static is None else static),
        spinup_steps=0,
        timestep="PT1D",
    )
    table = pd.DataFrame({"time": time, "dis": discharge, "stage": stage})
    return RunResult(case, table, {}, 0.0)


def _params(**overrides):
    params = {
        "tolerance": 0.05,
        "max_cv": 0.015,
        "max_relative_trend": 0.01,
        "steady_days": 90,
        "plateaus": list(LABELS),
    }
    params.update(overrides)
    return params


def test_three_exact_rectangular_normal_depths_pass():
    q, stage = _steady_series()
    result = get("uniform_flow_friction")(_run(q, stage), None, _params())
    assert result.status == PASS
    assert result.value < 1e-10
    assert result.diagnostics["steady_steps"] == 90
    assert set(result.diagnostics["plateaus"]) == set(LABELS)
    assert result.diagnostics["plateaus"]["high"]["max_froude_number"] < 1.0


def test_case_allows_five_slowest_store_time_constants_before_scoring():
    probe = registry.find_probe("momentum/uniform-flow-friction-consistency")
    case = build_case(probe, gate_seeds(probe.id, 1)[0])
    scored = case.forcing.iloc[case.spinup_steps:]
    counts = scored["_plateau"].value_counts()
    steady_days = int(probe.criteria[0].params["steady_days"])
    assert float(probe.criteria[0].params["max_cv"]) == 0.015
    assert probe.min_scored_fraction == 1.0
    minimum_settling = math.ceil(
        5.0 / float(case.static["baseflow_coefficient"])
    )
    assert counts.to_dict() == {label: 924 for label in LABELS}
    assert int(counts.min()) - steady_days >= minimum_settling


@pytest.mark.parametrize(
    ("model_name", "seeds"),
    [
        ("reference_bucket", [0, 22, 23, 27, 30, 39]),
        ("sacsma_snow17", [4, 13, 40, 56]),
    ],
)
def test_extended_plateaus_score_previously_nonsteady_must_pass_runs(
    model_name, seeds, tmp_path,
):
    probe = registry.find_probe("momentum/uniform-flow-friction-consistency")
    outcome = run_probe(
        registry.find_model(model_name),
        probe,
        seeds,
        workdir=tmp_path,
    )
    assert outcome.verdict == PROBE_PASS


def test_near_boundary_wrong_roughness_fails():
    q, stage = _steady_series(roughnesses=(0.035 * 1.03,) * 3)
    result = get("uniform_flow_friction")(_run(q, stage), None, _params())
    assert result.status == FAIL
    assert 0.05 < result.value < 0.07


def test_near_boundary_wrong_slope_fails():
    q, stage = _steady_series(slopes=(STATIC["slope"] * 1.06,) * 3)
    result = get("uniform_flow_friction")(_run(q, stage), None, _params())
    assert result.status == FAIL
    assert result.value == pytest.approx(0.06, abs=1e-10)


def test_opposite_signed_residuals_cannot_cancel():
    q, stage = _steady_series()
    # Put equal and opposite signed slope errors inside each scored block.
    for start in range(0, len(q), ROWS_PER_PLATEAU):
        discharge = q[start]
        stage[start + 30:start + 75] = (
            STATIC["bed_elevation_m"] + _depth(discharge, 0.5 * STATIC["slope"])
        )
        stage[start + 75:start + 120] = (
            STATIC["bed_elevation_m"] + _depth(discharge, 1.5 * STATIC["slope"])
        )
    result = get("uniform_flow_friction")(
        _run(q, stage), None, _params(max_cv=10.0, max_relative_trend=10.0)
    )
    assert result.status == FAIL
    assert result.value > 0.45


def test_a_slow_drift_skips_only_its_plateau_even_when_cv_is_loose():
    q, stage = _steady_series()
    q[-ROWS_PER_PLATEAU:] = np.linspace(22.0, 24.0, ROWS_PER_PLATEAU)
    stage[-ROWS_PER_PLATEAU:] = STATIC["bed_elevation_m"] + np.array([
        _depth(value, STATIC["slope"]) for value in q[-ROWS_PER_PLATEAU:]
    ])
    result = get("uniform_flow_friction")(
        _run(q, stage), None, _params(max_cv=1.0)
    )
    assert result.status == PASS
    assert result.diagnostics["steady_plateaus"] == ["low", "medium"]
    assert result.diagnostics["skipped_plateaus"] == ["high"]
    assert "skipped non-steady plateau(s): high" in result.message


def test_a_nonsteady_plateau_cannot_hide_steady_friction_failures():
    q, stage = _steady_series(slopes=(0.25 * STATIC["slope"],) * 3)
    high = slice(-ROWS_PER_PLATEAU, None)
    ripple = 1.0 + 0.03 * np.sin(np.linspace(0.0, 12.0 * np.pi, ROWS_PER_PLATEAU))
    q[high] *= ripple

    result = get("uniform_flow_friction")(_run(q, stage), None, _params())

    assert result.status == FAIL
    assert result.value == pytest.approx(0.75, abs=1e-10)
    assert result.diagnostics["steady_plateaus"] == ["low", "medium"]
    assert result.diagnostics["skipped_plateaus"] == ["high"]
    assert "low residual 75.00%" in result.message
    assert "skipped non-steady plateau(s): high" in result.message


def test_dry_or_nonfinite_answers_fail_instead_of_being_skipped():
    q, stage = _steady_series()
    dry = stage.copy()
    dry[-1] = STATIC["bed_elevation_m"]
    dry_result = get("uniform_flow_friction")(_run(q, dry), None, _params())
    assert dry_result.status == FAIL
    assert "at or below the bed" in dry_result.message

    nonfinite = stage.copy()
    nonfinite[-1] = np.nan
    nonfinite_result = get("uniform_flow_friction")(
        _run(q, nonfinite), None, _params()
    )
    assert nonfinite_result.status == FAIL
    assert "finite" in nonfinite_result.message


def test_malformed_case_inputs_raise_configuration_errors():
    q, stage = _steady_series()
    missing_width = dict(STATIC)
    del missing_width["width_m"]
    with pytest.raises(ValueError, match="width_m"):
        get("uniform_flow_friction")(
            _run(q, stage, static=missing_width), None, _params()
        )

    missing_high = np.repeat(("low", "medium", "medium"), ROWS_PER_PLATEAU)
    with pytest.raises(ValueError, match="high"):
        get("uniform_flow_friction")(
            _run(q, stage, labels=missing_high), None, _params()
        )

    short_q, short_stage = _steady_series(rows_per_plateau=60)
    with pytest.raises(ValueError, match="needs 90 rows"):
        get("uniform_flow_friction")(
            _run(short_q, short_stage), None, _params()
        )


@pytest.mark.parametrize("single_row", [None, "nan", "elevation"])
def test_depth_reported_in_the_stage_column_is_incompatible(
    monkeypatch, tmp_path, single_row,
):
    """A different datum convention is N/A, not a public physics violation."""
    from hydroturing import harness

    probe = registry.find_probe("momentum/uniform-flow-friction-consistency")
    model = registry.find_model("reference_uniform_flow")
    real_runner = get_runner(model)

    class DepthAsStageRunner:
        def run(self, model, probe, case, io_dir):
            result = real_runner.run(model, probe, case, io_dir)
            result.table["stage"] -= float(case.static["bed_elevation_m"])
            if single_row == "nan":
                result.table.loc[5, "stage"] = np.nan
            elif single_row == "elevation":
                result.table.loc[5, "stage"] += float(
                    case.static["bed_elevation_m"]
                )
            return result

    monkeypatch.setattr(harness, "get_runner", lambda _: DepthAsStageRunner())
    outcome = run_probe(
        model,
        probe,
        gate_seeds(probe.id, 1),
        workdir=tmp_path,
    )
    assert outcome.verdict == NOT_SCORED
    assert outcome.reason == INCOMPATIBLE
    assert "fixed vertical datum" in outcome.incompatible[0]


def test_nonsteady_output_is_incompatible_end_to_end(monkeypatch, tmp_path):
    from hydroturing import harness

    probe = registry.find_probe("momentum/uniform-flow-friction-consistency")
    model = registry.find_model("reference_uniform_flow")
    real_runner = get_runner(model)

    class DriftingRunner:
        def run(self, model, probe, case, io_dir):
            result = real_runner.run(model, probe, case, io_dir)
            for label in LABELS:
                plateau = np.flatnonzero(
                    case.forcing["_plateau"].to_numpy() == label
                )
                scored = plateau[-90:]
                result.table.loc[scored, "dis"] *= np.linspace(
                    0.8, 1.2, len(scored)
                )
            return result

    monkeypatch.setattr(harness, "get_runner", lambda _: DriftingRunner())
    outcome = run_probe(
        model,
        probe,
        gate_seeds(probe.id, 1),
        workdir=tmp_path,
    )
    assert outcome.verdict == NOT_SCORED
    assert outcome.reason == INCOMPATIBLE
    assert "uniform-flow precondition was not reached" in outcome.incompatible[0]
    assert "low residual" in outcome.incompatible[0]
    assert "medium residual" in outcome.incompatible[0]
    assert "high residual" in outcome.incompatible[0]


@pytest.mark.parametrize("negative_depth", [-1.0e-12, -0.2])
def test_small_negative_depth_is_scored_as_a_fault_not_a_datum_mismatch(
    negative_depth,
):
    q, _ = _steady_series()
    stage = np.full(len(q), STATIC["bed_elevation_m"] + negative_depth)
    derived = depth_series(_run(q, stage))
    assert np.allclose(derived, negative_depth)


def test_bed_zero_does_not_guess_between_depth_and_elevation_conventions():
    q, _ = _steady_series()
    static = dict(STATIC, bed_elevation_m=0.0)
    stage = np.full(len(q), -0.2)
    derived = depth_series(_run(q, stage, static=static))
    assert np.allclose(derived, -0.2)

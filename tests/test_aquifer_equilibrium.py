"""A prescribed head bounds storage independently of a closed water ledger."""

from __future__ import annotations

from dataclasses import replace
import importlib.util

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import CRITERIA, CriterionIncompatibleError
from hydroturing.protocol import Case, RunResult
from hydroturing.spec import REPO_ROOT, load_model


PROBE_ID = "mass/aquifer-equilibrium-bounds"
STATIC = {
    "area_km2": 1.0,
    "aquifer_storage_coefficient": 0.002,
    "aquifer_specific_yield": 0.2,
    "aquifer_initial_head_m": 95.0,
    "aquifer_top_m": 90.0,
    "aquifer_bottom_m": 80.0,
    "river_conductance_m2_per_day": 200.0,
    "river_bottom_offset_m": 4.0,
}


def _runs(raised=(0.3, 0.8, 1.2, 1.6, 1.8, 2.0),
          lowered=(-0.3, -0.8, -1.2, -1.6, -1.8, -2.0), *, offset=3000.0):
    """The literal capacity here is 2 mm for a one-metre head change."""
    n = len(raised)
    probe = replace(
        registry.find_probe("mass/groundwater-datum-invariance"),
        id=PROBE_ID, period_days=float(n), spinup_days=2,
        variants=("control", "raised", "lowered"),
    )
    runs = {}
    for variant, changes, stage in (
        ("control", np.zeros(n), 95.0),
        ("raised", raised, 96.0),
        ("lowered", lowered, 94.0),
    ):
        time = pd.date_range("2000-01-01", periods=n + 2).strftime("%Y-%m-%d")
        forcing = pd.DataFrame({
            "time": time, "gw_recharge": np.zeros(n + 2),
            "sw_stage_m": [95.0, 95.0, *([stage] * n)],
        })
        storage = offset + np.asarray([0.0, 0.0, *changes], dtype=float)
        q = np.diff(storage, prepend=offset)
        table = pd.DataFrame({
            "time": time, "gw": storage, "gw_sw_exchange": q,
            "gw_to_sw": np.minimum(q, 0.0), "sw_to_gw": np.maximum(q, 0.0),
        })
        case = Case(PROBE_ID, 7, forcing, dict(STATIC), 2)
        runs[variant] = RunResult(case, table, {}, 0.0)
    return runs, probe


def _score(runs, probe, **params):
    fn = CRITERIA.get("aquifer_equilibrium_bounds")
    assert fn is not None, "the paired hydraulic-envelope criterion is missing"
    return fn(runs, probe, params)


def _generator():
    path = REPO_ROOT / "probes/mass/aquifer-equilibrium-bounds/generate.py"
    assert path.exists(), "the deterministic confined-aquifer generator is missing"
    spec = importlib.util.spec_from_file_location("aquifer_equilibrium_generator", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_storage_within_both_prescribed_head_bounds_passes():
    runs, probe = _runs()
    result = _score(runs, probe)
    assert result.passed
    assert result.diagnostics["variants"]["raised"]["upper_bound_mm"] == 2.0
    assert result.diagnostics["variants"]["lowered"]["lower_bound_mm"] == -2.0


@pytest.mark.parametrize("variant,values", [
    ("raised", [0.2, 0.8, 1.5, 1.9, 2.0, 2.000201]),
    ("lowered", [-0.2, -0.8, -1.5, -1.9, -2.0, -2.000201]),
    ("raised", [-0.000201, 0.0, 0.2, 0.5, 0.7, 1.0]),
    ("lowered", [0.000201, 0.0, -0.2, -0.5, -0.7, -1.0]),
])
def test_storage_on_either_wrong_side_of_envelope_fails(variant, values):
    runs, probe = _runs(**{variant: values})
    result = _score(runs, probe)
    assert not result.passed
    assert result.value > 1.0
    assert result.diagnostics["variants"][variant]["allowed_mm"] == pytest.approx(0.0002)


def test_temporary_overshoot_cannot_be_cancelled_at_endpoint():
    runs, probe = _runs(raised=[0.5, 2.5, 1.8, 1.7, 1.6, 1.6])
    result = _score(runs, probe)
    assert not result.passed
    assert result.diagnostics["variants"]["raised"]["max_excess_mm"] == pytest.approx(0.5)
    assert result.diagnostics["variants"]["raised"]["worst_step"] == 1


@pytest.mark.parametrize("offset", [30.0, 3000.0, 1.0e9])
def test_large_absolute_inventory_cannot_buy_more_tolerance(offset):
    runs, probe = _runs(raised=[0.0, 0.0, 0.0, 0.0, 0.0, 2.2], offset=offset)
    result = _score(runs, probe)
    assert not result.passed
    assert result.diagnostics["variants"]["raised"]["allowed_mm"] == pytest.approx(0.0002)


def test_exact_binary_boundary_of_configured_allowance_passes():
    # B=2 mm, rel_tol=2^-13: the permitted excess is exactly 2^-12 mm.
    runs, probe = _runs(raised=[0.0, 0.0, 0.0, 0.0, 0.0, 2.000244140625])
    assert _score(runs, probe, rel_tol=0.0001220703125, abs_tol_mm=0.0001).passed


@pytest.mark.parametrize("excess,passes", [(0.00009, True), (0.00011, False)])
def test_absolute_precision_floor_accepts_roundoff_but_not_larger_excess(excess, passes):
    runs, probe = _runs(raised=[0.0, 0.0, 0.0, 0.0, 0.0, 2.0 + excess])
    assert _score(runs, probe, rel_tol=0.0).passed is passes


@pytest.mark.parametrize("column", ["gw", "gw_sw_exchange"])
def test_neutral_precision_guard_rejects_activity_above_a_tenth_micron(column):
    runs, probe = _runs()
    runs["control"].table.loc[3, column] += 0.00011
    assert not _score(runs, probe).passed


def test_common_spinup_inventory_mismatch_is_not_subtracted_away():
    runs, probe = _runs()
    runs["raised"].table["gw"] += 0.1
    result = _score(runs, probe)
    assert not result.passed
    assert "initial" in result.message


def test_small_spinup_inventory_mismatch_exceeds_the_precision_guard():
    runs, probe = _runs()
    runs["raised"].table["gw"] += 0.00011
    result = _score(runs, probe)
    assert not result.passed
    assert "initial" in result.message


def test_nonstationary_control_cannot_move_the_comparison_envelope():
    runs, probe = _runs()
    runs["control"].table.loc[3:, "gw"] += 1.0
    result = _score(runs, probe)
    assert not result.passed
    assert "control" in result.message


def test_neutral_flux_cannot_cancel_between_directions():
    runs, probe = _runs()
    runs["control"].table["gw_to_sw"] = -0.1
    runs["control"].table["sw_to_gw"] = 0.1
    result = _score(runs, probe)
    assert not result.passed
    assert "control" in result.message


def test_frozen_storage_is_reported_as_unexercised_not_failed():
    runs, probe = _runs(raised=[0.0] * 6, lowered=[0.0] * 6)
    result = _score(runs, probe)
    assert result.passed
    assert result.diagnostics["no_hydraulic_response_exercised"]


@pytest.mark.parametrize("variant", ["raised", "lowered"])
def test_near_frozen_arm_is_flagged_without_failing_a_necessary_bound(variant):
    values = [0.002] * 6 if variant == "raised" else [-0.002] * 6
    runs, probe = _runs(**{variant: values})
    result = _score(runs, probe)
    assert result.passed
    assert result.diagnostics["no_hydraulic_response_exercised"]
    detail = result.diagnostics["variants"][variant]
    assert detail["response_fraction"] == pytest.approx(0.001)
    assert detail["no_hydraulic_response_exercised"]
    assert variant in result.message
    other = "lowered" if variant == "raised" else "raised"
    assert not result.diagnostics["variants"][other]["no_hydraulic_response_exercised"]


def test_response_fraction_threshold_changes_only_the_descriptive_flag():
    runs, probe = _runs(raised=[0.03] * 6)
    ordinary = _score(runs, probe)
    labelled = _score(runs, probe, response_fraction_threshold=0.02)
    assert ordinary.passed and labelled.passed
    assert ordinary.value == labelled.value
    assert not ordinary.diagnostics["no_hydraulic_response_exercised"]
    assert labelled.diagnostics["no_hydraulic_response_exercised"]
    assert labelled.diagnostics["variants"]["raised"]["response_fraction"] == pytest.approx(0.015)


@pytest.mark.parametrize("threshold", [-0.01, 1.01, float("nan"), float("inf")])
def test_invalid_response_fraction_threshold_is_rejected(threshold):
    runs, probe = _runs()
    with pytest.raises(ValueError, match="response"):
        _score(runs, probe, response_fraction_threshold=threshold)


@pytest.mark.parametrize("n", [6, 90, 900])
def test_neutral_flux_precision_floor_does_not_shrink_with_record_length(n):
    runs, probe = _runs(raised=[0.1] * n, lowered=[-0.1] * n)
    # Alternating signed rounding residuals remain below 1e-4 mm each day.
    q = np.resize([1e-5, -1e-5], n + 2)
    runs["control"].table["gw_sw_exchange"] = q
    runs["control"].table["gw_to_sw"] = np.minimum(q, 0.0)
    runs["control"].table["sw_to_gw"] = np.maximum(q, 0.0)
    result = _score(runs, probe)
    assert result.passed
    assert result.diagnostics["control"]["max_interval_gross_exchange_mm"] == pytest.approx(1e-5)
    assert result.diagnostics["ratios"]["control exchange"] == pytest.approx(0.1)


def test_equilibrium_control_rounded_to_four_and_five_decimals_is_stationary():
    runs, probe = _runs(offset=3000.1234567)
    control = runs["control"].table
    control["gw"] = control["gw"].round(4)
    for column in ("gw_sw_exchange", "gw_to_sw", "sw_to_gw"):
        control[column] = control[column].round(5)
    assert _score(runs, probe).passed


def test_spinup_drift_is_not_hidden_by_a_later_stationary_control():
    runs, probe = _runs()
    # Every arm violates the prescribed equilibrium in the same first row.
    for run in runs.values():
        run.table.loc[0, "gw"] += 0.01
    result = _score(runs, probe)
    assert not result.passed
    assert "control storage" in result.message


def test_active_additional_boundary_is_incompatible_even_if_net_zero():
    runs, probe = _runs()
    runs["lowered"].table["gw_boundary"] = [0.0, 0.0, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0]
    with pytest.raises(CriterionIncompatibleError, match="boundary"):
        _score(runs, probe)


def test_declared_zero_additional_boundary_stays_eligible():
    runs, probe = _runs()
    for run in runs.values():
        run.table["gw_boundary"] = 0.0
    assert _score(runs, probe).passed


def test_missing_declared_boundary_is_a_contract_error():
    runs, probe = _runs()
    manifest = load_model(REPO_ROOT / "models/modflow6/model.yaml")
    runs["raised"].model = replace(manifest, emits_fluxes=manifest.emits_fluxes + ("gw_boundary",))
    with pytest.raises(ValueError, match="gw_boundary"):
        _score(runs, probe)


@pytest.mark.parametrize("column,value", [
    ("gw", np.nan), ("gw", np.inf), ("gw_sw_exchange", -np.inf),
    ("gw_boundary", np.nan),
])
def test_nonfinite_reported_values_never_pass(column, value):
    runs, probe = _runs()
    if column not in runs["raised"].table:
        runs["raised"].table[column] = 0.0
    runs["raised"].table.loc[4, column] = value
    assert not _score(runs, probe).passed


def test_subtracting_finite_extreme_states_cannot_overflow_into_a_pass():
    runs, probe = _runs()
    for run in runs.values():
        run.table["gw"] = 1.0e308
    runs["raised"].table.loc[3, "gw"] = -1.0e308
    assert not _score(runs, probe).passed


@pytest.mark.parametrize("mutation", ["time", "static", "recharge", "stage", "spinup"])
def test_invalid_paired_experiment_is_not_scored_as_physics(mutation):
    runs, probe = _runs()
    run = runs["raised"]
    if mutation == "time":
        run.case.forcing.loc[3, "time"] = "2000-02-01"
    elif mutation == "static":
        run.case.static["aquifer_storage_coefficient"] = 0.003
    elif mutation == "recharge":
        run.case.forcing.loc[4, "gw_recharge"] = 1.0
    elif mutation == "stage":
        run.case.forcing.loc[5, "sw_stage_m"] = 96.5
    else:
        run.case.forcing.loc[1, "sw_stage_m"] = 96.0
    with pytest.raises(ValueError):
        _score(runs, probe)


def test_generator_gives_identical_spinup_and_safe_confined_geometry():
    module = _generator()
    for seed in (0, 1, 7, 42, 2147483646):
        cases = {v: module.generate(seed, v) for v in ("control", "raised", "lowered")}
        control, static = cases["control"]
        assert len(control) == 92
        pd.testing.assert_frame_equal(control, module.generate(seed, "control")[0])
        stages = []
        for variant, (forcing, attrs) in cases.items():
            assert attrs == static
            assert len(forcing) == 92
            assert np.all(forcing["gw_recharge"] == 0.0)
            pd.testing.assert_frame_equal(forcing.iloc[:2], control.iloc[:2])
            assert forcing["time"].equals(control["time"])
            assert forcing["sw_stage_m"].iloc[2:].nunique() == 1
            stages.extend(forcing["sw_stage_m"])
        minimum_head = min(stages)
        maximum_bed = max(stages) - static["river_bottom_offset_m"]
        assert minimum_head > static["aquifer_top_m"]
        assert minimum_head - maximum_bed >= 1.0
        assert min(abs(cases[v][0]["sw_stage_m"].iloc[2] - static["aquifer_initial_head_m"])
                   for v in ("raised", "lowered")) * static["aquifer_storage_coefficient"] * 1000 >= 1.0


def test_generator_explicitly_rejects_insufficient_riverbed_clearance():
    module = _generator()
    with pytest.raises(ValueError, match="clearance"):
        module.check_geometry(dict(STATIC), 93.0, 97.0)


def test_generator_rejects_unknown_variant():
    module = _generator()
    with pytest.raises(ValueError, match="variant"):
        module.generate(7, "unexpected")

"""Unit tests for the paired routing-lag criteria.

The acceptance gate exercises real subprocess reference models.  These tests
keep the measurement itself pinned: Snyder geometry, baseline removal, event and
peak centroids, discrete-time tolerance, response guards, and malformed paired
cases each have a small counterexample here.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get, is_paired
from hydroturing.harness import build_case, compatibility_issues, run_probe
from hydroturing.protocol import Case, RunResult
from hydroturing.scoring import PASS
from hydroturing.seeds import gate_seeds
from hydroturing.spec import Criterion, ProbeSpec

VARIANTS = ("small", "medium", "large", "xlarge")
AREAS = (30.0, 300.0, 3000.0, 10000.0)
LENGTHS = (9.7965015486, 39.0005751285, 155.2640861436, 319.7409482737)
TARGET_LAGS = (1, 1, 2, 3)
SPINUP_DAYS = 5
PERIOD_DAYS = 90
EVENT = slice(SPINUP_DAYS + 20, SPINUP_DAYS + 21)

COMMON = {
    "event_column": "_event_pr",
    "runoff": "mrro",
    "baseline_days": 5,
    "min_response_fraction": 0.01,
    "min_pre_event_days": 10,
    "min_post_event_days": 30,
}

GEOMETRY_KEYS = {
    "area_km2",
    "main_channel_length_km",
    "centroid_channel_length_km",
}
EVALUATED_MODELS = (
    "cwatm",
    "dhbv2",
    "flex_topo",
    "google_flood_forecast",
    "lisflood",
    "sacsma_snow17",
    "summa",
    "wflow_sbm",
)


@pytest.fixture(scope="module")
def registered_probe() -> ProbeSpec:
    return registry.find_probe("momentum/routing-lag-consistency")


@pytest.fixture(scope="module")
def flex_adapter():
    path = registry.find_model("flex_lumped").path / "ht_adapter.py"
    spec = importlib.util.spec_from_file_location("routing_lag_flex_lumped", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def snyder_adapter():
    path = registry.find_model("reference_snyder_router").path / "ht_adapter.py"
    spec = importlib.util.spec_from_file_location("routing_lag_snyder_reference", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _probe() -> ProbeSpec:
    return ProbeSpec(
        id="momentum/routing-lag-consistency",
        title="Routing lag follows catchment geometry",
        law="momentum",
        track="synthetic",
        version=1,
        authors=({"name": "Test Author"},),
        citation="",
        requires_fluxes=("mrro",),
        requires_states=(),
        generator="generate.py",
        n_seeds=3,
        timestep="PT1D",
        period_years=PERIOD_DAYS / 365,
        spinup_days=SPINUP_DAYS,
        max_output_mb=2.0,
        max_runtime_s=120.0,
        variants=VARIANTS,
        criteria=(
            Criterion("lag_time_bounds", {}),
            Criterion("scaling_monotonicity", {}),
        ),
        must_pass=("reference_router",),
        must_fail={"reference_bad_router": "lag_time_bounds"},
        provenance="synthetic",
        path=Path("."),
        period_days=PERIOD_DAYS,
        min_window_days=PERIOD_DAYS,
    )


def _forcing() -> pd.DataFrame:
    n = SPINUP_DAYS + PERIOD_DAYS
    extra = np.zeros(n)
    extra[EVENT] = 50.0
    return pd.DataFrame(
        {
            "time": pd.date_range("2000-01-01", periods=n, freq="D").strftime(
                "%Y-%m-%d"
            ),
            "pr": extra.copy(),
            "tas": np.full(n, 15.0),
            "pet": np.zeros(n),
            "_event_pr": extra,
        }
    )


def _run(variant: str, peak_lag_days: int, *, tie: bool = False) -> RunResult:
    forcing = _forcing()
    baseline = 0.2
    runoff = np.full(len(forcing), baseline)
    # The one-row rain pulse has its centroid in that event row.
    centre = EVENT.start + peak_lag_days
    if tie:
        runoff[centre : centre + 2] += 5.0
    else:
        runoff[centre - 1 : centre + 2] += (2.0, 5.0, 2.0)
    index = VARIANTS.index(variant)
    static = {
        "area_km2": AREAS[index],
        "main_channel_length_km": LENGTHS[index],
        "centroid_channel_length_km": 0.5 * LENGTHS[index],
    }
    case = Case(
        probe_id=f"momentum/routing-lag-consistency@{variant}",
        seed=11,
        forcing=forcing,
        static=static,
        spinup_steps=SPINUP_DAYS,
        timestep="PT1D",
    )
    table = pd.DataFrame({"time": forcing["time"], "mrro": runoff})
    return RunResult(case, table, {"status": "ok"}, 0.0)


def _runs(lags=TARGET_LAGS) -> dict[str, RunResult]:
    return {name: _run(name, lag) for name, lag in zip(VARIANTS, lags)}


def test_criteria_are_registered_as_paired():
    assert is_paired("lag_time_bounds")
    assert is_paired("scaling_monotonicity")


def test_registered_probe_declares_every_input_its_verdict_depends_on(registered_probe):
    assert registered_probe.requires_forcing == ("pr",)
    assert registered_probe.requires_static == (
        "area_km2",
        "main_channel_length_km",
        "centroid_channel_length_km",
    )


def test_registered_variants_share_forcing_and_change_only_geometry(registered_probe):
    cases = {
        variant: build_case(registered_probe, 20260914, variant)
        for variant in registered_probe.variants
    }
    control = cases[registered_probe.variants[0]]

    for case in cases.values():
        pd.testing.assert_frame_equal(case.forcing, control.forcing, check_exact=True)
        assert set(case.static) == set(control.static)
        assert {
            key: value for key, value in case.static.items() if key not in GEOMETRY_KEYS
        } == {
            key: value for key, value in control.static.items() if key not in GEOMETRY_KEYS
        }

    geometries = [
        tuple(cases[variant].static[key] for key in GEOMETRY_KEYS)
        for variant in registered_probe.variants
    ]
    assert len(set(geometries)) == len(registered_probe.variants)
    for key in GEOMETRY_KEYS:
        values = [cases[variant].static[key] for variant in registered_probe.variants]
        assert values == sorted(values)
        assert len(set(values)) == len(registered_probe.variants)

    scored = control.after_spinup(control.forcing)
    assert (scored["_event_pr"] > 0.0).sum() == 1
    np.testing.assert_array_equal(scored["pr"], scored["_event_pr"])


@pytest.mark.parametrize(
    "model_name",
    [
        "reference_snyder_router",
        "flex_lumped",
        "reference_instant_router",
        "reference_inverse_router",
    ],
)
def test_routing_references_declare_compatible_inputs(registered_probe, model_name):
    case = build_case(registered_probe, 20260914, "small")
    assert compatibility_issues(
        registry.find_model(model_name), registered_probe, case
    ) == []


def test_flex_lumped_passes_all_registered_gate_seeds(registered_probe):
    outcome = run_probe(
        registry.find_model("flex_lumped"),
        registered_probe,
        gate_seeds(registered_probe.id, registered_probe.n_seeds),
    )
    assert outcome.verdict == PASS, outcome.reason


def test_flex_lumped_passes_a_seed_whose_rain_centroid_rounds_late(registered_probe):
    # Seed 887452083 puts a 49.939823 mm storm on scored day 21, whose volume
    # centroid computes as 21.500000000000004. FLEX's small catchment peaks in
    # the storm row, so its lag measured -8.5e-14 h against a lower bound of 0.
    outcome = run_probe(
        registry.find_model("flex_lumped"),
        registered_probe,
        [887452083],
    )
    assert outcome.verdict == PASS, outcome.reason


def test_a_peak_in_the_storm_row_is_not_rejected_by_centroid_rounding():
    runs = {}
    for name, run in _runs((0, 1, 2, 3)).items():
        forcing = run.case.forcing.copy()
        # This depth makes the one-row centroid 20.500000000000004 days.
        forcing.loc[EVENT.start, ["pr", "_event_pr"]] = 49.95412
        case = Case(
            probe_id=run.case.probe_id,
            seed=run.case.seed,
            forcing=forcing,
            static=run.case.static,
            spinup_steps=SPINUP_DAYS,
            timestep="PT1D",
        )
        runs[name] = RunResult(case, run.table, run.meta, run.wall_seconds)
    result = get("lag_time_bounds")(
        runs,
        _probe(),
        {**COMMON, "lower_ratio": 0.5, "upper_ratio": 2.0,
         "discretization_tolerance_days": 0.5},
    )
    small = result.diagnostics["variants"]["small"]
    assert small["rain_centroid_day"] != 20.5
    assert result.passed, result.message


def test_flex_lumped_maps_complete_geometry_to_its_native_lag(flex_adapter):
    static = {
        "area_km2": AREAS[2],
        "main_channel_length_km": LENGTHS[2],
        "centroid_channel_length_km": 0.5 * LENGTHS[2],
    }
    dt_days = 1.0
    base_days, routing = flex_adapter.routing_parameters(static, dt_days)
    expected_travel_days = static["centroid_channel_length_km"] / 86.4
    expected_mode_from_row_start = expected_travel_days + 0.5 * dt_days
    assert base_days == pytest.approx(2.0 * expected_mode_from_row_start)
    assert routing["flood_wave_celerity_m_s"] == pytest.approx(1.0)
    assert routing["celerity_reference"] == (
        "https://doi.org/10.5194/hess-24-2655-2020"
    )
    assert "not site-specific" in routing["celerity_assumption"]
    assert routing["travel_distance_km"] == pytest.approx(
        static["centroid_channel_length_km"]
    )
    assert "centroid-to-outlet" in routing["distance_definition"]
    assert routing["channel_travel_time_days"] == pytest.approx(
        expected_travel_days
    )
    assert routing["travel_time_origin"] == "generated-runoff interval centroid"
    assert not any("snyder" in key.lower() for key in routing)
    assert routing["source_interval_centroid_offset_days"] == pytest.approx(0.5)


@pytest.mark.parametrize("dt_days", [1.0, 1.0 / 24.0])
def test_flex_lag_kernel_represents_channel_travel_from_source_row_centroid(
    flex_adapter, dt_days
):
    static = {
        "area_km2": AREAS[2],
        "main_channel_length_km": LENGTHS[2],
        "centroid_channel_length_km": 0.5 * LENGTHS[2],
    }
    base_days, routing = flex_adapter.routing_parameters(static, dt_days)
    weights = np.asarray(flex_adapter.lag_weights(base_days / dt_days))
    peak = np.max(weights)
    tied = np.flatnonzero(np.isclose(weights, peak, rtol=1e-12, atol=1e-12))
    peak_centroid_days = float(np.mean((tied + 0.5) * dt_days))
    source_centroid_days = 0.5 * dt_days
    measured_lag_days = peak_centroid_days - source_centroid_days
    assert measured_lag_days == pytest.approx(
        routing["channel_travel_time_days"], abs=0.5 * dt_days
    )


def test_flex_lumped_keeps_wark_lag_when_only_area_is_available(flex_adapter):
    base_days, routing = flex_adapter.routing_parameters(
        {"area_km2": 2500.0}, 1.0
    )
    assert base_days == pytest.approx(1.1)
    assert routing["parameter_source"] == "calibrated Wark fallback"


@pytest.mark.parametrize(
    "static",
    [
        {"area_km2": None},
        {"area_km2": 300.0, "main_channel_length_km": 39.0},
        {
            "area_km2": 300.0,
            "main_channel_length_km": 39.0,
            "centroid_channel_length_km": 50.0,
        },
    ],
)
def test_flex_lumped_rejects_invalid_or_partial_geometry(flex_adapter, static):
    with pytest.raises(ValueError):
        flex_adapter.routing_parameters(static, 1.0)


@pytest.mark.parametrize(
    "seed", gate_seeds("momentum/routing-lag-consistency", 3)
)
def test_ct2_snyder_reference_uses_rainfall_centroid_origin(
    registered_probe, snyder_adapter, seed
):
    runs = {}
    routing = {}
    for variant in registered_probe.variants:
        case = build_case(registered_probe, seed, variant)
        staged_forcing = case.forcing.loc[
            :, [name for name in case.forcing.columns if not name.startswith("_")]
        ]
        rows, routing[variant] = snyder_adapter.simulate(
            staged_forcing.to_dict("records"), case.static, 1.0, ct=2.0
        )
        runs[variant] = RunResult(
            case,
            pd.DataFrame(rows),
            {"status": "ok", "routing": routing[variant]},
            0.0,
        )

    params = {item.name: item.params for item in registered_probe.criteria}
    bounds = get("lag_time_bounds")(
        runs, registered_probe, dict(params["lag_time_bounds"])
    )
    scaling = get("scaling_monotonicity")(
        runs, registered_probe, dict(params["scaling_monotonicity"])
    )
    assert bounds.passed, bounds.message
    assert scaling.passed, scaling.message
    observed = [
        bounds.diagnostics["variants"][variant]["observed_lag_days"]
        for variant in registered_probe.variants
    ]
    assert observed == [0.0, 1.0, 1.0, 2.0]
    for variant in registered_probe.variants:
        meta = routing[variant]
        assert meta["snyder_ct"] == pytest.approx(2.0)
        assert meta["rainfall_centroid_offset_hours"] == pytest.approx(12.0)
        assert meta["kernel_mode_from_storm_start_hours"] == pytest.approx(
            meta["lag_hours"] + 12.0
        )


def test_evaluated_models_are_not_judged_without_geometry_inputs(registered_probe):
    case = build_case(registered_probe, 20260914, "small")
    for model_name in EVALUATED_MODELS:
        model = registry.find_model(model_name)
        issues = compatibility_issues(model, registered_probe, case)
        message = "; ".join(issues)
        assert "main_channel_length_km" in message, model.name
        assert "centroid_channel_length_km" in message, model.name


def test_snyder_bounds_and_area_scaling_pass_for_physical_lags():
    probe = _probe()
    runs = _runs()
    bounds = get("lag_time_bounds")(
        runs,
        probe,
        {**COMMON, "lower_ratio": 0.5, "upper_ratio": 2.0,
         "discretization_tolerance_days": 0.5},
    )
    scaling = get("scaling_monotonicity")(
        runs,
        probe,
        {
            **COMMON,
            "reversal_tolerance_days": 0.5,
            "min_span_days": 2.0,
        },
    )
    assert bounds.passed, bounds.message
    assert scaling.passed, scaling.message
    expected_days = [
        bounds.diagnostics["variants"][name]["expected_lag_days"]
        for name in VARIANTS
    ]
    np.testing.assert_allclose(
        expected_days,
        [0.63110113, 1.12305225, 2.25004715, 3.33515212],
        rtol=1e-8,
    )
    assert scaling.diagnostics["increments_days"] == {
        "small->medium": 0.0,
        "medium->large": 1.0,
        "large->xlarge": 1.0,
    }
    assert scaling.diagnostics["span_days"] == 2.0


def test_instantaneous_runoff_fails_the_lower_lag_bound():
    result = get("lag_time_bounds")(
        _runs((0, 0, 0, 0)),
        _probe(),
        {**COMMON, "lower_ratio": 0.5, "upper_ratio": 2.0,
         "discretization_tolerance_days": 0.5},
    )
    assert not result.passed
    assert result.value > 0
    assert "outside" in result.message


def test_inverse_area_scaling_fails_monotonicity_without_failing_the_span():
    result = get("scaling_monotonicity")(
        _runs((1, 2, 1, 3)),
        _probe(),
        {
            **COMMON,
            "reversal_tolerance_days": 0.5,
            "min_span_days": 2.0,
        },
    )
    assert not result.passed
    assert result.value == pytest.approx(0.5)
    assert result.diagnostics["minimum_adjacent_increment_days"] == -1.0
    assert "medium to large" in result.message


def test_equal_runoff_maxima_use_their_temporal_centroid():
    runs = _runs()
    runs["small"] = _run("small", 1, tie=True)
    result = get("lag_time_bounds")(
        runs,
        _probe(),
        {**COMMON, "lower_ratio": 0.5, "upper_ratio": 2.0,
         "discretization_tolerance_days": 0.5},
    )
    small = result.diagnostics["variants"]["small"]
    assert small["tied_peak_steps"] == 2
    assert small["observed_lag_days"] == pytest.approx(1.5)


def test_response_baseline_uses_only_rows_immediately_before_the_event():
    runs = _runs()
    run = runs["small"]
    table = run.table.copy()
    table.loc[: EVENT.start - 6, "mrro"] = 100.0
    table.loc[EVENT.start - 5 : EVENT.start - 1, "mrro"] = 0.2
    runs["small"] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get("lag_time_bounds")(
        runs,
        _probe(),
        {
            **COMMON,
            "lower_ratio": 0.5,
            "upper_ratio": 2.0,
            "discretization_tolerance_days": 0.5,
        },
    )
    assert result.passed, result.message
    assert result.diagnostics["variants"]["small"][
        "baseline_runoff"
    ] == pytest.approx(0.2)


@pytest.mark.parametrize("criterion_name", ["lag_time_bounds", "scaling_monotonicity"])
def test_zero_response_fails_as_a_model_answer(criterion_name):
    runs = _runs()
    for name, run in list(runs.items()):
        table = run.table.copy()
        table["mrro"] = 0.2
        runs[name] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get(criterion_name)(runs, _probe(), dict(COMMON))
    assert not result.passed
    assert "no measurable storm response" in result.message
    assert result.diagnostics["response_fraction"] == pytest.approx(0.0)


def test_direct_criterion_call_defensively_fails_nonfinite_runoff():
    runs = _runs()
    table = runs["medium"].table.copy()
    table.loc[SPINUP_DAYS + 25, "mrro"] = np.nan
    run = runs["medium"]
    runs["medium"] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get("lag_time_bounds")(runs, _probe(), dict(COMMON))
    assert not result.passed
    assert result.diagnostics == {"variant": "medium", "nonfinite_count": 1}


def test_a_peak_on_the_final_row_has_no_measurable_lag():
    runs = _runs()
    run = runs["xlarge"]
    table = run.table.copy()
    table["mrro"] = 0.2
    table.loc[EVENT.start:, "mrro"] = np.linspace(
        0.21, 5.0, len(table) - EVENT.start
    )
    runs["xlarge"] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get("lag_time_bounds")(runs, _probe(), dict(COMMON))
    assert not result.passed
    assert "final row" in result.message


def test_forcing_must_be_exactly_the_same_between_catchments():
    runs = _runs()
    run = runs["large"]
    forcing = run.case.forcing.copy()
    forcing.loc[0, "tas"] += 0.01
    runs["large"] = RunResult(
        Case(
            probe_id=run.case.probe_id,
            seed=run.case.seed,
            forcing=forcing,
            static=run.case.static,
            spinup_steps=run.case.spinup_steps,
            timestep=run.case.timestep,
        ),
        run.table,
        run.meta,
        run.wall_seconds,
    )
    with pytest.raises(ValueError, match="exactly the same forcing"):
        get("lag_time_bounds")(runs, _probe(), dict(COMMON))


def test_event_annotation_must_be_contiguous_and_keep_response_tail():
    runs = _runs()
    for name, run in list(runs.items()):
        forcing = run.case.forcing.copy()
        forcing["pr"] = 0.0
        forcing["_event_pr"] = 0.0
        forcing.loc[len(forcing) - 3 :, ["pr", "_event_pr"]] = 10.0
        case = Case(
            probe_id=run.case.probe_id,
            seed=run.case.seed,
            forcing=forcing,
            static=run.case.static,
            spinup_steps=run.case.spinup_steps,
            timestep=run.case.timestep,
        )
        runs[name] = RunResult(case, run.table, run.meta, run.wall_seconds)
    with pytest.raises(ValueError, match="truncates the routing experiment"):
        get("lag_time_bounds")(runs, _probe(), dict(COMMON))


def test_unannotated_rain_in_the_scored_record_is_a_case_error():
    runs = _runs()
    for name, run in list(runs.items()):
        forcing = run.case.forcing.copy()
        forcing.loc[SPINUP_DAYS + 60, "pr"] = 1.0
        runs[name] = RunResult(
            Case(
                probe_id=run.case.probe_id,
                seed=run.case.seed,
                forcing=forcing,
                static=run.case.static,
                spinup_steps=run.case.spinup_steps,
                timestep=run.case.timestep,
            ),
            run.table,
            run.meta,
            run.wall_seconds,
        )
    with pytest.raises(ValueError, match="must be dry outside"):
        get("lag_time_bounds")(runs, _probe(), dict(COMMON))


def test_daily_quantisation_may_flatten_one_adjacent_lag_step():
    result = get("scaling_monotonicity")(
        _runs((0, 0, 1, 2)),
        _probe(),
        {
            **COMMON,
            "reversal_tolerance_days": 0.5,
            "min_span_days": 2.0,
        },
    )
    assert result.passed, result.message


def test_a_fixed_lag_fails_the_full_ladder_span():
    result = get("scaling_monotonicity")(
        _runs((1, 1, 1, 1)),
        _probe(),
        {
            **COMMON,
            "reversal_tolerance_days": 0.5,
            "min_span_days": 2.0,
        },
    )
    assert not result.passed
    assert result.value > 0.0
    assert result.threshold == 0.0
    assert result.diagnostics["span_shortfall_days"] == pytest.approx(
        result.value
    )
    assert "lag span" in result.message


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("main_channel_length_km", 0.0, "positive area"),
        ("centroid_channel_length_km", np.nan, "non-finite static field"),
        ("centroid_channel_length_km", 1000.0, "longer than its main channel"),
        ("area_km2", 20.0, "outside Snyder's published"),
    ],
)
def test_invalid_static_geometry_is_a_case_error(field, value, message):
    runs = _runs()
    run = runs["small"]
    static = dict(run.case.static)
    static[field] = value
    runs["small"] = RunResult(
        Case(
            probe_id=run.case.probe_id,
            seed=run.case.seed,
            forcing=run.case.forcing,
            static=static,
            spinup_steps=run.case.spinup_steps,
            timestep=run.case.timestep,
        ),
        run.table,
        run.meta,
        run.wall_seconds,
    )
    with pytest.raises(ValueError, match=message):
        get("lag_time_bounds")(runs, _probe(), dict(COMMON))

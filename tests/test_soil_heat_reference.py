"""Reference physics and the hourly experiment reaching the model intact."""

from __future__ import annotations

import importlib.util
import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from hydroturing import harness, registry
from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.harness import (
    build_case,
    load_generator,
    resolve_window_days,
    run_probe,
    select_window,
    window_case,
)
from hydroturing.protocol import FORCING_FILE, RunResult, stage
from hydroturing.runner import get_runner
from hydroturing.scoring import INCOMPATIBLE, NOT_SCORED


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("energy/soil-heat-storage-consistency")


@pytest.fixture(scope="module")
def simulate():
    path = registry.find_model("reference_soil_heat").path / "ht_adapter.py"
    spec = importlib.util.spec_from_file_location("soil_heat_reference", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.simulate


def test_generator_reproduces_seed_and_preserves_the_declared_experiment(probe):
    generate = load_generator(probe).generate
    first, static = generate(11)
    repeated, repeated_static = generate(11)
    other, other_static = generate(29)
    assert first.to_csv(index=False) == repeated.to_csv(index=False)
    assert static == repeated_static
    assert not first.rsds.equals(other.rsds)
    assert static != other_static
    assert len(first) == 96
    assert first._phase.tolist() == ["spinup"] * 24 + ["heating"] * 12 + ["recovery"] * 60
    assert first.rsds.iloc[:24].eq(0).all()
    assert first.rsds.iloc[24:36].gt(0).all()
    assert first.rsds.iloc[36:].eq(0).all()
    assert "rn" not in first
    assert {"rlds", "sfcWind", "huss", "ps"}.issubset(first.columns)
    assert first.pr.eq(0).all() and first.pet.eq(0).all()
    assert first.tas.gt(0).all()
    assert pd.to_datetime(first.time).diff().iloc[1:].eq(pd.Timedelta(hours=1)).all()
    assert 0.1 <= static["soil_layer_depth_m"] <= 0.35
    assert 1.8e6 <= static["soil_solid_heat_capacity"] <= 2.2e6
    assert 0 < static["soil_water_content_initial"] < static["soil_porosity"] < 1
    assert not {"soil_sensible_exchange", "soil_top_conductance", "soil_bottom_conductance"}.intersection(static)


def test_short_window_keeps_heating_recovery_and_staging_hides_labels(probe, tmp_path):
    model = registry.find_model("reference_soil_heat")
    submitted = replace(model, name="submitted_thermal_model", window_days=1)
    assert resolve_window_days(submitted, probe, override=1) == 3
    case = build_case(probe, 11)
    cut = window_case(case, select_window(case, probe, resolve_window_days(submitted, probe)))
    assert cut.n_steps == 96 and cut.spinup_steps == 24
    assert cut.window["rows"] == 72
    pd.testing.assert_frame_equal(cut.forcing, case.forcing)
    request = json.loads(stage(tmp_path, cut, probe, model).read_text())
    staged = pd.read_csv(tmp_path / FORCING_FILE)
    assert "_phase" not in staged
    assert request["n_steps"] == 96
    assert request["timestep"] == "PT1H"
    assert request["request"]["diagnostics"] == ["tsoil_layer"]
    assert request["request"]["states"] == []
    pd.testing.assert_frame_equal(staged, case.forcing.drop(columns="_phase"))


@pytest.mark.parametrize("seed", [11, 29])
def test_reference_equilibrium_heating_recovery_and_both_budgets(probe, simulate, seed):
    case = build_case(probe, seed)
    table = pd.DataFrame(simulate(case.forcing.to_dict("records"), case.static))
    initial = case.static["soil_temperature_initial"]
    assert np.max(np.abs(table.tsoil_layer.iloc[:24] - initial)) < 1e-10
    assert np.diff(table.tsoil_layer.iloc[23:36]).min() > 0
    assert np.diff(table.tsoil_layer.iloc[35:]).max() < 0
    assert table.tsoil_layer.iloc[-1] > initial
    residual = (table.hfg - table.hfg_bottom
                - case.static["soil_heat_capacity_areal"]
                * np.diff(np.r_[initial, table.tsoil_layer]) / 3600)
    assert np.max(np.abs(residual)) < 1e-8
    # Net radiation is calculated by the reference from incoming radiation
    # and its skin temperature; it is not prescribed in the shared forcing.
    assert np.max(np.abs(table.rn - table.hfss - table.hfls - table.hfg)) < 1e-10
    assert table.rn.iloc[24:36].gt(0).all()
    assert table.rn.iloc[36:].lt(0).all()
    run = RunResult(case, table, {}, 0.0)
    assert get("soil_heat_storage")(run, probe, {}).status == PASS


def test_radiative_reference_resolves_transient_at_sixty_seconds(probe, simulate):
    case = build_case(probe, 11)
    forcing = case.forcing.to_dict("records")
    regular = pd.DataFrame(simulate(forcing, case.static, substep_seconds=60.0))
    finer = pd.DataFrame(simulate(forcing, case.static, substep_seconds=30.0))
    assert np.max(np.abs(regular.tsoil_layer - finer.tsoil_layer)) < 1e-6
    for column in ("rn", "hfss", "hfg", "hfg_bottom"):
        assert np.max(np.abs(regular[column] - finer[column])) < 1e-5


def test_rounded_reference_keeps_margin_at_largest_declared_capacity(probe, simulate):
    case = build_case(probe, 11)
    # Largest layer and solid capacity in the declared dry-material range.
    porosity = case.static["soil_porosity"]
    water = case.static["soil_water_content_initial"]
    capacity = 0.35 * ((1 - porosity) * 2.2e6
                       + water * case.static["soil_water_heat_capacity"]
                       + (porosity - water) * case.static["soil_pore_air_heat_capacity"])
    case.static.update(soil_layer_depth_m=0.35, soil_solid_heat_capacity=2.2e6,
                       soil_heat_capacity_areal=capacity)
    table = pd.DataFrame(simulate(case.forcing.to_dict("records"), case.static))
    table["tsoil_layer"] = table.tsoil_layer.round(3)
    table[["hfg", "hfg_bottom"]] = table[["hfg", "hfg_bottom"]].round(2)
    result = get("soil_heat_storage")(RunResult(case, table, {}, 0.0), probe, {})
    assert result.status == PASS
    # A 0.001 K difference-rounding error and two 0.005 W m-2 flux errors.
    assert result.diagnostics["max_abs_residual_w_m2"] <= capacity * 0.001 / 3600 + 0.01


@pytest.mark.parametrize("name", ["reference_soil_heat", "reference_frozen_soil", "reference_half_soil"])
def test_prescribed_net_radiation_case_is_incompatible_before_adapter_execution(tmp_path, name, monkeypatch):
    # These references consume incoming radiation and a configured layer.
    # The older net-radiation surface case cannot supply those inputs.
    surface_probe = registry.find_probe("energy/surface-energy-closure")
    model = registry.find_model(name)
    monkeypatch.setattr(type(get_runner(model)), "run", lambda *args: pytest.fail("incompatible adapter ran"))
    outcome = run_probe(model, surface_probe, [11], workdir=tmp_path)
    assert (outcome.verdict, outcome.reason) == (NOT_SCORED, INCOMPATIBLE)
    assert any("rsds, rlds" in issue for issue in outcome.incompatible)
    assert any("soil_heat_capacity_areal" in issue for issue in outcome.incompatible)


@pytest.mark.parametrize("kind,key", [
    ("forcing", "rsds"), ("forcing", "rlds"),
    ("static", "soil_heat_capacity_areal"),
    ("static", "soil_layer_depth_m"),
    ("static", "soil_temperature_initial"),
])
def test_missing_case_input_is_incompatible_instead_of_adapter_error(probe, tmp_path, monkeypatch, kind, key):
    case = build_case(probe, 11)
    if kind == "forcing":
        case.forcing = case.forcing.drop(columns=key)
    else:
        del case.static[key]
    model = registry.find_model("reference_soil_heat")
    monkeypatch.setattr(harness, "build_case", lambda *args: case)
    monkeypatch.setattr(type(get_runner(model)), "run", lambda *args: pytest.fail("incompatible adapter ran"))
    outcome = run_probe(model, probe, [11], workdir=tmp_path)
    assert (outcome.verdict, outcome.reason) == (NOT_SCORED, INCOMPATIBLE)
    assert f"{kind} does not provide {key}" in outcome.incompatible


@pytest.mark.parametrize("kind,key", [
    ("forcing", "rsds"), ("forcing", "rlds"),
    ("static", "soil_heat_capacity_areal"),
    ("static", "soil_layer_depth_m"),
    ("static", "soil_temperature_initial"),
])
def test_outputs_alone_do_not_qualify_a_model_without_declared_case_consumption(probe, tmp_path, kind, key):
    model = registry.find_model("reference_soil_heat")
    model = replace(model, **{f"needs_{kind}": tuple(
        item for item in getattr(model, f"needs_{kind}") if item != key
    )})
    assert model.missing_for(probe) == []
    outcome = run_probe(model, probe, [11], workdir=tmp_path)
    assert (outcome.verdict, outcome.reason) == (NOT_SCORED, INCOMPATIBLE)
    assert f"model does not declare that it consumes {kind} {key}" in outcome.incompatible


def test_configured_depth_and_capacity_define_the_actual_lumped_layer(probe, tmp_path):
    model = registry.find_model("reference_soil_heat")
    case = build_case(probe, 11)
    original = get_runner(model).run(model, probe, case, tmp_path / "original")
    # Depth defines the reporting control volume. At unchanged areal
    # capacity a deeper homogeneous layer has lower Cv, not extra storage.
    deeper = replace(case, static=dict(case.static, soil_layer_depth_m=2 * case.static["soil_layer_depth_m"]))
    same_capacity = get_runner(model).run(model, probe, deeper, tmp_path / "same-capacity")
    pd.testing.assert_frame_equal(original.table, same_capacity.table)
    layer = same_capacity.meta["thermal_layer"]
    assert layer["bottom_depth_m"] == 2 * original.meta["thermal_layer"]["bottom_depth_m"]
    assert layer["heat_capacity_volumetric_j_m3_k"] == 0.5 * original.meta["thermal_layer"]["heat_capacity_volumetric_j_m3_k"]
    # With the same material Cv, doubling depth doubles C_A. The actual
    # simulated response must then change, while the layer budget closes.
    deeper.static["soil_heat_capacity_areal"] *= 2
    same_material = get_runner(model).run(model, probe, deeper, tmp_path / "same-material")
    assert same_material.table.tsoil_layer.iloc[35] < original.table.tsoil_layer.iloc[35]
    assert get("soil_heat_storage")(same_material, probe, {}).status == PASS
    assert same_material.meta["thermal_layer"]["heat_capacity_volumetric_j_m3_k"] == original.meta["thermal_layer"]["heat_capacity_volumetric_j_m3_k"]


def test_initial_temperature_is_used_even_when_it_differs_from_air(probe, tmp_path):
    case = build_case(probe, 11)
    air = case.forcing.tas.iloc[0] + 273.15
    case.static["soil_temperature_initial"] = air + 3.0
    model = registry.find_model("reference_soil_heat")
    run = get_runner(model).run(model, probe, case, tmp_path)
    initial = case.static["soil_temperature_initial"]
    # A warmer initialized layer cools toward the air, rather than silently
    # starting at air temperature. The first interval uses that initial T.
    assert air < run.table.tsoil_layer.iloc[0] < initial
    storage = case.static["soil_heat_capacity_areal"] * (run.table.tsoil_layer.iloc[0] - initial) / 3600
    assert run.table.hfg.iloc[0] - run.table.hfg_bottom.iloc[0] == pytest.approx(storage, abs=1e-8)


@pytest.mark.parametrize("name,scale", [
    ("reference_soil_heat", 1.0),
    ("reference_frozen_soil", 0.0),
    ("reference_half_soil", 0.5),
])
def test_actual_adapters_keep_fluxes_while_controls_change_temperature(probe, simulate, tmp_path, name, scale):
    case = build_case(probe, 11)
    model = registry.find_model(name)
    run = get_runner(model).run(model, probe, case, tmp_path)
    correct = pd.DataFrame(simulate(case.forcing.to_dict("records"), case.static))
    pd.testing.assert_frame_equal(run.table.drop(columns="tsoil_layer"), correct.drop(columns="tsoil_layer"))
    initial = case.static["soil_temperature_initial"]
    np.testing.assert_allclose(run.table.tsoil_layer, initial + scale * (correct.tsoil_layer - initial), atol=1e-12)
    assert run.meta["model"]["name"] == name
    assert run.meta["radiation_input"] == "incoming rsds and rlds"
    assert run.meta["thermal_layer"] == {
        "top_depth_m": 0.0,
        "bottom_depth_m": case.static["soil_layer_depth_m"],
        "heat_capacity_areal_j_m2_k": case.static["soil_heat_capacity_areal"],
        "heat_capacity_volumetric_j_m3_k": case.static["soil_heat_capacity_areal"] / case.static["soil_layer_depth_m"],
        "initial_temperature_k": initial,
    }
    assert run.meta["time_convention"]["time"] == "interval start"
    assert "interval-end" in run.meta["time_convention"]["tsoil_layer"]
    # Use the reference's own net-radiation output for this supplemental
    # surface identity. The shared storage case supplies incoming SW/LW,
    # not the prescribed Rn expected by the old surface-energy criteria.
    assert np.max(np.abs(run.table.rn - run.table.hfss - run.table.hfls - run.table.hfg)) < 1e-10
    assert get("soil_heat_storage")(run, probe, {}).status == (PASS if scale == 1 else FAIL)

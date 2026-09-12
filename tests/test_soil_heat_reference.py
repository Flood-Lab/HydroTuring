"""Reference physics and the hourly experiment reaching the model intact."""

from __future__ import annotations

import importlib.util
import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.harness import (
    build_case,
    load_generator,
    resolve_window_days,
    select_window,
    window_case,
)
from hydroturing.protocol import FORCING_FILE, RunResult, stage
from hydroturing.runner import get_runner


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


def test_existing_net_radiation_fluxes_match_integrated_native_conduction(probe, simulate):
    case = build_case(probe, 11)
    # This is the older prescribed-net-radiation input format, still used by
    # the existing surface probes. These coefficients belong to the model.
    static = dict(case.static, soil_sensible_exchange=12.0,
                  soil_top_conductance=6.0, soil_bottom_conductance=1.0)
    forcing = case.forcing.assign(rn=case.forcing.rsds).drop(columns=["rsds", "rlds"])
    output = simulate(forcing.to_dict("records"), static)[24]
    a, ks, kb = 12.0, 6.0, 1.0
    # Resolve the transient at 1-second intervals and integrate the physical
    # flux laws numerically, rather than use the adapter's analytic mean.
    rn = forcing.rn.iloc[24]
    exchange = a * ks / (a + ks)
    equilibrium = ks * rn / ((a + ks) * (exchange + kb))
    seconds = np.arange(3601)
    layer = equilibrium * (1 - np.exp(-(exchange + kb) * seconds / static["soil_heat_capacity_areal"]))
    skin = (rn + ks * layer) / (a + ks)
    for column, values in (("hfg", ks * (skin - layer)), ("hfg_bottom", kb * layer),
                           ("hfss", a * skin)):
        mean = (0.5 * (values[:-1] + values[1:]) * np.diff(seconds)).sum() / 3600
        assert output[column] == pytest.approx(mean, abs=1e-6)
    assert output["tsoil_layer"] == pytest.approx(static["soil_temperature_initial"] + layer[-1])
    assert abs(output["hfg_bottom"] - kb * layer[-1]) > 0.1


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
def test_reference_runs_existing_surface_case_without_thermal_settings(tmp_path, name):
    # `ht run --model reference_soil_heat` also selects this existing probe.
    # Its forcing specifies no thermal layer, so fixed model defaults apply.
    surface_probe = registry.find_probe("energy/surface-energy-closure")
    case = build_case(surface_probe, 11)
    assert "soil_heat_capacity_areal" not in case.static
    model = registry.find_model(name)
    run = get_runner(model).run(model, surface_probe, case, tmp_path)
    assert len(run.table) == case.n_steps
    assert np.isfinite(run.table.tsoil_layer).all()
    assert get("energy_closure")(run, surface_probe, {}).status == PASS
    assert get("energy_closure_by_phase")(run, surface_probe, {}).status == PASS


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
    # Use the reference's own net-radiation output for this supplemental
    # surface identity. The shared storage case supplies incoming SW/LW,
    # not the prescribed Rn expected by the old surface-energy criteria.
    assert np.max(np.abs(run.table.rn - run.table.hfss - run.table.hfls - run.table.hfg)) < 1e-10
    assert get("soil_heat_storage")(run, probe, {}).status == (PASS if scale == 1 else FAIL)

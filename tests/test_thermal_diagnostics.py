"""Temperature crosses the adapter contract without becoming water storage."""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest
import yaml

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.criteria.base import make_window, reported_states
from hydroturing.harness import run_probe, verify_adapter_contract
from hydroturing.protocol import Case, ProtocolError, RunResult, read_result, stage
from hydroturing.scoring import INCOMPLETE, NOT_SCORED
from hydroturing.spec import SpecError, load_model, load_probe


@pytest.fixture
def mass_probe():
    return registry.find_probe("mass/catchment-closure")


@pytest.fixture
def thermal_probe(mass_probe):
    return replace(
        mass_probe,
        requires_fluxes=("hfg", "hfg_bottom"),
        requires_states=(),
        requires_diagnostics=("tsoil_layer",),
    )


@pytest.fixture
def case():
    return Case(
        probe_id="energy/soil-heat-storage-consistency",
        seed=11,
        forcing=pd.DataFrame({
            "time": pd.date_range("2000-01-01", periods=3, freq="D"),
            "pr": [0.0, 1.0, 1.0],
        }),
        static={},
        spinup_steps=1,
    )


def _write_result(path, case, **columns):
    output = path / "output"
    output.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame({"time": case.forcing["time"], **columns})
    table.to_csv(output / "result.csv", index=False)
    (output / "run.json").write_text(json.dumps({"status": "ok"}))


def test_optional_diagnostic_manifests_load_and_old_defaults_remain_empty(mass_probe, tmp_path):
    model = registry.find_model("reference_bucket")
    assert model.emits_diagnostics == ()
    assert mass_probe.requires_diagnostics == ()

    probe_dir = tmp_path / "mass" / "catchment-closure"
    probe_dir.mkdir(parents=True)
    raw = yaml.safe_load((mass_probe.path / "probe.yaml").read_text())
    raw["requires"]["diagnostics"] = ["tsoil_layer"]
    (probe_dir / "generate.py").write_text("# Only manifest loading is tested.\n")
    (probe_dir / "probe.yaml").write_text(yaml.safe_dump(raw))
    loaded_probe = load_probe(probe_dir)

    model_dir = tmp_path / "thermal_model"
    model_dir.mkdir()
    (model_dir / "model.yaml").write_text(yaml.safe_dump({
        "name": "thermal_model", "version": "1", "entrypoint": ["python", "adapter.py"],
        "timestep": "PT1D",
        "emits": {"fluxes": ["hfg", "hfg_bottom"], "states": [],
                  "diagnostics": ["tsoil_layer"]},
    }))
    loaded_model = load_model(model_dir)
    assert loaded_probe.requires_diagnostics == loaded_model.emits_diagnostics == ("tsoil_layer",)
    assert "tsoil_layer" in loaded_model.emitted
    assert "tsoil_layer" not in loaded_model.emits_states


def test_temperature_cannot_be_declared_as_a_water_requirement(mass_probe, tmp_path):
    probe_dir = tmp_path / "mass" / "catchment-closure"
    probe_dir.mkdir(parents=True)
    raw = yaml.safe_load((mass_probe.path / "probe.yaml").read_text())
    raw["requires"]["states"].append("tsoil_layer")
    (probe_dir / "generate.py").write_text("# Only manifest loading is tested.\n")
    (probe_dir / "probe.yaml").write_text(yaml.safe_dump(raw))
    with pytest.raises(SpecError, match="requires.diagnostics"):
        load_probe(probe_dir)


def test_diagnostic_request_and_result_roundtrip(case, thermal_probe, mass_probe, tmp_path):
    model = replace(
        registry.find_model("reference_bucket"),
        emits_fluxes=("hfg", "hfg_bottom"), emits_states=(),
        emits_diagnostics=("tsoil_layer",),
    )
    # The model sees its declared outputs even when the current probe only
    # needs water. Changing the probe must not reveal which check will run.
    request = json.loads(stage(tmp_path, case, mass_probe, model).read_text())
    assert request["request"]["diagnostics"] == ["tsoil_layer"]
    assert request["request"]["states"] == []
    assert request["units"]["tsoil_layer"] == "K"
    assert request["units"]["hfg_bottom"] == "W m-2"

    _write_result(tmp_path, case, hfg=[0, 1, 1], hfg_bottom=[0, 0, 0],
                  tsoil_layer=[280.0, 280.1, 280.2])
    result = read_result(tmp_path, case, thermal_probe, 0.0)
    window = make_window(result, thermal_probe)
    assert window.state0["tsoil_layer"] == 280.0
    assert list(window.table["tsoil_layer"]) == [280.1, 280.2]


def test_unreported_required_diagnostic_is_incomplete_before_execution(thermal_probe):
    model = replace(registry.find_model("reference_bucket"),
                    emits_fluxes=("hfg", "hfg_bottom"), emits_states=())
    outcome = run_probe(model, thermal_probe, [11])
    assert outcome.verdict == NOT_SCORED
    assert outcome.reason == INCOMPLETE
    assert outcome.missing == ["tsoil_layer"]


@pytest.mark.parametrize("temperature", [None, [280, np.nan, 281], [280, np.inf, 281],
                                         [280, "not-a-temperature", 281]])
def test_missing_or_nonfinite_required_temperature_is_rejected(case, thermal_probe,
                                                              tmp_path, temperature):
    columns = {"hfg": [0, 1, 1], "hfg_bottom": [0, 0, 0]}
    if temperature is not None:
        columns["tsoil_layer"] = temperature
    _write_result(tmp_path, case, **columns)
    with pytest.raises(ProtocolError, match="tsoil_layer"):
        read_result(tmp_path, case, thermal_probe, 0.0)


def test_layer_warming_does_not_change_water_closure(case, mass_probe):
    table = pd.DataFrame({
        "time": case.forcing["time"], "pr": case.forcing["pr"],
        "evspsbl": 0.0, "mrro": 0.0,
        "mrso": [0.0, 1.0, 2.0], "snw": 0.0, "canopy": 0.0,
        "tsoil_layer": [280.0, 290.0, 300.0],
    })
    probe = replace(mass_probe, requires_diagnostics=("tsoil_layer",))
    run = RunResult(case, table, {}, 0.0)
    window = make_window(run, probe)
    assert "tsoil_layer" not in reported_states(window, probe)
    np.testing.assert_array_equal(window.storage(reported_states(window, probe)), [1.0, 2.0])
    params = next(c.params for c in probe.criteria if c.name == "closure")
    result = get("closure")(run, probe, params)
    assert result.passed, result.message
    assert result.value == 0.0


def test_verify_adapter_checks_emitted_diagnostics_even_on_water_probe(mass_probe, tmp_path):
    # reference_bucket does not write temperature; declaring it makes the
    # adapter contract false, even though the selected water probe needs none.
    model = replace(registry.find_model("reference_bucket"),
                    emits_diagnostics=("tsoil_layer",))
    with pytest.raises(ProtocolError, match="tsoil_layer"):
        verify_adapter_contract(model, mass_probe, 11, workdir=tmp_path)


def test_verify_adapter_does_not_require_an_undeclared_diagnostic(thermal_probe, tmp_path):
    model = registry.find_model("reference_bucket")
    result = verify_adapter_contract(model, thermal_probe, 11, workdir=tmp_path)
    assert "tsoil_layer" not in result.table
    assert len(result.table) == result.case.n_steps

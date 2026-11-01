"""Independent native diagnostics must validate assumptions, not only a ledger."""

from __future__ import annotations

import importlib.util

import numpy as np
import pandas as pd
import pytest

from hydroturing.spec import REPO_ROOT
from test_aquifer_equilibrium import STATIC
from hydroturing import registry
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult


def _audit():
    path = REPO_ROOT / "scripts/validate_aquifer_equilibrium.py"
    assert path.exists(), "the reproducible native-head audit is missing"
    spec = importlib.util.spec_from_file_location("validate_equilibrium", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _inputs():
    # A one-cell literal case: b = 2 mm/m; head gains 0.1 and 0.2 m.
    forcing = pd.DataFrame({"sw_stage_m": [96.0, 96.0]})
    table = pd.DataFrame({"gw": [3000.2, 3000.4]})
    heads = np.array([[[[95.1]]], [[[95.2]]]])
    sy = np.zeros((2, 1, 1, 1))
    return forcing, table, heads, sy


def test_native_heads_independently_match_confined_storage_change():
    result = _audit().native_metrics(*_inputs(), STATIC)
    assert result["minimum_head_m"] == 95.1
    assert result["minimum_head_above_top_m"] == pytest.approx(5.1)
    assert result["minimum_riverbed_clearance_m"] == pytest.approx(3.1)
    assert result["storage_head_mapping_max_error_mm"] < 1.0e-9
    assert result["sto_sy_gross_mm"] == 0.0


@pytest.mark.parametrize("fault", ["unconfined", "disconnected", "sy", "mapping", "nonfinite", "overflow"])
def test_native_audit_rejects_invalid_physical_evidence(fault):
    forcing, table, heads, sy = _inputs()
    static = dict(STATIC)
    if fault == "unconfined":
        heads[0, 0, 0, 0] = 89.0
    elif fault == "disconnected":
        static["river_bottom_offset_m"] = 0.1
    elif fault == "sy":
        sy[0, 0, 0, 0] = 0.01
    elif fault == "mapping":
        table.loc[0, "gw"] += 0.05
    elif fault == "nonfinite":
        heads[0, 0, 0, 0] = np.nan
    else:
        sy[:] = 1.0e308
    with pytest.raises(ValueError):
        _audit().native_metrics(forcing, table, heads, sy, static)


def test_native_audit_rejects_a_head_outside_prescribed_envelope():
    forcing, table, heads, sy = _inputs()
    heads[1, 0, 0, 0] = 96.2
    table.loc[1, "gw"] = 3002.4
    with pytest.raises(ValueError, match="head envelope"):
        _audit().native_metrics(forcing, table, heads, sy, STATIC)


def test_diagnostic_capture_must_preserve_standard_result_bytes(tmp_path):
    first = tmp_path / "first.csv"
    second = tmp_path / "second.csv"
    first.write_bytes(b"gw\n3000.2\n")
    second.write_bytes(first.read_bytes())
    assert _audit().compare_standard_output(first, second)
    second.write_bytes(b"gw\n3000.3\n")
    with pytest.raises(ValueError, match="changes"):
        _audit().compare_standard_output(first, second)


def test_default_audit_really_launches_an_independent_standard_run(tmp_path, monkeypatch):
    module = _audit()
    case = build_case(registry.find_probe(module.PROBE_ID), 7, "control")
    initial = (case.static["aquifer_specific_yield"]
               * (case.static["aquifer_initial_head_m"] - case.static["aquifer_bottom_m"]) * 1000)
    n = len(case.forcing)
    table = pd.DataFrame({
        "time": case.forcing["time"], "gw": np.full(n, initial),
        "gw_sw_exchange": np.zeros(n), "gw_to_sw": np.zeros(n), "sw_to_gw": np.zeros(n),
    })

    class Container:
        def __init__(self):
            self.invocations = []

        def run(self, model, probe, generated_case, directory):
            self.invocations.append(directory)
            (directory / "output").mkdir(parents=True, exist_ok=True)
            table.to_csv(directory / "output/result.csv", index=False)
            return RunResult(generated_case, table, {}, 0.0, model=model)

    # Only the external container and native binary I/O are replaced. The
    # audit's directory selection, comparison and physical checks are real.
    heads = np.full((n, 1, 1, 1), case.static["aquifer_initial_head_m"])
    monkeypatch.setattr(module, "read_native", lambda directory: (heads, np.zeros_like(heads)))
    monkeypatch.setattr(module, "read_native_inputs", lambda directory, forcing, static: {})
    runners = {"capture": Container(), "standard": Container()}
    result = module.native_case(7, "control", "standard", tmp_path, None, runners)
    assert result["standard_output_byte_identical"]
    assert len(runners["standard"].invocations) == 1
    assert runners["capture"].invocations[0] != runners["standard"].invocations[0]


def test_diagnostic_runner_extracts_native_files_before_contract_validation(tmp_path, monkeypatch):
    module = _audit()
    probe = registry.find_probe(module.PROBE_ID)
    model = registry.find_model("modflow6")
    case = build_case(probe, 7, "control")
    initial = (case.static["aquifer_specific_yield"]
               * (case.static["aquifer_initial_head_m"] - case.static["aquifer_bottom_m"]) * 1000)

    def native_invocation(self, model, probe, io_dir, request_path):
        output = io_dir / "output"
        (output / "native").mkdir()
        (output / "native/ht_gw.hds").write_bytes(b"native-head-evidence")
        pd.DataFrame({
            "time": case.forcing["time"], "gw": initial,
            "gw_sw_exchange": 0.0, "gw_to_sw": 0.0, "sw_to_gw": 0.0,
        }).to_csv(output / "result.csv", index=False)
        (output / "run.json").write_text('{"status": "ok"}')

    monkeypatch.setattr(module.DockerRunner, "invoke", native_invocation)
    captured = tmp_path / "captured-standard"
    run = module.NativeCaptureRunner().run(model, probe, case, captured)
    assert len(run.table) == len(case.forcing)
    assert not (captured / "output/native").exists()
    assert (tmp_path / "captured-standard-native/ht_gw.hds").read_bytes() == b"native-head-evidence"


def test_native_envelope_uses_verified_serialized_boundary_precision():
    forcing, table, heads, sy = _inputs()
    forcing["sw_stage_m"] = 95.0
    heads[:] = 95.0000002
    table["gw"] = 3000.0000004
    result = _audit().native_metrics(
        forcing, table, heads, sy, STATIC,
        native_initial_head=95.0,
        native_stage=np.full(2, 95.0000004),
        native_bed=np.full(2, 91.0000004),
    )
    assert result["native_head_envelope_excess_m"] == 0.0
    assert result["requested_head_envelope_excess_m"] == pytest.approx(2e-7)
    assert result["input_head_roundoff_m"] == pytest.approx(4e-7)


def test_native_input_reader_checks_actual_values_and_rejects_changed_physics(tmp_path):
    module = _audit()
    forcing = pd.DataFrame({"sw_stage_m": [195.123456789, 196.123456789]})
    static = dict(STATIC, aquifer_initial_head_m=195.123456789)
    (tmp_path / "ht_gw.ic").write_text('BEGIN griddata\nstrt\nCONSTANT 195.12345679\nEND griddata\n')
    (tmp_path / "ht_gw.sto").write_text('BEGIN options\nSTORAGECOEFFICIENT\nEND options\nBEGIN griddata\nss\nCONSTANT 0.00200000\nEND griddata\n')
    river = '\n'.join(f'BEGIN period {i+1}\n1 6 6 {stage:.8E} 200.0 {stage-4:.8E}\nEND period {i+1}'
                      for i, stage in enumerate(forcing["sw_stage_m"]))
    (tmp_path / "ht_gw.riv").write_text(river)
    result = module.read_native_inputs(tmp_path, forcing, static)
    assert result["native_stage"][0] == 195.123457
    (tmp_path / "ht_gw.riv").write_text(river.replace('1.95123457E+02', '1.97123457E+02'))
    with pytest.raises(ValueError, match="serialized"):
        module.read_native_inputs(tmp_path, forcing, static)


def test_native_storage_mapping_uses_verified_serialized_storativity():
    forcing, table, heads, sy = _inputs()
    heads[:] = 95.9
    table["gw"] = 3001.8000036
    result = _audit().native_metrics(
        forcing, table, heads, sy, STATIC, native_storativity=0.002000004,
    )
    assert result["storage_head_mapping_max_error_mm"] < 1e-9
    assert result["requested_storage_head_mapping_max_error_mm"] == pytest.approx(3.6e-6)

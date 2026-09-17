"""Counterexamples for soil storage accounting and its time/phase boundaries."""

from __future__ import annotations

import json
from dataclasses import asdict

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.harness import build_case
from hydroturing.protocol import Case, RunResult, read_result
from hydroturing.spec import TIMESTEP_DAYS


def build(net, error=0.0, phases=None, spinup=1, timestep="PT1H", bottom=7.0):
    net = np.asarray(net, dtype=float)
    n = len(net)
    error = np.broadcast_to(error, (n,))
    dt = TIMESTEP_DAYS[timestep] * 86400
    capacity, initial = 3600.0, 290.0
    temperatures = initial + np.cumsum((net - error) * dt / capacity)
    times = pd.date_range("2001-07-01", periods=n, freq=pd.Timedelta(seconds=dt))
    if phases is None:
        phases = ["heating"] * n
    case = Case(
        probe_id="energy/soil-heat-storage-consistency", seed=1,
        forcing=pd.DataFrame({"time": times, "_phase": phases}),
        static={"soil_heat_capacity_areal": capacity, "soil_temperature_initial": initial},
        spinup_steps=spinup, timestep=timestep,
    )
    return RunResult(
        case, pd.DataFrame({"time": times, "hfg": net + bottom,
                            "hfg_bottom": bottom, "tsoil_layer": temperatures}),
        {}, 0.0,
    )


def score(run, **params):
    return get("soil_heat_storage")(run, None, params)


def generated_csv_run(tmp_path, top=0.0, bottom=0.0, temperature_jump=0.0):
    """Exercise the generated case and normal result reader, including spinup."""
    probe = registry.find_probe("energy/soil-heat-storage-consistency")
    case = build_case(probe, 11)
    temperatures = np.full(case.n_steps, case.static["soil_temperature_initial"])
    temperatures[case.spinup_steps:] += temperature_jump
    table = pd.DataFrame({
        "time": case.forcing["time"], "hfg": top, "hfg_bottom": bottom,
        "tsoil_layer": temperatures,
    })
    output = tmp_path / "output"
    output.mkdir()
    table.to_csv(output / "result.csv", index=False)
    (output / "run.json").write_text(json.dumps({"status": "ok", "n_steps": case.n_steps}))
    return read_result(tmp_path, case, probe, 0.0), probe


@pytest.mark.parametrize("timestep", ["PT1H", "PT15M"])
def test_correct_heating_and_cooling_use_case_dt_and_previous_endpoint(timestep):
    run = build([30, 20, 10, -10, -5], phases=["spinup", "heating", "heating", "recovery", "recovery"], timestep=timestep)
    result = score(run)
    assert result.status == PASS
    assert result.diagnostics["max_abs_residual_w_m2"] < 1e-10
    assert [(b["start"], b["stop"]) for b in result.diagnostics["blocks"]] == [(0, 2), (2, 4)]


def test_no_spinup_uses_initial_temperature_not_first_output():
    run = build([10, 10], spinup=0)
    assert score(run).status == PASS
    run.table.loc[0, "tsoil_layer"] = run.case.static["soil_temperature_initial"]
    assert score(run).status == FAIL


def test_last_spinup_temperature_is_used_even_when_static_initial_differs():
    run = build([100, 10, 10])
    assert run.table.loc[0, "tsoil_layer"] != run.case.static["soil_temperature_initial"]
    assert score(run).status == PASS
    run.table.loc[0, "tsoil_layer"] += 30
    assert score(run).status == FAIL


def test_opposite_storage_errors_cannot_cancel():
    run = build([0, 20, 20, 20, 20], error=[0, 4, -4, 4, -4])
    result = score(run)
    assert result.status == FAIL
    assert result.diagnostics["worst_block"]["mean_abs_w_m2"] == pytest.approx(4)


def test_bad_heating_not_diluted_by_long_recovery():
    run = build([0, 20, 20] + [0] * 60, error=[0, 4, 4] + [0] * 60,
                phases=["spinup", "heating", "heating"] + ["recovery"] * 60)
    result = score(run)
    assert result.status == FAIL
    assert result.diagnostics["failed_blocks"] == 1
    assert result.diagnostics["blocks"][1]["passed"]


@pytest.mark.parametrize("net,error,expected", [
    (100, 5, PASS), (100, 6, FAIL), (10, 1, PASS), (10, 2, FAIL),
    (0, 0, PASS), (-100, 5, PASS), (-100, 6, FAIL),
])
def test_phase_allowance_uses_absolute_net_input_and_floor(net, error, expected):
    assert score(build([0, net, net], error=[0, error, error])).status == expected


def test_large_throughflow_does_not_inflate_net_input_denominator():
    result = score(build([0, 10, 10], error=[0, 2, 2], bottom=990))
    assert result.status == FAIL
    assert result.diagnostics["worst_block"]["allowance_w_m2"] == 1


@pytest.mark.parametrize("column,row", [("tsoil_layer", 0), ("tsoil_layer", 2), ("hfg_bottom", 2)])
def test_nonfinite_spinup_endpoint_or_scored_values_fail(column, row):
    run = build([0, 20, 20])
    run.table.loc[row, column] = np.nan
    assert score(run).status == FAIL


def test_missing_heat_capacity_or_no_spinup_initial_is_not_assumed():
    run = build([10, 10], spinup=0)
    del run.case.static["soil_temperature_initial"]
    with pytest.raises(ValueError, match="soil_temperature_initial"):
        score(run)
    del run.case.static["soil_heat_capacity_areal"]
    with pytest.raises(ValueError, match="soil_heat_capacity_areal"):
        score(run)


@pytest.mark.parametrize("top,bottom,temperature_jump,stage", [
    pytest.param(1e308, -1e308, 0.0, "derived", id="boundary-subtraction"),
    pytest.param(1e308, 0.0, 0.0, "phase", id="phase-reduction"),
    pytest.param(0.0, 0.0, 1e308, "derived", id="storage-product"),
])
def test_finite_csv_values_cannot_pass_after_calculation_overflow(
    tmp_path, top, bottom, temperature_jump, stage,
):
    run, probe = generated_csv_run(tmp_path, top, bottom, temperature_jump)
    assert len(run.table) == 96 and run.case.spinup_steps == 24
    assert np.isfinite(run.table[["hfg", "hfg_bottom", "tsoil_layer"]]).all().all()
    result = get("soil_heat_storage")(run, probe, {})
    assert result.status == FAIL
    assert "non-finite" in result.message
    assert result.diagnostics["non_finite_stage"] == stage
    assert result.value is None
    # The diagnostic itself must remain serializable without NaN or Infinity.
    json.dumps(asdict(result), allow_nan=False)


@pytest.mark.parametrize("throughflow", [0.0, 25.0])
def test_inert_case_is_identified_without_confusing_steady_throughflow(tmp_path, throughflow):
    run, probe = generated_csv_run(tmp_path, top=throughflow, bottom=throughflow)
    result = get("soil_heat_storage")(run, probe, {})
    assert result.status == PASS  # Budget consistency alone is a necessary condition.
    assert result.diagnostics["zero_temperature_change"] is True
    assert result.diagnostics["zero_boundary_fluxes"] is (throughflow == 0)
    assert ("no thermal response exercised" in result.message) is (throughflow == 0)
    assert all(block["mean_abs_w_m2"] == 0 for block in result.diagnostics["blocks"])

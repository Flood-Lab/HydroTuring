"""Counterexamples for soil storage accounting and its time/phase boundaries."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.protocol import Case, RunResult
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

"""Hand-built checks of the radiation identity, tolerance and invalid inputs."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.protocol import Case, RunResult

# CODATA 2018, typed here rather than imported so that the module under test
# cannot vouch for its own constant.
SIGMA = 5.670374419e-8


def upward(ts, rlds, eps):
    """The exact identity, the only physics these tests need."""
    return eps * SIGMA * np.asarray(ts, dtype=float) ** 4 + (1.0 - eps) * np.asarray(rlds, dtype=float)


def build(ts, rlus=None, rlds=300.0, eps=0.98, spinup=0, static=None):
    """A run reporting `ts` and `rlus` under a sky of `rlds`; exact unless told otherwise."""
    ts = np.asarray(ts, dtype=float)
    n = len(ts)
    rlds = np.broadcast_to(np.asarray(rlds, dtype=float), (n,))
    if rlus is None:
        rlus = upward(ts, rlds, eps)
    times = pd.date_range("2001-07-01T06:00", periods=n, freq="h")
    forcing = pd.DataFrame({"time": times, "rlds": rlds})
    table = pd.DataFrame({"time": times, "ts": ts, "rlus": np.asarray(rlus, dtype=float)})
    case = Case(
        probe_id="t", seed=1, forcing=forcing,
        static={"eps": eps} if static is None else static,
        spinup_steps=spinup, timestep="PT1H",
    )
    return RunResult(case=case, table=table, meta={}, wall_seconds=0.0)


def score(run, **params):
    return get("radiative_identity")(run, None, params)


def diurnal(n=48, mean=290.0, amplitude=8.0):
    hours = 6.0 + np.arange(n)
    return mean + amplitude * np.cos(2.0 * np.pi * (hours - 14.0) / 24.0)


def test_exact_identity_passes_at_every_step():
    run = build(diurnal(), rlds=300.0 + 40.0 * np.sin(np.arange(48) / 5.0))
    result = score(run)
    assert result.status == PASS
    assert result.diagnostics["violating_steps"] == 0
    assert result.diagnostics["scored_steps"] == 48
    assert result.diagnostics["max_abs_residual_w_m2"] < 1e-9
    assert result.value < 1e-9


@pytest.mark.parametrize("ts,error,expected,floor_steps", [
    # At 290 K the exact flux is 399.0 W m-2, so 0.5 percent is about 2.0.
    (290.0, 1.9, PASS, 0),
    (290.0, -1.9, PASS, 0),
    (290.0, 2.1, FAIL, 0),
    (290.0, -2.1, FAIL, 0),
    # At 180 K the flux is about 64 W m-2 and 0.5 percent of it is below the
    # 0.5 W m-2 floor, so the floor is the bound on every step.
    (180.0, 0.45, PASS, 24),
    (180.0, 0.55, FAIL, 24),
])
def test_tolerance_is_relative_to_the_reported_flux_with_a_floor(ts, error, expected, floor_steps):
    ts = np.full(24, ts)
    run = build(ts, rlus=upward(ts, 300.0, 0.98) + error)
    result = score(run)
    assert result.status == expected
    assert result.diagnostics["floor_steps"] == floor_steps


def test_the_relative_bound_is_taken_against_the_reported_flux():
    ts = np.full(12, 290.0)
    exact = upward(ts, 300.0, 0.98)
    # Over-reporting by 0.5 percent of the expected value passes, because the
    # bound is 0.5 percent of the larger, reported, flux.
    assert score(build(ts, rlus=exact * 1.005)).status == PASS
    # Under-reporting by the same amount fails: the bound shrank with the flux.
    assert score(build(ts, rlus=exact * 0.995)).status == FAIL


@pytest.mark.parametrize("eps", [0.95, 0.98, 0.99])
def test_reporting_emission_alone_as_upward_longwave_fails(eps):
    """Omit reflected sky across the probe's emissivity range. At eps 0.99,
    the fixed sky contributes 3 W m-2 of missing reflection."""
    ts = diurnal()
    run = build(ts, rlus=eps * SIGMA * ts**4, rlds=300.0, eps=eps)
    result = score(run)
    assert result.status == FAIL
    assert result.diagnostics["violating_steps"] == 48
    assert result.diagnostics["worst_step"]["residual_w_m2"] == pytest.approx(-(1.0 - eps) * 300.0)


@pytest.mark.parametrize("offset_k,expected", [
    # Around 290 K under a 300 W m-2 sky, an offset below about 0.4 K sits
    # inside 0.5 percent of the flux: the criterion's resolution under these
    # conditions, stated rather than hidden. It is coarser at a warmer surface.
    (0.3, PASS),
    (1.0, FAIL),
    (5.0, FAIL),
    (-5.0, FAIL),
])
def test_emitting_at_a_different_temperature_than_reported_fails(offset_k, expected):
    """The proposal's first negative control: rlus from the air, ts from the surface."""
    ts = diurnal()
    run = build(ts, rlus=upward(ts + offset_k, 300.0, 0.98))
    result = score(run)
    assert result.status == expected
    if expected == FAIL:
        assert result.diagnostics["violating_steps"] > 0


@pytest.mark.parametrize("amplitude_k,expected", [(5.0, PASS), (10.0, FAIL)])
def test_a_linearised_stefan_boltzmann_law_is_seen_only_on_large_excursions(amplitude_k, expected):
    """The proposal's optional third control, reported rather than gated: the
    second-order term eps*sigma*6*T0**2*dT**2 is 0.7 W m-2 at 5 K and 2.8 at
    10 K around 288 K, against a bound near 2 W m-2."""
    t0, eps = 288.0, 0.98
    ts = diurnal(mean=t0, amplitude=amplitude_k)
    linear = eps * SIGMA * (t0**4 + 4.0 * t0**3 * (ts - t0)) + (1.0 - eps) * 300.0
    result = score(build(ts, rlus=linear, eps=eps))
    assert result.status == expected
    second_order = eps * SIGMA * 6.0 * t0**2 * amplitude_k**2
    assert result.diagnostics["max_abs_residual_w_m2"] == pytest.approx(second_order, rel=0.05)


def test_spinup_is_not_scored_and_indices_start_at_the_window():
    ts = diurnal(n=72)
    rlus = upward(ts, 300.0, 0.98)
    rlus[:24] += 50.0
    rlus[24 + 7] += 3.0
    result = score(build(ts, rlus=rlus, spinup=24))
    assert result.status == FAIL
    assert result.diagnostics["scored_steps"] == 48
    assert result.diagnostics["violating_steps"] == 1
    assert result.diagnostics["worst_step"]["index"] == 7
    assert result.diagnostics["worst_step"]["time"].startswith("2001-07-02 13:00")
    assert result.diagnostics["worst_step"]["residual_w_m2"] == pytest.approx(3.0)


@pytest.mark.parametrize("column", ["ts", "rlus"])
@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf])
def test_non_finite_model_output_is_a_failure(column, invalid):
    run = build(diurnal(n=24))
    run.table.loc[12, column] = invalid
    result = score(run)
    assert result.status == FAIL
    assert result.diagnostics["non_finite_steps"] == 1


def test_finite_temperature_that_overflows_the_identity_fails():
    run = build([1e100], rlus=[500.0])
    result = score(run)
    assert result.status == FAIL
    assert result.diagnostics["non_finite_calculation_steps"] == 1


def test_large_finite_residuals_keep_the_report_finite():
    run = build([290.0, 290.0], rlus=[1e308, 1e308])
    with np.errstate(over="raise", invalid="raise"):
        result = score(run)
    assert result.status == FAIL
    assert result.diagnostics["mean_abs_residual_w_m2"] == pytest.approx(1e308)
    json.dumps(result.diagnostics, allow_nan=False)


@pytest.mark.parametrize("temperature", [0.0, -290.0])
def test_a_temperature_at_or_below_zero_kelvin_is_named_as_a_unit_mistake(temperature):
    # Even a matching fourth-power emission must not validate a non-positive K.
    ts = np.full(24, temperature)
    result = score(build(ts))
    assert result.status == FAIL
    assert "kelvin" in result.message
    assert result.diagnostics["non_positive_kelvin_steps"] == 24


def test_a_warm_surface_reported_in_celsius_fails_on_the_residual():
    ts = np.full(24, 17.0)
    result = score(build(ts, rlus=upward(ts + 273.15, 300.0, 0.98)))
    assert result.status == FAIL
    assert result.diagnostics["violating_steps"] == 24
    assert "17.00 K" in result.message


@pytest.mark.parametrize("setting", [
    {"sigma": 0.0}, {"sigma": -1.0}, {"sigma": np.nan}, {"sigma": np.inf},
    {"rel_tol": -0.01}, {"rel_tol": np.nan}, {"rel_tol": np.inf},
    {"rel_tol": 1e308}, {"sigma": 1e308, "rel_tol": 1e308},
    {"abs_floor": 0.0}, {"abs_floor": -1.0}, {"abs_floor": np.nan}, {"abs_floor": np.inf},
], ids=lambda s: "=".join(f"{k}{v}" for k, v in s.items()))
def test_an_invalid_setting_is_an_error_rather_than_a_pass(setting):
    # About 100 W m-2 of residual, which any sane setting rejects. A NaN or an
    # infinite bound would count zero violations and pass it.
    ts = np.full(24, 290.0)
    run = build(ts, rlus=np.full(24, 500.0))
    assert score(run).status == FAIL
    with pytest.raises(ValueError, match="finite"):
        score(run, **setting)


def test_missing_inputs_and_a_bad_emissivity_are_errors_rather_than_verdicts():
    run = build(diurnal(n=24))

    run.case.forcing = run.case.forcing.drop(columns="rlds")
    with pytest.raises(ValueError, match="forcing column 'rlds'"):
        score(run)

    run = build(diurnal(n=24), static={})
    with pytest.raises(ValueError, match="'eps' in the case's static"):
        score(run)

    for bad in (0.0, 1.2, -0.5, np.nan, np.inf):
        with pytest.raises(ValueError, match="outside"):
            score(build(diurnal(n=24), static={"eps": bad}))

    run = build(diurnal(n=24))
    run.case.forcing.loc[3, "rlds"] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        score(run)

    run = build(diurnal(n=24))
    run.table = run.table.drop(columns="ts")
    with pytest.raises(ValueError, match="needs 'ts'"):
        score(run)


def test_column_names_emissivity_key_and_constants_are_configurable():
    ts = diurnal(n=24)
    sigma = 5.67e-8
    run = build(ts, rlus=0.97 * sigma * ts**4 + 0.03 * 280.0, rlds=280.0, static={"emissivity": 0.97})
    run.case.forcing = run.case.forcing.rename(columns={"rlds": "lw_down"})
    run.table = run.table.rename(columns={"ts": "skin", "rlus": "lw_up"})
    result = score(
        run, upward="lw_up", temperature="skin", downward="lw_down",
        emissivity="emissivity", sigma=sigma,
    )
    assert result.status == PASS
    assert result.diagnostics["emissivity"] == 0.97
    assert result.diagnostics["sigma_w_m2_k4"] == sigma
    # The two constants differ by 6.6e-5 relative, far inside the tolerance,
    # so only an exact residual shows the parameter was used in the arithmetic.
    assert result.diagnostics["max_abs_residual_w_m2"] < 1e-9

    # Loosening either bound changes the verdict where that bound governs,
    # and the report says which values were in force.
    warm = build(ts, rlus=upward(ts, 300.0, 0.98) + 10.0)
    assert score(warm).status == FAIL
    loose = score(warm, rel_tol=0.05)
    assert loose.status == PASS
    assert loose.diagnostics["rel_tol"] == 0.05
    cold = np.full(24, 180.0)
    off = build(cold, rlus=upward(cold, 300.0, 0.98) + 3.0)
    assert score(off).status == FAIL
    lifted = score(off, abs_floor=5.0)
    assert lifted.status == PASS
    assert lifted.diagnostics["abs_floor_w_m2"] == 5.0
    assert lifted.diagnostics["floor_steps"] == 24

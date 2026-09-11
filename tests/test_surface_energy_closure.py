"""Counterexamples for temporal surface energy closure, independent of the gate."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.protocol import Case, RunResult
from hydroturing.spec import TIMESTEP_DAYS


def build(residual, rn=100.0, phases=None, spinup=0, timestep="PT1H"):
    """Give H a known accounting error; no model physics is needed to hide it."""
    residual = np.asarray(residual, dtype=float)
    n = len(residual)
    rn = np.broadcast_to(rn, (n,))
    if phases is None:
        phases = np.where(np.arange(n) % 24 < 12, "day", "night")
    times = pd.date_range(
        "2001-07-01T06:00", periods=n,
        freq=pd.Timedelta(days=TIMESTEP_DAYS[timestep]),
    )
    forcing = pd.DataFrame({"time": times, "rn": rn, "_regime": phases})
    table = pd.DataFrame({
        "time": times, "hfls": 30.0, "hfss": rn - 40.0 - residual, "hfg": 10.0,
    })
    case = Case(
        probe_id="t", seed=1, forcing=forcing, static={},
        spinup_steps=spinup, timestep=timestep,
    )
    return RunResult(case=case, table=table, meta={}, wall_seconds=0.0)


def score(run, **params):
    return get("energy_closure_by_phase")(run, None, params)


@pytest.mark.parametrize("residual", [
    np.repeat([20.0, -20.0], 12),
    np.tile([20.0, -20.0], 12),
], ids=["day-night-cancellation", "within-phase-cancellation"])
def test_opposite_errors_that_pass_cumulative_closure_fail(residual):
    run = build(residual, rn=np.repeat([100.0, -20.0], 12))
    assert get("energy_closure")(run, None, {}).status == PASS
    result = score(run)
    assert result.status == FAIL
    assert result.diagnostics["failed_blocks"] == 2
    assert result.diagnostics["worst_block"]["mean_abs_w_m2"] == pytest.approx(20.0)


@pytest.mark.parametrize("rn,error,expected", [
    (100.0, 4.0, PASS),
    (100.0, 5.0, PASS),
    (100.0, 6.0, FAIL),
    (2.0, 1.0, PASS),
    (2.0, 2.0, PASS),
    (2.0, 3.0, FAIL),
    (0.0, 1.0, PASS),
], ids=[
    "relative-pass", "relative-boundary", "relative-fail",
    "floor-pass", "floor-boundary", "floor-fail", "zero-radiation",
])
def test_relative_allowance_and_absolute_floor(rn, error, expected):
    result = score(build(np.full(24, error), rn=rn))
    assert result.status == expected, result.message


def test_exact_closure_allows_negative_turbulent_fluxes():
    run = build(np.zeros(24), rn=np.repeat([100.0, -20.0], 12))
    run.table["hfls"] = np.repeat([-10.0, 10.0], 12)
    run.table["hfss"] = run.case.forcing["rn"] - run.table["hfls"] - run.table["hfg"]
    assert (run.table["hfls"] < 0).any() and (run.table["hfss"] < 0).any()
    assert score(run).status == PASS


@pytest.mark.parametrize("amplitude_k", [0.5, 2.0])
def test_plate_depth_flux_needs_storage_correction(amplitude_k):
    # Known heat content of a uniform 5 cm layer above a plate. Temperature
    # is its layer mean at interval endpoints, with a peak at 14:00.
    run = build(np.zeros(48), rn=np.tile(np.repeat([100.0, -30.0], 12), 2))
    endpoint_hours = 6.0 + np.arange(49)
    temperature = 20.0 + amplitude_k * np.cos(
        2.0 * np.pi * (endpoint_hours - 14.0) / 24.0
    )
    heat_content = 2.0e6 * 0.05 * temperature  # J m-2
    storage_rate = np.diff(heat_content) / 3600.0  # W m-2
    surface_flux = run.table["hfg"].copy()
    run.table["hfg"] = surface_flux - storage_rate

    # Daily storage changes cancel, hiding the wrong boundary cumulatively.
    assert get("energy_closure")(run, None, {}).status == PASS
    uncorrected = score(run)
    assert uncorrected.status == FAIL
    assert "uncorrected deeper-boundary flux" in uncorrected.message

    # The adapter supplies the independently known layer storage, not a
    # correction diagnosed from the surface-budget residual.
    run.table["hfg"] += storage_rate
    assert score(run).status == PASS


def test_spinup_is_not_scored_and_block_indices_start_at_the_scored_window():
    run = build(np.r_[np.full(24, 500.0), np.zeros(24)], spinup=24)
    result = score(run)
    assert result.status == PASS
    blocks = result.diagnostics["blocks"]
    assert [(b["start"], b["stop"]) for b in blocks] == [(0, 12), (12, 24)]


def test_a_bad_day_is_not_diluted_by_other_days_with_the_same_label():
    # The first day has 8 W m-2 of error; the next has none. Grouping all
    # 'day' rows would produce 4 W m-2 and incorrectly pass the 5 W m-2 bound.
    run = build(np.r_[np.full(12, 8.0), np.zeros(36)])
    result = score(run)
    assert result.status == FAIL
    assert result.diagnostics["failed_blocks"] == 1
    assert len(result.diagnostics["blocks"]) == 4
    assert result.diagnostics["worst_block"]["start"] == 0


def test_duration_uses_the_case_timestep_and_configured_segment_column():
    run = build(np.full(8, 3.0), phases=["day"] * 8, timestep="PT15M")
    run.case.forcing = run.case.forcing.rename(columns={"_regime": "_phase"})
    result = score(run, segment_column="_phase", threshold=0.01, floor=4.0)
    assert result.status == PASS
    block = result.diagnostics["blocks"][0]
    assert block["duration_hours"] == pytest.approx(2.0)
    assert block["mean_abs_w_m2"] == pytest.approx(3.0)
    assert block["allowance_w_m2"] == pytest.approx(4.0)


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf])
def test_non_finite_model_flux_is_a_failure(invalid):
    run = build(np.zeros(24))
    run.table.loc[12, "hfss"] = invalid
    result = score(run)
    assert result.status == FAIL
    assert result.diagnostics["non_finite_steps"] == 1

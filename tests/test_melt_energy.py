"""Counterexamples for the melt-fusion identity, independent of the gate.

The gate shows that `melt_energy` separates the reference models. These show
*why*: what it charges for, what it refuses to charge for, and where its
boundaries sit, on tables built by hand with no model physics behind them.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.protocol import Case, RunResult

LAMBDA_F = 3.337e5
SECONDS_PER_DAY = 86400.0

RN, HFLS, HFG = 80.0, 20.0, 8.0


def demand(melt_mm, n_melt=60, dt_days=1.0):
    """The mean flux a given melt over a given block demands, W m-2."""
    return LAMBDA_F * melt_mm / (n_melt * dt_days * SECONDS_PER_DAY)


def build(
    *, melt_mm=600.0, residual=None, peak=600.0, snowfall=600.0,
    n_accum=100, n_melt=60, sbl_rate=0.0, spinup=0,
):
    """A pack that builds and then goes, with the energy budget dialled by hand.

    `residual` is what `rn - hfls - hfss - hfg` comes to over the melt block;
    it is imposed through the sensible flux, because no model physics is
    needed to state that a budget did or did not pay for something.
    """
    if residual is None:
        residual = demand(melt_mm, n_melt)
    n = n_accum + n_melt

    sbl_total = sbl_rate * n_melt
    loss = melt_mm + sbl_total
    snw = np.concatenate([
        np.linspace(peak / n_accum, peak, n_accum),
        np.linspace(peak - loss / n_melt, peak - loss, n_melt),
    ])
    pr = np.concatenate([
        np.full(n_accum, snowfall / n_accum), np.zeros(n_melt),
    ])
    sbl = np.concatenate([np.zeros(n_accum), np.full(n_melt, sbl_rate)])
    regime = np.array(["accumulation"] * n_accum + ["melt"] * n_melt, dtype=object)

    # Zero residual outside the melt block, `residual` inside it.
    imposed = np.where(regime == "melt", residual, 0.0)
    times = pd.date_range("2000-10-01", periods=n, freq="D")

    forcing = pd.DataFrame({"time": times, "pr": pr, "rn": RN, "_regime": regime})
    table = pd.DataFrame({
        "time": times, "snw": snw, "sbl": sbl,
        "hfls": HFLS, "hfg": HFG, "hfss": RN - HFLS - HFG - imposed,
    })
    case = Case(
        probe_id="t", seed=1, forcing=forcing, static={},
        spinup_steps=spinup, timestep="PT1D",
    )
    return RunResult(case=case, table=table, meta={}, wall_seconds=0.0)


def score(run, **params):
    return get("melt_energy")(run, None, params)


def test_a_budget_that_pays_for_its_melt_passes():
    assert score(build()).status == PASS


def test_a_degree_day_model_that_pays_nothing_fails():
    """The probe's reason for existing: both budgets close, fusion is free."""
    result = score(build(residual=0.0))
    assert result.status == FAIL
    assert result.diagnostics["blocks"][0]["demanded_w_m2"] == pytest.approx(38.62, abs=0.01)
    assert result.diagnostics["blocks"][0]["available_w_m2"] == pytest.approx(0.0)


@pytest.mark.parametrize("share,expected", [
    (0.00, FAIL), (0.94, FAIL), (0.95, PASS), (1.00, PASS),
    (1.05, PASS), (1.06, FAIL), (1.20, FAIL),
], ids=["pays-nothing", "-6%", "-5%", "on-target", "+5%", "+6%", "+20%"])
def test_the_tolerance_is_where_it_says_it_is(share, expected):
    """600 mm over 60 days demands 38.62 W m-2, so 5 percent of it is 1.93.

    The 2 W m-2 floor is the larger of the two here and is therefore what
    binds, giving a band of +/-5.2 percent. That is the normal case for this
    probe rather than a corner: the floor is the operative bound for any melt
    demanding under 40 W m-2, and the relative rule takes over above it.
    """
    assert score(build(residual=share * demand(600.0))).status == expected


def test_the_floor_binds_when_there_is_almost_nothing_to_melt():
    """A thin pack demands little, and a percentage of little is not a bound."""
    thin = build(melt_mm=5.0, peak=5.0, snowfall=8.0, residual=1.9)
    assert demand(5.0) < 1.0
    assert score(thin).status == PASS
    assert score(build(melt_mm=5.0, peak=5.0, snowfall=8.0, residual=4.0)).status == FAIL


def test_a_model_reporting_no_pack_fails_rather_than_passing_on_two_zeroes():
    run = build(melt_mm=0.0, peak=0.0, snowfall=600.0, residual=0.0)
    result = score(run)
    assert result.status == FAIL
    assert "no pack to melt" in result.message


def test_a_pack_that_only_partly_melts_passes_if_it_paid_for_what_it_melted():
    """Melting slowly is a calibration choice, not a conservation violation.

    This criterion is entitled to judge whether the two ledgers agree about
    the ice that moved, and nothing else. A model that melts a tenth of its
    pack and pays exactly the fusion energy for that tenth is coupled
    correctly; failing it for the other nine tenths would be scoring a
    parameter. Catching a model that had the energy and declined to melt
    needs the surface temperature, which the contract does not carry.
    """
    result = score(build(melt_mm=60.0, peak=600.0, residual=demand(60.0)))
    assert result.status == PASS
    assert result.diagnostics["blocks"][0]["share_of_peak_left"] == pytest.approx(0.9)


def test_a_pack_that_only_partly_melts_still_fails_if_it_paid_nothing():
    """The identity is what does the work once the melt-out guard is gone."""
    assert score(build(melt_mm=60.0, peak=600.0, residual=0.0)).status == FAIL


def test_sublimation_is_not_charged_as_melt():
    """Water that left as vapour never underwent fusion and must not be billed.

    The pack loses 700 mm, 100 of it to sublimation. A criterion reading the
    whole loss as melt would demand the energy for 700 and fail a budget that
    correctly paid for 600.
    """
    run = build(melt_mm=600.0, sbl_rate=100.0 / 60.0, residual=demand(600.0))
    result = score(run)
    assert result.status == PASS
    assert result.diagnostics["blocks"][0]["melt_mm"] == pytest.approx(600.0, abs=1e-6)


def test_melt_is_measured_from_the_step_before_the_block():
    """The classic off-by-one: starting at the block's own first row loses a day.

    With the window opening exactly on the melt block, the pack it is compared
    against is the one carried in from spinup, not the block's first row.
    """
    run = build(spinup=100)
    result = score(run)
    assert result.status == PASS
    assert result.diagnostics["blocks"][0]["melt_mm"] == pytest.approx(600.0, abs=1e-6)


def test_an_unlabelled_record_is_an_error_not_a_pass():
    run = build()
    run.case.forcing = run.case.forcing.drop(columns=["_regime"])
    with pytest.raises(ValueError, match="_regime"):
        score(run)


def test_a_nan_in_sublimation_is_reported_as_a_contract_violation():
    """Not as "melt of nan mm demands nan W/m2", which reads as bad physics."""
    run = build()
    run.table.loc[len(run.table) - 5, "sbl"] = float("nan")
    result = score(run)
    assert result.status == FAIL
    assert "non-finite" in result.message and "sublimation" in result.message


def test_a_scored_block_that_is_not_dry_is_refused():
    """Melt is inferred from the pack, which only holds where nothing is added.

    A probe reusing this criterion over a block with snowfall in it would have
    that snowfall netted out of the pack change, shrinking the melt silently.
    """
    run = build()
    run.case.forcing.loc[len(run.case.forcing) - 3, "pr"] = 4.0
    with pytest.raises(ValueError, match="has to be dry"):
        score(run)


def test_paying_for_melt_that_never_happened_does_not_read_as_short():
    """The failing branch fires both ways; the message must say which way."""
    over = score(build(residual=2.0 * demand(600.0)))
    assert over.status == FAIL
    assert "over by" in over.message and "short by" not in over.message
    under = score(build(residual=0.0))
    assert "short by" in under.message

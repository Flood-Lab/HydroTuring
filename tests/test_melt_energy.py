"""Counterexamples for the melt-fusion identity, independent of the gate.

The gate shows that `melt_energy` separates the reference models. These show
*why*: what it charges for, what it refuses to charge for, and where its
boundaries sit, on tables built by hand with no model physics behind them.

The balance under test has two terms, because a snowpack can absorb energy in
two ways and only one of them is a phase change:

    mean(rn - hfls - hfss - hfg) == [ lambda_f * M + dU ] / (N * dt * 86400)

    M  = -d(snw - lwsnl) - sum(sbl dt)   the ice that changed phase
    dU = -d(csnow)                       what warming a cold pack cost
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.protocol import Case, RunResult

LAMBDA_F = 3.337e5
C_ICE = 2100.0
T0 = 273.15
SECONDS_PER_DAY = 86400.0

RN, HFLS, HFG = 80.0, 20.0, 8.0


def demand(melt_mm, n_melt=60, d_internal_j=0.0, dt_days=1.0):
    """The mean flux a given melt and warming demand over a block, W m-2."""
    return (LAMBDA_F * melt_mm + d_internal_j) / (n_melt * dt_days * SECONDS_PER_DAY)


def build(
    *, melt_mm=600.0, residual=None, peak=600.0, snowfall=600.0,
    n_accum=100, n_melt=60, sbl_rate=0.0, spinup=0,
    liquid_mm=0.0, cold_start_c=0.0, cold_end_c=0.0, air_c=-20.0,
):
    """A pack that builds and then goes, with every term dialled by hand.

    `liquid_mm` is water held in the pore space and counted inside `snw`, so a
    criterion reading `snw` alone would mistake it for ice. `cold_start_c` and
    `cold_end_c` say how far below freezing the pack sits at the two ends of
    the scored block; `csnow` is built from them as `c_ice * ice * dT`, which
    is the identity the contract asks an adapter to satisfy.
    """
    n = n_accum + n_melt
    sbl_total = sbl_rate * n_melt
    ice_loss = melt_mm + sbl_total

    ice = np.concatenate([
        np.linspace(peak / n_accum, peak, n_accum),
        np.linspace(peak - ice_loss / n_melt, peak - ice_loss, n_melt),
    ])
    # Ramped rather than stepped, and taken out of the ice below, because in a
    # real pack the retained liquid *is* melted ice: stepping it at the block
    # boundary would invent pack water in a block that carries no precipitation.
    liquid = np.concatenate([np.zeros(n_accum), np.linspace(liquid_mm / n_melt,
                                                            liquid_mm, n_melt)])
    below = np.full(n, cold_start_c)
    below[n_accum:] = np.linspace(cold_start_c, cold_end_c, n_melt)
    csnow = C_ICE * ice * np.maximum(below, 0.0)

    if residual is None:
        residual = demand(melt_mm, n_melt, -(csnow[-1] - csnow[n_accum - 1]))

    regime = np.array(["accumulation"] * n_accum + ["melt"] * n_melt, dtype=object)
    imposed = np.where(regime == "melt", residual, 0.0)
    times = pd.date_range("2000-10-01", periods=n, freq="D")

    forcing = pd.DataFrame({
        "time": times,
        "pr": np.concatenate([np.full(n_accum, snowfall / n_accum), np.zeros(n_melt)]),
        "rn": RN, "tas": air_c, "_regime": regime,
    })
    table = pd.DataFrame({
        "time": times, "snw": ice + liquid, "lwsnl": liquid, "csnow": csnow,
        "sbl": np.concatenate([np.zeros(n_accum), np.full(n_melt, sbl_rate)]),
        "hfls": HFLS, "hfg": HFG, "hfss": RN - HFLS - HFG - imposed,
    })
    case = Case(probe_id="t", seed=1, forcing=forcing, static={},
                spinup_steps=spinup, timestep="PT1D")
    return RunResult(case=case, table=table, meta={}, wall_seconds=0.0)


def score(run, **params):
    return get("melt_energy")(run, None, params)


# --- the identity ----------------------------------------------------------

def test_a_budget_that_pays_for_its_melt_passes():
    assert score(build()).status == PASS


def test_a_degree_day_model_that_pays_nothing_fails():
    """The probe's reason for existing: both budgets close, fusion is free."""
    result = score(build(residual=0.0))
    assert result.status == FAIL
    assert result.diagnostics["blocks"][0]["demanded_w_m2"] == pytest.approx(38.62, abs=0.01)


@pytest.mark.parametrize("share,expected", [
    (0.00, FAIL), (0.94, FAIL), (0.95, PASS), (1.00, PASS),
    (1.05, PASS), (1.06, FAIL), (1.20, FAIL),
], ids=["pays-nothing", "-6%", "-5%", "on-target", "+5%", "+6%", "+20%"])
def test_the_tolerance_is_where_it_says_it_is(share, expected):
    """600 mm over 60 days demands 38.62 W m-2, so 5 percent of it is 1.93.

    The 2 W m-2 floor is the larger of the two here and is therefore what
    binds, giving a band of +/-5.2 percent. That is the normal case for this
    probe rather than a corner: the floor is the operative bound for any
    demand under 40 W m-2, and the relative rule takes over above it.
    """
    assert score(build(residual=share * demand(600.0))).status == expected


def test_the_floor_binds_when_there_is_almost_nothing_to_melt():
    """A thin pack demands little, and a percentage of little is not a bound."""
    assert demand(5.0) < 1.0
    assert score(build(melt_mm=5.0, peak=5.0, snowfall=8.0, residual=1.9)).status == PASS
    assert score(build(melt_mm=5.0, peak=5.0, snowfall=8.0, residual=4.0)).status == FAIL


# --- what `snw` alone cannot tell you --------------------------------------

def test_liquid_held_in_the_pack_is_not_mistaken_for_ice():
    """`snw` is total pack water in Snow-17 and SUMMA, ice plus liquid.

    Here 100 mm of the melt stays in the pore space, so `snw` falls by 100 mm
    less than the ice did. A criterion differencing `snw` would charge 500 mm
    of fusion where 600 mm happened, and fail a budget that paid correctly.
    """
    run = build(melt_mm=600.0, liquid_mm=100.0)
    result = score(run)
    assert result.status == PASS
    assert result.diagnostics["blocks"][0]["melt_mm"] == pytest.approx(600.0, abs=1e-6)
    snw_only = -(float(run.table["snw"].iloc[-1]) - float(run.table["snw"].iloc[99]))
    assert snw_only == pytest.approx(500.0, abs=1e-6)


def test_sublimation_is_not_charged_as_melt():
    """Water that left as vapour never underwent fusion and must not be billed.

    700 mm leaves a 700 mm pack, 100 of it as vapour, so the fusion is 600.
    """
    result = score(build(melt_mm=600.0, peak=700.0, snowfall=700.0,
                         sbl_rate=100.0 / 60.0))
    assert result.status == PASS
    assert result.diagnostics["blocks"][0]["melt_mm"] == pytest.approx(600.0, abs=1e-6)


# --- cold content ----------------------------------------------------------

def test_warming_a_cold_pack_is_charged_as_well_as_melting_it():
    """The term that makes ripeness measurable instead of assumed."""
    run = build(melt_mm=400.0, peak=600.0, cold_start_c=6.0, cold_end_c=0.0)
    b = score(run).diagnostics["blocks"][0]
    assert b["cold_content_w_m2"] > 1.0, "the case does not exercise cold content"
    assert b["demanded_w_m2"] == pytest.approx(
        b["fusion_w_m2"] + b["cold_content_w_m2"], abs=1e-9
    )


def test_a_pack_warmed_for_nothing_fails_by_exactly_the_warming():
    """Melt paid for, the climb from -12 C to 0 C not: the gap is the warming.

    -12 C rather than -6 C because 600 mm of ice warming 6 K demands 1.46 W m-2,
    under the 2 W m-2 floor: the probe would pass such a model and the test
    would prove nothing. At 12 K it is 2.92 W m-2 and the term has to be paid.
    """
    run = build(melt_mm=400.0, peak=600.0, cold_start_c=12.0, cold_end_c=0.0)
    fusion_only = score(run).diagnostics["blocks"][0]["fusion_w_m2"]
    result = score(build(melt_mm=400.0, peak=600.0, cold_start_c=12.0,
                         cold_end_c=0.0, residual=fusion_only))
    assert result.status == FAIL
    b = result.diagnostics["blocks"][0]
    assert b["gap_w_m2"] == pytest.approx(b["cold_content_w_m2"], abs=1e-6)


# --- guards and contract ---------------------------------------------------

def test_a_model_reporting_no_pack_fails_rather_than_passing_on_two_zeroes():
    result = score(build(melt_mm=0.0, peak=0.0, snowfall=600.0, residual=0.0))
    assert result.status == FAIL
    assert "no pack to melt" in result.message


def test_a_pack_that_only_partly_melts_passes_if_it_paid_for_what_it_melted():
    """Melting slowly is a calibration choice, not a conservation violation."""
    assert score(build(melt_mm=60.0, peak=600.0)).status == PASS


def test_a_pack_that_only_partly_melts_still_fails_if_it_paid_nothing():
    assert score(build(melt_mm=60.0, peak=600.0, residual=0.0)).status == FAIL


def test_melt_is_measured_from_the_step_before_the_block():
    """The classic off-by-one: starting at the block's own first row loses a day."""
    result = score(build(spinup=100))
    assert result.status == PASS
    assert result.diagnostics["blocks"][0]["melt_mm"] == pytest.approx(600.0, abs=1e-6)


def test_an_unlabelled_record_is_an_error_not_a_pass():
    run = build()
    run.case.forcing = run.case.forcing.drop(columns=["_regime"])
    with pytest.raises(ValueError, match="_regime"):
        score(run)


def test_a_nan_in_sublimation_is_reported_as_a_contract_violation():
    run = build()
    run.table.loc[len(run.table) - 5, "sbl"] = float("nan")
    result = score(run)
    assert result.status == FAIL
    assert "non-finite" in result.message and "sublimation" in result.message


def test_a_scored_block_that_is_not_dry_is_refused():
    """Ice arriving during the block cannot be told from ice melting."""
    run = build()
    run.case.forcing.loc[len(run.case.forcing) - 3, "pr"] = 4.0
    with pytest.raises(ValueError, match="has to be dry"):
        score(run)


def test_paying_for_melt_that_never_happened_does_not_read_as_short():
    over = score(build(residual=2.0 * demand(600.0)))
    assert over.status == FAIL
    assert "over by" in over.message and "short by" not in over.message
    assert "short by" in score(build(residual=0.0)).message


def test_an_overflowing_calculation_cannot_report_pass():
    """Finite inputs, infinite arithmetic, and `inf <= inf` is True."""
    huge = 1e308
    run = build()
    n = len(run.table)
    run.table["snw"] = np.concatenate([
        np.linspace(huge / 100, huge, 100), np.linspace(huge, 0.0, n - 100),
    ])
    run.case.forcing.loc[:99, "pr"] = huge / 100
    result = score(run)
    assert result.status == FAIL
    assert "overflowed" in result.message


# --- the case's own invariants ---------------------------------------------

_GEN = Path(__file__).resolve().parents[1] / "probes/energy/snowmelt-energy-water/generate.py"
_spec = importlib.util.spec_from_file_location("snowmelt_generate", _GEN)
snowmelt_generate = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(snowmelt_generate)

SWEEP = range(0, 8000, 23)          # 348 seeds


def test_a_negative_stored_quantity_is_named_as_a_contract_error():
    """The one-row route: negative liquid on the step before the block.

    A model that pays for 600 mm and melts 300 mm could otherwise report
    `lwsnl = -300` on that row and have the difference read as ice that left.
    """
    for column in ("snw", "lwsnl", "csnow"):
        run = build()
        run.table.loc[run.table.index[99], column] = -300.0
        result = score(run)
        assert result.status == FAIL, column
        assert "contract:" in result.message and column in result.message


def test_cold_content_without_ice_is_named_as_a_contract_error():
    """Cold content is the energy to bring ice to 0 C, and there is no ice."""
    run = build(melt_mm=600.0, peak=600.0)
    run.table.loc[run.table.index[-1], "snw"] = 0.0
    run.table.loc[run.table.index[-1], "lwsnl"] = 0.0
    run.table.loc[run.table.index[-1], "csnow"] = 5.0e6
    result = score(run)
    assert result.status == FAIL
    assert "contract:" in result.message and "no ice" in result.message


CAP_J = C_ICE * 600.0 * 30.0     # 600 mm of ice, 10 K below air held at -20 C


@pytest.mark.parametrize("factor,refused", [
    (0.98, False), (1.02, True),
], ids=["just-under-the-cap", "just-over-the-cap"])
def test_the_cold_content_cap_binds_where_it_says_it_does(factor, refused):
    """Bounded by what the reported ice could hold ten kelvin below the coldest
    air in the record, which is a number the case supplies and a model cannot
    move. Raising the opening row is the route that buys room on the balance.

    The assertion is on whether the *contract* check fires, not on the verdict:
    just under the cap the declaration is allowed through and then costs the
    model on the identity, because inflating the opening cold content raises
    the demand by at least as much as it buys. That is the point of the bound
    rather than a side effect of it.
    """
    run = build(melt_mm=400.0, peak=600.0, air_c=-20.0)
    run.table.loc[run.table.index[99], "csnow"] = factor * CAP_J
    message = score(run).message
    assert ("physically plausible on the rows" in message) is refused


def test_the_cap_is_not_applied_to_rows_the_identity_never_reads():
    """Snow-17 bounds its deficit as a mass fraction, not as a temperature.

    `NEGHS <= 0.33 * WE` is worth 52 K, and on a thin pack early in
    accumulation that is legitimately colder than a temperature cap allows.
    Those rows are never read by the balance, so applying the bound to them
    would refuse an honest Snow-17-type pack for being colder than the air.
    """
    run = build(melt_mm=400.0, peak=600.0, air_c=-20.0)
    run.table.loc[run.table.index[5], "csnow"] = 50.0 * CAP_J   # mid-accumulation
    assert score(run).status == PASS


def test_the_opening_ice_guard_looks_at_the_step_before_the_block():
    """A degree-day model that reports its pack as liquid on that one row.

    Measured on the peak anywhere in the window this passes, because the pack
    was real earlier; measured where the block opens, it is caught.
    """
    run = build()
    run.table.loc[run.table.index[99], "lwsnl"] = float(run.table["snw"].iloc[99])
    result = score(run)
    assert result.status == FAIL
    # Caught by the opening-liquid check, which is the stronger statement: this
    # case opens frozen, so liquid declared there is refused before the share
    # guard is reached at all.
    assert "opens its block on a frozen pack" in result.message




def test_every_seed_gives_a_dry_scored_block_that_opens_cold_and_ends_warm():
    """Dryness is what the identity rests on; the cold opening is what gives
    the cold-content term something to measure."""
    for seed in SWEEP:
        frame, _ = snowmelt_generate.generate(seed)
        assert len(frame) == snowmelt_generate.N_STEPS, seed
        melt = frame[frame["_regime"] == "melt"]
        assert melt["pr"].sum() == 0.0, f"seed {seed}: scored block is not dry"
        opening = float(melt["tas"].iloc[0])
        assert opening <= -14.0, (
            f"seed {seed}: block opens at {opening:.1f} C, warmer than the -14 C "
            "the generator holds. Two things rest on that opening: a pack less "
            "than 8 K below freezing carries a cold content under five percent "
            "of its own fusion, which the relative rule cannot see; and the "
            "opening-liquid check is only valid because an honest model cannot "
            "hold liquid at this temperature"
        )
        assert float(melt["tas"].iloc[-1]) > 3.0, (
            f"seed {seed}: block never warms enough to melt anything"
        )


def test_the_generator_is_byte_identical_for_a_repeated_seed():
    a, sa = snowmelt_generate.generate(4242)
    b, sb = snowmelt_generate.generate(4242)
    assert a.to_csv(index=False) == b.to_csv(index=False) and sa == sb


def test_more_liquid_than_pack_is_named_as_a_contract_error():
    """`lwsnl` is held inside `snw`; more of it than all of it is impossible."""
    run = build()
    run.table["lwsnl"] = run.table["snw"] * 1.5
    result = score(run)
    assert result.status == FAIL
    assert "contract:" in result.message and "cannot be more than all of it" in result.message


def test_an_empty_pack_with_no_cold_content_reported_is_accepted():
    """The companion to the refusal above: zero on a snow-free row is fine.

    Reporting cold content *with* no ice is refused, and that is asserted in
    `test_cold_content_without_ice_is_named_as_a_contract_error`. This one
    pins the other side, that an empty pack reporting zero passes.
    """
    run = build(spinup=0)
    run.table.loc[run.table.index[:2], ["snw", "lwsnl", "csnow"]] = 0.0
    assert score(run).status == PASS


# --- the relative rule, which the floor hides at this case's melt rates -----
#
# Every test above demands under 40 W/m2, where `max(0.05 * demand, 2.0)` is
# the floor and the 5 percent rule is never the operative bound: delete the
# relative rule and they all still pass. These two demand about 80 W/m2, where
# 5 percent is 4 W/m2 and the floor is not reached, so they fail if the
# relative rule is removed and the pair pins both halves of the tolerance.

HEAVY_MM = 1250.0     # 1250 mm over 60 days demands 80.5 W/m2


@pytest.mark.parametrize("share,expected", [
    (1.04, PASS), (1.06, FAIL),
], ids=["+4%-inside", "+6%-outside"])
def test_the_five_percent_rule_binds_where_the_floor_does_not(share, expected):
    d = demand(HEAVY_MM)
    assert d > 40.0 and 0.05 * d > 2.0, "this case is not in the relative regime"
    run = build(melt_mm=HEAVY_MM, peak=HEAVY_MM, snowfall=HEAVY_MM,
                residual=share * d)
    assert score(run).status == expected


# --- a guard that can silently vanish is not a guard ------------------------

@pytest.mark.parametrize("column,needle", [
    ("pr", "is dry"),
    ("tas", "to bound"),
], ids=["dry-block-check", "cold-content-cap"])
def test_a_missing_forcing_column_refuses_rather_than_dropping_a_guard(column, needle):
    """Both checks read a forcing column, and both would otherwise skip.

    The dry-block check is the precondition the melt inference rests on; the
    cap is what bounds a diagnostic the model reports about itself. A probe
    reusing this criterion without either column would lose the corresponding
    protection with no error, which is the one failure mode a conservation
    gate must not have.
    """
    run = build()
    run.case.forcing = run.case.forcing.drop(columns=[column])
    with pytest.raises(ValueError, match=needle):
        score(run)


# --- the partial flip, which the share guard alone does not close -----------

def flipped(run, extra_mm=1.0):
    """Reduce the ice the block opens with, on that one row and nowhere else.

    `M = ice_before - ice_after - sublimation`, and `ice_before` is a
    self-report on a single row that no other term constrains. Declaring part
    of it liquid shrinks the melt attributed to the model without touching any
    other reported quantity. `min_peak_share_of_snowfall` is evaluated on that
    same row, so it moves with the perturbation instead of opposing it: the
    declared ice can stay above half the snowfall while the melt goes to
    almost nothing.
    """
    tbl = run.table
    b0, b1 = 100, len(tbl) - 1
    ice_end = float(tbl["snw"].iloc[b1] - tbl["lwsnl"].iloc[b1])
    subl = float(tbl["sbl"].iloc[b0:b1 + 1].sum())
    tbl.loc[tbl.index[b0 - 1], "lwsnl"] = float(tbl["snw"].iloc[b0 - 1]) - (
        ice_end + subl + extra_mm
    )
    tbl.loc[tbl.index[b0 - 1], "csnow"] = float(tbl["csnow"].iloc[b1])
    return run


FLIP = dict(melt_mm=250.0, peak=600.0, snowfall=600.0,
            cold_start_c=12.0, cold_end_c=6.0, residual=0.0)
# The block has to leave ice behind, or the declared opening ice falls under
# the share guard and the flip is caught for the wrong reason; and the pack has
# to end with some cold content, or the specific-cold-content check has nothing
# to compare. Both hold for `reference_degree_day` on the probe's own seeds.


def test_the_partial_flip_passes_with_both_new_checks_off():
    """The loophole is real: without them a model that paid nothing passes."""
    run = flipped(build(**FLIP))
    loose = score(run, max_opening_liquid_mm=1e9,
                  specific_cold_content_cannot_rise=False)
    assert loose.status == PASS
    assert loose.diagnostics["blocks"][0]["melt_mm"] == pytest.approx(1.0, abs=1e-6)
    assert float(run.table["snw"].iloc[99] - run.table["lwsnl"].iloc[99]) > 0.5 * 600.0, (
        "the declared ice must still clear the share guard, or the flip proves nothing"
    )


def test_the_opening_liquid_check_closes_the_partial_flip_on_its_own():
    """(a) is what carries this route, and it holds for both constructions."""
    result = score(flipped(build(**FLIP)), specific_cold_content_cannot_rise=False)
    assert result.status == FAIL
    assert "opens its block on a frozen pack" in result.message


def test_the_specific_cold_content_check_closes_only_one_construction():
    """(b) is not independent, and the tests should not imply that it is.

    It sees the flip only when `csnow` is copied from the block's last row,
    which leaves the specific cold content visibly raised. Scale `csnow` to the
    declared ice instead and the excess is round-off: (b) is silent, and only
    (a) stands between that construction and a pass.
    """
    copied = flipped(build(**FLIP))
    assert score(copied, max_opening_liquid_mm=1e9).status == FAIL

    scaled = build(**FLIP)
    tbl, b0 = scaled.table, 100
    ice_b = float(tbl["snw"].iloc[b0 - 1] - tbl["lwsnl"].iloc[b0 - 1])
    flipped(scaled)
    declared = float(tbl["snw"].iloc[b0 - 1] - tbl["lwsnl"].iloc[b0 - 1])
    tbl.loc[tbl.index[b0 - 1], "csnow"] = (
        float(build(**FLIP).table["csnow"].iloc[b0 - 1]) * declared / ice_b
    )
    assert score(scaled, max_opening_liquid_mm=1e9).status == PASS, (
        "(b) is expected to be silent on the scaled construction"
    )
    assert score(scaled).status == FAIL, "(a) must still catch it"


def test_an_honest_table_sits_far_inside_both_new_checks():
    """A synthetic honest table, with the reference margins recorded for scale.

    Measured separately on the three reference models over the probe's gate
    seeds, and quoted here rather than re-run: the opening row carries exactly
    0 mm of liquid; the largest excess over
    `csnow_after <= ice_after * csnow_before / ice_before` is +2.2e-08 J/m2
    over a 20-seed sweep, +1.1e-08 over the five gate seeds;
    and the largest per-step pack gain is +2.5e-13 mm. The tolerances in use
    are 1 J/m2 and 0.5 mm, which clear those by eight and twelve orders of
    magnitude. This test asserts only that a hand-built honest table is inside
    the same checks -- it does not run the references.
    """
    run = build(melt_mm=400.0, peak=600.0, cold_start_c=12.0, cold_end_c=0.0)
    assert float(run.table["lwsnl"].iloc[99]) == 0.0
    ice_b = float(run.table["snw"].iloc[99] - run.table["lwsnl"].iloc[99])
    ice_a = float(run.table["snw"].iloc[-1] - run.table["lwsnl"].iloc[-1])
    excess = float(run.table["csnow"].iloc[-1]) - ice_a * float(run.table["csnow"].iloc[99]) / ice_b
    assert excess <= 1.0, excess
    assert score(run).status == PASS


# --- guards whose removal must break something ------------------------------

def test_a_negative_pack_is_refused_even_with_no_liquid_reported():
    """`snw < 0` on its own, which the ice check cannot see: ice is snw - lwsnl
    and stays non-negative when both are negative together."""
    run = build()
    run.table.loc[run.table.index[99], ["snw", "lwsnl"]] = [-5.0, 0.0]
    result = score(run)
    assert result.status == FAIL
    # The specific message, not just "snw": the `lwsnl > snw` check also names
    # that column, so a looser assertion would pass with the sign check deleted.
    assert "A stored mass and a cold content are both" in result.message
    assert result.diagnostics.get("column") == "snw"


def test_a_block_that_gains_ice_is_refused_rather_than_credited():
    """`M <= 0` is refused rather than credited.

    Built the only physical way a dry block can gain ice: liquid already in the
    pack refreezes, so the pack water never grows and the pack-gain check has
    nothing to say. That requires liquid where the block opens, which the
    opening-liquid check refuses first, so it is relaxed here to reach this
    guard in isolation.

    Note what that combination means for the probe as configured: with an
    opening pack that is frozen and a block that cannot gain water, `M > 0` is
    already implied. The guard is kept as the last line rather than as the
    first, and this test is what fails if it is deleted.
    """
    run = build(melt_mm=0.0, peak=400.0, snowfall=400.0)
    tbl, n_a = run.table, 100
    tbl["snw"] = float(tbl["snw"].iloc[n_a])
    tbl["lwsnl"] = np.concatenate([
        np.full(n_a, 80.0), np.linspace(80.0, 30.0, len(tbl) - n_a),
    ])
    result = score(run, max_opening_liquid_mm=1.0e9)
    assert result.status == FAIL
    assert "did not melt" in result.message


# --- the mirror of the flip, through `snw` on the same row -------------------

def dipped(run, extra_mm=1.0):
    """Lower the *pack* on the opening row instead of declaring part of it liquid.

    `ice_before = snw - lwsnl`, so the same understatement is reachable from
    either side. No liquid is declared, so the opening-liquid check is silent;
    `csnow` is scaled to the declared ice, so the specific-cold-content check is
    silent too; and `closure` is cumulative, so a dip that returns nets to
    nothing and never reaches it.
    """
    tbl = run.table
    b0, b1 = 100, len(tbl) - 1
    ice_end = float(tbl["snw"].iloc[b1] - tbl["lwsnl"].iloc[b1])
    subl = float(tbl["sbl"].iloc[b0:b1 + 1].sum())
    ice_before = float(tbl["snw"].iloc[b0 - 1] - tbl["lwsnl"].iloc[b0 - 1])
    declared = ice_end + subl + extra_mm
    tbl.loc[tbl.index[b0 - 1], "csnow"] = (
        float(tbl["csnow"].iloc[b0 - 1]) * declared / ice_before
    )
    tbl.loc[tbl.index[b0 - 1], "snw"] = float(tbl["lwsnl"].iloc[b0 - 1]) + declared
    return run


def test_the_snw_dip_passes_without_the_pack_gain_check():
    loose = score(dipped(build(**FLIP)), max_pack_gain_mm_per_step=1e9)
    assert loose.status == PASS
    assert loose.diagnostics["blocks"][0]["melt_mm"] == pytest.approx(1.0, abs=1e-6)


def test_the_pack_gain_check_closes_the_snw_dip():
    result = score(dipped(build(**FLIP)))
    assert result.status == FAIL
    assert "gains" in result.message and "in one step" in result.message


def test_deposition_a_model_reports_is_netted_rather_than_refused():
    """The tolerance does not have to cover honest deposition.

    A pack that gains water by deposition and says so in `sbl` nets to zero in
    the check, so it passes at any tolerance. What the tolerance covers is
    *unreported* gain.

    The deposition has to exceed the melt, or the pack never grows and the test
    would pass with the `sbl` netting deleted.
    """
    melt_per_step = 400.0 / 60.0
    deposit = 3.0 * melt_per_step          # the pack grows on every scored step
    run = build(melt_mm=400.0, peak=600.0)
    tbl = run.table
    grow = np.concatenate([np.zeros(100), np.full(len(tbl) - 100, deposit)])
    tbl["snw"] = tbl["snw"].to_numpy() + np.cumsum(grow)
    tbl["sbl"] = tbl["sbl"].to_numpy() - grow          # reported as deposition
    gain = np.diff(tbl["snw"].to_numpy()[99:]).max()
    assert gain > 0.5, "the pack must actually grow, or the netting is untested"
    assert score(run).status == PASS
    # and without the netting it would be refused
    assert "gains" in score(run, sublimation="not_a_column").message


# --- tolerances that round-off cannot trip ----------------------------------

@pytest.mark.parametrize("column,value", [
    ("csnow", -6.0e-8), ("snw", -1.0e-8), ("lwsnl", -1.0e-8),
], ids=["csnow-roundoff", "snw-roundoff", "lwsnl-roundoff"])
def test_round_off_does_not_fail_an_honest_model(column, value):
    """Applied as a uniform shift, which is what round-off looks like.

    Editing one row instead would put a real discontinuity in the series and
    the pack-gain check would be right to object; the sign checks are what is
    under test here, and they must not fire on values this small.
    """
    run = build()
    run.table[column] = run.table[column].to_numpy() + value
    assert "contract:" not in score(run).message


def test_the_cap_reads_state0_when_a_block_opens_the_window():
    """The opening row can be outside the window, and is still bounded."""
    run = build(spinup=100, melt_mm=400.0, peak=600.0, air_c=-20.0)
    run.table.loc[run.table.index[99], "csnow"] = 50.0 * CAP_J
    result = score(run)
    assert result.status == FAIL
    assert "physically plausible on the rows" in result.message
    assert result.diagnostics.get("row") == "the step before"


def test_the_last_row_cap_binds_with_the_specific_cold_content_check_off():
    """Raising only the last row's cold content is (b)-invisible, cap-visible."""
    run = build(melt_mm=400.0, peak=600.0, air_c=-20.0)
    run.table.loc[run.table.index[-1], "csnow"] = 50.0 * CAP_J
    result = score(run, specific_cold_content_cannot_rise=False)
    assert result.status == FAIL
    assert result.diagnostics.get("row") == "the last step of"


def test_a_quoted_false_switches_a_check_off():
    """`bool("false")` is True; a spec author writing it means the opposite."""
    run = flipped(build(**FLIP))
    assert score(run, max_opening_liquid_mm=1e9,
                 specific_cold_content_cannot_rise="false").status == PASS
    assert score(run, max_opening_liquid_mm=1e9,
                 specific_cold_content_cannot_rise="true").status == FAIL


# --- the pack-gain check, at its edges --------------------------------------

def test_the_canopy_cannot_carry_the_dip():
    """`snw` alone, because that is the store the identity differences.

    Summing the canopy with it would leave the same dip reachable by moving
    the water into `canopy` on the opening row: the sum does not move, and
    nothing else looks at either store.
    """
    run = build(**FLIP)
    tbl, b0 = run.table, 100
    before = float(tbl["snw"].iloc[b0 - 1])
    dipped(run)
    moved = before - float(tbl["snw"].iloc[b0 - 1])
    tbl["canopy"] = 0.0                               # the column must exist first
    tbl.loc[tbl.index[b0 - 1], "canopy"] = moved      # hidden, not removed
    assert moved > 100.0
    result = score(run)
    assert result.status == FAIL
    assert "gains" in result.message


@pytest.mark.parametrize("gain,refused", [
    (0.49, False), (0.51, True),
], ids=["inside-the-tolerance", "outside-it"])
def test_the_pack_gain_tolerance_binds_where_it_says(gain, refused):
    """A flat pack, so the added step *is* the per-step gain.

    Against a melting pack the same addition nets a decrease and tests
    nothing: at 400 mm over 60 steps the pack falls 6.67 mm a step, and
    0.51 mm of growth on top of that is still a 6.16 mm fall.
    """
    run = build(melt_mm=0.0, peak=600.0, snowfall=600.0)
    tbl = run.table
    step = np.concatenate([np.zeros(100), np.full(len(tbl) - 100, gain)])
    tbl["snw"] = tbl["snw"].to_numpy() + np.cumsum(step)
    observed = np.diff(tbl["snw"].to_numpy()[99:]).max()
    assert observed == pytest.approx(gain, abs=1e-9)
    assert ("gains" in score(run).message) is refused


def test_sublimation_on_the_row_before_the_block_is_not_counted_as_a_gain():
    """The off-by-one: pairing the opening row with its own `sbl`.

    An honest model that sublimates 0.6 mm on the last accumulation day and
    reports it consistently would otherwise be refused on every seed, because
    the difference starts at that row while the flux is taken from it too.
    """
    run = build(melt_mm=400.0, peak=600.0)
    tbl = run.table
    tbl.loc[tbl.index[99]:, "snw"] = tbl.loc[tbl.index[99]:, "snw"] - 0.6
    tbl.loc[tbl.index[99], "sbl"] = 0.6
    assert "contract:" not in score(run).message


@pytest.mark.parametrize("consistent,refused", [
    (True, False), (False, True),
], ids=["snw-falls-with-it", "snw-stays-flat"])
def test_sublimation_varying_inside_the_block_is_netted_step_by_step(consistent, refused):
    """Alternating 0 and 1.2 mm a day, on a pack that would otherwise be flat.

    Swinging by more than the 0.5 mm tolerance is what makes this test say
    anything: against a melting pack `sbl` moves by hundredths while the pack
    falls 6.67 mm a step, and the check never comes close either way.

    Reported consistently -- `snw` falling by what `sbl` says left -- the
    netting cancels it step by step and nothing is refused. Reported with the
    pack held flat, the same `sbl` says water left while `snw` says it did
    not, and that inconsistency is exactly what the check exists to catch.
    """
    run = build(melt_mm=0.0, peak=600.0, snowfall=600.0)
    tbl = run.table
    n = len(tbl) - 100
    swing = np.concatenate([np.zeros(100),
                            np.tile([0.0, 1.2], (n + 1) // 2)[:n]])
    assert swing.max() > 0.5, "the swing must exceed the tolerance"
    if consistent:
        tbl["snw"] = tbl["snw"].to_numpy() - np.cumsum(swing)
    tbl["sbl"] = tbl["sbl"].to_numpy() + swing
    assert ("gains" in score(run).message) is refused


def test_a_spike_on_the_last_row_is_caught():
    run = build(melt_mm=400.0, peak=600.0)
    run.table.loc[run.table.index[-1], "snw"] = float(run.table["snw"].iloc[-1]) + 50.0
    result = score(run)
    assert result.status == FAIL
    assert "gains" in result.message


@pytest.mark.parametrize("liquid,expected", [
    (1.0e-3, PASS), (2.0, FAIL),
], ids=["trace-liquid-allowed", "a-millimetre-too-much"])
def test_the_opening_liquid_tolerance_binds_at_its_default(liquid, expected):
    """1 mm, not a hair: a smooth freezing curve holds trace liquid."""
    run = build(melt_mm=400.0, peak=600.0)
    tbl = run.table
    tbl.loc[tbl.index[99], "lwsnl"] = liquid
    tbl.loc[tbl.index[99], "snw"] = float(tbl["snw"].iloc[99]) + liquid
    refused = "opens its block on a frozen pack" in score(run).message
    assert refused is (expected is FAIL)


def test_the_probe_pins_both_tolerances_rather_than_taking_the_default():
    """Raising either in `probe.yaml` must be a visible change, not a silent one."""
    import yaml
    spec = yaml.safe_load(
        (Path(__file__).resolve().parents[1]
         / "probes/energy/snowmelt-energy-water/probe.yaml").read_text()
    )
    params = next(c["melt_energy"] for c in spec["criteria"] if "melt_energy" in c)
    assert params["max_pack_gain_mm_per_step"] == 0.5
    assert params["max_opening_liquid_mm"] == 1.0


@pytest.mark.parametrize("value", [None, ""], ids=["yaml-null", "empty-string"])
def test_an_empty_specific_cold_content_setting_keeps_the_check_on(value):
    """An empty YAML value is None, and `bool("")` is False; neither means off."""
    run = flipped(build(**FLIP))
    result = score(run, max_opening_liquid_mm=1e9,
                   specific_cold_content_cannot_rise=value)
    assert result.status == FAIL
    assert "per unit ice rises" in result.message

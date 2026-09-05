"""The latent heat a model reports and the water it evaporates have to
describe one evaporation. These test `flux_identity` on tables written by
hand, because the gate alone cannot say which of its branches did the work,
and because the branch that used to guess the sublimating share needed a
counterexample rather than a green gate to be found."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing.criteria import get
from hydroturing.criteria.base import PASS
from hydroturing.protocol import Case, RunResult

LAMBDA_A, LAMBDA_B = 2.501e6, -2361.0
LAMBDA_F = 3.337e5
SECONDS = 86400.0
SPINUP = 5
N = 200


def lambda_v(tas):
    return LAMBDA_A + LAMBDA_B * tas


def build(evspsbl, tas, hfls, snw, pr, sbl=None):
    """One run of a made-up model, as the criterion receives it.

    `snw=None` is a model that reports no snow state at all."""
    columns = {"time": np.arange(N), "evspsbl": evspsbl, "hfls": hfls}
    if snw is not None:
        columns["snw"] = snw
    if sbl is not None:
        columns["sbl"] = sbl
    forcing = pd.DataFrame({"time": np.arange(N), "pr": pr, "tas": tas})
    case = Case(probe_id="t", seed=1, forcing=forcing, static={}, spinup_steps=SPINUP)
    return RunResult(case=case, table=pd.DataFrame(columns), meta={}, wall_seconds=0.0)


@pytest.fixture
def winter():
    """Half the record below freezing with a pack, half hot and snow-free.

    The warm half is hot and evaporating hard on purpose: the error a constant
    latent heat makes is 3.4 percent at most, so a test that wants to see it
    has to look where lambda_v is furthest from the constant and where there
    is enough water moving for the absolute floor not to swallow it."""
    tas = np.concatenate([np.full(N // 2, -6.0), np.full(N // 2, 32.0)])
    pr = np.zeros(N)
    pr[10:20] = 4.0                      # snowfall, since it is below freezing
    snw = np.where(np.arange(N) < N // 2, 50.0, 0.0)
    evspsbl = np.where(np.arange(N) < N // 2, 1.5, 5.0)
    return tas, pr, snw, evspsbl


def run(**kwargs):
    return get("flux_identity")(build(**kwargs), None, {})


def test_the_exact_identity_passes(winter):
    tas, pr, snw, evspsbl = winter
    sbl = np.where(snw > 0, evspsbl, 0.0)
    hfls = (lambda_v(tas) * (evspsbl - sbl) + (LAMBDA_A + LAMBDA_F) * sbl) / SECONDS
    result = run(evspsbl=evspsbl, tas=tas, hfls=hfls, snw=snw, pr=pr, sbl=sbl)
    assert result.status == PASS, result.message
    assert result.diagnostics["reported_split"] is True


def test_a_constant_latent_heat_is_caught(winter):
    tas, pr, snw, evspsbl = winter
    sbl = np.where(snw > 0, evspsbl, 0.0)
    hfls = 2.45e6 * evspsbl / SECONDS
    assert run(evspsbl=evspsbl, tas=tas, hfls=hfls, snw=snw, pr=pr, sbl=sbl).status != PASS


def test_vaporisation_applied_to_sublimation_is_caught(winter):
    """The reported split says these kilograms left as ice; charging them at
    the latent heat of vaporisation is 13.3 percent short."""
    tas, pr, snw, evspsbl = winter
    sbl = np.where(snw > 0, evspsbl, 0.0)
    hfls = lambda_v(tas) * evspsbl / SECONDS
    result = run(evspsbl=evspsbl, tas=tas, hfls=hfls, snw=snw, pr=pr, sbl=sbl)
    assert result.status != PASS
    assert result.diagnostics["violating_steps"] > 0


def test_a_pack_that_loses_water_at_its_base_is_not_a_violation(winter):
    """Snow-17's DAYGM moves water from the pack into the soil without any
    latent heat. The pack's mass loss is therefore not the sublimated mass,
    and a criterion that reads it as such fails an honest model. This is the
    counterexample the first version of the criterion did not survive."""
    tas, pr, snw, evspsbl = winter
    # The pack drains 0.3 mm/day at its base on top of the 1.5 mm/day that
    # sublimates, so the store falls faster than the evaporation alone would
    # explain. Under the old inference the extra 0.3 was charged the latent
    # heat of sublimation and the model was told it had violated physics.
    step = np.arange(N)
    snw = np.where(step < N // 2, np.clip(400.0 - 1.8 * step, 0.0, None), 0.0)
    sbl = np.where(snw > 0, evspsbl, 0.0)
    hfls = (lambda_v(tas) * (evspsbl - sbl) + (LAMBDA_A + LAMBDA_F) * sbl) / SECONDS
    assert run(evspsbl=evspsbl, tas=tas, hfls=hfls, snw=snw, pr=pr, sbl=sbl).status == PASS


def test_a_model_that_reports_no_split_is_bounded_not_guessed_at(winter):
    """Without `sbl` the criterion asserts only the interval the two latent
    heats span wherever a pack is or could be present. Both extremes pass;
    a value outside them does not."""
    tas, pr, snw, evspsbl = winter
    inside_low = lambda_v(tas) * evspsbl / SECONDS
    inside_high = np.where(snw > 0, (LAMBDA_A + LAMBDA_F) * evspsbl / SECONDS, inside_low)
    outside = np.where(snw > 0, 3.2e6 * evspsbl / SECONDS, inside_low)

    for hfls in (inside_low, inside_high):
        r = run(evspsbl=evspsbl, tas=tas, hfls=hfls, snw=snw, pr=pr)
        assert r.status == PASS, r.message
        assert r.diagnostics["reported_split"] is False
    assert run(evspsbl=evspsbl, tas=tas, hfls=outside, snw=snw, pr=pr).status != PASS


def test_snow_free_steps_are_still_an_equality_without_a_split(winter):
    """The interval is only for steps where a pack is or could be present.
    Where the model reports none and none can fall, lambda_v is exact, which
    is what keeps a constant-lambda model catchable when it reports no split."""
    tas, pr, snw, evspsbl = winter
    hfls = np.where(snw > 0, lambda_v(tas), 2.45e6) * evspsbl / SECONDS
    assert run(evspsbl=evspsbl, tas=tas, hfls=hfls, snw=snw, pr=pr).status != PASS


def test_a_fictitious_split_is_rejected(winter):
    """A split reported where the model itself says there is no snow would let
    a warm-season model bend its effective lambda upwards."""
    tas, pr, snw, evspsbl = winter
    sbl = np.full(N, 1.5)                       # claimed even where snw is zero
    hfls = (LAMBDA_A + LAMBDA_F) * evspsbl / SECONDS
    result = run(evspsbl=evspsbl, tas=tas, hfls=hfls, snw=snw, pr=pr, sbl=sbl)
    assert result.status != PASS
    assert "no snow" in result.message


def test_a_split_larger_than_the_evaporation_is_rejected(winter):
    tas, pr, snw, evspsbl = winter
    sbl = np.where(snw > 0, evspsbl * 2.0, 0.0)
    hfls = (LAMBDA_A + LAMBDA_F) * evspsbl / SECONDS
    result = run(evspsbl=evspsbl, tas=tas, hfls=hfls, snw=snw, pr=pr, sbl=sbl)
    assert result.status != PASS
    assert "share of" in result.message


def test_a_model_with_no_snow_state_is_bounded_everywhere(winter):
    """A model that reports no snow store has told us nothing about phase, so
    every step is one where a pack could be present and the interval is all
    the criterion may assert. Holding such a model to the liquid equality
    would be an implicit assumption that it never sublimates."""
    tas, pr, _, evspsbl = winter
    sublimating = (LAMBDA_A + LAMBDA_F) * evspsbl / SECONDS
    liquid = lambda_v(tas) * evspsbl / SECONDS
    outside = 3.2e6 * evspsbl / SECONDS

    for hfls in (liquid, sublimating):
        r = run(evspsbl=evspsbl, tas=tas, hfls=hfls, snw=None, pr=pr)
        assert r.status == PASS, r.message
        assert r.diagnostics["bounded_steps"] == len(tas) - SPINUP
    assert run(evspsbl=evspsbl, tas=tas, hfls=outside, snw=None, pr=pr).status != PASS


# --------------------------------------------------------------- partition_shift

W_START, W_LEN = 40, 60


def paired(d_latent, d_sensible, d_ground, rn=120.0, drift=0.0, label=True):
    """Two runs differing only inside the marked window, by construction.

    The control is flat and the perturbed run differs from it by exactly the
    per-step changes asked for, so the identity the criterion asserts is
    whatever the caller chose to make it."""
    window = np.zeros(N)
    window[W_START:W_START + W_LEN] = 1.0
    columns = {"time": np.arange(N), "pr": np.zeros(N), "tas": np.full(N, 15.0),
               "rn": np.full(N, rn), "_perturbed": window}
    if label:
        onset = np.zeros(N)
        onset[W_START + 30:W_START + W_LEN] = 1.0
        columns["_intensification"] = onset
    control_forcing = pd.DataFrame(columns)
    perturbed_forcing = control_forcing.copy()
    perturbed_forcing["rn"] = control_forcing["rn"] + drift * window

    def table(dle, dh, dg):
        return pd.DataFrame({
            "time": np.arange(N),
            "hfls": 80.0 + dle * window,
            "hfss": 30.0 + dh * window,
            "hfg": 10.0 + dg * window,
        })

    def result(forcing, dle, dh, dg):
        case = Case(probe_id="t", seed=1, forcing=forcing, static={}, spinup_steps=SPINUP)
        return RunResult(case=case, table=table(dle, dh, dg), meta={}, wall_seconds=0.0)

    return {
        "control": result(control_forcing, 0.0, 0.0, 0.0),
        "drought": result(perturbed_forcing, d_latent, d_sensible, d_ground),
    }


def shift(**kwargs):
    params = {"control": "control", "perturbed": "drought"}
    return get("partition_shift")(paired(**kwargs), None, params)


def test_the_partition_identity_passes():
    """Latent falls, sensible takes exactly what it gave up, ground stays put."""
    r = shift(d_latent=-40.0, d_sensible=+40.0, d_ground=0.0)
    assert r.status == PASS, r.message
    assert r.diagnostics["d_evaporative_fraction"] < 0
    assert r.diagnostics["onset"]["steps"] == 30


def test_a_model_that_does_not_respond_is_caught():
    """reference_two_head in miniature: its energy side never reads its water
    side, so nothing moves at all."""
    r = shift(d_latent=0.0, d_sensible=0.0, d_ground=0.0)
    assert r.status != PASS
    assert "barely moved" in r.message


def test_energy_absorbed_by_the_ground_is_caught():
    """reference_ground_dodge in miniature: the sum is satisfied and the
    energy went into the soil instead of the air."""
    r = shift(d_latent=-40.0, d_sensible=0.0, d_ground=+40.0)
    assert r.status != PASS
    assert "absorbed" in r.message


def test_energy_that_does_not_add_up_is_caught():
    r = shift(d_latent=-40.0, d_sensible=+10.0, d_ground=0.0)
    assert r.status != PASS
    assert "not zero" in r.message


def test_a_driver_that_drifts_is_a_generator_fault_not_a_verdict():
    """The identity holds only while net radiation is held fixed. A generator
    that lets it drift should be diagnosed as a generator fault rather than
    reported as a model failure."""
    with pytest.raises(ValueError, match="held fixed"):
        shift(d_latent=-40.0, d_sensible=+40.0, d_ground=0.0, drift=5.0)


def test_the_intensification_stretch_is_reported_and_never_gated():
    """Without the label the criterion still scores; with it, it also says how
    much of the shift landed in the marked stretch."""
    unlabelled = shift(d_latent=-40.0, d_sensible=+40.0, d_ground=0.0, label=False)
    assert unlabelled.status == PASS
    assert unlabelled.diagnostics["onset"] == {}


def test_variants_that_disagree_about_the_window_are_a_generator_fault():
    """The criterion indexes both variants with one mask, so a generator that
    marks a different stretch in each, or produces records of different
    lengths, has to be told so rather than crashing inside numpy or silently
    comparing two unrelated stretches of record."""
    runs = paired(d_latent=-40.0, d_sensible=+40.0, d_ground=0.0)
    params = {"control": "control", "perturbed": "drought"}

    moved = runs["drought"].case.forcing.copy()
    shifted = np.zeros(N)
    shifted[W_START + 5:W_START + W_LEN + 5] = 1.0
    moved["_perturbed"] = shifted
    runs["drought"].case.forcing = moved
    with pytest.raises(ValueError, match="different steps"):
        get("partition_shift")(runs, None, params)

    dropped = moved.drop(columns=["_perturbed"])
    runs["drought"].case.forcing = dropped
    with pytest.raises(ValueError, match="perturbed variant has none"):
        get("partition_shift")(runs, None, params)


def test_variants_of_different_lengths_are_a_generator_fault():
    runs = paired(d_latent=-40.0, d_sensible=+40.0, d_ground=0.0)
    params = {"control": "control", "perturbed": "drought"}
    runs["drought"].table = runs["drought"].table.iloc[:-10].reset_index(drop=True)
    runs["drought"].case.forcing = runs["drought"].case.forcing.iloc[:-10].reset_index(drop=True)
    with pytest.raises(ValueError, match="different scored lengths"):
        get("partition_shift")(runs, None, params)

"""counterfactual_response on one perturbed variant or several.

The probe mass/precipitation-counterfactual pairs three perturbed variants
with one control. These tests pin the single-name behaviour, score every
variant when a list is given, isolate a failure to the variant that caused
it, fail a term that moves against the change without calling a share a
hair below zero one, and refuse an empty list or a variant that changes
nothing.
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult
from hydroturing.runner import get_runner
from hydroturing.seeds import gate_seeds


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/precipitation-counterfactual")


@pytest.fixture(scope="module")
def runs(probe, tmp_path_factory):
    model = registry.find_model("reference_bucket")
    seed = gate_seeds(probe.id, 1)[0]
    root = tmp_path_factory.mktemp("cf")
    return {
        variant: get_runner(model).run(model, probe, build_case(probe, seed, variant), root / variant)
        for variant in probe.variants
    }


@pytest.fixture(scope="module")
def params(probe):
    return next(c.params for c in probe.criteria if c.name == "counterfactual_response")


def score(runs, probe, params, **override):
    return get("counterfactual_response")(runs, probe, {**params, **override})


def test_a_single_name_keeps_the_original_shape(runs, probe, params):
    result = score(runs, probe, params, perturbed="wetter20")
    assert result.passed
    spinup = runs["control"].case.spinup_steps
    scored_rain = runs["control"].case.forcing["pr"].iloc[spinup:].sum()
    assert result.diagnostics["added_mm"] == pytest.approx(0.2 * scored_rain, rel=1e-6)
    assert set(result.diagnostics) == {"added_mm", "shares", "largest_share", "accounted"}
    assert result.value == pytest.approx(sum(result.diagnostics["shares"].values()))
    assert "wetter20" not in result.message


def test_a_list_scores_every_variant(runs, probe, params):
    result = score(runs, probe, params)
    assert result.passed
    variants = result.diagnostics["variants"]
    assert set(variants) == {"wetter20", "wetter10", "drier20"}
    assert variants["drier20"]["added_mm"] < 0 < variants["wetter10"]["added_mm"] < variants["wetter20"]["added_mm"]
    for name, v in variants.items():
        assert 0.5 < v["shares"]["mrro"] < 0.95, name
        assert 0.03 < v["shares"]["evspsbl"] < 0.5, name
        assert v["accounted"] == pytest.approx(1.0, abs=1e-6), name
    assert "drier20" in result.message


def test_an_unresponsive_variant_fails_alone(runs, probe, params):
    unresponsive = dict(runs)
    d = runs["drier20"]
    unresponsive["drier20"] = RunResult(d.case, runs["control"].table, d.meta, 0.0)
    result = score(unresponsive, probe, params)
    assert not result.passed
    assert "drier20" in result.message
    assert "wetter20" not in result.message and "wetter10" not in result.message


def with_evaporation_scaled_from_control(runs, name, factor):
    """`runs` with the scored evaporation of `name` set to `factor` times the control's.

    What evaporation gives up or takes on is put in the final soil store, so
    the accounted sum stays where it was and only the minimum can object.
    """
    w, c = runs[name], runs["control"]
    spinup = w.case.spinup_steps
    table = w.table.copy()
    moved = (table["evspsbl"] - factor * c.table["evspsbl"]).iloc[spinup:]
    table.loc[spinup:, "evspsbl"] -= moved
    table.loc[table.index[-1], "mrso"] += moved.sum()
    return {**runs, name: RunResult(w.case, table, w.meta, 0.0)}


@pytest.mark.parametrize("name, factor", [("wetter20", 0.9999), ("drier20", 1.0001)])
def test_a_share_just_below_zero_barely_responds(runs, probe, params, name, factor):
    # A term that all but ignores the change can land a hair below zero, as
    # sacsma_snow17's evaporation does on a 30-day window. That fails the
    # minimum as a term that barely responds, not as one moving the wrong way.
    against = with_evaporation_scaled_from_control(runs, name, factor)
    result = score(against, probe, params, perturbed=name)
    assert not result.passed
    assert -params["min_share"] < result.diagnostics["shares"]["evspsbl"] < 0
    assert result.message.startswith("evspsbl barely responds (-0.000") and ";" not in result.message


@pytest.mark.parametrize("name, factor", [("wetter20", 0.95), ("drier20", 1.05)])
def test_a_term_moving_against_the_change_fails(runs, probe, params, name, factor):
    # Evaporation ends 5% below the control's when rain is added, or 5% above
    # it when rain is removed. The water it trades sits in storage, so only
    # the signed minimum has anything to object to, alone or in the list.
    against = with_evaporation_scaled_from_control(runs, name, factor)
    alone = score(against, probe, params, perturbed=name)
    assert not alone.passed
    assert alone.diagnostics["shares"]["evspsbl"] < 0
    assert alone.message.startswith("evspsbl moves the wrong way") and ";" not in alone.message
    listed = score(against, probe, params)
    assert listed.message.startswith(f"{name}: evspsbl moves the wrong way")
    assert ";" not in listed.message


def test_an_empty_list_of_variants_raises(runs, probe, params):
    with pytest.raises(ValueError, match="at least one variant"):
        score(runs, probe, params, perturbed=[])


def test_a_variant_that_changes_nothing_raises(runs, probe, params):
    same = dict(runs)
    c = runs["control"]
    same["same"] = RunResult(c.case, c.table, c.meta, 0.0)
    with pytest.raises(ValueError, match="does not change pr"):
        score(same, probe, params, perturbed=["wetter20", "same"])


def test_a_window_of_another_length_raises(runs, probe, params):
    d = runs["drier20"]
    short_case = replace(d.case, forcing=d.case.forcing.iloc[:-10].reset_index(drop=True))
    short = dict(runs)
    short["drier20"] = RunResult(short_case, d.table.iloc[:-10].reset_index(drop=True), d.meta, 0.0)
    with pytest.raises(ValueError, match="different length"):
        score(short, probe, params)

"""counterfactual_response on one perturbed variant or several.

The probe mass/precipitation-counterfactual pairs three perturbed variants
with one control. These tests pin the single-name behaviour, score every
variant when a list is given, isolate a failure to the variant that caused
it, and refuse a variant that changes nothing.
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

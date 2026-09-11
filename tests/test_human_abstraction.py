"""Tests for the human-abstraction probe: a prescribed withdrawal must leave the budget.

The same seed is run twice with byte-identical weather, the second adding a
net withdrawal in the visible forcing column `abstr`. These check that only the
withdrawal changes between the variants, that the exact and the four physical
baselines pass while a model blind to the column fails by the whole abstracted
volume, that a difference which does not account for the withdrawal is caught,
that a weather mismatch is refused rather than attributed, and that the
absolute floor rescues a negligible withdrawal.
"""

from __future__ import annotations

import tempfile
from dataclasses import replace
from pathlib import Path

import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.harness import build_case, run_probe
from hydroturing.protocol import RunResult
from hydroturing.runner import get_runner
from hydroturing.scoring import FAIL, PASS
from hydroturing.seeds import gate_seeds


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/human-abstraction")


def _params(probe):
    return dict(next(c.params for c in probe.criteria if c.name == "human_abstraction"))


def _runs(probe, seed):
    model = registry.find_model("reference_bucket")
    runs = {}
    for variant in probe.variants:
        case = build_case(probe, seed, variant)
        runs[variant] = get_runner(model).run(model, probe, case, Path(tempfile.mkdtemp()))
    return runs


def test_only_the_withdrawal_changes(probe):
    control = build_case(probe, 7, "natural")
    irrigated = build_case(probe, 7, "irrigated")
    for column in ("pr", "tas", "pet"):
        assert control.forcing[column].equals(irrigated.forcing[column])
    assert control.forcing["time"].equals(irrigated.forcing["time"])
    assert (control.forcing["abstr"] == 0.0).all()
    assert irrigated.forcing["abstr"].sum() > 0.0


def test_exact_model_passes(probe):
    outcome = run_probe(registry.find_model("reference_bucket"), probe, gate_seeds(probe.id, 2))
    assert outcome.verdict == PASS, outcome.failing
    result = next(c for c in outcome.criteria if c.name == "human_abstraction")
    assert result.diagnostics["abstracted_mm"] > 0.0
    assert abs(result.diagnostics["residual_mm"]) <= 5.0
    assert result.value <= 0.05


@pytest.mark.parametrize(
    "name",
    ["reference_bucket", "flex_lumped", "flex_topo", "sacsma_snow17"],
)
def test_every_must_pass_baseline_passes(probe, name):
    outcome = run_probe(registry.find_model(name), probe, gate_seeds(probe.id, 1))
    assert outcome.verdict == PASS, (name, outcome.failing)


def test_abstraction_blind_is_caught(probe):
    outcome = run_probe(registry.find_model("reference_abstraction_blind"), probe, gate_seeds(probe.id, 1))
    assert outcome.verdict == FAIL
    assert "human_abstraction" in outcome.failing
    result = next(c for c in outcome.criteria if c.name == "human_abstraction")
    # The whole abstracted volume is unaccounted for.
    assert result.value == pytest.approx(1.0, abs=0.02)


def test_a_difference_that_does_not_account_is_caught(probe):
    seed = gate_seeds(probe.id, 1)[0]
    runs = _runs(probe, seed)
    # The irrigated run reports exactly the control's fluxes: the withdrawal
    # left no trace in the reported budget, though both runs still close.
    runs["irrigated"] = RunResult(
        runs["irrigated"].case, runs["natural"].table.copy(), runs["irrigated"].meta, 0.0
    )
    result = get("human_abstraction")(runs, probe, _params(probe))
    assert not result.passed
    assert result.value == pytest.approx(1.0, abs=0.02)


def test_weather_must_match_between_variants(probe):
    seed = gate_seeds(probe.id, 1)[0]
    runs = _runs(probe, seed)
    forcing = runs["irrigated"].case.forcing.copy()
    forcing["pr"] = forcing["pr"] + 1.0
    runs["irrigated"] = RunResult(
        replace(runs["irrigated"].case, forcing=forcing),
        runs["irrigated"].table,
        runs["irrigated"].meta,
        0.0,
    )
    with pytest.raises(ValueError):
        get("human_abstraction")(runs, probe, _params(probe))


def test_the_absolute_floor_rescues_a_tiny_withdrawal(probe):
    seed = gate_seeds(probe.id, 1)[0]
    runs = _runs(probe, seed)
    forcing = runs["irrigated"].case.forcing.copy()
    forcing["abstr"] = 0.001  # ~3.65 mm over the scored decade, under the 5 mm floor
    runs["irrigated"] = RunResult(
        replace(runs["irrigated"].case, forcing=forcing),
        runs["natural"].table.copy(),
        runs["irrigated"].meta,
        0.0,
    )
    params = _params(probe)
    params["floor_mm"] = 5.0
    result = get("human_abstraction")(runs, probe, params)
    assert result.passed
    assert result.diagnostics["abstracted_mm"] < 5.0
    # One effective limit, so a pass never reads as value above threshold.
    assert result.value <= result.threshold


def test_an_absent_or_zero_withdrawal_raises(probe):
    seed = gate_seeds(probe.id, 1)[0]
    runs = _runs(probe, seed)
    forcing = runs["irrigated"].case.forcing.copy()
    forcing["abstr"] = 0.0
    runs["irrigated"] = RunResult(
        replace(runs["irrigated"].case, forcing=forcing),
        runs["irrigated"].table,
        runs["irrigated"].meta,
        0.0,
    )
    with pytest.raises(ValueError):
        get("human_abstraction")(runs, probe, _params(probe))

"""mass/exchange-response: the criterion fires, and fires for the right reasons.

Two things make a test worth having here beyond the gate. The criterion is
conditional on what the model declared, and `RunResult.model` defaults to
None, so a hand-built RunResult — the repository's usual way to test a
criterion — has to be shown to judge rather than silently switch the gate off.
And the gate is a paired counterfactual whose timing is load-bearing: the
shift has to begin with the scored record, or a fast aquifer answers it during
spinup and is never seen answering.
"""

from __future__ import annotations

import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.criteria.exchange import (reversal_fraction, significant_reversals,
                                           turning_points)
from hydroturing.harness import build_case, run_probe
from hydroturing.runner import get_runner
from hydroturing.scoring import FAIL, NOT_SCORED, PASS
from hydroturing.seeds import gate_seeds


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/exchange-response")


def _params(probe):
    return dict(next(c.params for c in probe.criteria if c.name == "exchange_response"))


def _runs(probe, model_name, seed):
    model = registry.find_model(model_name)
    runs = {}
    for variant in probe.variants:
        case = build_case(probe, seed, variant)
        runs[variant] = get_runner(model).run(model, probe, case, Path(tempfile.mkdtemp()))
    return runs


# --- the case ---------------------------------------------------------------

def test_only_the_head_changes_and_only_after_spinup(probe):
    control = build_case(probe, 7, "control")
    raised = build_case(probe, 7, "raised")
    lowered = build_case(probe, 7, "lowered")
    for column in ("time", "pr", "tas", "pet", "_regime"):
        assert control.forcing[column].equals(raised.forcing[column])
        assert control.forcing[column].equals(lowered.forcing[column])
    spin = control.spinup_steps
    h0 = control.forcing["gwh"].to_numpy(dtype=float)
    up = raised.forcing["gwh"].to_numpy(dtype=float)
    down = lowered.forcing["gwh"].to_numpy(dtype=float)
    # identical spinup, so every variant enters the scored record in one state
    assert np.array_equal(h0[:spin], up[:spin])
    assert np.array_equal(h0[:spin], down[:spin])
    # a constant shift over the whole scored record
    assert np.allclose(up[spin:] - h0[spin:], +0.5)
    assert np.allclose(down[spin:] - h0[spin:], -0.5)


def test_the_probe_requires_the_driver(probe):
    """A model that does not declare gwh is INCOMPATIBLE, not passed."""
    assert "gwh" in probe.requires_forcing
    outcome = run_probe(registry.find_model("reference_bucket"), probe, gate_seeds(probe.id, 1))
    assert outcome.verdict == NOT_SCORED
    assert outcome.reason == "INCOMPATIBLE"


# --- the controls -----------------------------------------------------------

@pytest.mark.parametrize("name", ["reference_driven_exchange", "reference_evolving_exchange"])
def test_positive_controls_pass(probe, name):
    outcome = run_probe(registry.find_model(name), probe, gate_seeds(probe.id, 1))
    assert outcome.verdict == PASS, outcome.failing
    result = next(c for c in outcome.criteria if c.name == "exchange_response")
    assert result.passed
    assert result.value is not None and result.value >= result.threshold


def test_evolving_control_answers_with_the_water_that_moves_its_head(probe):
    """S dh/dt = C (H - h): a 0.5 m shift admits about S * 0.5 mm, and no more."""
    runs = _runs(probe, "reference_evolving_exchange", gate_seeds(probe.id, 1)[0])
    result = get("exchange_response")(runs, probe, _params(probe))
    assert result.passed, result.message
    up = result.diagnostics["response_mm"]["raised"]
    down = result.diagnostics["response_mm"]["lowered"]
    # storativity 10 mm/m, shift 0.5 m -> 5 mm; a little more on the raised
    # side as the head keeps tracking the moving driver, a little less on the
    # lowered side where an empty soil limits the loss
    assert 3.0 < up < 8.0, up
    assert -8.0 < down < -3.0, down


@pytest.mark.parametrize("name", ["reference_noise_sink", "reference_token_exchange"])
def test_negative_controls_are_caught_by_the_headline(probe, name):
    outcome = run_probe(registry.find_model(name), probe, gate_seeds(probe.id, 1))
    assert outcome.verdict == FAIL
    assert outcome.failing == ["exchange_response"], outcome.failing


def test_token_term_fails_on_share_not_on_sign(probe):
    """The token cheat answers the head with the right sign; it fails because
    the answer is a millionth of the exchange it declared."""
    runs = _runs(probe, "reference_token_exchange", gate_seeds(probe.id, 1)[0])
    result = get("exchange_response")(runs, probe, _params(probe))
    assert not result.passed
    resp = result.diagnostics["response_mm"]
    assert resp["raised"] > 0 and resp["lowered"] < 0          # right sign
    assert result.diagnostics["response_share_of_gross"] < 1e-5  # wrong size


# --- the conditional gate ----------------------------------------------------

def test_a_hand_built_run_without_a_manifest_is_still_judged(probe):
    """RunResult.model defaults to None. That must mean 'judge', not 'skip':
    the repository's usual hand-built RunResult would otherwise pass every
    cheat by turning the gate off."""
    seed = gate_seeds(probe.id, 1)[0]
    bad = {k: replace(r, model=None) for k, r in _runs(probe, "reference_noise_sink", seed).items()}
    assert not get("exchange_response")(bad, probe, _params(probe)).passed
    good = {k: replace(r, model=None) for k, r in _runs(probe, "reference_driven_exchange", seed).items()}
    assert get("exchange_response")(good, probe, _params(probe)).passed


def test_a_manifest_that_does_not_declare_the_driver_is_a_contract_error(probe):
    """The harness makes such a model INCOMPATIBLE before any run; reaching
    the criterion anyway is an error, not an excuse."""
    seed = gate_seeds(probe.id, 1)[0]
    runs = _runs(probe, "reference_driven_exchange", seed)
    stripped = {k: replace(r, model=replace(r.model, uses_forcing=(), needs_forcing=("pr", "tas", "pet")))
                for k, r in runs.items()}
    with pytest.raises(ValueError, match="does not declare"):
        get("exchange_response")(stripped, probe, _params(probe))


def test_a_shift_during_spinup_is_refused(probe):
    """The timing is load-bearing, so a generator that shifted the head from
    the start of spinup is refused rather than scored."""
    seed = gate_seeds(probe.id, 1)[0]
    runs = _runs(probe, "reference_driven_exchange", seed)
    forcing = runs["raised"].case.forcing.copy()
    forcing.loc[: runs["raised"].case.spinup_steps - 1, "gwh"] += 0.5
    runs["raised"] = replace(runs["raised"], case=replace(runs["raised"].case, forcing=forcing))
    with pytest.raises(ValueError, match="during spinup"):
        get("exchange_response")(runs, probe, _params(probe))


# --- the reversal diagnostics ------------------------------------------------

def test_reversal_fraction_charges_flips_not_slow_crossings():
    one_signed = [-1.0] * 100
    alternating = [1.0, -1.0] * 50
    slow = np.sin(2 * np.pi * np.arange(730) / 365)
    assert reversal_fraction(one_signed)[0] == 0.0
    assert reversal_fraction(alternating)[0] > 0.98
    assert reversal_fraction(slow)[0] < 1e-3
    assert reversal_fraction([3.0, -1.0] * 50)[0] > 0.45      # the +3,-1 loophole


def test_significant_reversals_ignore_rounding_jitter():
    jitter = np.array([1.0, 1e-9, -1e-9, 1e-9, 1.0])
    assert significant_reversals(jitter, deadband=0.05) == 0
    assert significant_reversals(np.array([1.0, -1.0, 1.0]), deadband=0.05) == 2


def test_turning_points_count_extrema():
    h = np.sin(2 * np.pi * np.arange(90) / 20.0)
    assert 7 <= turning_points(h, deadband=0.05) <= 9

"""Interception response: the paired static counterfactual."""

import pytest

from hydroturing import registry
from hydroturing.criteria import is_paired
from hydroturing.harness import run_probe
from hydroturing.seeds import gate_seeds

PROBE = "mass/interception-response"
POSITIVE = ("reference_bucket", "flex_lumped", "flex_topo")
NEGATIVE = (
    "reference_interception_blind",
    "reference_interception_wetter",
    "reference_interception_drier",
    "reference_interception_leaky",
)


def _probe():
    return registry.find_probe(PROBE)


def _seed():
    return gate_seeds(PROBE, 1)[0]


def test_the_criterion_is_paired():
    assert is_paired("interception_response")


@pytest.mark.parametrize("model", POSITIVE)
def test_must_pass_models_respond_and_close(model):
    outcome = run_probe(registry.find_model(model), _probe(), [_seed()])
    assert outcome.verdict == "PASS", outcome.failing


@pytest.mark.parametrize("model", NEGATIVE)
def test_each_negative_control_fails_the_paired_criterion(model):
    outcome = run_probe(registry.find_model(model), _probe(), [_seed()])
    assert outcome.verdict == "FAIL", outcome.failing
    assert "interception_response" in outcome.failing


def test_the_probe_declares_the_pair_and_the_static():
    probe = _probe()
    assert probe.variants == ("capacity_low", "capacity_high")
    assert probe.requires_static == ("canopy_capacity_mm",)
    assert [c.name for c in probe.criteria] == ["interception_response"]

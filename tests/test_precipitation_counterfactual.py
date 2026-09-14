"""mass/precipitation-counterfactual: four rainfalls, one weather."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from hydroturing import registry
from hydroturing.harness import (
    build_case,
    load_generator,
    resolve_window_days,
    run_probe,
)
from hydroturing.scoring import FAIL, PASS
from hydroturing.seeds import gate_seeds

VALIDATION_SEED = 20260911


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/precipitation-counterfactual")


@pytest.mark.parametrize("variant, factor", [("wetter20", 1.20), ("wetter10", 1.10), ("drier20", 0.80)])
def test_only_the_scored_rain_changes(probe, variant, factor):
    control = build_case(probe, 7, "control")
    other = build_case(probe, 7, variant)
    assert control.n_steps == other.n_steps == 4015
    assert control.static == other.static
    for column in ("time", "tas", "pet"):
        assert control.forcing[column].equals(other.forcing[column]), column
    spinup = control.spinup_steps
    assert control.forcing["pr"][:spinup].equals(other.forcing["pr"][:spinup])
    # Each value is rounded to 1e-6 mm after scaling, so compare depths, not ratios.
    np.testing.assert_allclose(
        other.forcing["pr"][spinup:], factor * control.forcing["pr"][spinup:], atol=2e-6
    )


def test_generation_is_deterministic_and_variants_are_checked(probe):
    generator = load_generator(probe)
    first = generator.generate(7, "wetter20")[0]
    assert first.to_csv(index=False) == generator.generate(7, "wetter20")[0].to_csv(index=False)
    assert not first["pr"].equals(generator.generate(8, "wetter20")[0]["pr"])
    with pytest.raises(ValueError, match="unknown variant"):
        generator.generate(7, "wetter")


def test_submitted_models_are_scored_on_the_whole_record(probe):
    submitted = replace(registry.find_model("reference_bucket"), name="submitted_model")
    assert resolve_window_days(submitted, probe) == 3650
    assert resolve_window_days(replace(submitted, window_days=30), probe) == 3650


def test_baselines_hold_on_a_seed_outside_the_gate(probe):
    assert VALIDATION_SEED not in gate_seeds(probe.id, probe.n_seeds)
    for name in probe.must_pass:
        outcome = run_probe(registry.find_model(name), probe, [VALIDATION_SEED])
        assert outcome.verdict == PASS, (name, outcome.error or outcome.failing)
    for name, criterion in probe.must_fail.items():
        outcome = run_probe(registry.find_model(name), probe, [VALIDATION_SEED])
        assert outcome.verdict == FAIL and criterion in outcome.failing, (name, outcome.failing)

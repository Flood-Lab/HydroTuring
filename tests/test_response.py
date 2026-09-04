"""Tests for the driver-ablation probe: same rain, warmer and cooler air.

The expectation is a sign, not a closure: more evaporative demand over the
same water must leave less runoff, and less demand more. These check that
the generator changes only the driver, in both directions, that the
criterion reads the sign and the size of the response per unit of demand
and scores each direction on its own, that the exact model responds the
way physics says while a temperature-blind or one-way model does not, and
that a submitted model's window is widened to the year the expectation
needs.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from hydroturing import registry
from hydroturing.harness import build_case, resolve_window_days, run_probe
from hydroturing.scoring import FAIL, PASS
from hydroturing.seeds import gate_seeds


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/warming-response")


@pytest.mark.parametrize("variant, shift", [("warmer", 3.0), ("cooler", -3.0)])
def test_only_the_driver_changes(probe, variant, shift):
    control = build_case(probe, 7, "control")
    other = build_case(probe, 7, variant)
    assert control.forcing["pr"].equals(other.forcing["pr"])
    assert control.forcing["time"].equals(other.forcing["time"])

    spinup = control.spinup_steps
    assert control.forcing["tas"][:spinup].equals(other.forcing["tas"][:spinup])
    delta = (other.forcing["tas"] - control.forcing["tas"])[spinup:]
    np.testing.assert_allclose(delta, shift, atol=1e-5)
    pet_delta = (other.forcing["pet"] - control.forcing["pet"])[spinup:]
    assert (np.sign(pet_delta) == np.sign(shift)).mean() > 0.95
    # A rain-dominated catchment: snow timing must not confound the sign.
    assert (other.forcing["tas"] < 0).mean() < 0.02


def test_exact_model_responds_the_right_way_in_both_directions(probe):
    outcome = run_probe(registry.find_model("reference_bucket"), probe, gate_seeds(probe.id, 2))
    assert outcome.verdict == PASS, outcome.failing
    response = next(c for c in outcome.criteria if c.name == "response_sign")
    variants = response.diagnostics["variants"]
    assert set(variants) == {"warmer", "cooler"}
    assert variants["warmer"]["driver_change_mm"] > 1000
    assert variants["cooler"]["driver_change_mm"] < -1000
    for name, v in variants.items():
        shares = v["shares"]
        # Per unit of demand the derivative has the same sign both ways.
        assert shares["mrro"] < -0.1, name
        assert shares["evspsbl"] > 0.1, name
        # Over ten years storage change is small, so what evaporates is
        # what stopped running off, and vice versa.
        assert abs(shares["mrro"] + shares["evspsbl"]) < 0.1, name
    assert "both directions" in response.message


def test_temperature_blind_models_are_caught(probe):
    """A model that never reads the temperature cannot respond to it, and a
    model that evaporates everything responds by exactly nothing."""
    for name in ("reference_streamflow_only", "reference_degenerate"):
        outcome = run_probe(registry.find_model(name), probe, gate_seeds(probe.id, 1))
        assert outcome.verdict == FAIL
        assert "response_sign" in outcome.failing, (name, outcome.failing)
        response = next(c for c in outcome.criteria if c.name == "response_sign")
        assert abs(response.value) < 1e-6, name


def test_a_one_way_model_is_caught(probe):
    """A model that responds to warming but not to cooling has learned a
    one-way rule, and the cooling variant alone must fail it."""
    from hydroturing.criteria import get
    from hydroturing.protocol import RunResult
    from hydroturing.runner import get_runner
    import tempfile
    from pathlib import Path

    bucket = registry.find_model("reference_bucket")
    seed = gate_seeds(probe.id, 1)[0]
    runs = {}
    for variant in probe.variants:
        case = build_case(probe, seed, variant)
        runs[variant] = get_runner(bucket).run(bucket, probe, case, Path(tempfile.mkdtemp()))
    params = next(c.params for c in probe.criteria if c.name == "response_sign")
    assert get("response_sign")(runs, probe, dict(params)).passed

    # Cooling changes nothing: the cooler run reports the control's answer.
    runs["cooler"] = RunResult(runs["cooler"].case, runs["control"].table.copy(), runs["cooler"].meta, 0.0)
    result = get("response_sign")(runs, probe, dict(params))
    assert not result.passed
    assert "cooler" in result.message and "warmer" not in result.message
    assert result.value == pytest.approx(0.0)


def test_a_response_larger_than_the_demand_is_a_failure(probe):
    """Losing more runoff than the warming asked for means the model invented
    a loss. Built by scaling the exact model's response."""
    from hydroturing.criteria import get
    from hydroturing.protocol import RunResult

    from hydroturing.runner import get_runner
    import tempfile
    from pathlib import Path

    bucket = registry.find_model("reference_bucket")
    seed = gate_seeds(probe.id, 1)[0]
    runs = {}
    for variant in probe.variants:
        case = build_case(probe, seed, variant)
        runs[variant] = get_runner(bucket).run(bucket, probe, case, Path(tempfile.mkdtemp()))

    params = next(c.params for c in probe.criteria if c.name == "response_sign")
    assert get("response_sign")(runs, probe, dict(params)).passed

    # Exaggerate: the warmer run loses eight times what the bucket lost,
    # which is more water than the warming asked for.
    control_q = runs["control"].table["mrro"].to_numpy()
    warmer_q = runs["warmer"].table["mrro"].to_numpy()
    exaggerated = runs["warmer"].table.copy()
    exaggerated["mrro"] = np.maximum(control_q + 8.0 * (warmer_q - control_q), 0.0)
    runs["warmer"] = RunResult(runs["warmer"].case, exaggerated, runs["warmer"].meta, 0.0)
    result = get("response_sign")(runs, probe, dict(params))
    assert not result.passed
    assert "more than" in result.message


def test_submitted_models_get_at_least_a_year(probe):
    """Inside a month the sign can be dominated by storage. The probe asks
    for a year and a submitted model's window is widened to it."""
    submitted = replace(registry.find_model("reference_bucket"), name="submitted_model")
    assert resolve_window_days(submitted, probe) == 365
    assert resolve_window_days(replace(submitted, window_days=30), probe) == 365
    assert resolve_window_days(replace(submitted, window_days=730), probe) == 730
    assert resolve_window_days(submitted, probe, override="full") is None

    outcome = run_probe(submitted, probe, [5])
    assert outcome.window_days == 365
    assert outcome.windows[0]["rows"] == 365
    assert outcome.verdict == PASS, outcome.failing

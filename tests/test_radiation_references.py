"""Adapter parity and radiation margins on two gate seeds and emissivity bounds."""

from __future__ import annotations

import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.harness import build_case
from hydroturing.runner import get_runner
from hydroturing.seeds import gate_seeds

KELVIN = 273.15
# The design margin: each negative control's worst step must exceed the
# bound by this factor on every seed, so that the verdict does not rest on
# a step that sits at the tolerance. The criterion's own threshold stays 1.
MARGIN = 1.3


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("energy/radiation-consistency")


def run(name, probe, case):
    # Read back what the model declares, as verify-adapter does, so that the
    # coupled reference can be run on a probe that asks for more than it has.
    model = registry.find_model(name)
    asks = replace(
        probe, requires_fluxes=model.emits_fluxes, requires_states=model.emits_states,
        requires_diagnostics=model.emits_diagnostics,
    )
    with tempfile.TemporaryDirectory(prefix="ht-radiation-test-") as workdir:
        return get_runner(model).run(model, asks, case, Path(workdir))


def score(result, probe):
    return get("radiative_identity")(result, probe, dict(probe.criteria[0].params))


def cases(probe, eps=None):
    for seed in gate_seeds(probe.id, 2):
        case = build_case(probe, seed)
        if eps is not None:
            case.static["eps"] = eps
        yield case


def test_positive_control_is_the_coupled_reference_with_a_skin(probe):
    for case in cases(probe):
        skin = run("reference_radiative", probe, case).table
        coupled = run("reference_coupled", probe, case).table
        assert set(skin.columns) - set(coupled.columns) == {"rlus", "ts"}
        for column in coupled.columns:
            assert np.array_equal(skin[column].to_numpy(), coupled[column].to_numpy()), column
        ts = skin["ts"].to_numpy()
        assert np.isfinite(ts).all()
        assert ts.min() > KELVIN, "the snow-free boundary needs a skin that never freezes"
        assert ts.max() < KELVIN + 60.0


@pytest.mark.parametrize("negative", ["reference_air_emitter", "reference_no_reflection"])
def test_each_negative_control_differs_in_upward_longwave_alone(probe, negative):
    for case in cases(probe):
        positive = run("reference_radiative", probe, case).table
        broken = run(negative, probe, case).table
        assert list(broken.columns) == list(positive.columns)
        for column in positive.columns:
            same = np.array_equal(positive[column].to_numpy(), broken[column].to_numpy())
            assert same == (column != "rlus"), column


@pytest.mark.parametrize("eps", [None, 0.95, 0.99], ids=["drawn", "eps=0.95", "eps=0.99"])
def test_the_probe_separates_the_controls_with_margin(probe, eps):
    """At the drawn emissivity and at both ends of the range the case can
    draw. The reflected term the second control drops shrinks with (1 - eps),
    so 0.99 is where its margin is thinnest."""
    for case in cases(probe, eps):
        result = score(run("reference_radiative", probe, case), probe)
        assert result.status == PASS
        assert result.diagnostics["worst_slack"] < 1e-9
        for negative in ("reference_air_emitter", "reference_no_reflection"):
            result = score(run(negative, probe, case), probe)
            assert result.status == FAIL, negative
            assert result.diagnostics["worst_slack"] >= MARGIN, (negative, case.static["eps"])
        # Dropping the reflected sky is wrong by (1 - eps) * rlds on every
        # step, so the second control must violate at every scored step.
        assert result.diagnostics["violating_steps"] == result.diagnostics["scored_steps"]

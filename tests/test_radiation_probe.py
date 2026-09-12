"""Radiation forcing, reproducibility, staging and evaluation-window checks."""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.harness import (
    build_case,
    compatibility_issues,
    load_generator,
    resolve_window_days,
    run_probe,
    select_window,
    verify_adapter_contract,
    window_case,
)
from hydroturing.protocol import FORCING_FILE, STATIC_FILE, stage
from hydroturing.scoring import INCOMPATIBLE, NOT_SCORED
from hydroturing.seeds import gate_seeds

SIGMA = 5.670374419e-8
COLUMNS = ["time", "pr", "tas", "pet", "rn", "rlds"]


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("energy/radiation-consistency")


def test_spec_asks_for_the_identity_and_nothing_else(probe):
    # Every extra required variable excludes models from being testable at
    # all, so the requirement is pinned to the two outputs the identity needs.
    assert probe.requires_fluxes == ("rlus",)
    assert probe.requires_states == ()
    assert probe.requires_diagnostics == ("ts",)
    assert [c.name for c in probe.criteria] == ["radiative_identity"]
    assert probe.criteria[0].params["emissivity"] == "eps"
    # The verdict rests on the case's sky and emissivity, so a model must
    # declare that it consumes them.
    assert probe.requires_forcing == ("rlds",)
    assert probe.requires_static == ("eps",)


def test_adapter_verification_asks_only_what_the_adapter_needs(probe):
    """The coupled reference reads neither rlds nor eps, so it cannot be
    judged here, but its adapter can still be smoke-tested on this case."""
    result = verify_adapter_contract(
        registry.find_model("reference_coupled"), probe, gate_seeds(probe.id, 1)[0]
    )
    assert len(result.table) == 768


@pytest.mark.parametrize("missing", [("rlds",), ("eps",), ("rlds", "eps")])
def test_a_model_that_does_not_consume_the_sky_or_the_emissivity_is_not_judged(probe, missing):
    """Emitting ts and rlus is not enough: a model with its own downward
    longwave or its own emissivity would be scored against values it never
    read, so it is INCOMPATIBLE rather than a candidate for VIOLATION."""
    reference = registry.find_model("reference_radiative")
    assert compatibility_issues(reference, probe, build_case(probe, 0)) == []
    own_inputs = replace(
        reference,
        uses_forcing=() if "rlds" in missing else ("rlds",),
        uses_static=() if "eps" in missing else ("eps",),
    )
    outcome = run_probe(own_inputs, probe, gate_seeds(probe.id, 1))
    assert (outcome.verdict, outcome.reason) == (NOT_SCORED, INCOMPATIBLE)
    expected = []
    if "rlds" in missing:
        expected.append("model does not declare that it consumes forcing rlds")
    if "eps" in missing:
        expected.append("model does not declare that it consumes static eps")
    assert outcome.incompatible == expected


@pytest.mark.parametrize("forcing_kind", ["needs", "uses"])
@pytest.mark.parametrize("static_kind", ["needs", "uses"])
def test_required_inputs_accept_mandatory_or_optional_declarations(probe, forcing_kind, static_kind):
    model = registry.find_model("reference_radiative")
    inputs = dict(needs_forcing=model.needs_forcing, needs_static=(), uses_forcing=(), uses_static=())
    inputs[f"{forcing_kind}_forcing"] += ("rlds",)
    inputs[f"{static_kind}_static"] += ("eps",)
    model = replace(model, **inputs)
    assert compatibility_issues(model, probe, build_case(probe, 0)) == []

    other = registry.find_probe("energy/surface-energy-closure")
    expected = []
    if forcing_kind == "needs":
        expected.append("forcing does not provide rlds")
    if static_kind == "needs":
        expected.append("static does not provide eps")
    assert compatibility_issues(model, other, build_case(other, 4242)) == expected


def test_generator_reproduces_bytes_and_changes_the_case_with_seed(probe):
    generate = load_generator(probe).generate
    first, static = generate(0)
    repeated, repeated_static = generate(0)
    other, other_static = generate(1)
    assert first.to_csv(index=False) == repeated.to_csv(index=False)
    assert static == repeated_static
    for column in ("pr", "tas", "rn", "rlds"):
        assert not first[column].equals(other[column]), column
    assert static["eps"] != other_static["eps"]


@pytest.mark.parametrize("seed", [0, 19])
def test_case_is_hourly_and_starts_at_dawn(probe, seed):
    case = build_case(probe, seed)
    assert case.n_steps == 768
    assert case.spinup_steps == 48
    assert case.timestep == "PT1H"
    assert list(case.forcing.columns) == COLUMNS
    time = pd.to_datetime(case.forcing["time"])
    assert time.diff().iloc[1:].eq(pd.Timedelta(hours=1)).all()
    assert time.iloc[0].hour == 6
    assert len(case.after_spinup(case.forcing)) == 720


def test_every_seed_stays_in_the_warm_gray_surface_regime(probe):
    """Check input bounds on the gate seeds and twenty additional seeds."""
    generate = load_generator(probe).generate
    for seed in [*gate_seeds(probe.id, probe.n_seeds), *range(20)]:
        frame, static = generate(seed)
        where = f"seed {seed}"
        assert 0.95 <= static["eps"] <= 0.99, where
        assert static["canopy_capacity_mm"] == 0.0, where
        assert np.isfinite(frame[COLUMNS[1:]].to_numpy()).all(), where

        # Warm enough that nothing freezes, which the snow-free boundary needs.
        assert frame["tas"].between(15.0, 31.0).all(), where
        assert (frame["pr"] >= 0.0).all() and (frame["pet"] > 0.0).all(), where
        assert frame["rn"].between(-45.0, 560.0).all(), where

        # The sky radiates as a gray body at the air temperature, between a
        # clear and an overcast emissivity, and varies hour to hour.
        air_k = frame["tas"].to_numpy() + 273.15
        sky = frame["rlds"].to_numpy() / (SIGMA * air_k**4)
        assert np.all((sky > 0.70) & (sky < 0.97)), where
        assert frame["rlds"].between(280.0, 470.0).all(), where
        assert frame["rlds"].std() > 20.0, where

        # The reflected term the second negative control drops is well above
        # the criterion's 0.5 W m-2 floor at every hour, at any emissivity
        # the case can draw. Whether it also clears the relative bound is a
        # question for the reference models' actual fluxes, not for the input.
        assert ((1.0 - static["eps"]) * frame["rlds"] >= 2.5).all(), where


def test_the_same_cloud_dims_the_sun_and_brightens_the_sky(probe):
    """At noon, lower net radiation must go with a more emissive sky."""
    frame = build_case(probe, 3).forcing
    noon = frame[pd.to_datetime(frame["time"]).dt.hour == 12]
    air_k = noon["tas"].to_numpy() + 273.15
    sky = noon["rlds"].to_numpy() / (SIGMA * air_k**4)
    assert len(noon) == 32
    assert np.corrcoef(noon["rn"].to_numpy(), sky)[0, 1] < -0.5


def test_staging_hands_the_model_the_sky_and_the_emissivity(probe, tmp_path):
    case = build_case(probe, 0)
    model = registry.find_model("reference_bucket")
    request = json.loads(stage(tmp_path, case, probe, model).read_text())
    staged = pd.read_csv(tmp_path / FORCING_FILE)
    pd.testing.assert_frame_equal(staged, case.forcing)
    static = json.loads((tmp_path / STATIC_FILE).read_text())
    assert static["eps"] == case.static["eps"]
    assert request["timestep"] == "PT1H"
    assert request["n_steps"] == 768
    assert "ts and rlus are instantaneous values at row i's time" in request["notes"]
    assert "the same instant as row i's rlds" in request["notes"]


def test_submitted_models_see_every_scored_day(probe):
    # A submitted hourly model normally gets seven days; this probe widens
    # that to the whole scored record, even when a caller asks for less.
    model = replace(registry.find_model("reference_bucket"), name="submitted_skin")
    assert resolve_window_days(model, probe) == 30
    assert resolve_window_days(model, probe, override=7) == 30
    case = build_case(probe, 0)
    cut = window_case(case, select_window(case, probe, 30))
    assert cut.spinup_steps == 48
    assert cut.window["rows"] == 720
    pd.testing.assert_frame_equal(cut.forcing, case.forcing)

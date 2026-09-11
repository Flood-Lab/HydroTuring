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
    load_generator,
    resolve_window_days,
    select_window,
    window_case,
)
from hydroturing.protocol import FORCING_FILE, STATIC_FILE, stage
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

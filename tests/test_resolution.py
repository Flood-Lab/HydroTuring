"""Tests for running a probe at more than one step.

A resolution transform is the same weather at two steps. These check that a
probe can declare a step per variant, that a model declares the steps it can
run at and is INCOMPATIBLE otherwise, that the record length and the window
are counted in time rather than rows, and that the reference bucket is the
same model at every step while the fixed-step one is not.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from hydroturing import registry
from hydroturing.harness import (
    build_case,
    compatibility_issues,
    run_probe,
    select_window,
    window_case,
)
from hydroturing.scoring import FAIL, INCOMPATIBLE, PASS
from hydroturing.seeds import gate_seeds
from hydroturing.spec import TIMESTEP_DAYS, SpecError, load_model


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/resolution-invariance")


def test_probe_declares_a_step_per_variant(probe):
    assert probe.variants == ("minute", "hourly")
    assert probe.timestep == "PT1M"
    assert probe.timestep_for("minute") == "PT1M"
    assert probe.timestep_for("hourly") == "PT1H"
    assert probe.timesteps == ("PT1M", "PT1H")
    assert probe.period_days == 30
    # The same days, counted at each step.
    assert probe.n_steps_for("minute") == 40 * 1440
    assert probe.n_steps_for("hourly") == 40 * 24
    assert probe.spinup_steps_for("minute") == 10 * 1440
    assert probe.spinup_steps_for("hourly") == 10 * 24


def test_variants_carry_the_same_water(probe):
    """Draw once, then aggregate: every hour holds exactly its sixty minutes."""
    fine = build_case(probe, 7, "minute")
    coarse = build_case(probe, 7, "hourly")
    assert fine.timestep == "PT1M" and coarse.timestep == "PT1H"
    assert fine.n_steps == 60 * coarse.n_steps

    minutes = fine.forcing["pr"].to_numpy().reshape(-1, 60)
    hours = coarse.forcing["pr"].to_numpy()
    np.testing.assert_allclose(minutes.mean(axis=1), hours, atol=2e-6)
    total_fine = fine.forcing["pr"].sum() * fine.dt_days
    total_coarse = coarse.forcing["pr"].sum() * coarse.dt_days
    assert abs(total_fine - total_coarse) < 1e-3 * total_fine
    # And the storms have structure inside the hour, or there is nothing to test.
    assert fine.forcing["pr"].max() > 2 * coarse.forcing["pr"].max()


def test_manifest_timestep_may_be_a_list(tmp_path):
    import shutil

    target = tmp_path / "submitted_model"
    shutil.copytree(registry.MODELS_DIR / "_template", target)
    manifest = target / "model.yaml"
    text = manifest.read_text().replace("name: _template", "name: submitted_model")

    manifest.write_text(text.replace("timestep: PT1D", "timestep: [PT1H, PT1M]"))
    model = load_model(target)
    assert model.timesteps == ("PT1H", "PT1M")
    assert model.timestep == "PT1H"
    assert model.supports_timestep("PT1M") and not model.supports_timestep("PT1D")

    manifest.write_text(text.replace("timestep: PT1D", "timestep: [PT1H, PT1H]"))
    with pytest.raises(SpecError):
        load_model(target)


def test_single_step_model_is_incompatible_not_wrong(probe):
    """A model that only runs at one resolution says so and is not run."""
    daily = registry.find_model("reference_leaky")
    assert daily.timesteps == ("PT1D",)
    outcome = run_probe(daily, probe, [11])
    assert outcome.reason == INCOMPATIBLE
    assert "PT1M" in outcome.incompatible[0]

    hourly_only = replace(registry.find_model("reference_bucket"), timesteps=("PT1H",))
    issues = compatibility_issues(hourly_only, probe)
    assert issues and "PT1M" in issues[0] and "PT1H" not in issues[0].split("cover")[1]


def test_reference_bucket_is_the_same_model_at_every_step(probe):
    seeds = gate_seeds(probe.id, 1)
    outcome = run_probe(registry.find_model("reference_bucket"), probe, seeds)
    assert outcome.verdict == PASS, outcome.failing
    invariance = next(c for c in outcome.criteria if c.name == "resolution_invariance")
    assert invariance.value < 0.01
    assert invariance.diagnostics["fine_step"] == "PT1M"
    assert invariance.diagnostics["coarse_step"] == "PT1H"


def test_fixed_step_model_is_caught(probe):
    """The failure the probe exists for: arithmetic that assumes a step."""
    outcome = run_probe(registry.find_model("reference_fixed_step"), probe, gate_seeds(probe.id, 1))
    assert outcome.verdict == FAIL
    assert outcome.failing == ["resolution_invariance"]
    invariance = next(c for c in outcome.criteria if c.name == "resolution_invariance")
    # Far past the 5 percent limit on every seed; the smallest gate seed is 0.198.
    assert invariance.value > 0.1


def test_window_is_the_same_stretch_of_time_at_both_steps(probe):
    """A submitted model gets a window; both variants must be cut to the
    same days even though one has sixty times the rows of the other."""
    fine = build_case(probe, 5, "minute")
    coarse = build_case(probe, 5, "hourly")
    bounds = select_window(fine, probe, 7)
    cut_fine = window_case(fine, bounds)
    cut_coarse = window_case(coarse, bounds)

    assert cut_fine.window["days"] == cut_coarse.window["days"] == 7
    assert cut_fine.window["rows"] == 7 * 1440
    assert cut_coarse.window["rows"] == 7 * 24
    assert cut_fine.window["start"][:13] == cut_coarse.window["start"][:13]
    assert cut_fine.n_steps == 17 * 1440 and cut_coarse.n_steps == 17 * 24

    submitted = replace(registry.find_model("reference_bucket"), name="submitted_model")
    outcome = run_probe(submitted, probe, [5])
    assert outcome.window_days == 7
    assert outcome.verdict == PASS, outcome.failing


def test_step_lengths_are_consistent():
    assert TIMESTEP_DAYS["PT1H"] * 24 == pytest.approx(1.0)
    assert TIMESTEP_DAYS["PT1M"] * 1440 == pytest.approx(1.0)
    assert TIMESTEP_DAYS["PT5M"] * 288 == pytest.approx(1.0)
    assert TIMESTEP_DAYS["PT15M"] * 96 == pytest.approx(1.0)

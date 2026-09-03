"""Tests for running a probe at more than one step.

A resolution transform is the same weather at two or more steps. These
check that a probe can declare a step per variant, that each model is run
at its own step and the next finer one whatever its manifest claims, that
record length and window are counted in time rather than rows, and that the
reference bucket is the same model at every step while the fixed-step one
is not, by a measured amount.
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
    select_variants,
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
    assert probe.variants == ("minute", "hourly", "daily")
    assert probe.timestep == "PT1M"
    assert probe.timestep_for("hourly") == "PT1H"
    assert probe.timestep_for("daily") == "PT1D"
    assert probe.timesteps == ("PT1M", "PT1H", "PT1D")
    assert probe.period_days == 30
    # The same days, counted at each step.
    assert probe.n_steps_for("minute") == 40 * 1440
    assert probe.n_steps_for("hourly") == 40 * 24
    assert probe.n_steps_for("daily") == 40
    assert probe.spinup_steps_for("hourly") == 10 * 24


def test_variants_carry_the_same_water(probe):
    """Draw once, then aggregate: every hour holds exactly its sixty minutes
    and every day its twenty-four hours."""
    fine = build_case(probe, 7, "minute")
    hourly = build_case(probe, 7, "hourly")
    daily = build_case(probe, 7, "daily")
    assert (fine.timestep, hourly.timestep, daily.timestep) == ("PT1M", "PT1H", "PT1D")
    assert fine.n_steps == 60 * hourly.n_steps == 1440 * daily.n_steps

    np.testing.assert_allclose(
        fine.forcing["pr"].to_numpy().reshape(-1, 60).mean(axis=1),
        hourly.forcing["pr"].to_numpy(), atol=2e-6,
    )
    np.testing.assert_allclose(
        hourly.forcing["pr"].to_numpy().reshape(-1, 24).mean(axis=1),
        daily.forcing["pr"].to_numpy(), atol=2e-6,
    )
    # And the storms have structure inside the hour, or there is nothing to test.
    assert fine.forcing["pr"].max() > 2 * hourly.forcing["pr"].max()


def test_each_model_is_run_at_its_own_step_and_the_next_finer(probe):
    """The manifest says where a model lives; the probe decides what to run
    and does not consult the manifest's list of other steps at all."""
    daily = registry.find_model("reference_leaky")
    assert daily.timesteps == ("PT1D",)
    assert select_variants(daily, probe) == ["daily", "hourly"]

    hourly = replace(daily, timesteps=("PT1H",))
    assert select_variants(hourly, probe) == ["hourly", "minute"]

    # The finest step has nothing finer and pairs with the next coarser.
    minute = replace(daily, timesteps=("PT1M",))
    assert select_variants(minute, probe) == ["minute", "hourly"]

    # A step the probe does not serve is rounded to the nearest coarser one.
    quarter = replace(daily, timesteps=("PT15M",))
    assert select_variants(quarter, probe) == ["hourly", "minute"]

    # No step gating on a multi-step probe, whatever the manifest declares.
    assert not compatibility_issues(daily, probe)
    assert not compatibility_issues(daily, probe, build_case(probe, 7, "hourly"))

    # An ordinary paired probe still runs every variant, control first.
    closure = registry.find_probe("mass/catchment-closure")
    assert select_variants(daily, closure) == [None]


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

    manifest.write_text(text.replace("timestep: PT1D", "timestep: [PT1H, PT1H]"))
    with pytest.raises(SpecError):
        load_model(target)


def test_reference_bucket_is_the_same_model_at_every_step(probe):
    seeds = gate_seeds(probe.id, 1)
    outcome = run_probe(registry.find_model("reference_bucket"), probe, seeds)
    assert outcome.verdict == PASS, outcome.failing
    invariance = next(c for c in outcome.criteria if c.name == "resolution_invariance")
    # A lumped model integrated at a daily step differs from itself at an
    # hourly step by the discretisation of its own equations, about two
    # percent of the rain here. That is the floor the 5 percent rule sits on.
    assert invariance.value < 0.03
    # The bucket lives at a daily step, so it is compared daily against hourly.
    assert invariance.diagnostics["steps"] == ["PT1H", "PT1D"]

    # At its finest pair the same model agrees far more closely.
    minute_native = replace(registry.find_model("reference_bucket"), timesteps=("PT1M", "PT1H", "PT1D"))
    outcome = run_probe(minute_native, probe, seeds)
    invariance = next(c for c in outcome.criteria if c.name == "resolution_invariance")
    assert invariance.diagnostics["steps"] == ["PT1M", "PT1H"]
    assert invariance.value < 0.005


def test_fixed_step_model_is_measured_not_declared(probe):
    """The failure the probe exists for, as a number: arithmetic that assumes
    a daily step, run on hourly rows, moves the month's runoff by a large
    share of the rain."""
    outcome = run_probe(registry.find_model("reference_fixed_step"), probe, gate_seeds(probe.id, 1))
    assert outcome.verdict == FAIL
    assert "resolution_invariance" in outcome.failing
    invariance = next(c for c in outcome.criteria if c.name == "resolution_invariance")
    assert invariance.value > 0.1
    assert "PT1H" in invariance.message and "PT1D" in invariance.message


def test_a_model_that_declares_one_step_is_run_anyway(probe):
    """The old bucket knows nothing about steps and says it runs daily. It is
    run at hourly regardless and the disagreement is measured."""
    outcome = run_probe(registry.find_model("reference_leaky"), probe, [11])
    assert outcome.verdict == FAIL
    assert "resolution_invariance" in outcome.failing
    invariance = next(c for c in outcome.criteria if c.name == "resolution_invariance")
    assert invariance.value > 0.05
    assert not outcome.incompatible


def test_streamflow_only_models_are_quantified(probe):
    """Runoff is enough to be scored here, so a model that reports nothing
    else is measured rather than stopped at INCOMPLETE."""
    outcome = run_probe(registry.find_model("reference_streamflow_only"), probe, [11])
    assert not outcome.missing
    assert outcome.criteria, "the model was not run"
    invariance = next(c for c in outcome.criteria if c.name == "resolution_invariance")
    assert invariance.value is not None
    assert set(invariance.diagnostics["deviations"]) == {"mrro PT1D"}


def test_step_mismatch_on_a_single_step_probe_is_still_incompatible():
    """Closure at a daily step measures closure. A model that cannot be fed
    daily rows has not violated closure; the probe cannot say."""
    closure = registry.find_probe("mass/catchment-closure")
    hourly_only = replace(registry.find_model("reference_bucket"), timesteps=("PT1H",))
    outcome = run_probe(hourly_only, closure, [11])
    assert outcome.reason == INCOMPATIBLE
    assert not outcome.criteria


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
    assert cut_fine.window["start"][:10] == cut_coarse.window["start"][:10]
    assert cut_fine.n_steps == 17 * 1440 and cut_coarse.n_steps == 17 * 24

    # A submitted daily model gets the daily default window, not the minute
    # one, so its scored stretch has days in it rather than a week of rows.
    submitted = replace(registry.find_model("reference_bucket"), name="submitted_model")
    outcome = run_probe(submitted, probe, [5])
    assert outcome.window_days == 30
    assert outcome.verdict == PASS, outcome.failing

    hourly_native = replace(submitted, timesteps=("PT1H",))
    outcome = run_probe(hourly_native, probe, [5])
    assert outcome.window_days == 7
    assert outcome.windows[0]["rows"] == 7 * 24
    assert outcome.verdict == PASS, outcome.failing


def test_step_lengths_are_consistent():
    assert TIMESTEP_DAYS["PT1H"] * 24 == pytest.approx(1.0)
    assert TIMESTEP_DAYS["PT1M"] * 1440 == pytest.approx(1.0)
    assert TIMESTEP_DAYS["PT5M"] * 288 == pytest.approx(1.0)
    assert TIMESTEP_DAYS["PT15M"] * 96 == pytest.approx(1.0)

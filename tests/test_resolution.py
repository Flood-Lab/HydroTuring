"""Tests for running a probe at more than one step.

A resolution transform is the same weather at two or more steps. These
check that a probe can declare a step per variant, that each model is run
at its own step and the next finer one whatever its manifest claims, that
record length and window are counted in time rather than rows, that the
reference bucket is the same model at every step while the fixed-step one
is not, by a measured amount, and that a model whose runoff depends on rain
intensity is not failed for that physics.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.harness import (
    build_case,
    compatibility_issues,
    load_generator,
    run_probe,
    select_variants,
    select_window,
    window_case,
)
from hydroturing.protocol import Case, RunResult
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


def test_variants_carry_the_same_water_at_the_same_intensity(probe):
    """One daily record held flat: every hour and every minute of a day
    carries the day's rates, so the steps differ in the length of the row
    and in nothing else."""
    fine = build_case(probe, 7, "minute")
    hourly = build_case(probe, 7, "hourly")
    daily = build_case(probe, 7, "daily")
    assert (fine.timestep, hourly.timestep, daily.timestep) == ("PT1M", "PT1H", "PT1D")
    assert fine.n_steps == 60 * hourly.n_steps == 1440 * daily.n_steps

    for name in ("pr", "tas", "pet"):
        day = daily.forcing[name].to_numpy()
        np.testing.assert_array_equal(hourly.forcing[name].to_numpy(), np.repeat(day, 24))
        np.testing.assert_array_equal(fine.forcing[name].to_numpy(), np.repeat(day, 1440))
    # Rain falls on some days and not on others, or there is nothing to test.
    assert (daily.forcing["pr"] > 0).sum() >= 5 and (daily.forcing["pr"] == 0).sum() >= 5


def _horton(case: Case, capacity_mm_per_hour: float = 2.0) -> RunResult:
    """A step-aware bucket whose surface runoff is infiltration excess: rain
    arriving faster than the soil takes it runs off at once, and the rest
    fills a 120 mm store that spills when full and drains slowly."""
    dt = case.dt_days
    pr = case.forcing["pr"].to_numpy(dtype=float)
    pet = case.forcing["pet"].to_numpy(dtype=float)
    capacity = capacity_mm_per_hour * 24.0
    store, store_max, recession = 60.0, 120.0, 0.006
    runoff = np.empty(len(pr))
    for i, (rain, demand) in enumerate(zip(pr, pet)):
        infiltration = min(rain, capacity) * dt
        surface = rain * dt - infiltration
        store += infiltration
        spill = max(0.0, store - store_max)
        store -= spill
        store -= min(store, demand * dt * store / store_max)
        drain = (1.0 - (1.0 - recession) ** dt) * store
        store -= drain
        runoff[i] = (surface + spill + drain) / dt
    table = pd.DataFrame({"time": case.forcing["time"], "pr": pr, "mrro": runoff})
    return RunResult(case=case, table=table, meta={}, wall_seconds=0.0)


def _as_drawn(probe, case: Case) -> Case:
    """The case with version 2's weather: the minutes as drawn, storms in
    five-minute bursts, averaged over the case's rows."""
    pr, tas, pet = load_generator(probe)._draw_minutes(np.random.default_rng(case.seed))
    per_row = round(case.dt_days * 1440)
    forcing = case.forcing.copy()
    for name, values in (("pr", pr), ("tas", tas), ("pet", pet)):
        forcing[name] = values.reshape(-1, per_row).mean(axis=1)
    return replace(case, forcing=forcing)


@pytest.mark.parametrize("seed", [1, 3, 11])
def test_infiltration_excess_is_not_failed_for_its_physics(probe, seed):
    """Averaging rain lowers its intensity, and infiltration excess depends
    on intensity, so a correct Horton model partitions the same water
    differently at each step when the storms keep their structure. Version
    2 served them that way and failed such a model for its physics: on
    these seeds a daily Horton bucket with a 2 mm/h capacity moves 21 to 35
    percent of the rain between hours and days. Held flat, the intensity is
    the same at every step and the model agrees with itself to about a
    percent, the discretisation of its own daily step. Over 200 seeds the
    version 2 weather failed it on 87 and this version on none."""
    criterion = get("resolution_invariance")
    params = dict(next(c.params for c in probe.criteria if c.name == "resolution_invariance"))
    cases = {v: build_case(probe, seed, v) for v in ("daily", "hourly")}

    held = criterion({v: _horton(c) for v, c in cases.items()}, probe, params)
    assert held.passed, held.message
    assert held.value < 0.02

    drawn = criterion({v: _horton(_as_drawn(probe, c)) for v, c in cases.items()}, probe, params)
    assert not drawn.passed
    assert drawn.value > 0.2


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
    # hourly step by the discretisation of its own equations, about one
    # percent of the rain here. That is the floor the 10 percent limit sits on.
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

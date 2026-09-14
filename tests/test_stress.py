"""Tests for the stress probes: causality, dry-down, steady state, extreme rain,
antecedent memory.

Each states a limit any correct physical model satisfies exactly. These
check that the generators build the case they claim, that the exact bucket
passes every one, that the broken baseline built for each is caught by the
criterion named for it, that the criteria read edge cases the way they
should, and that paired variants hand a stochastic model the same seed.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.harness import build_case, run_probe
from hydroturing.protocol import Case, RunResult, _opaque_case_metadata
from hydroturing.runner import get_runner
from hydroturing.scoring import FAIL, PASS
from hydroturing.seeds import gate_seeds


def _runs(probe, model_name, seed):
    model = registry.find_model(model_name)
    runs = {}
    for variant in probe.variants or (None,):
        case = build_case(probe, seed, variant)
        runs[variant or "_"] = get_runner(model).run(model, probe, case, Path(tempfile.mkdtemp()))
    return runs


def _params(probe, name):
    return dict(next(c.params for c in probe.criteria if c.name == name))


# --- causality -----------------------------------------------------------------


def test_causality_variants_differ_on_one_day():
    probe = registry.find_probe("mass/causality")
    control = build_case(probe, 3, "control")
    pulse = build_case(probe, 3, "pulse")
    diff = pulse.forcing["pr"].to_numpy() - control.forcing["pr"].to_numpy()
    (where,) = np.nonzero(diff)
    assert len(where) == 1 and diff[where[0]] == pytest.approx(40.0)
    assert where[0] > control.spinup_steps + 300
    assert control.forcing["tas"].equals(pulse.forcing["tas"])


def test_paired_variants_share_the_model_seed():
    """A stochastic model must draw the same numbers in both runs, or the
    comparison sees its sampling noise instead of the perturbation."""
    probe = registry.find_probe("mass/causality")
    control = build_case(probe, 3, "control")
    pulse = build_case(probe, 3, "pulse")
    assert _opaque_case_metadata(control) == _opaque_case_metadata(pulse)
    assert _opaque_case_metadata(control) != _opaque_case_metadata(build_case(probe, 4, "control"))


def test_exact_model_is_causal_and_a_smoother_is_not():
    probe = registry.find_probe("mass/causality")
    seeds = gate_seeds(probe.id, 1)
    good = run_probe(registry.find_model("reference_bucket"), probe, seeds)
    assert good.verdict == PASS, good.failing
    causal = next(c for c in good.criteria if c.name == "causality")
    assert causal.value == 0.0
    assert causal.diagnostics["response_share"] > 0.02

    bad = run_probe(registry.find_model("reference_anticipating"), probe, seeds)
    assert "causality" in bad.failing
    causal = next(c for c in bad.criteria if c.name == "causality")
    assert causal.value > 1e-3
    assert "before the storm" in causal.message


def test_causality_needs_a_response_after_the_cause():
    """Ignoring the rain altogether is not causal behaviour, it is no
    behaviour: the storm has to show up afterwards."""
    probe = registry.find_probe("mass/causality")
    runs = _runs(probe, "reference_bucket", gate_seeds(probe.id, 1)[0])
    runs["pulse"] = RunResult(runs["pulse"].case, runs["control"].table.copy(), runs["pulse"].meta, 0.0)
    result = get("causality")(runs, probe, _params(probe, "causality"))
    assert not result.passed and "at least" in result.message


# --- dry-down ------------------------------------------------------------------


def test_dry_down_record_is_rainless_after_spinup():
    probe = registry.find_probe("mass/dry-down")
    case = build_case(probe, 3)
    assert case.forcing["pr"][: case.spinup_steps].sum() > 500
    assert case.forcing["pr"][case.spinup_steps :].sum() == 0
    assert case.forcing["pet"][case.spinup_steps :].sum() > 1000


def test_exact_model_drains_and_a_climatology_keeps_flowing():
    probe = registry.find_probe("mass/dry-down")
    seeds = gate_seeds(probe.id, 1)
    good = run_probe(registry.find_model("reference_bucket"), probe, seeds)
    assert good.verdict == PASS, good.failing
    drain = next(c for c in good.criteria if c.name == "dry_down")
    assert 0 < drain.value < 1
    assert all(r <= 0 for r in drain.diagnostics["worst_block_rise"].values())

    bad = run_probe(registry.find_model("reference_climatology"), probe, seeds)
    assert "dry_down" in bad.failing
    drain = next(c for c in bad.criteria if c.name == "dry_down")
    assert drain.value > 1
    assert "more than" in drain.message and "rose" in drain.message


def test_dry_down_refuses_a_record_with_rain():
    closure = registry.find_probe("mass/catchment-closure")
    runs = _runs(closure, "reference_bucket", 3)
    with pytest.raises(ValueError, match="rainless"):
        get("dry_down")(runs["_"], closure, {})


# --- steady state ---------------------------------------------------------------


def test_steady_state_forcing_is_constant():
    probe = registry.find_probe("mass/steady-state")
    case = build_case(probe, 3)
    scored = case.forcing.iloc[case.spinup_steps :]
    for column in ("pr", "tas", "pet"):
        assert scored[column].nunique() == 1
    assert case.forcing["pr"][: case.spinup_steps].nunique() > 10


def test_exact_model_settles_and_a_restless_one_does_not():
    probe = registry.find_probe("mass/steady-state")
    seeds = gate_seeds(probe.id, 2)
    good = run_probe(registry.find_model("reference_bucket"), probe, seeds)
    assert good.verdict == PASS, good.failing
    settled = next(c for c in good.criteria if c.name == "steady_state")
    means = settled.diagnostics["means"]
    assert means["mrro"] <= 2.5 * 1.01
    assert abs(settled.diagnostics["budget_residual"]) < 0.025
    assert settled.value < 0.01

    bad = run_probe(registry.find_model("reference_restless"), probe, seeds)
    assert "steady_state" in bad.failing
    assert "still varies" in next(c for c in bad.criteria if c.name == "steady_state").message


def test_steady_state_catches_runoff_above_the_rain():
    probe = registry.find_probe("mass/steady-state")
    runs = _runs(probe, "reference_bucket", gate_seeds(probe.id, 1)[0])
    table = runs["_"].table.copy()
    table["mrro"] = 3.0  # more than the 2.5 mm/day falling, perfectly steady
    result = get("steady_state")(RunResult(runs["_"].case, table, {}, 0.0), probe, _params(probe, "steady_state"))
    assert not result.passed and "runs off" in result.message


# --- extreme rain ----------------------------------------------------------------


def test_extreme_rain_scales_one_storm():
    probe = registry.find_probe("mass/extreme-rain")
    control = build_case(probe, 3, "control")
    top = build_case(probe, 3, "x10")
    diff = top.forcing["pr"].to_numpy() - control.forcing["pr"].to_numpy()
    (where,) = np.nonzero(diff)
    assert len(where) == 1
    assert top.forcing["pr"].iloc[where[0]] == pytest.approx(10 * control.forcing["pr"].iloc[where[0]])
    assert top.forcing["pr"].max() > 500


def test_exact_model_runs_off_the_extra_and_a_capped_one_does_not():
    probe = registry.find_probe("mass/extreme-rain")
    seeds = gate_seeds(probe.id, 1)
    good = run_probe(registry.find_model("reference_bucket"), probe, seeds)
    assert good.verdict == PASS, good.failing
    ladder = next(c for c in good.criteria if c.name == "monotone_response")
    assert ladder.value > 0.5
    assert all(r["share"] >= 0 for r in ladder.diagnostics["rungs"])

    bad = run_probe(registry.find_model("reference_saturating"), probe, seeds)
    assert "monotone_response" in bad.failing
    ladder = next(c for c in bad.criteria if c.name == "monotone_response")
    assert ladder.value < 0.1


def test_extreme_rain_catches_amplification():
    """Runoff that grows faster than the rain added is water from nowhere."""
    probe = registry.find_probe("mass/extreme-rain")
    runs = _runs(probe, "reference_bucket", gate_seeds(probe.id, 1)[0])
    table = runs["x10"].table.copy()
    table["mrro"] = table["mrro"] * 3.0
    runs["x10"] = RunResult(runs["x10"].case, table, runs["x10"].meta, 0.0)
    result = get("monotone_response")(runs, probe, _params(probe, "monotone_response"))
    assert not result.passed and "more than" in result.message


# --- antecedent monotonicity ------------------------------------------------------


def _antecedent_pair(gap, memory, next_day=0.0):
    """A dry and a wet record around one 60 mm storm that falls `gap` rainless
    days after the wet record's 120 mm of antecedent rain. The wet run drains
    3 mm/day of that rain until the storm arrives, and then answers the storm
    with `memory` mm/day more than the dry run for five days. `next_day` mm on
    the day straight after the storm is part of it; a shower on the third
    day, past a dry one, is weather."""
    probe = registry.find_probe("mass/antecedent-monotonicity")
    n, first = 90, 10
    storm = first + 10 + gap
    time = pd.date_range("2001-06-01", periods=n, freq="D").strftime("%Y-%m-%d")
    runs = {}
    for variant in ("dry", "wet"):
        pr, mrro = np.zeros(n), np.zeros(n)
        pr[2], mrro[2] = 8.0, 1.0                      # both: rain long before
        pr[storm], mrro[storm:storm + 5] = 60.0, 6.0   # both: the same storm
        pr[storm + 1] = next_day                       # both: its second day, if any
        pr[storm + 3] = 5.0                            # both: a later shower
        if variant == "wet":
            pr[first:first + 10] = 12.0
            mrro[first:storm] = 3.0
            mrro[storm:storm + 5] += memory
        case = Case(probe.id, 0, pd.DataFrame({"time": time, "pr": pr}), {}, spinup_steps=0)
        runs[variant] = RunResult(case, pd.DataFrame({"time": time, "mrro": mrro}), {}, 0.0)
    return probe, runs, time[storm]


@pytest.mark.parametrize("gap", [0, 1, 10, 25])
def test_antecedent_window_opens_on_the_storm_not_when_the_antecedent_rain_stops(gap):
    """The wet run answers the storm exactly as the dry one does, so it has no
    memory at the storm. A window opened when the antecedent rain stopped
    would count its recession through the quiet days as storm runoff and
    pass it."""
    probe, runs, when = _antecedent_pair(gap, memory=0.0)
    result = get("antecedent_monotonicity")(runs, probe, _params(probe, "antecedent_monotonicity"))
    assert not result.passed and "at least" in result.message
    assert result.diagnostics["runoff_wet_mm"] == pytest.approx(result.diagnostics["runoff_dry_mm"])
    assert result.diagnostics["storm_time"] == when


@pytest.mark.parametrize("gap", [0, 1, 10, 25])
def test_antecedent_counts_only_runoff_from_the_storm_on(gap):
    probe, runs, when = _antecedent_pair(gap, memory=1.0)
    result = get("antecedent_monotonicity")(runs, probe, _params(probe, "antecedent_monotonicity"))
    assert result.passed, result.message
    assert result.diagnostics["storm_mm"] == pytest.approx(60.0)
    assert result.diagnostics["runoff_dry_mm"] == pytest.approx(30.0)
    assert result.diagnostics["runoff_wet_mm"] == pytest.approx(35.0)
    assert result.diagnostics["storm_time"] == when


def test_antecedent_storm_runs_to_its_first_dry_step():
    probe, runs, when = _antecedent_pair(10, memory=1.0, next_day=12.0)
    result = get("antecedent_monotonicity")(runs, probe, _params(probe, "antecedent_monotonicity"))
    assert result.diagnostics["storm_mm"] == pytest.approx(72.0)
    assert result.diagnostics["storm_time"] == when


def test_antecedent_window_is_never_shorter_than_one_step():
    probe, runs, _ = _antecedent_pair(10, memory=1.0)
    params = {**_params(probe, "antecedent_monotonicity"), "window_days": 0.25}
    result = get("antecedent_monotonicity")(runs, probe, params)
    assert result.diagnostics["storm_mm"] == pytest.approx(60.0)
    assert result.diagnostics["runoff_wet_mm"] - result.diagnostics["runoff_dry_mm"] == pytest.approx(1.0)


def test_antecedent_needs_a_storm_after_the_antecedent_rain():
    probe, runs, _ = _antecedent_pair(10, memory=0.0)
    for run in runs.values():
        run.case.forcing.loc[20:, "pr"] = 0.0
    with pytest.raises(ValueError, match="no storm follows"):
        get("antecedent_monotonicity")(runs, probe, _params(probe, "antecedent_monotonicity"))


def test_antecedent_window_opens_on_the_generators_storm():
    """On the probe's own forcing the criterion finds the storm the generator
    placed, without being told how long the quiet stretch before it is."""
    import importlib.util

    probe = registry.find_probe("mass/antecedent-monotonicity")
    spec = importlib.util.spec_from_file_location("antecedent_generator", probe.path / probe.generator)
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    for seed in gate_seeds(probe.id, probe.n_seeds):
        runs = {}
        for variant in probe.variants:
            case = build_case(probe, seed, variant)
            runs[variant] = RunResult(case, pd.DataFrame({"mrro": np.zeros(case.n_steps)}), {}, 0.0)
        result = get("antecedent_monotonicity")(runs, probe, _params(probe, "antecedent_monotonicity"))
        assert result.diagnostics["storm_step"] == generator.STORM_DAY - runs["dry"].case.spinup_steps


def test_streamflow_only_models_are_scored_on_every_stress_probe():
    """Runoff is enough for all four, so the model that reports nothing else
    gets a verdict on each rather than INCOMPLETE."""
    model = registry.find_model("reference_streamflow_only")
    for pid in ("mass/causality", "mass/dry-down", "mass/steady-state", "mass/extreme-rain"):
        probe = registry.find_probe(pid)
        outcome = run_probe(model, probe, gate_seeds(probe.id, 1))
        assert not outcome.missing, pid
        assert outcome.criteria, pid
        assert outcome.verdict in (PASS, FAIL)

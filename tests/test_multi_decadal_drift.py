"""A short conservative-looking record must not hide impossible long storage."""

from __future__ import annotations

from dataclasses import replace
import runpy

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.base import get, make_window
from hydroturing.harness import build_case, resolve_window_days, select_window
from hydroturing.protocol import RunResult
from hydroturing.runner import get_runner
from hydroturing.seeds import gate_seeds
from hydroturing.spec import REPO_ROOT


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/multi-decadal-drift")


def evaluate(run, probe):
    return {c.name: get(c.name)(run, probe, c.params) for c in probe.criteria}


def test_repeated_forcing_contract(probe):
    a, b, c = [build_case(probe, seed) for seed in (7, 7, 8)]
    assert a.forcing.to_csv(index=False) == b.forcing.to_csv(index=False)
    assert a.static == b.static
    assert not a.forcing.pr.equals(c.forcing.pr)
    assert a.n_steps == 20075
    assert a.spinup_steps == 1825
    assert a.forcing.tas.min() > a.static["snow_threshold_degC"]
    assert a.forcing.pr.sum() > 0
    times = pd.to_datetime(a.forcing.time)
    assert times.diff().iloc[1:].eq(pd.Timedelta(days=1)).all()
    weather = a.forcing[["pr", "tas", "pet"]].to_numpy()
    for start in range(1825, a.n_steps, 1825):
        np.testing.assert_array_equal(weather[:1825], weather[start:start + 1825])


@pytest.mark.parametrize("choice", [None, 1, 30, 365])
def test_full_horizon_cannot_be_replaced_by_event_window(probe, choice):
    model = registry.find_model("google_flood_forecast")
    days = resolve_window_days(model, probe, choice)
    assert days == 18250
    bounds = select_window(build_case(probe, 7), probe, days)
    assert bounds.offset_days == 0
    assert bounds.days == 18250


@pytest.mark.parametrize("seed", [7, 99, 12345])
def test_drift_closes_but_only_long_record_exceeds_bounds(probe, seed, tmp_path):
    model = registry.find_model("reference_slow_drift")
    case = build_case(probe, seed)
    runner = get_runner(model)
    long = runner.run(model, probe, case, tmp_path / "long")
    short_case = replace(case, forcing=case.forcing.iloc[:1825 + 3650].copy())
    short = runner.run(model, probe, short_case, tmp_path / "short")
    pd.testing.assert_frame_equal(short.table, long.table.iloc[:len(short.table)])
    assert long.table.mrro.min() >= 0
    short_scores, long_scores = evaluate(short, probe), evaluate(long, probe)
    assert all(c.passed for c in short_scores.values()), short_scores
    assert {name for name, c in long_scores.items() if not c.passed} == {"state_bounds"}
    assert long_scores["closure"].value < 1e-10
    window = make_window(long, probe)
    assert window.state0.equals(long.table.iloc[1824])


def test_reporter_changes_only_the_two_declared_terms(probe):
    exact = runpy.run_path(str(REPO_ROOT / "models/reference_bucket/ht_adapter.py"))["simulate"]
    drift = runpy.run_path(str(REPO_ROOT / "models/reference_slow_drift/ht_adapter.py"))["simulate"]
    case = build_case(probe, 19)
    forcing = case.forcing.to_dict("records")
    baseline = pd.DataFrame(exact(forcing, case.static))
    biased = pd.DataFrame(drift(forcing, case.static))
    np.testing.assert_allclose(biased.mrro, baseline.mrro - 0.012, atol=1e-12)
    np.testing.assert_allclose(biased.mrso, baseline.mrso + 0.012 * np.arange(1, len(biased) + 1))
    pd.testing.assert_frame_equal(biased.drop(columns=["mrro", "mrso"]),
                                  baseline.drop(columns=["mrro", "mrso"]))


@pytest.mark.parametrize("name", ["reference_bucket", "flex_lumped", "flex_topo", "sacsma_snow17"])
def test_physical_baseline_and_longer_spinup(probe, name, tmp_path):
    model = registry.find_model(name)
    case = build_case(probe, gate_seeds(probe.id, 1)[0])
    # Append one full forcing block: 10-year spinup followed by 50 scored years.
    extra = case.forcing.iloc[:1825].copy()
    forcing = pd.concat([case.forcing, extra], ignore_index=True)
    forcing["time"] = pd.date_range("2000-01-01", periods=len(forcing), freq="D").strftime("%Y-%m-%d")
    longer_case = replace(case, forcing=forcing, spinup_steps=3650)
    longer_probe = replace(probe, spinup_days=3650)
    runner = get_runner(model)
    normal = runner.run(model, probe, case, tmp_path / "normal")
    longer = runner.run(model, longer_probe, longer_case, tmp_path / "longer")
    for result, spec in ((normal, probe), (longer, longer_probe)):
        scores = evaluate(result, spec)
        assert all(c.passed for c in scores.values()), scores
    # Compare late repeated-cycle means, allowing slow initial equilibration.
    for state in ("mrso", "snw", "canopy", "gw", "channel"):
        if state in normal.table:
            a = normal.table[state].iloc[-1825:].mean()
            b = longer.table[state].iloc[-1825:].mean()
            assert abs(a - b) < 0.01, (name, state, a, b)


def test_longer_spinup_preserves_discrimination(probe):
    simulate = runpy.run_path(str(REPO_ROOT / "models/reference_slow_drift/ht_adapter.py"))["simulate"]
    case = build_case(probe, 42)
    forcing = pd.concat([case.forcing, case.forcing.iloc[:1825]], ignore_index=True)
    forcing["time"] = pd.date_range("2000-01-01", periods=len(forcing), freq="D").strftime("%Y-%m-%d")
    case = replace(case, forcing=forcing, spinup_steps=3650)
    spec = replace(probe, spinup_days=3650)
    run = RunResult(case, pd.DataFrame(simulate(forcing.to_dict("records"), case.static)), {}, 0)
    scores = evaluate(run, spec)
    assert scores["closure"].passed and not scores["state_bounds"].passed
    short = replace(run, case=replace(case, forcing=forcing.iloc[:7300]), table=run.table.iloc[:7300])
    assert all(c.passed for c in evaluate(short, spec).values())

"""The hourly probe must reach models intact, with complete scored phases.

Criterion arithmetic is tested separately. These tests cover mistakes in
generation and staging that would otherwise change the experiment itself.
"""

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
from hydroturing.protocol import FORCING_FILE, stage


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("energy/surface-energy-closure")


def test_generator_reproduces_bytes_and_changes_the_weather_with_seed(probe):
    generate = load_generator(probe).generate
    first, static = generate(0)
    repeated, repeated_static = generate(0)
    other, _ = generate(1)
    assert first.to_csv(index=False) == repeated.to_csv(index=False)
    assert static == repeated_static
    assert not first["rn"].equals(other["rn"])


@pytest.mark.parametrize("seed", [0, 19])
def test_case_has_warm_bare_conditions_and_complete_hourly_phases(probe, seed):
    case = build_case(probe, seed)
    assert case.n_steps == 384
    assert case.spinup_steps == 48
    assert case.timestep == "PT1H"
    time = pd.to_datetime(case.forcing["time"])
    assert time.diff().iloc[1:].eq(pd.Timedelta(hours=1)).all()
    assert time.iloc[0].hour == 6
    assert case.forcing["tas"].gt(0.0).all()
    assert case.static["canopy_capacity_mm"] == 0.0

    scored = case.after_spinup(case.forcing)
    assert len(scored) == 336
    labels = scored["_regime"]
    groups = labels.ne(labels.shift()).cumsum()
    sizes = labels.groupby(groups).size()
    assert len(sizes) == 28
    assert sizes.eq(12).all()
    assert labels.groupby(groups).first().tolist() == ["day", "night"] * 14
    hours = pd.to_datetime(scored["time"]).dt.hour
    assert np.array_equal(labels.eq("day"), (hours >= 6) & (hours < 18))
    # Weak radiation exercises the absolute allowance in the actual forcing.
    assert scored["rn"].abs().lt(40.0).any()


def test_staging_strips_phase_annotation_but_preserves_hourly_forcing(probe, tmp_path):
    case = build_case(probe, 0)
    model = registry.find_model("reference_coupled")
    request_path = stage(tmp_path, case, probe, model)
    request = json.loads(request_path.read_text())
    staged = pd.read_csv(tmp_path / FORCING_FILE)
    assert "_regime" not in staged.columns
    assert "_regime" in case.forcing.columns
    assert request["timestep"] == "PT1H"
    assert request["n_steps"] == 384
    pd.testing.assert_frame_equal(staged, case.forcing.drop(columns="_regime"))


def test_default_and_short_requested_windows_keep_all_scored_phases(probe):
    # A submitted model normally gets seven hourly days. Here the minimum
    # must retain fourteen, including when a caller explicitly asks for less.
    model = replace(registry.find_model("reference_coupled"), name="submitted_energy")
    assert resolve_window_days(model, probe) == 14
    assert resolve_window_days(model, probe, override=1) == 14
    case = build_case(probe, 0)
    bounds = select_window(case, probe, resolve_window_days(model, probe))
    cut = window_case(case, bounds)
    assert cut.spinup_steps == 48
    assert cut.window["rows"] == 336
    pd.testing.assert_frame_equal(cut.forcing, case.forcing)

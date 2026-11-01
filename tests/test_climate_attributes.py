"""The attributes the network-driven adapters read from a catchment's own forcing.

The adapters import torch, which the harness does not depend on, and these
functions do not use it, so a stand-in lets them run without it.
"""

import importlib.util
import math
import sys
from unittest import mock

import numpy as np
import pytest

from hydroturing.spec import REPO_ROOT


def load(monkeypatch, name):
    monkeypatch.setitem(sys.modules, "torch", mock.MagicMock())
    spec = importlib.util.spec_from_file_location(
        f"{name}_adapter", REPO_ROOT / "models" / name / "ht_adapter.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rain_only_after_the_last_whole_day():
    """Twenty-five hourly rows: the last is the only wet one, and no whole day holds it."""
    pr = np.zeros(25)
    pr[-1] = 24.0
    return pr, np.full(25, -5.0), np.full(25, 2.0)


def rain_in_cold_and_warm_days():
    """Two days of hourly rows, cold and wet then warm and twice as wet."""
    pr = np.r_[np.full(24, 1.0), np.full(24, 2.0)]
    tas = np.r_[np.full(24, -5.0), np.full(24, 5.0)]
    return pr, tas, np.full(48, 2.0)


def test_dhbv2_snowfall_fraction_has_a_value_when_no_whole_day_holds_rain(monkeypatch):
    adapter = load(monkeypatch, "dhbv2")
    pr, tas, pet = rain_only_after_the_last_whole_day()
    attributes = adapter.climate_attributes(pr, tas, pet, 0.0, dt_days=1 / 24)
    assert attributes["snowfall_fraction"] == 0.0
    assert all(math.isfinite(value) for value in attributes.values())


def test_dhbv2_snowfall_fraction_is_taken_over_days(monkeypatch):
    adapter = load(monkeypatch, "dhbv2")
    pr, tas, pet = rain_in_cold_and_warm_days()
    attributes = adapter.climate_attributes(pr, tas, pet, 0.0, dt_days=1 / 24)
    assert attributes["snowfall_fraction"] == pytest.approx(1 / 3)


def test_gff_frac_snow_has_a_value_when_no_whole_day_holds_rain(monkeypatch):
    adapter = load(monkeypatch, "google_flood_forecast")
    pr, tas, pet = rain_only_after_the_last_whole_day()
    attributes = adapter.climate_attributes(pr, tas, pet, np.ones(25, dtype=int), dt_days=1 / 24)
    assert attributes["frac_snow"] == 0.0
    assert all(math.isfinite(value) for value in attributes.values())


def test_gff_frac_snow_is_taken_over_days(monkeypatch):
    adapter = load(monkeypatch, "google_flood_forecast")
    pr, tas, pet = rain_in_cold_and_warm_days()
    attributes = adapter.climate_attributes(pr, tas, pet, np.ones(48, dtype=int), dt_days=1 / 24)
    assert attributes["frac_snow"] == pytest.approx(1 / 3)

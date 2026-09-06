"""Forcing generator for energy/pet-consistency: the closure probe's ten years
of continental weather, so the soil spends long stretches both near
capacity (spring) and near empty (late summer). One variant."""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 10
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 40.0,
}

def _weather(rng: np.random.Generator, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Continental daily weather with a real winter: about 830 mm/yr, three
    months below freezing, so a third of the precipitation falls as snow."""
    day = np.arange(n)
    doy = day % 365
    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(n) < p_wet
    pr = np.where(wet, rng.gamma(shape=0.7, scale=13.0, size=n), 0.0)
    seasonal = 9.0 - 12.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(n)
    innovation = rng.normal(0.0, 2.6, n)
    for t in range(1, n):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    tas = seasonal + noise
    daylength = 1.0 + 0.35 * np.cos(2 * np.pi * (doy - 172) / 365)
    pet = np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength
    return pr, tas, pet

def _frame(pr, tas, pet) -> pd.DataFrame:
    time = pd.date_range("2000-01-01", periods=len(pr), freq="D")
    return pd.DataFrame({
        "time": time.strftime("%Y-%m-%d"),
        "pr": np.round(pr, 6),
        "tas": np.round(tas, 6),
        "pet": np.round(pet, 6),
    })


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    pr, tas, pet = _weather(rng, N_STEPS)
    return _frame(pr, tas, pet), dict(STATIC)

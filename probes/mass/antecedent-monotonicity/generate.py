"""Forcing generator for mass/antecedent-monotonicity: the same storm after a
dry month and after a wet one. The `dry` variant has a rainless month
before a 60 mm storm; the `wet` variant puts 120 mm into that month, spread
over its first 20 days, and ends it with the same ten dry days so the storm
falls on a catchment that is wetter but not raining."""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 2
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS
VARIANTS = ("dry", "wet")
STORM_DAY = SPINUP_DAYS + 200          # mid-summer of the first scored year
STORM_MM = 60.0
ANTECEDENT_DAYS = 30
ANTECEDENT_MM = 120.0
QUIET_DAYS = 10                        # rainless days just before the storm in both variants

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 38.0,
}

def _weather(rng: np.random.Generator, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Temperate, rain-dominated daily weather: about 880 mm/yr, mean air
    temperature near 15 degC and almost never below freezing, so snow timing
    cannot confound the expectation."""
    day = np.arange(n)
    doy = day % 365
    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(n) < p_wet
    pr = np.where(wet, rng.gamma(shape=0.7, scale=13.0, size=n), 0.0)
    seasonal = 15.0 - 7.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(n)
    innovation = rng.normal(0.0, 2.0, n)
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


def generate(seed: int, variant: str = "dry") -> tuple[pd.DataFrame, dict]:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")
    rng = np.random.default_rng(seed)
    pr, tas, pet = _weather(rng, N_STEPS)
    pr = pr.copy()
    start = STORM_DAY - ANTECEDENT_DAYS
    pr[start:STORM_DAY] = 0.0                       # both: a dry month before the storm
    pr[STORM_DAY] = STORM_MM                        # both: the same storm
    if variant == "wet":
        wet_days = ANTECEDENT_DAYS - QUIET_DAYS
        pr[start:start + wet_days] = ANTECEDENT_MM / wet_days
    return _frame(pr, tas, pet), dict(STATIC)

"""Forcing generator for momentum/channel-routing-mass: three years of
temperate weather punctuated by a handful of large storms, each followed by a
long dry spell, so that every event leaves the channel a recession to drain
through.

The recessions are the point. The probe scores what the reach does when no
water is entering it: a store that is only releasing what it holds must fall,
step after step, until it is empty or the next storm arrives. The dry spells
are therefore made long — several times the length of the storm that precedes
them — so that a model which leaks or invents water inside its routing has to
apply that error on a great many consecutive steps where nothing is coming in
to hide it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 3
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

# Each event is (first day past spinup, storm days, mm/day). Every storm is
# followed by a dry spell to the next one, so the record alternates a hard
# fill with a long, uninterrupted drain.
STORMS = [
    (300, 6, 42.0),
    (520, 5, 35.0),
    (760, 7, 48.0),
    (1000, 5, 30.0),
]

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


def _recessions(n: int) -> np.ndarray:
    """A mask that is one on the dry spell that follows each named storm.

    The spell runs from the day after the storm ends to the day before the
    next storm begins, so it is as long as the gap between events. Those are
    the steps the criterion scores, and the generator guarantees there are
    many of them rather than leaving it to chance.
    """
    mask = np.zeros(n, dtype=bool)
    for i, (start, storm_days, _rate) in enumerate(STORMS):
        lo = SPINUP_DAYS + start + storm_days
        if i + 1 < len(STORMS):
            hi = SPINUP_DAYS + STORMS[i + 1][0]
        else:
            hi = n
        mask[lo:hi] = True
    return mask


def _frame(pr, tas, pet) -> pd.DataFrame:
    time = pd.date_range("2000-01-01", periods=len(pr), freq="D")
    return pd.DataFrame({
        "time": time.strftime("%Y-%m-%d"),
        "pr": np.round(pr, 6),
        "tas": np.round(tas, 6),
        "pet": np.round(pet, 6),
    })


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    """Deterministic: the seed fully determines the record, and nothing is
    read from disk."""
    rng = np.random.default_rng(seed)
    pr, tas, pet = _weather(rng, N_STEPS)
    pr = pr.copy()
    for start, days, rate in STORMS:
        lo = SPINUP_DAYS + start
        hi = min(lo + days, N_STEPS)
        pr[lo:hi] += rate
    return _frame(pr, tas, pet), dict(STATIC)

"""Forcing generator for momentum/routing-conservation: three years of
temperate weather punctuated by four large storms, each followed by a long,
uninterrupted dry spell.

The dry spells are the point, and they are what the bound needs. The criterion
allows a channel at most `max_lag_days` times the largest runoff of the
preceding month. Over a rainless spell the runoff falls to the baseflow and then
to nothing, that peak decays out of the window, and the ceiling comes down with
it. A reach that keeps water it should have released is therefore measured
against a shrinking allowance rather than against a storm-time peak, which is
what makes a kernel that holds back a small share of every day's runoff visible
at all: the residue is small, but the allowance it is compared against is
smaller.

So each storm is made to fill the store hard, and the gap that follows it is
several times the length of the storm, so that the drain is long enough for the
peak to leave the window and for the store to have to keep falling.
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
    next storm begins, so it is as long as the gap between events: `generate`
    zeroes the rain over it, which is what makes the drain the storm leaves
    behind rather than the weather's own chance.
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
    read from disk.

    The weather is drawn first and the storms are added on top of it; the mask
    `_recessions` returns then zeroes the rain over every post-storm spell.
    Without that last step the spells would still carry the weather's `p_wet`
    of about 0.25, no stretch of the record would be rainless for longer than a
    month, the peak runoff would never leave the bound's window, and the
    ceiling would stay at a storm-time level all year.
    """
    rng = np.random.default_rng(seed)
    pr, tas, pet = _weather(rng, N_STEPS)
    pr = pr.copy()
    for start, days, rate in STORMS:
        lo = SPINUP_DAYS + start
        hi = min(lo + days, N_STEPS)
        pr[lo:hi] += rate
    pr[_recessions(N_STEPS)] = 0.0
    return _frame(pr, tas, pet), dict(STATIC)

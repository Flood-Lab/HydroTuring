"""Forcing generator for momentum/stage-discharge-monotonic: four years of
temperate weather with a sequence of multi-day storms, so a model that routes
its runoff traces a sequence of distinct flood waves — each with a rising
limb and a longer recession — and the rating curve it reports is sampled
across the full range of flow from baseflow to peak.

The floods are built from ordinary rainfall: a multi-day wet spell that fills
the catchment, then a dry spell that lets the channel drain. That is the
shape a real hydrograph has and the shape `rating_loop` needs — a rising limb
long enough to bin, and a falling limb longer still, so the two can be paired
at a shared store value.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 4
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

# Each flood is (first day past spinup, wet days, drainage days, mm/day).
# The wet spell is short and the drainage long, which is what makes the
# recession the more gradual of the two limbs.
FLOODS = [
    (400, 5, 14, 26.0),
    (520, 6, 18, 32.0),
    (700, 5, 12, 22.0),
    (830, 7, 20, 38.0),
    (1010, 6, 16, 28.0),
]

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 38.0,
    # Reach geometry: one channel for every model to be judged in, so an
    # adapter that builds a hydraulic stage has a bank to build it in. The
    # cross-section is a rectangular bed of `width_m`; `bankfull_depth_m` is
    # where it first spills, and the reach is long enough that a flood wave
    # is still in transit while the next one is forming, which is what lets
    # a rating loop at all.
    "width_m": 18.0,
    "bankfull_depth_m": 2.5,
    "slope": 0.0015,
    "manning_n": 0.035,
    "reach_length_m": 4500.0,
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


def _storm(rng: np.random.Generator, wet_days: int, drain_days: int, rate: float) -> np.ndarray:
    """One event as a rainfall-rate series: a steady multi-day wet spell, then
    a dry tail longer than the spell itself.

    The tail is not rainfall — it is the days after the storm, left dry so
    that whatever the model stored has to drain on its own, which is what
    makes the recession the slow limb of the event. The spell is given a
    little texture rather than a flat rate, so the rising limb is not a
    single step.
    """
    spell = np.full(wet_days, rate)
    spell = spell * (1.0 + 0.08 * rng.normal(0.0, 1.0, wet_days))
    tail = np.zeros(drain_days)
    return np.clip(np.concatenate([spell, tail]), 0.0, None)


def _storms(rng: np.random.Generator, n: int) -> np.ndarray:
    """Rainfall for every named event, summed onto the record."""
    series = np.zeros(n)
    for start, wet_days, drain_days, rate in FLOODS:
        storm = _storm(rng, wet_days, drain_days, rate)
        lo = SPINUP_DAYS + start
        hi = min(lo + len(storm), n)
        series[lo:hi] += storm[: hi - lo]
    return series


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
    pr = pr + _storms(rng, N_STEPS)
    return _frame(pr, tas, pet), dict(STATIC)

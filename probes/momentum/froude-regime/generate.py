"""Forcing generator for momentum/froude-regime: four years of temperate
weather with a sequence of multi-day storms, in one reach on a mild slope.

The reach is the point. Froude is a statement about a cross-section, so the
case fixes the section — width, bed and slope are static attributes every
model is judged in — and then varies only the water going through it. Mild
slope is declared, not inferred: on a steep reach Fr > 1 is the expected
regime and the criterion would be wrong to call it a violation.

What the storms have to do is carry the reach from baseflow to bankfull and
back. A Froude number is `Q / (w * d**1.5 * sqrt(g))`, so it is only
discriminating where the flow actually changes: at low flow every plausible
cross-section is subcritical whatever the model does with its gauge, and
the separation between a gauge that reads the flow and one that does not
appears on the rising limb, where the discharge climbs and a pinned depth
has to absorb all of it in the velocity.

So each event is a short wet spell that fills the catchment followed by a
long dry tail, which is also what makes the falling limb long enough to
score. Between events the weather is ordinary rain, so the reach spends most
of the record at low flow and the scored high-flow steps are the ones the
storms made — which is where the finding is.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 4
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

# Each flood is (first day past spinup, wet days, drainage days, mm/day).
# The wet spell is short and the drainage long, so the rising limb is steep
# — a lot of water arriving in a few steps — and the falling limb is the
# long one.
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
    # Reach geometry. One channel for every model to be judged in, so the
    # width in the Froude denominator is not one the model chose. The slope
    # is mild by construction: at 0.0015 the reach is a lowland channel, and
    # supercritical flow in it is the violation the probe is looking for
    # rather than the expected regime.
    "width_m": 18.0,
    "bankfull_depth_m": 2.5,
    "slope": 0.0015,
    "manning_n": 0.035,
    "reach_length_m": 4500.0,
    # The bed the stage is measured from. Zero: every adapter in the suite
    # reports its stage as a depth above the bed, from a normal-depth
    # relation, so there is nothing to subtract. Declared rather than assumed
    # so a model that reports an absolute level can be scored by subtracting
    # it, instead of being read as a reach two metres deep.
    "bed_elevation_m": 0.0,
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

    The spell is given a little texture rather than a flat rate so the rising
    limb is several steps and not one, and the tail is left dry so that
    whatever the reach stored has to drain on its own.
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

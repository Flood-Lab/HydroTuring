"""Forcing generator for mass/extreme-rain: the same weather, the biggest
storm scaled up a ladder.

Deterministic given a seed and a variant. All variants carry the same
weather, drawn once; the largest scored day of rain is found once, and each
variant multiplies that one day by its factor. Ten times the largest day of
a temperate record is several hundred millimetres, outside anything the
model will have been fitted to, which is the point.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 3
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

FACTORS = {"control": 1.0, "x2": 2.0, "x5": 5.0, "x10": 10.0}
VARIANTS = tuple(FACTORS)

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
    """Temperate, rain-dominated daily weather: about 880 mm/yr of rain, mean
    air temperature near 15 degC and almost never below freezing, so snow
    timing cannot confound any of the expectations built on it."""
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


def generate(seed: int, variant: str = "control") -> tuple[pd.DataFrame, dict]:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")
    rng = np.random.default_rng(seed)
    pr, tas, pet = _weather(rng, N_STEPS)
    storm = SPINUP_DAYS + int(np.argmax(pr[SPINUP_DAYS:]))
    pr = pr.copy()
    pr[storm] *= FACTORS[variant]
    return _frame(pr, tas, pet), dict(STATIC)


if __name__ == "__main__":
    control, _ = generate(20260903, "control")
    storm = SPINUP_DAYS + int(np.argmax(control["pr"].to_numpy()[SPINUP_DAYS:]))
    print(f"{N_STEPS} steps; largest storm {control['pr'].iloc[storm]:.0f} mm on {control['time'].iloc[storm]}")
    for name in VARIANTS:
        frame, _ = generate(20260903, name)
        print(f"  {name:8s} storm day {frame['pr'].iloc[storm]:.0f} mm, record total {frame['pr'].sum():.0f} mm")

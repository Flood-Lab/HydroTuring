"""Forcing generator for mass/steady-state: a year of weather, then the same
day three years running.

Deterministic given a seed. The seed only shapes the spinup year, so the
three seeds start the constant stretch from three different states and
have to settle to the same place. The constant stretch carries 2.5 mm/day
of rain, 12 degC and 2.0 mm/day of potential evaporation, every day.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 3
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

RAIN_MM_PER_DAY = 2.5
TAS_DEGC = 12.0
PET_MM_PER_DAY = 2.0

# A faster recession than the closure probe's catchment, so that the stores
# settle well inside the three years.
STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.02,
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


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    pr, tas, pet = _weather(rng, N_STEPS)
    pr, tas, pet = pr.copy(), tas.copy(), pet.copy()
    pr[SPINUP_DAYS:] = RAIN_MM_PER_DAY
    tas[SPINUP_DAYS:] = TAS_DEGC
    pet[SPINUP_DAYS:] = PET_MM_PER_DAY
    return _frame(pr, tas, pet), dict(STATIC)


if __name__ == "__main__":
    frame, _ = generate(20260903)
    scored = frame.iloc[SPINUP_DAYS:]
    print(f"{N_STEPS} steps; scored rain {scored['pr'].min():.2f} to {scored['pr'].max():.2f} mm/day, "
          f"tas {scored['tas'].min():.1f} to {scored['tas'].max():.1f}, pet {scored['pet'].min():.2f} to {scored['pet'].max():.2f}")

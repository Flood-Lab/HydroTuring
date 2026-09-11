"""Identical weather and catchment under two seasonally aligned calendars."""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 10
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS
ORIGINS = {"control": "2000-01-01", "shifted": "1972-01-01"}

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 38.0,
}


def generate(seed: int, variant: str = "control") -> tuple[pd.DataFrame, dict]:
    if variant not in ORIGINS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {tuple(ORIGINS)}")

    rng = np.random.default_rng(seed)
    control_time = pd.date_range(ORIGINS["control"], periods=N_STEPS, freq="D")
    # Use actual calendar phase, including leap years. Both origins have the
    # same month, day, day of year, weekday and leap-day positions in this
    # fixed window, so seasonal features remain valid under the transform.
    phase = (control_time.dayofyear.to_numpy() - 1) / np.where(
        control_time.is_leap_year, 366.0, 365.0
    )
    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (phase - 30 / 365))
    wet = rng.random(N_STEPS) < p_wet
    pr = np.where(wet, rng.gamma(shape=0.7, scale=13.0, size=N_STEPS), 0.0)

    seasonal = 15.0 - 7.0 * np.cos(2 * np.pi * (phase - 15 / 365))
    noise = np.zeros(N_STEPS)
    innovation = rng.normal(0.0, 2.0, N_STEPS)
    for t in range(1, N_STEPS):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    # A rain-only regime keeps the rainfall-response guard meaningful. Snow
    # storage is still compared; a legitimately empty store must be allowed.
    tas = np.maximum(2.0, seasonal + noise)
    daylength = 1.0 + 0.35 * np.cos(2 * np.pi * (phase - 172 / 365))
    pet = np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength

    time = pd.date_range(ORIGINS[variant], periods=N_STEPS, freq="D")
    frame = pd.DataFrame({
        "time": time.strftime("%Y-%m-%d"),
        "pr": np.round(pr, 6),
        "tas": np.round(tas, 6),
        "pet": np.round(pet, 6),
    })
    return frame, dict(STATIC)

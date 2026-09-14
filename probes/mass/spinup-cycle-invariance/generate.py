"""One deterministic annual cycle on a calendar-aligned repeated record."""

from __future__ import annotations

import numpy as np
import pandas as pd

CYCLE_DAYS = 365
SHORT_CYCLES = 5
EXTRA_CYCLES = 4
TOTAL_CYCLES = SHORT_CYCLES + 1 + EXTRA_CYCLES
SPINUP_DAYS = 365
VARIANTS = ("short", "long")

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 38.0,
}


def _cycle(seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    day = np.arange(CYCLE_DAYS)
    p_wet = 0.24 + 0.13 * np.cos(2.0 * np.pi * (day - 30) / CYCLE_DAYS)
    wet = rng.random(CYCLE_DAYS) < p_wet
    pr = np.where(wet, rng.gamma(shape=0.8, scale=12.0, size=CYCLE_DAYS), 0.0)
    seasonal = 15.0 - 6.0 * np.cos(2.0 * np.pi * (day - 15) / CYCLE_DAYS)
    noise = np.zeros(CYCLE_DAYS)
    innovation = rng.normal(0.0, 1.5, CYCLE_DAYS)
    for index in range(1, CYCLE_DAYS):
        noise[index] = 0.65 * noise[index - 1] + innovation[index]
    tas = np.maximum(2.0, seasonal + noise)
    daylength = 1.0 + 0.35 * np.cos(2.0 * np.pi * (day - 172) / CYCLE_DAYS)
    pet = np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength
    return pr, tas, pet


def _cycle_day(times: pd.DatetimeIndex) -> np.ndarray:
    """Map calendar dates to the 365 values in the synthetic annual cycle.

    A real Gregorian axis contains leap days, while the synthetic weather
    cycle intentionally has 365 values.  Reusing February 28 on February 29
    explicitly, then shifting later dates back by one, keeps every January 1
    and every month/day aligned across repetitions.  The scored comparison is
    four cycles apart so both years have the same leap-day phase.
    """
    day = times.dayofyear.to_numpy() - 1
    feb_29 = times.is_leap_year & (times.month == 2) & (times.day == 29)
    leap_after_feb = times.is_leap_year & (times.month > 2)
    day = np.where(feb_29, 58, day)
    day = day - leap_after_feb.astype(int)
    return np.clip(day, 0, CYCLE_DAYS - 1)


def generate(seed: int, variant: str = "short") -> tuple[pd.DataFrame, dict]:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")
    pr, tas, pet = _cycle(seed)
    period_start = pd.Timestamp("2002-01-01")
    period_end = period_start + pd.DateOffset(years=TOTAL_CYCLES)
    period_times = pd.date_range(
        period_start, period_end - pd.Timedelta(days=1), freq="D"
    )
    spinup_times = pd.date_range("2001-01-01", periods=SPINUP_DAYS, freq="D")
    times = spinup_times.append(period_times)
    values = _cycle_day(times)
    forcing = pd.DataFrame({
        "time": times.strftime("%Y-%m-%d"),
        "pr": np.round(pr[values], 6),
        "tas": np.round(tas[values], 6),
        "pet": np.round(pet[values], 6),
    })
    evaluation_year = period_start.year + (
        SHORT_CYCLES if variant == "short" else SHORT_CYCLES + EXTRA_CYCLES
    )
    forcing["_phase"] = np.where(
        forcing["time"].str[:4].astype(int).eq(evaluation_year),
        "evaluation",
        np.where(
            forcing["time"].str[:4].astype(int).lt(period_start.year),
            "spinup",
            "spinup_or_tail",
        ),
    )
    return forcing, dict(STATIC)

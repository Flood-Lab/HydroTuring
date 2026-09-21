"""One deterministic annual cycle on records with a common evaluation year."""

from __future__ import annotations

import numpy as np
import pandas as pd

CYCLE_DAYS = 365
SHORT_CYCLES = 5
PLUS3_CYCLES = SHORT_CYCLES + 3
LONG_CYCLES = SHORT_CYCLES + 4
TOTAL_CYCLES = LONG_CYCLES + 1
PERIOD_DAYS = 3652
SPINUP_DAYS = 365
VARIANTS = ("short", "plus3", "long")
EVALUATION_CYCLES = {
    "short": SHORT_CYCLES,
    "plus3": PLUS3_CYCLES,
    "long": LONG_CYCLES,
}
EVALUATION_YEAR = 2007
# Each variant reaches the same calendar evaluation year after a different
# number of prior cycles. The differing starts avoid comparing leap and
# non-leap evaluation years while retaining the coprime history offsets.
START_YEARS = {
    "short": 2001,
    "plus3": 1998,
    "long": 1997,
}

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
    and every month/day aligned across repetitions.  All variants score the
    non-leap calendar year 2007, so the paired comparison is not confounded by
    different leap-day phases.
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
    start_year = START_YEARS[variant]
    period_start = pd.Timestamp(f"{start_year + 1}-01-01")
    period_times = pd.date_range(
        period_start, periods=PERIOD_DAYS, freq="D"
    )
    spinup_times = pd.date_range(
        f"{start_year}-01-01", periods=SPINUP_DAYS, freq="D"
    )
    times = spinup_times.append(period_times)
    values = _cycle_day(times)
    forcing = pd.DataFrame({
        "time": times.strftime("%Y-%m-%d"),
        "pr": np.round(pr[values], 6),
        "tas": np.round(tas[values], 6),
        "pet": np.round(pet[values], 6),
    })
    evaluation_year = period_start.year + EVALUATION_CYCLES[variant]
    if evaluation_year != EVALUATION_YEAR:  # pragma: no cover - constants guard
        raise AssertionError(
            f"{variant} evaluates in {evaluation_year}, expected {EVALUATION_YEAR}"
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

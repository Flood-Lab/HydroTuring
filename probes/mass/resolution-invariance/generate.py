"""Seeded daily weather for the resolution-invariance probe, served at three steps.

One month per seed, written at the day, the hour and the minute. The month
is drawn minute by minute and averaged over each day, which gives the same
daily record version 2 of this probe served. The `hourly` and `minute`
variants are that daily record held flat: every hour and every minute of a
day carries the day's rates. Each step therefore delivers the same water at
the same intensity, and only the length of the row changes.

That is the point of holding it flat. Version 2 served the minutes as
drawn, with storms in five-minute bursts, and averaged them to hours and
days. Averaging lowers the intensity a model sees, and a model whose
runoff depends on intensity, infiltration excess above all, then partitions
the same water differently at each step. That is physics, not an error, and
it moved the integrated runoff of a correct Horton bucket by up to 18
percent of the rain between minutes and hours. With the intensity held,
the physics is the same at every step, and what is left to move the volumes
is the step itself: the discretisation of a model's own equations, which is
small, and arithmetic that assumes one step, which is not. A model built at
a day and handed hours still treats each hour as a day and receives
twenty-four times the rain.

Forcing is in the contract's units at every step: precipitation and
potential evapotranspiration are rates in mm per day, temperature in degC.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_DAYS = 30
SPINUP_DAYS = 10
MINUTES_PER_DAY = 1440
N_MINUTES = (PERIOD_DAYS + SPINUP_DAYS) * MINUTES_PER_DAY
N_HOURS = N_MINUTES // 60
N_DAYS = PERIOD_DAYS + SPINUP_DAYS

VARIANTS = ("minute", "hourly", "daily")
MINUTES_PER = {"minute": 1, "hourly": 60, "daily": MINUTES_PER_DAY}
FREQ = {"minute": "min", "hourly": "h", "daily": "D"}
START = "2000-04-01"

STORMS_PER_DAY = 0.6
BURST_MINUTES = 5

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 120.0,  # shallow, so storms produce saturation excess
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 40.0,
}


def _draw_minutes(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precipitation, temperature and potential ET, one value per minute."""
    minute = np.arange(N_MINUTES)
    hour_of_day = (minute % MINUTES_PER_DAY) / 60.0
    day = minute // MINUTES_PER_DAY
    n_days = PERIOD_DAYS + SPINUP_DAYS

    # Storms: Poisson arrivals, lognormal durations, gamma depths, and the
    # depth spread over five-minute bursts with lognormal weights so the
    # intensity varies inside the hour.
    depth_mm = np.zeros(N_MINUTES)
    n_storms = int(rng.poisson(STORMS_PER_DAY * n_days))
    starts = np.sort(rng.integers(0, N_MINUTES, n_storms))
    durations = np.clip(rng.lognormal(np.log(90.0), 0.6, n_storms), 10, 600).astype(int)
    depths = rng.gamma(shape=1.2, scale=9.0, size=n_storms)
    for start, duration, depth in zip(starts, durations, depths):
        stop = min(start + duration, N_MINUTES)
        n_bursts = max(1, (stop - start) // BURST_MINUTES)
        weights = rng.lognormal(0.0, 0.8, n_bursts)
        weights /= weights.sum()
        for burst, weight in enumerate(weights):
            b0 = start + burst * BURST_MINUTES
            b1 = min(b0 + BURST_MINUTES, stop)
            if b1 > b0:
                depth_mm[b0:b1] += depth * weight / (b1 - b0)
    pr = depth_mm * MINUTES_PER_DAY  # mm per minute -> mm per day

    # Temperature: a diurnal cycle on a slowly wandering daily anomaly, cooled
    # a little while it rains. Mostly above freezing, with cold nights.
    anomaly = np.zeros(n_days)
    innovation = rng.normal(0.0, 2.5, n_days)
    for d in range(1, n_days):
        anomaly[d] = 0.7 * anomaly[d - 1] + innovation[d]
    tas = (
        9.0
        + 5.5 * np.sin(2.0 * np.pi * (hour_of_day - 9.0) / 24.0)
        + anomaly[day]
        - 2.0 * (depth_mm > 0)
    )

    # Potential ET: a daily amount tied to the day's mean temperature, spread
    # over daylight with a half-sine so the rate is zero at night and the
    # daily mean equals the daily amount.
    daily_mean_t = np.bincount(day, tas) / MINUTES_PER_DAY
    daily_pet = np.maximum(0.0, 0.13 * (daily_mean_t + 5.0))
    daylight = np.maximum(0.0, np.sin(np.pi * (hour_of_day - 6.0) / 12.0))
    pet = daily_pet[day] * daylight * np.pi
    return pr, tas, pet


def generate(seed: int, variant: str = "minute") -> tuple[pd.DataFrame, dict]:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")

    rng = np.random.default_rng(seed)
    pr, tas, pet = _draw_minutes(rng)

    # Average to days, then hold each day flat over the rows of the step
    # served, so the intensity a model sees does not depend on the step.
    rows_per_day = MINUTES_PER_DAY // MINUTES_PER[variant]
    pr, tas, pet = (
        np.repeat(x.reshape(-1, MINUTES_PER_DAY).mean(axis=1), rows_per_day)
        for x in (pr, tas, pet)
    )
    time = pd.date_range(START, periods=N_DAYS * rows_per_day, freq=FREQ[variant])

    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d" if variant == "daily" else "%Y-%m-%d %H:%M"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    fine, _ = generate(20260903, "minute")
    coarse, _ = generate(20260903, "hourly")
    daily, _ = generate(20260903, "daily")
    total_fine = fine["pr"].sum() / MINUTES_PER_DAY
    total_coarse = coarse["pr"].sum() / 24.0
    total_daily = daily["pr"].sum()
    wet_days = int((daily["pr"] > 0).sum())
    print(f"{len(fine)} minutes, {len(coarse)} hours, {len(daily)} days")
    print(f"precipitation {total_fine:.2f} mm at the minute step, {total_coarse:.2f} mm hourly, "
          f"{total_daily:.2f} mm daily")
    print(f"{wet_days} wet days, peak {daily['pr'].max():.0f} mm/day at every step "
          f"(minute {fine['pr'].max():.0f}, hourly {coarse['pr'].max():.0f})")
    print(f"temperature {daily['tas'].mean():.1f} degC mean, {int((daily['tas'] < 0).sum())} sub-zero days")

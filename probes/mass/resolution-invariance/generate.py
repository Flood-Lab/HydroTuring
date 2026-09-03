"""Seeded minute-by-minute weather for the resolution-invariance probe.

One sequence per seed, written at two steps. The `minute` variant is the
record as drawn. The `hourly` variant is the same record averaged over each
hour, so that every hour carries exactly the water its sixty minutes did.
Draw once, then aggregate: a second draw at the coarse step would be a
different month, and the comparison would no longer isolate the step.

Storms have structure inside the hour on purpose. A storm that is uniform
over its hour aggregates to itself and tests nothing; one delivered in
five-minute bursts is what a fixed-step model gets wrong.

Forcing is in the contract's units at every step: precipitation and
potential evapotranspiration are rates in mm per day, temperature in degC.
A one-minute burst of one millimetre is therefore a rate of 1440 mm/day.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_DAYS = 30
SPINUP_DAYS = 10
MINUTES_PER_DAY = 1440
N_MINUTES = (PERIOD_DAYS + SPINUP_DAYS) * MINUTES_PER_DAY
N_HOURS = N_MINUTES // 60

VARIANTS = ("minute", "hourly")
START = "2000-04-01"

STORMS_PER_DAY = 0.6
BURST_MINUTES = 5

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
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

    if variant == "hourly":
        pr = pr.reshape(N_HOURS, 60).mean(axis=1)
        tas = tas.reshape(N_HOURS, 60).mean(axis=1)
        pet = pet.reshape(N_HOURS, 60).mean(axis=1)
        time = pd.date_range(START, periods=N_HOURS, freq="h")
    else:
        time = pd.date_range(START, periods=N_MINUTES, freq="min")

    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d %H:%M"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    fine, _ = generate(20260903, "minute")
    coarse, _ = generate(20260903, "hourly")
    total_fine = fine["pr"].sum() / MINUTES_PER_DAY
    total_coarse = coarse["pr"].sum() / 24.0
    wet_minutes = int((fine["pr"] > 0).sum())
    print(f"{len(fine)} minutes, {len(coarse)} hours, {PERIOD_DAYS + SPINUP_DAYS} days")
    print(f"precipitation {total_fine:.2f} mm at the minute step, {total_coarse:.2f} mm hourly")
    print(f"{wet_minutes} wet minutes, peak {fine['pr'].max():.0f} mm/day; "
          f"hourly peak {coarse['pr'].max():.0f} mm/day")
    print(f"temperature {fine['tas'].mean():.1f} degC mean, {int((fine['tas'] < 0).sum())} sub-zero minutes")

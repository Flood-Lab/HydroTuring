"""Seeded hourly forcing for a warm, snow-free, bare surface energy budget.

Timestamps are interval starts in local solar time. The imposed net radiation
is an interval mean, not a radiative-temperature prediction. The harness
retains the day/night annotation for scoring and strips it before staging.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_DAYS = 14
SPINUP_DAYS = 2
STEPS_PER_DAY = 24
N_STEPS = (PERIOD_DAYS + SPINUP_DAYS) * STEPS_PER_DAY

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 180.0,
    "canopy_capacity_mm": 0.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 0.0,
}


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    time = pd.date_range("2000-03-20T06:00:00", periods=N_STEPS, freq="h")
    hour = time.hour.to_numpy()
    daytime = (hour >= 6) & (hour < 18)
    solar_day = np.arange(N_STEPS) // STEPS_PER_DAY
    n_days = PERIOD_DAYS + SPINUP_DAYS

    # Exact hourly mean of the positive sine between sunrise at 06:00 and
    # sunset at 18:00: integral sin(omega * (h - 6)) dh over [h, h + 1).
    # Cloud attenuation and longwave cooling are constant within each hour,
    # so their combination below is also an interval-mean forcing.
    omega = np.pi / 12.0
    daylight = np.where(
        daytime,
        (np.cos(omega * (hour - 6)) - np.cos(omega * (hour - 5))) / omega,
        0.0,
    )
    peak_shortwave = rng.uniform(420.0, 540.0, n_days)[solar_day]
    cloud = rng.uniform(0.45, 1.0, N_STEPS)
    longwave_cooling = rng.uniform(25.0, 40.0, n_days)[solar_day]
    rn = peak_shortwave * cloud * daylight - longwave_cooling

    # Warm air excludes snow and freezing. A daily weather offset and a
    # smooth diurnal cycle are prescribed interval means; they are not a
    # surface temperature or an energy storage state.
    weather = rng.uniform(-2.0, 2.0, n_days)[solar_day]
    midpoint = hour + 0.5
    tas = 22.0 + weather + 4.0 * np.cos(2 * np.pi * (midpoint - 14.0) / 24.0)

    # Water-side inputs for the accounting reference, in mm/day even at an
    # hourly step. Small nonzero nighttime demand is allowed: this probe
    # imposes no sign restriction on latent or sensible heat.
    pet = 0.12 + 3.2 * daylight * cloud
    wet = rng.random(N_STEPS) < 0.025
    pr = np.where(wet, rng.uniform(12.0, 48.0, N_STEPS), 0.0)

    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
            "rn": np.round(rn, 6),
            # Phase follows the solar clock, never the sign of net radiation.
            "_regime": np.where(daytime, "day", "night"),
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    frame, _ = generate(20260910)
    print(f"{len(frame)} hourly rows: {SPINUP_DAYS} spinup + {PERIOD_DAYS} scored days")
    print(f"Net radiation: {frame['rn'].min():.1f} to {frame['rn'].max():.1f} W/m2")

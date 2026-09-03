"""Seeded synthetic forcing for the catchment water balance probe.

Deterministic given a seed, and never committed as data. Generating the case
at run time is what stops a model from memorising it, and shipping the
generator instead of a NetCDF file is what lets a reviewer see exactly what
the model will be given.

Produces a lumped daily record: precipitation, air temperature and potential
evapotranspiration, plus the catchment attributes every model receives.
Targets roughly 830 mm of precipitation and 700 mm of potential ET per year,
with a real seasonal snowpack, so the water-limited and energy-limited
regimes both appear in a ten year record.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 10
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 40.0,
}


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    day = np.arange(N_STEPS)
    doy = day % 365

    # Precipitation: seasonal occurrence, gamma depths on wet days.
    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(N_STEPS) < p_wet
    depth = rng.gamma(shape=0.7, scale=13.0, size=N_STEPS)
    pr = np.where(wet, depth, 0.0)

    # Temperature: annual cycle plus a persistent (AR1) weather anomaly.
    seasonal = 9.0 - 12.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(N_STEPS)
    innovation = rng.normal(0.0, 2.6, N_STEPS)
    for t in range(1, N_STEPS):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    tas = seasonal + noise

    # Potential evapotranspiration: temperature driven with a daylength factor.
    daylength = 1.0 + 0.35 * np.cos(2 * np.pi * (doy - 172) / 365)
    pet = np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    frame, static = generate(20260903)
    annual_p = frame["pr"].sum() / (N_STEPS / 365)
    annual_pet = frame["pet"].sum() / (N_STEPS / 365)
    snow_days = int((frame["tas"] < 0).sum())
    print(f"{len(frame)} steps, {annual_p:.0f} mm/yr precipitation, "
          f"{annual_pet:.0f} mm/yr potential ET, {snow_days} sub-zero days")

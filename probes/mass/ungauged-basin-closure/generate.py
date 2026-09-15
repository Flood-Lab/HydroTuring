"""One independent catchment per seed; capacities stay fixed within each run."""
from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 10
SPINUP_DAYS = 730
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS
ATTRIBUTES = (
    ("reference", (320., 320.), (2., 2.)),
    ("reference_medium", (230., 280.), (2., 2.)),
    ("soil_small", (80., 140.), (2., 2.)),
    ("soil_large", (500., 700.), (2., 2.)),
    ("canopy_small", (320., 320.), (.5, .8)),
    ("canopy_large", (320., 320.), (2.2, 2.5)),
    ("joint_small", (80., 140.), (.5, .8)),
    ("joint_large", (500., 700.), (2.2, 2.5)),
)
MAX_DAILY_RAIN_MM = 35.0

# (name, wet-day intercept, seasonal wet amplitude, mean T, T amplitude,
#  mean PET, PET amplitude, latitude). These are controlled climate classes,
# not fitted frequencies or claims of a representative global sample.
CLIMATES = (
    ("temperate_winter_wet", .28, .10, 15., 9., 2.2, 1.3, 38.),
    ("warm_humid", .40, .03, 24., 3., 3.0, .5, 20.),
    ("summer_rainfall", .30, -.20, 23., 5., 3.0, 1.0, 25.),
    ("seasonally_water_limited", .18, .12, 19., 8., 2.6, 1.3, 35.),
    ("cool_maritime_rain", .32, .06, 10., 6., 1.6, .8, 48.),
)


def climate_for_seed(seed: int):
    return CLIMATES[int(seed) % len(CLIMATES)]


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    attr_seed, weather_seed = np.random.SeedSequence(seed).spawn(2)
    attrs = np.random.default_rng(attr_seed)
    _, soil, canopy = ATTRIBUTES[int(seed) % len(ATTRIBUTES)]
    static = {"area_km2": 250.0,
              "soil_capacity_mm": float(attrs.uniform(*soil)),
              "canopy_capacity_mm": float(attrs.uniform(*canopy)),
              "degree_day_factor_mm_per_C_day": 3.2,
              "baseflow_coefficient": 0.006, "snow_threshold_degC": 0.0}

    rng = np.random.default_rng(weather_seed)
    _, wet_base, wet_amplitude, mean_t, amplitude_t, mean_pet, amplitude_pet, latitude = climate_for_seed(seed)
    static["latitude_deg"] = latitude
    day = np.arange(N_STEPS)
    phase = 2 * np.pi * ((day % 365) - 20) / 365
    seasonal = np.cos(phase)
    # Wet-day persistence and climate-specific rainfall seasonality, with strictly bounded amounts.
    # A beta amount is drawn directly: no post-hoc clipping or seed rejection.
    uniforms = rng.random(N_STEPS)
    wet = np.zeros(N_STEPS, dtype=bool)
    for i in range(N_STEPS):
        probability = wet_base + wet_amplitude * seasonal[i] + 0.22 * (wet[i - 1] if i else 0)
        wet[i] = uniforms[i] < probability
    pr = np.where(wet, MAX_DAILY_RAIN_MM * rng.beta(0.9, 3.2, N_STEPS), 0.0)
    # Positive temperatures deliberately isolate soil/canopy attributes from
    # differences in snow-process support between physical models.
    tas = mean_t - amplitude_t * seasonal + rng.uniform(-3.0, 3.0, N_STEPS)
    pet = (mean_pet - amplitude_pet * seasonal) * rng.uniform(0.85, 1.15, N_STEPS)
    forcing = pd.DataFrame({
        "time": pd.date_range("2000-01-01", periods=N_STEPS, freq="D").strftime("%Y-%m-%d"),
        "pr": np.round(pr, 6), "tas": np.round(tas, 6), "pet": np.round(pet, 6),
    })
    return forcing, static

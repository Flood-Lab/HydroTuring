"""Daily synthetic catchments for momentum/routing-lag-consistency.

All four variants receive byte-identical weather and the same isolated,
one-day design storm. Only catchment area and the two channel lengths derived
from it change. The hidden ``_event_pr`` column marks the design-storm
component for the criterion and is removed before forcing reaches a model.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


SPINUP_DAYS = 365
PERIOD_DAYS = 90
N_STEPS = SPINUP_DAYS + PERIOD_DAYS
SQ_KM_PER_SQ_MI = 2.589988110336
KM_PER_MI = 1.609344


def _hack_geometry(area_km2: float) -> dict[str, float]:
    """Derive the public channel geometry from area with Hack's relation."""
    area_mi2 = area_km2 / SQ_KM_PER_SQ_MI
    length_mi = 1.4 * area_mi2**0.6
    length_km = length_mi * KM_PER_MI
    return {
        "area_km2": area_km2,
        "main_channel_length_km": length_km,
        "centroid_channel_length_km": 0.5 * length_km,
    }


VARIANTS = {
    "small": _hack_geometry(30.0),
    "medium": _hack_geometry(300.0),
    "large": _hack_geometry(3000.0),
    "xlarge": _hack_geometry(10000.0),
}

# Standard synthetic-catchment attributes are held fixed across the ladder.
# No expected lag is handed to the model in static.json.
COMMON_STATIC = {
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 38.0,
    "elevation_m": 500.0,
}


def _spinup_weather(
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Temperate rain-only weather followed by an isolated scored event."""
    day = np.arange(N_STEPS)
    doy = day % 365

    wet_probability = 0.25 + 0.10 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(N_STEPS) < wet_probability
    pr = np.where(wet, rng.gamma(shape=0.8, scale=11.0, size=N_STEPS), 0.0)

    seasonal_temperature = 15.0 - 6.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    temperature_noise = np.zeros(N_STEPS)
    innovations = rng.normal(0.0, 1.5, N_STEPS)
    for i in range(1, N_STEPS):
        temperature_noise[i] = 0.70 * temperature_noise[i - 1] + innovations[i]
    tas = np.maximum(3.0, seasonal_temperature + temperature_noise)

    daylength = 1.0 + 0.30 * np.cos(2 * np.pi * (doy - 172) / 365)
    pet = np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength

    # The scored record contains no other rain. Its first 18--22 days let
    # pre-event drainage settle, while the remaining tail is long enough for
    # the largest synthetic catchment's response.
    pr[SPINUP_DAYS:] = 0.0
    return pr, tas, pet


def _frame(seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    pr, tas, pet = _spinup_weather(rng)

    centre = SPINUP_DAYS + 20 + int(rng.integers(-2, 3))
    total_mm = float(rng.uniform(48.0, 52.0))
    event_pr = np.zeros(N_STEPS)
    event_pr[centre] = total_mm
    pr = pr + event_pr

    time = pd.date_range("2001-01-01", periods=N_STEPS, freq="D")
    return pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
            "_event_pr": np.round(event_pr, 6),
        }
    )


def generate(seed: int, variant: str = "small") -> tuple[pd.DataFrame, dict]:
    """Return one member of the four-catchment geometry ladder."""
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {tuple(VARIANTS)}")

    static = dict(COMMON_STATIC)
    static.update(VARIANTS[variant])
    return _frame(seed), static

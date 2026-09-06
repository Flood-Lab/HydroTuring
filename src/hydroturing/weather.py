"""Synthetic weather for a real catchment.

A probe's catchment is real: its attributes are a HydroATLAS row for a gauged
basin, loaded from `catchments/<id>.json`. Only the weather is generated, and
it is generated to match that catchment's own monthly climatology, so the
attributes a model reads and the forcing it receives describe the same place.
The sequence is still unseen: which days are wet, how much falls, how warm
the anomaly runs are all drawn from the seed. What the seed cannot change is
the climate.

Precipitation: each month has a wet-day frequency and a gamma depth chosen so
that the expected monthly total is the HydroATLAS total for that month; wet
days cluster with a first-order chain so storms last more than a day.
Temperature: the HydroATLAS monthly means, interpolated to the day, plus a
persistent AR(1) anomaly. Potential evaporation: the HydroATLAS monthly
totals, interpolated to the day, scaled with the temperature anomaly so a
warm spell raises demand as it does in any temperature-based PET.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

CATCHMENTS_DIR = Path(__file__).resolve().parents[2] / "catchments"
MONTH_DAYS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=float)
MONTH_MID = np.cumsum(MONTH_DAYS) - MONTH_DAYS / 2.0  # day of year at mid-month

# Wet-day frequency and storm persistence are not in HydroATLAS. These are
# the ranges the Caravan climate indices span for the eight basins shipped
# here; a catchment can override them from its own `caravan` block when it
# carries high_prec_freq / low_prec_freq.
DEFAULT_WET_FRACTION = 0.35
WET_AFTER_WET = 0.65        # persistence of wet days
GAMMA_SHAPE = 0.7           # depth distribution on wet days
TEMPERATURE_AR1 = 0.72
TEMPERATURE_SIGMA = 2.6     # degC, innovation of the daily anomaly
PET_PER_DEGREE = 0.06       # relative change in demand per degC of anomaly


def load_catchment(name: str) -> dict:
    path = CATCHMENTS_DIR / f"{name}.json"
    if not path.exists():
        raise FileNotFoundError(f"no catchment '{name}' under {CATCHMENTS_DIR}")
    return json.loads(path.read_text())


def monthly(catchment: dict, prefix: str, scale: float = 1.0) -> np.ndarray:
    """Twelve HydroATLAS monthly values, e.g. pre_mm_s01..s12."""
    atlas = catchment["hydroatlas"]
    return np.array([float(atlas[f"{prefix}{m:02d}"]) for m in range(1, 13)]) * scale


def _daily_from_monthly(values: np.ndarray, doy: np.ndarray) -> np.ndarray:
    """Interpolate mid-month values around the year, periodically."""
    x = np.concatenate([MONTH_MID - 365.0, MONTH_MID, MONTH_MID + 365.0])
    y = np.tile(values, 3)
    return np.interp(doy + 0.5, x, y)


def generate_weather(
    catchment: dict, n_days: int, seed: int, wet_fraction: float | None = None
) -> pd.DataFrame:
    """Daily pr, tas, pet (mm/day, degC, mm/day) for a real catchment."""
    rng = np.random.default_rng(seed)
    day = np.arange(n_days)
    doy = day % 365
    month = np.searchsorted(np.cumsum(MONTH_DAYS), doy, side="right")  # 0..11

    pre_month = monthly(catchment, "pre_mm_s")            # mm per month
    tmp_month = monthly(catchment, "tmp_dc_s", 0.1)        # degC
    pet_month = monthly(catchment, "pet_mm_s")             # mm per month

    # Wet-day chain, then depths whose expectation hits the monthly total.
    wet_frac = DEFAULT_WET_FRACTION if wet_fraction is None else float(wet_fraction)
    wet_after_dry = wet_frac * (1.0 - WET_AFTER_WET) / max(1.0 - wet_frac, 1e-9)
    wet = np.zeros(n_days, dtype=bool)
    wet[0] = rng.random() < wet_frac
    u = rng.random(n_days)
    for t in range(1, n_days):
        wet[t] = u[t] < (WET_AFTER_WET if wet[t - 1] else wet_after_dry)
    mean_depth = pre_month / (MONTH_DAYS * wet_frac)       # mm per wet day, by month
    depth = rng.gamma(GAMMA_SHAPE, mean_depth[month] / GAMMA_SHAPE)
    pr = np.where(wet, depth, 0.0)

    # Temperature: climatology plus a persistent anomaly.
    innovation = rng.normal(0.0, TEMPERATURE_SIGMA, n_days)
    anomaly = np.zeros(n_days)
    for t in range(1, n_days):
        anomaly[t] = TEMPERATURE_AR1 * anomaly[t - 1] + innovation[t]
    tas = _daily_from_monthly(tmp_month, doy) + anomaly

    # Demand: climatology, warmer days ask for more.
    pet = _daily_from_monthly(pet_month / MONTH_DAYS, doy) * np.clip(
        1.0 + PET_PER_DEGREE * anomaly, 0.0, None
    )

    time = pd.date_range("2000-01-01", periods=n_days, freq="D")
    return pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(np.maximum(pet, 0.0), 6),
        }
    )


def static_attributes(catchment: dict) -> dict:
    """What a model receives in static.json: the real row, flattened.

    The HydroATLAS and Caravan attributes keep their published names so an
    adapter maps them rather than guesses. Area and latitude sit at the top
    because every adapter needs them.
    """
    return {
        "catchment": catchment["id"],
        "area_km2": catchment["area_km2"],
        "latitude_deg": catchment["latitude_deg"],
        "longitude_deg": catchment["longitude_deg"],
        **{f"caravan.{k}": v for k, v in catchment["caravan"].items()},
        **{f"hydroatlas.{k}": v for k, v in catchment["hydroatlas"].items()},
    }

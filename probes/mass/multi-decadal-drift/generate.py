"""Repeated warm five-year weather blocks for long-horizon storage bounds."""

from __future__ import annotations

import math
import random

import pandas as pd

PERIOD_YEARS = 50
SPINUP_DAYS = 5 * 365
BLOCK_DAYS = 5 * 365
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
    rng = random.Random(seed)
    block = []
    for day in range(BLOCK_DAYS):
        season = math.cos(2 * math.pi * (day % 365) / 365)
        block.append({
            "pr": round((5.4 + 1.2 * season) * rng.uniform(0.8, 1.2), 6)
            if day % 3 == 0 else 0.0,
            "tas": round(15.0 + 5.0 * season, 6),
            "pet": round(1.0 - 0.3 * season, 6),
        })
    forcing = pd.DataFrame([block[i % BLOCK_DAYS] for i in range(N_STEPS)])
    # Model years have 365 steps; dates remain continuous across leap days.
    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing.insert(0, "time", time.strftime("%Y-%m-%d"))
    return forcing, dict(STATIC)

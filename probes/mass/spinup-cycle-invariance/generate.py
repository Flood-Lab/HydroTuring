"""One deterministic 365-day weather cycle, repeated for two spin-up lengths."""

from __future__ import annotations

import numpy as np
import pandas as pd

CYCLE_DAYS = 365
SHORT_CYCLES = 5
EXTRA_CYCLES = 3
TOTAL_CYCLES = SHORT_CYCLES + 1 + EXTRA_CYCLES
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


def generate(seed: int, variant: str = "short") -> tuple[pd.DataFrame, dict]:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")
    pr, tas, pet = _cycle(seed)
    forcing = pd.DataFrame({
        "time": pd.date_range("2001-01-01", periods=TOTAL_CYCLES * CYCLE_DAYS, freq="D").strftime("%Y-%m-%d"),
        "pr": np.round(np.tile(pr, TOTAL_CYCLES), 6),
        "tas": np.round(np.tile(tas, TOTAL_CYCLES), 6),
        "pet": np.round(np.tile(pet, TOTAL_CYCLES), 6),
    })
    evaluation_cycle = SHORT_CYCLES if variant == "short" else SHORT_CYCLES + EXTRA_CYCLES
    cycle_index = np.arange(TOTAL_CYCLES * CYCLE_DAYS) // CYCLE_DAYS
    forcing["_phase"] = np.where(cycle_index == evaluation_cycle, "evaluation", "spinup_or_tail")
    return forcing, dict(STATIC)

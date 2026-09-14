"""Seeded recharge and river-stage forcing for the groundwater probe."""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 2
SPINUP_DAYS = 90
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

STATIC = {
    "area_km2": 1.0,
    "aquifer_specific_yield": 0.20,
    "aquifer_storage_coefficient": 0.002,
    "aquifer_initial_head_m": 10.0,
    "river_stage_reference_m": 10.0,
    "river_conductance_m2_per_day": 150.0,
    "river_bottom_offset_m": 1.0,
}


def _ar1(rng: np.random.Generator, n: int, scale: float, persistence: float) -> np.ndarray:
    noise = np.zeros(n)
    innovation = rng.normal(0.0, scale, n)
    for index in range(1, n):
        noise[index] = persistence * noise[index - 1] + innovation[index]
    return noise


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    day = np.arange(N_STEPS)
    doy = day % 365

    recharge = np.maximum(
        0.02
        + 0.018 * np.sin(2.0 * np.pi * (doy - 80) / 365.0)
        + _ar1(rng, N_STEPS, scale=0.004, persistence=0.85),
        0.001,
    )
    stage = (
        10.0
        + 1.8 * np.sin(2.0 * np.pi * (day - SPINUP_DAYS) / 90.0)
        + _ar1(rng, N_STEPS, scale=0.025, persistence=0.90)
    )

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "gw_recharge": np.round(recharge, 6),
            "sw_stage_m": np.round(stage, 6),
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    frame, _ = generate(20260912)
    print(f"{len(frame)} daily rows: {SPINUP_DAYS} spinup + {PERIOD_YEARS} scored years")
    print(f"recharge: {frame['gw_recharge'].sum() / (N_STEPS / 365):.1f} mm/year")

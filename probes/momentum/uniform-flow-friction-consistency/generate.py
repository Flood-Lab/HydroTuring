"""Constant forcing for a steady uniform-flow momentum check."""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 3
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    """Return one wet, temperate, deterministic steady-flow case per seed."""
    rng = np.random.default_rng(seed)
    pr = float(rng.uniform(4.5, 7.0))
    pet = float(rng.uniform(1.2, min(2.8, pr - 1.0)))
    tas = float(rng.uniform(12.0, 19.0))

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame({
        "time": time.strftime("%Y-%m-%d"),
        "pr": np.full(N_STEPS, pr),
        "tas": np.full(N_STEPS, tas),
        "pet": np.full(N_STEPS, pet),
    })

    area_km2 = float(rng.uniform(80.0, 320.0))
    width_m = float(rng.uniform(90.0, 160.0))
    static = {
        # Hydrologic attributes used by the independent physical baselines.
        # area_km2 is deliberately not a required probe input: it helps those
        # models produce `dis`, but the residual consumes `dis` directly.
        "area_km2": area_km2,
        "soil_capacity_mm": float(rng.uniform(220.0, 380.0)),
        "canopy_capacity_mm": float(rng.uniform(1.5, 3.0)),
        "degree_day_factor_mm_per_C_day": 3.2,
        "baseflow_coefficient": 0.006,
        "snow_threshold_degC": 0.0,
        "latitude_deg": 38.0,
        # The depth produced by these ranges is small relative to width. That
        # makes the wide-channel ratings in the physical baselines an honest
        # approximation while the criterion retains the exact rectangular R_h.
        "width_m": width_m,
        "bed_elevation_m": float(rng.uniform(25.0, 180.0)),
        "slope": float(rng.uniform(5.0e-4, 3.0e-3)),
        "manning_n": float(rng.uniform(0.025, 0.045)),
        "reach_length_m": float(rng.uniform(3000.0, 8000.0)),
    }
    return forcing, static

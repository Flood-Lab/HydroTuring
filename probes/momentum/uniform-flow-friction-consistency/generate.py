"""Three constant-flow plateaus for a steady momentum-balance check."""

from __future__ import annotations

import numpy as np
import pandas as pd

SPINUP_DAYS = 365
BASEFLOW_COEFFICIENT_PER_DAY = 0.006
STEADY_DAYS = 90
SETTLING_TIME_CONSTANTS = 5.0
# The slowest must-pass store has tau ~= 1 / k. Give it five time constants
# before the final scored block starts, rather than relying on the five gate
# seeds to happen to settle inside a one-year plateau.
SETTLING_DAYS = int(np.ceil(
    SETTLING_TIME_CONSTANTS / BASEFLOW_COEFFICIENT_PER_DAY
))
PLATEAU_DAYS = SETTLING_DAYS + STEADY_DAYS
PERIOD_DAYS = 3 * PLATEAU_DAYS
N_STEPS = PERIOD_DAYS + SPINUP_DAYS


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    """Return low, medium and high steady-flow plateaus for one reach."""
    rng = np.random.default_rng(seed)
    pet = float(rng.uniform(1.2, 2.4))
    tas = float(rng.uniform(12.0, 19.0))

    # The spinup holds the low-flow forcing for a full year.  The scored
    # record then visits three separately labelled plateaus.  The
    # labels are stripped before the model runs and are visible only to the
    # criterion, so a model cannot identify which rows are scored.
    effective_low = float(rng.uniform(1.0, 1.8))
    effective = {
        "low": effective_low,
        "medium": 2.0 * effective_low,
        "high": 4.0 * effective_low,
    }
    labels = np.concatenate([
        np.full(SPINUP_DAYS, "spinup", dtype=object),
        np.full(PLATEAU_DAYS, "low", dtype=object),
        np.full(PLATEAU_DAYS, "medium", dtype=object),
        np.full(PLATEAU_DAYS, "high", dtype=object),
    ])
    pr = np.concatenate([
        np.full(SPINUP_DAYS, pet + effective["low"]),
        np.full(PLATEAU_DAYS, pet + effective["low"]),
        np.full(PLATEAU_DAYS, pet + effective["medium"]),
        np.full(PLATEAU_DAYS, pet + effective["high"]),
    ])

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame({
        "time": time.strftime("%Y-%m-%d"),
        "pr": pr,
        "tas": np.full(N_STEPS, tas),
        "pet": np.full(N_STEPS, pet),
        "_plateau": labels,
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
        "baseflow_coefficient": BASEFLOW_COEFFICIENT_PER_DAY,
        "snow_threshold_degC": 0.0,
        "latitude_deg": 38.0,
        # The depth produced by these ranges is small relative to width. That
        # makes the wide-channel ratings in the physical baselines an honest
        # approximation while the criterion retains the exact rectangular R_h.
        "width_m": width_m,
        "cross_section_shape": "rectangular",
        "bed_elevation_m": float(rng.uniform(25.0, 180.0)),
        "slope": float(rng.uniform(5.0e-4, 3.0e-3)),
        "manning_n": float(rng.uniform(0.025, 0.045)),
        "reach_length_m": float(rng.uniform(3000.0, 8000.0)),
    }
    return forcing, static

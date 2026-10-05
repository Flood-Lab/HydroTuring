from __future__ import annotations
import numpy as np
import pandas as pd

PERIOD_DAYS = 20
SPINUP_DAYS = 2
N_STEPS = PERIOD_DAYS + SPINUP_DAYS
VARIANTS = ("control", "datum_up", "datum_down")

STATIC = {
    "area_km2": 1.0,
    "aquifer_specific_yield": 0.20,
    "aquifer_storage_coefficient": 0.002,
    "aquifer_initial_head_m": 100.0,
    "river_conductance_m2_per_day": 1000.0,
    "river_bottom_offset_m": 1.0,
    "aquifer_top_m": 101.0,
    "aquifer_bottom_m": 90.0,
}

def generate(seed: int, variant: str = "control") -> tuple[pd.DataFrame, dict]:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected {VARIANTS}")
    rng = np.random.default_rng(seed)
    # All physical forcing is drawn before the datum branch.  The small jitter
    # makes the case deterministic per seed without changing the pulse logic.
    jitter = rng.normal(0.0, 0.005, N_STEPS)
    stage = np.full(N_STEPS, 100.0)
    stage[2:6] = 100.75 + jitter[2:6]
    stage[6:10] = 100.0 + jitter[6:10]
    stage[10:14] = 99.25 + jitter[10:14]
    stage[14:] = 100.0 + jitter[14:]
    delta = 0.0 if variant == "control" else (250.0 if variant == "datum_up" else -250.0)
    forcing = pd.DataFrame({
        "time": pd.date_range("2000-01-01", periods=N_STEPS, freq="D").strftime("%Y-%m-%d"),
        "gw_recharge": np.zeros(N_STEPS),
        "sw_stage_m": np.round(stage + delta, 6),
    })
    static = dict(STATIC)
    static["aquifer_initial_head_m"] += delta
    static["aquifer_top_m"] += delta
    static["aquifer_bottom_m"] += delta
    return forcing, static

"""Source-free confined aquifer, with a stationary or shifted river head."""

from __future__ import annotations

import numpy as np
import pandas as pd

SPINUP_DAYS = 2
PERIOD_DAYS = 90
VARIANTS = ("control", "raised", "lowered")
MINIMUM_CLEARANCE_M = 1.0


def check_geometry(static: dict, minimum_stage: float, maximum_stage: float) -> None:
    """Keep all heads confined and above the adapter's moving riverbed."""
    initial = float(static["aquifer_initial_head_m"])
    top = float(static["aquifer_top_m"])
    bottom = float(static["aquifer_bottom_m"])
    offset = float(static["river_bottom_offset_m"])
    if not np.isfinite((initial, top, bottom, offset, minimum_stage, maximum_stage)).all() or not bottom < top:
        raise ValueError("aquifer geometry must be finite with bottom below top")
    minimum_head = min(initial, minimum_stage)
    maximum_bed = maximum_stage - offset
    if minimum_head <= top:
        raise ValueError("the prescribed head envelope must remain confined")
    if minimum_head - maximum_bed < MINIMUM_CLEARANCE_M:
        raise ValueError("riverbed-to-head clearance must be at least 1 m")


def generate(seed: int, variant: str = "control") -> tuple[pd.DataFrame, dict]:
    """Draw every attribute before selecting the sole intervention: river head."""
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected {VARIANTS}")
    rng = np.random.default_rng(seed)
    area = float(rng.uniform(0.5, 2.0))
    storativity = float(rng.uniform(0.002, 0.01))
    head_change = float(rng.uniform(0.5, 1.5))
    timescale_days = float(rng.uniform(2.0, 10.0))
    top = float(rng.uniform(50.0, 200.0))
    thickness = float(rng.uniform(10.0, 40.0))
    specific_yield = float(rng.uniform(0.1, 0.3))
    initial = top + 5.0
    static = {
        "area_km2": area,
        "aquifer_storage_coefficient": storativity,
        "aquifer_specific_yield": specific_yield,
        "aquifer_initial_head_m": initial,
        "aquifer_top_m": top,
        "aquifer_bottom_m": top - thickness,
        "river_conductance_m2_per_day": area * 1.0e6 * storativity / timescale_days,
        "river_bottom_offset_m": 4.0,
    }
    # Check even the highest moving bed against the lowest admissible head.
    check_geometry(static, initial - head_change, initial + head_change)
    stage = np.full(SPINUP_DAYS + PERIOD_DAYS, initial)
    sign = {"control": 0.0, "raised": 1.0, "lowered": -1.0}[variant]
    stage[SPINUP_DAYS:] += sign * head_change
    forcing = pd.DataFrame({
        "time": pd.date_range("2000-01-01", periods=len(stage), freq="D").strftime("%Y-%m-%d"),
        "gw_recharge": np.zeros(len(stage)),
        "sw_stage_m": stage,
    })
    return forcing, static

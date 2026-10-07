"""A passive aquifer under a fixed, seeded river stage, with and without a
direct recharge pulse on the first two scored days.

The river conductance is the one mass/groundwater-datum-invariance uses. At
100 m2/day the pulse raised the exchange by about 0.08 mm over the record,
under the 0.1 mm tolerance, so no reversal of it could be seen; at 1000 m2/day
a reversed increment is 0.8 mm, eight times the tolerance (#152).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

N_STEPS = 20
PULSE_ROWS = slice(2, 4)
PULSE_MM_PER_DAY = 5.0
VARIANTS = ("control", "recharge_added")
STATIC = {
    "area_km2": 1.0,
    "aquifer_specific_yield": 0.20,
    "aquifer_storage_coefficient": 0.002,
    "aquifer_initial_head_m": 100.0,
    "aquifer_top_m": 110.0,
    "aquifer_bottom_m": 90.0,
    "river_conductance_m2_per_day": 1000.0,
    "river_bottom_offset_m": 1.0,
}


def generate(seed: int, variant: str = "control"):
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected {VARIANTS}")
    rng = np.random.default_rng(seed)
    # Drawn before the variant matters, so both runs see the same river.
    stage = 100.0 + 0.01 * rng.normal(size=N_STEPS)
    recharge = np.zeros(N_STEPS)
    if variant == "recharge_added":
        recharge[PULSE_ROWS] = PULSE_MM_PER_DAY
    forcing = pd.DataFrame(
        {
            "time": pd.date_range("2000-01-01", periods=N_STEPS, freq="D").strftime("%Y-%m-%d"),
            "gw_recharge": recharge,
            "sw_stage_m": np.round(stage, 6),
        }
    )
    return forcing, dict(STATIC)

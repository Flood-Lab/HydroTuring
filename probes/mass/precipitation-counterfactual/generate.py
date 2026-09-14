"""Forcing generator for mass/precipitation-counterfactual: one seed, four rainfalls.

The weather is drawn once per seed; a variant only scales the scored
precipitation, so `tas`, `pet`, wet-day occurrence, event order and the
spinup are shared exactly between variants.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 10
SPINUP_DAYS = 365
N_STEPS = int(PERIOD_YEARS * 365) + SPINUP_DAYS

# Multiplier on the scored precipitation (README, amplitude table).
PERTURBATION = {"control": 1.00, "wetter20": 1.20, "wetter10": 1.10, "drier20": 0.80}
VARIANTS = tuple(PERTURBATION)

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 40.0,
}


def generate(seed: int, variant: str = "control") -> tuple[pd.DataFrame, dict]:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")

    rng = np.random.default_rng(seed)
    day = np.arange(N_STEPS)
    doy = day % 365

    # Every draw happens before the variant is applied.
    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(N_STEPS) < p_wet
    pr = np.where(wet, rng.gamma(shape=0.7, scale=13.0, size=N_STEPS), 0.0)

    seasonal = 9.0 - 12.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(N_STEPS)
    innovation = rng.normal(0.0, 2.6, N_STEPS)
    for t in range(1, N_STEPS):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    tas = seasonal + noise

    daylength = 1.0 + 0.35 * np.cos(2 * np.pi * (doy - 172) / 365)
    pet = np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength

    # Scale depths in the scored window only: every variant enters it from
    # the same state, and occurrence and timing are held fixed.
    pr[SPINUP_DAYS:] *= PERTURBATION[variant]

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    control = generate(20260903)[0]["pr"][SPINUP_DAYS:].sum()
    for variant in VARIANTS:
        scored = generate(20260903, variant)[0]["pr"][SPINUP_DAYS:].sum()
        print(f"{variant:9s} {scored / PERIOD_YEARS:6.0f} mm/yr  ({scored - control:+.0f} mm)")

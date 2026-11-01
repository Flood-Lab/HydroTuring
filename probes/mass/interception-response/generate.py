"""Forcing generator for mass/interception-response.

Two variants of the same three-year temperate record, differing only in the
declared canopy capacity: 2 mm (a sparse canopy) and 10 mm (a closed one).
The weather is drawn once, before the variant branch, so both runs see
byte-identical rain, temperature and potential evaporation. The record ends
with a dry spell long enough that both canopies drain, which means the extra
water a larger canopy intercepts has to show up in the event peaks rather
than in a final store.

10 mm is chosen as the upper end because it is representable by every model
that currently declares `canopy_capacity_mm`; LISFLOOD's conversion tops out
near 11.7 mm. The change is still a factor of five.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 3
SPINUP_DAYS = 365
DRY_DAYS = 90
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 38.0,
}

CAPACITY_LOW_MM = 2.0
CAPACITY_HIGH_MM = 10.0


def _weather(rng: np.random.Generator, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Temperate, rain-dominated daily weather with a terminal dry spell."""
    day = np.arange(n)
    doy = day % 365
    p_wet = 0.27 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(n) < p_wet
    pr = np.where(wet, rng.gamma(shape=0.7, scale=14.0, size=n), 0.0)

    seasonal = 15.0 - 7.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(n)
    innovation = rng.normal(0.0, 2.0, n)
    for t in range(1, n):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    tas = seasonal + noise
    daylength = 1.0 + 0.35 * np.cos(2 * np.pi * (doy - 172) / 365)
    pet = np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength

    # The last dry spell is scored: both canopies have to empty before the
    # window ends, so the paired response cannot hide in a final store.
    pr = pr.copy()
    pr[-DRY_DAYS:] = 0.0
    return pr, tas, pet


def _frame(pr, tas, pet) -> pd.DataFrame:
    time = pd.date_range("2000-01-01", periods=len(pr), freq="D")
    return pd.DataFrame({
        "time": time.strftime("%Y-%m-%d"),
        "pr": np.round(pr, 6),
        "tas": np.round(tas, 6),
        "pet": np.round(pet, 6),
    })


def generate(seed: int, variant: str | None = None) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    pr, tas, pet = _weather(rng, N_STEPS)
    static = dict(STATIC)
    if variant == "capacity_low":
        static["canopy_capacity_mm"] = CAPACITY_LOW_MM
    elif variant == "capacity_high":
        static["canopy_capacity_mm"] = CAPACITY_HIGH_MM
    elif variant is not None:
        raise ValueError(f"unknown variant {variant!r}")
    return _frame(pr, tas, pet), static


if __name__ == "__main__":
    frame, static = generate(20261021)
    print(f"{N_STEPS} steps; {frame['pr'][:SPINUP_DAYS].sum():.0f} mm rain in spinup, "
          f"{frame['pr'][SPINUP_DAYS:].sum():.0f} mm in the scored years; "
          f"canopy capacity {static['canopy_capacity_mm']} mm")

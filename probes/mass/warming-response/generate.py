"""Forcing generator for mass/warming-response: the same rain, warmer or
cooler air.

Deterministic given a seed and a variant. All variants carry the same
precipitation, drawn once from the seed before anything branches on the
variant. The `warmer` variant adds three degrees to the air temperature
over the scored years and the `cooler` variant removes three, each
recomputing potential evaporation from its own air, because that is how
the demand for water changes when the air warms or cools. The spinup is
left alone, so every variant enters the window from the same state and the
difference between the runs is the temperature shift and nothing else.

The catchment is temperate and rain-dominated: mean air temperature around
15 degC and rarely below freezing, so that snow timing cannot confound the
sign of the response. Precipitation targets roughly 830 mm a year, potential
evaporation roughly 950, which puts the catchment on the water-limited side
where evaporation is free to rise with demand.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 10
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

VARIANTS = ("control", "warmer", "cooler")
SHIFT_DEGC = {"control": 0.0, "warmer": 3.0, "cooler": -3.0}

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 38.0,
}


def potential_evaporation(tas: np.ndarray, doy: np.ndarray) -> np.ndarray:
    """Temperature-driven demand with a daylength factor, as the other probes use."""
    daylength = 1.0 + 0.35 * np.cos(2 * np.pi * (doy - 172) / 365)
    return np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength


def generate(seed: int, variant: str = "control") -> tuple[pd.DataFrame, dict]:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")

    rng = np.random.default_rng(seed)
    day = np.arange(N_STEPS)
    doy = day % 365

    # Everything drawn from the seed happens here, before the branch.
    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(N_STEPS) < p_wet
    pr = np.where(wet, rng.gamma(shape=0.7, scale=13.0, size=N_STEPS), 0.0)

    seasonal = 15.0 - 7.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(N_STEPS)
    innovation = rng.normal(0.0, 2.0, N_STEPS)
    for t in range(1, N_STEPS):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    tas = seasonal + noise

    # --- the perturbation ---------------------------------------------------
    # Both ways. A model can carry "hotter means drier" as a one-way rule,
    # saturate in one direction, or clip at zero; pushing the driver both
    # ways and scoring each against the control is what tells.
    if SHIFT_DEGC[variant]:
        tas = tas.copy()
        tas[SPINUP_DAYS:] += SHIFT_DEGC[variant]
    pet = potential_evaporation(tas, doy)

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
    control, _ = generate(20260903, "control")
    years = PERIOD_YEARS
    scored = slice(SPINUP_DAYS, None)
    print(f"{N_STEPS} steps; precipitation {control['pr'][scored].sum() / years:.0f} mm/yr in every variant")
    print(f"control: tas {control['tas'][scored].mean():.1f} degC, "
          f"pet {control['pet'][scored].sum() / years:.0f} mm/yr, "
          f"{int((control['tas'] < 0).sum())} sub-zero days")
    for name in ("warmer", "cooler"):
        other, _ = generate(20260903, name)
        delta = (other["pet"][scored].sum() - control["pet"][scored].sum()) / years
        print(f"{name}:  tas {other['tas'][scored].mean():.1f} degC, "
              f"pet {other['pet'][scored].sum() / years:.0f} mm/yr ({delta:+.0f} mm/yr of demand), "
              f"{int((other['tas'] < 0).sum())} sub-zero days")
        assert control["pr"].equals(other["pr"]), "the rain must be identical"
        assert control["tas"][:SPINUP_DAYS].equals(other["tas"][:SPINUP_DAYS]), "spinup must be identical"
    print("rain identical, spinup identical, as they should be")

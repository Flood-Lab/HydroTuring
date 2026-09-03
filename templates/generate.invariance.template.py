"""Forcing generator for {id}: one weather sequence, written two ways.

Deterministic given a seed and a variant. Both variants carry the same
weather; the transform changes something the physics does not depend on. The
paired criterion then asks whether the model's answer moved.

The transform here is the time origin: the same series dated from a different
year. A catchment does not know what year it is, so nothing may change. Other
transforms and the trap to avoid are documented in the probe draft.

Requirements the harness enforces:
  - generate(seed, variant) returns (DataFrame, dict)
  - the frame has a 'time' column
  - it has exactly PERIOD_YEARS * 365 + SPINUP_DAYS rows
  - the same (seed, variant) produces byte-identical output

Run this file directly to confirm the weather really is identical.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = {period_years}
SPINUP_DAYS = {spinup_days}
N_STEPS = int(PERIOD_YEARS * 365) + SPINUP_DAYS

VARIANTS = ("control", "shifted")

# Whole years, so that day-of-year and therefore the seasonal cycle line up.
# A shift that is not a multiple of 365 moves the weather relative to the
# season and is a different case, not the same one relabelled.
ORIGINS = {{"control": "2000-01-01", "shifted": "1960-01-01"}}

STATIC = {{
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 40.0,
}}


def generate(seed: int, variant: str = "control") -> tuple[pd.DataFrame, dict]:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {{variant!r}}; expected one of {{VARIANTS}}")

    rng = np.random.default_rng(seed)
    day = np.arange(N_STEPS)
    doy = day % 365

    # Drawn once from the seed and shared by both variants. Nothing below this
    # point may depend on `variant` except the dates and, for a scaling
    # transform, the attribute being scaled.
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

    # --- the transform ------------------------------------------------------
    # TODO: this is the probe. The time origin is the cheapest defensible
    # symmetry and a good first one. For a scaling transform instead, scale
    # the attribute here and declare the expected factor under `scaled` in
    # probe.yaml — for example doubling `area_km2`, under which every depth in
    # millimetres is unchanged and `dis` in m3 s-1 doubles.
    static = dict(STATIC)
    time = pd.date_range(ORIGINS[variant], periods=N_STEPS, freq="D")

    forcing = pd.DataFrame(
        {{
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
        }}
    )
    return forcing, static


if __name__ == "__main__":
    control, static_c = generate(20260903, "control")
    shifted, static_s = generate(20260903, "shifted")

    print(f"control starts {{control['time'].iloc[0]}}, "
          f"shifted starts {{shifted['time'].iloc[0]}}")

    same = all(
        control[column].equals(shifted[column]) for column in ("pr", "tas", "pet")
    )
    print("weather identical in both variants" if same else
          "WARNING: the weather differs between variants; the transform is "
          "supposed to change nothing the physics depends on")

    if control["time"].equals(shifted["time"]):
        print("WARNING: the dates are identical, so no transform was applied "
              "and the criterion will pass for any model at all.")

    # Leap days are the classic way a date shift stops being a pure relabel.
    for frame, name in ((control, "control"), (shifted, "shifted")):
        stamps = pd.to_datetime(frame["time"])
        if stamps.dt.strftime("%m-%d").eq("02-29").any():
            print(f"note: {{name}} contains a 29 February; the 365-day record "
                  "means day-of-year drifts against the calendar")

"""Forcing generator for {id}: the same seed, perturbed and unperturbed.

Deterministic given a seed and a variant. The harness runs the model once per
variant and hands the paired criterion both results, so what is being compared
is one catchment under two rainfalls rather than two different draws.

The rule that makes the comparison mean anything: **draw the weather once,
then perturb it**. Re-drawing under a different seed changes the weather as
well as the rainfall, and the difference between the runs no longer isolates
the thing you perturbed.

Requirements the harness enforces:
  - generate(seed, variant) returns (DataFrame, dict)
  - the frame has a 'time' column
  - it has exactly PERIOD_YEARS * 365 + SPINUP_DAYS rows
  - the same (seed, variant) produces byte-identical output

Run this file directly to check that the variants differ only where you meant.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = {period_years}
SPINUP_DAYS = {spinup_days}
N_STEPS = int(PERIOD_YEARS * 365) + SPINUP_DAYS

VARIANTS = ("control", "wetter")

# How much more rain the perturbed variant gets. Large enough that the
# response is well clear of a model's own numerical noise, small enough that
# the catchment is not driven into a different regime entirely — at which
# point you are testing extrapolation, not response.
PERTURBATION = 1.20

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

    # Everything drawn from the seed happens here, before the branch, so the
    # two variants share one weather sequence exactly.
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

    # --- the perturbation ---------------------------------------------------
    # TODO: this is the probe. Scaling wet-day depths is the simplest
    # counterfactual and a reasonable default. Others worth considering, each
    # of which asks a different question: more wet days at the same intensity,
    # the same total delivered in fewer larger events, or a warmer year at
    # unchanged precipitation, which perturbs demand rather than supply.
    #
    # Perturb the scored window only. Leaving spinup alone means both variants
    # enter the window from the same state, so the difference between them is
    # the perturbation and not a different starting storage.
    if variant == "wetter":
        pr[SPINUP_DAYS:] = pr[SPINUP_DAYS:] * PERTURBATION

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {{
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
        }}
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    control, _ = generate(20260903, "control")
    wetter, _ = generate(20260903, "wetter")
    years = N_STEPS / 365

    added = wetter["pr"].sum() - control["pr"].sum()
    print(f"control: {{control['pr'].sum() / years:.0f}} mm/yr")
    print(f"wetter:  {{wetter['pr'].sum() / years:.0f}} mm/yr  (+{{added:.0f}} mm total)")

    for column in ("tas", "pet"):
        if not control[column].equals(wetter[column]):
            print(f"WARNING: '{{column}}' differs between the variants. Only the "
                  "perturbed variable may differ, or the comparison is confounded.")
    if control["pr"][:SPINUP_DAYS].equals(wetter["pr"][:SPINUP_DAYS]):
        print("spinup is identical in both variants, as it should be")

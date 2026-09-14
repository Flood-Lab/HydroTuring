"""Forcing generator for mass/human-abstraction: the same weather, with and
without a prescribed irrigation withdrawal.

Deterministic given a seed and a variant. The harness runs the model once per
variant and hands the paired criterion both results, so what is being compared
is one catchment with and without the human term, not two different draws.

The rule that makes the comparison mean anything: **draw the weather and fix
the abstraction schedule once, then branch**. Re-drawing under a different
seed changes the weather as well as the withdrawal, and the difference
between the runs no longer isolates the human term.

The withdrawal is written as a visible `abstr` column (mm/day, NET of return
flow) alongside time/pr/tas/pet, so the model can read what it is expected to
honour and the criterion can see the same prescribed value. The natural
variant carries the column as zeros: the driver exists in both runs and only
its value changes.

Requirements the harness enforces:
  - generate(seed, variant) returns (DataFrame, dict)
  - the frame has a 'time' column
  - it has exactly PERIOD_YEARS * 365 + SPINUP_DAYS rows
  - the same (seed, variant) produces byte-identical output

Run this file directly to check that the variants differ only where intended.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 10
SPINUP_DAYS = 365
N_STEPS = int(PERIOD_YEARS * 365) + SPINUP_DAYS

VARIANTS = ("natural", "irrigated")

# The abstraction schedule. Withdrawal is concentrated in the growing season
# (roughly May through September), shaped as a smooth bump rather than a
# step, because real irrigation demand follows the crop calendar. The peak
# rate and the net-of-return fraction are the probe's physics to defend; the
# values below give a catchment-mean seasonal total of ~38 mm/yr against
# ~830 mm/yr of rain — a modest irrigated fraction, not an irrigation
# district. Two constraints set the amplitude from opposite sides. From
# below: a blind model must fail by a wide margin, and at ~380 mm over the
# record the residual is twenty times the 5 percent tolerance. From above:
# the withdrawal must stay honourable — a catchment that cannot supply the
# prescription forces an honest model to under-remove and fail its own probe.
# At 2.5 mm/day peak the soil column empties every summer and even the
# reference model under-removes ~27-30% of the prescription; at 0.5 mm/day the
# honest model's residual is ~0.0% on the scored seeds (below ~1% on arbitrary
# seeds), far under the tolerance.
SEASON_START_DOY = 121   # May 1
SEASON_END_DOY = 273     # Sep 30
PEAK_ABSTR_MM_DAY = 0.5  # net peak rate at midsummer

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 40.0,
}


def abstraction_schedule(doy: np.ndarray) -> np.ndarray:
    """Net seasonal withdrawal, mm/day: a sin^2 bump over the growing season,
    zero outside it. sin^2 keeps the schedule and its slope continuous at the
    season boundaries, so no model can trip on a discontinuity in the driver."""
    in_season = (doy >= SEASON_START_DOY) & (doy <= SEASON_END_DOY)
    phase = (doy - SEASON_START_DOY) / (SEASON_END_DOY - SEASON_START_DOY)
    bump = np.sin(np.pi * np.clip(phase, 0.0, 1.0)) ** 2
    return np.where(in_season, PEAK_ABSTR_MM_DAY * bump, 0.0)


def generate(seed: int, variant: str = "natural") -> tuple[pd.DataFrame, dict]:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")

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

    # The schedule is deterministic given the calendar, so it is fixed before
    # the branch as well; only whether the driver is *applied* depends on the
    # variant. The spinup carries the schedule too: both variants then enter
    # the scored window from states that differ only by the human term's own
    # history, which is the physical situation — an irrigated catchment does
    # not start the decade with a natural catchment's stores.
    abstr = abstraction_schedule(doy) if variant == "irrigated" else np.zeros(N_STEPS)

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
            "abstr": np.round(abstr, 6),
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    natural, _ = generate(20260911, "natural")
    irrigated, _ = generate(20260911, "irrigated")
    years = N_STEPS / 365

    print(f"natural:   {natural['pr'].sum() / years:.0f} mm/yr rain, "
          f"{natural['abstr'].sum():.0f} mm abstracted over the record")
    print(f"irrigated: {irrigated['pr'].sum() / years:.0f} mm/yr rain, "
          f"{irrigated['abstr'].sum():.0f} mm abstracted over the record "
          f"(~{irrigated['abstr'].sum() / years:.0f} mm/yr in season)")

    for column in ("pr", "tas", "pet"):
        if not natural[column].equals(irrigated[column]):
            print(f"WARNING: '{column}' differs between the variants. Only the "
                  "abstraction may differ, or the comparison is confounded.")
    if natural["abstr"].sum() == 0 and irrigated["abstr"].sum() > 0:
        print("the variants differ only in the abstraction, as they should")

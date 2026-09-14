"""Paired forcing for energy/evaporative-partition: one summer without rain.

Deterministic given a seed and a variant. Both variants are drawn once, so
net radiation, air temperature and evaporative demand are byte-identical
between them. The `drought` variant sets precipitation to zero across a
120-day window in the second scored summer and changes nothing else.

A note on what this case is and is not. It is an idealised **seasonal
drydown**, not a flash drought. The flash drought literature keys its
definition on the rate of intensification rather than on duration (Otkin et
al. 2018, BAMS 99, 911-919), with onset criteria on the order of one to six
pentads and duration caps that exist precisely to separate these events from
seasonal drought (Christian et al. 2024, WIREs Water 11, e1714). A 120-day
rainfall-free window is on the wrong side of that line and is not claimed to
be on the right one. The repartitioning being tested is the same one that
governs flash drought onset; the case is simply longer than a flash drought,
because the identity needs a signal large enough to score.

Holding net radiation fixed is the whole design. It turns what would
otherwise be an expectation about how a surface ought to respond into an
identity: the energy a drying surface stops putting into evaporation has
nowhere to go but the sensible and ground fluxes, and the three changes must
sum to zero because the radiation driving them did not change. There is
nothing left to tune and no atmospheric feedback to model. What is being
tested is the model's internal partition operator, not the atmosphere's
response to it.

Real flash droughts do not come with radiation held fixed, and the coupled
warming would make the shift larger rather than smaller. Removing it from the
case removes the confounder, not the phenomenon.

The window is placed and sized so that the catchment crosses out of the
energy-limited regime inside it. Ninety days turned out not to be enough: the
reference model's latent flux fell by as little as 8 percent of the window's
net radiation on some seeds, because the soil never dropped below the stress
threshold in either variant and the rain that was removed left through runoff
and storage instead of through evaporation. At a hundred and twenty days and
a soil capacity of 120 mm the fall is between 29 and 43 percent across seeds,
and the drydown is a regime transition rather than a perturbation. The
capacity is smaller than the one energy/latent-heat-et-consistency uses, and
deliberately so: a bucket that cannot empty cannot become water-limited.

The `_perturbed` column marks the scoring window and `_intensification` the
stretch inside it where the drying actually bites. That stretch is not the
start of the window, and finding out why is worth recording: the control run
keeps receiving rain, so the two runs do not diverge until the drought run
has spent its soil buffer. Across seeds the reference model has moved less
than 1 percent of the shift by day 60 and half of it only after day 105.
Marking days 76 to 120, which carry about 98 percent of it, is what lets a
report state an intensification rate instead of a four-month total.

The marked stretch is reported and never gated. Its placement comes from the
reference model, which is how this suite already locates evaluation windows,
and it is stable across seeds; but a threshold whose position depends on the
reference model is a threshold worth not having.

Columns beginning with an underscore are stripped before the model sees the
forcing, so the model is never told which stretch it is being judged on.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 3
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

VARIANTS = ("control", "drought")
# The second scored summer: 1 June of year two.
DROUGHT_START = SPINUP_DAYS + 365 + 151
DROUGHT_DAYS = 120
# The stretch in which the reference model's shift actually accrues, marked
# separately and reported rather than scored. Measured, not assumed: across
# seeds the reference model has moved less than 1 percent of the shift by day
# 60, 10 percent by day 82 to 85, and half of it only after day 105. Days 76
# to 120 carry about 98 percent of it. Marking that stretch is what lets a
# report state an intensification rate rather than only a four-month total.
INTENSIFY_START = 75
INTENSIFY_DAYS = 45

PT_ALPHA = 1.26
GAMMA_KPA_PER_C = 0.0665
GROUND_SHARE = 0.10

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 120.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 40.0,
}


def _saturation_slope(tas: np.ndarray) -> np.ndarray:
    es = 0.6108 * np.exp(17.27 * tas / (tas + 237.3))
    return 4098.0 * es / np.square(tas + 237.3)


def _weather(rng: np.random.Generator, n: int):
    day = np.arange(n)
    doy = day % 365

    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(n) < p_wet
    pr = np.where(wet, rng.gamma(shape=0.7, scale=13.0, size=n), 0.0)

    seasonal = 9.0 - 12.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(n)
    innovation = rng.normal(0.0, 2.6, n)
    for t in range(1, n):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    tas = seasonal + noise

    cloud = np.zeros(n)
    shock = rng.normal(0.0, 0.45, n)
    for t in range(1, n):
        cloud[t] = 0.55 * cloud[t - 1] + shock[t]
    clearness = np.clip(0.72 + 0.18 * cloud, 0.35, 1.0)

    rsds = (245.0 - 135.0 * np.cos(2 * np.pi * (doy - 172) / 365)) * clearness
    albedo = np.where(tas < STATIC["snow_threshold_degC"], 0.60, 0.23)
    longwave_out = np.clip(32.0 + 1.2 * tas, 15.0, 95.0) * (0.55 + 0.45 * clearness)
    rn = (1.0 - albedo) * rsds - longwave_out

    slope = _saturation_slope(tas)
    lam = 2.501e6 - 2361.0 * tas
    available = np.maximum(rn * (1.0 - GROUND_SHARE), 0.0)
    pet = PT_ALPHA * (slope / (slope + GAMMA_KPA_PER_C)) * available * 86400.0 / lam
    return pr, tas, pet, rn


def generate(seed: int, variant: str = "control") -> tuple[pd.DataFrame, dict]:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")

    rng = np.random.default_rng(seed)
    pr, tas, pet, rn = _weather(rng, N_STEPS)

    window = np.zeros(N_STEPS)
    window[DROUGHT_START : DROUGHT_START + DROUGHT_DAYS] = 1.0
    intensification = np.zeros(N_STEPS)
    first = DROUGHT_START + INTENSIFY_START
    intensification[first : first + INTENSIFY_DAYS] = 1.0

    if variant == "drought":
        pr = pr.copy()
        pr[DROUGHT_START : DROUGHT_START + DROUGHT_DAYS] = 0.0

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
            "rn": np.round(rn, 6),
            "_perturbed": window,
            "_intensification": intensification,
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    control, _ = generate(20260903, "control")
    drought, _ = generate(20260903, "drought")
    w = control["_perturbed"].to_numpy() > 0
    for name in ("tas", "pet", "rn"):
        assert (control[name].to_numpy() == drought[name].to_numpy()).all(), name
    removed = float(control.loc[w, "pr"].sum())
    print(
        f"{N_STEPS} steps; tas, pet and rn identical between variants; "
        f"the drydown window is {int(w.sum())} days "
        f"(intensification stretch "
        f"{int((control['_intensification'].to_numpy() > 0).sum())} days) from "
        f"{control.loc[w, 'time'].iloc[0]} and removes {removed:.0f} mm of rain; "
        f"it receives {control.loc[w, 'rn'].sum():.0f} W m-2 day of net radiation"
    )

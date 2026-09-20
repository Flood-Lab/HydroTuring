"""Seeded synthetic forcing for the declared-exchange plausibility probe.

Deterministic given a seed and a variant, and never committed as data. The
harness runs the model once per variant and hands the paired criterion every
result, so what is compared is the same model under the same weather with one
thing changed: the external head it was given.

Contract:
  - generate(seed, variant) returns (DataFrame, dict)
  - the first variant is the control; the others shift `gwh` by a constant,
    **from the first scored step onward**: the spinup is byte-identical across
    variants, so every model enters the scored record in the same state and
    the response to the shift falls inside the window that is scored
  - the same (seed, variant) produces byte-identical output, and every column
    other than `gwh` is byte-identical across variants

Ten years of daily weather on a catchment that could plausibly be losing or
gaining water across its boundary, so that a model with a regional groundwater
term has every reason to use it. Two features matter for what this probe asks:

*   **A long rainless stretch in each of the last two years, labelled
    `_regime = dry`.** The criterion scores the whole record and gates no
    reversal statistic; these windows are where its reversal *diagnostics* are
    computed, so a reviewer can see how an exchange behaves with direct
    precipitation forcing absent. Reversal is never judged, here or anywhere:
    a gaining/losing reach changes state seasonally, bank storage reverses on
    the limbs of one flood, and a MODFLOW aquifer behind an oscillating
    boundary reverses with no rain at all.
*   **A prescribed external head, `gwh`, visible to the model.** The regime an
    exchange is judged against is specified rather than assumed: a slow annual
    cycle on which a faster 20-day oscillation rides, so the head crosses its
    own mean several times inside each rainless window. A model that declares
    it consumes this column is held to it; one that does not is not judged on
    reversal at all, because local rainlessness constrains nothing about a
    boundary whose driver the case did not supply.
*   **Strong day-to-day variability in the rain.** A residual sink's errors
    change sign with the weather, so the weather has to be able to. Gamma
    depths on a seasonal occurrence probability give the same spiky series the
    closure probe uses, and the drying tail of each rainless window still
    carries the baseflow an honest exchange draws on.

Climate targets are the closure probe's — about 843 mm of precipitation and
759 mm of potential ET a year, with a real seasonal snowpack — so a model that
passes that probe is not being asked to work in a new climate as well as to
account for its exchange.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 10
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

# The control head, and the same head shifted up and down by a constant. The
# shift is of the order of the head's own swing, so a genuinely head-driven
# exchange answers it unmistakably; everything else is identical. The shift
# starts with the scored record, not with the spinup: an aquifer that adjusts
# quickly would otherwise have finished answering before scoring began, and
# the paired difference would read as no response at all.
VARIANTS = ("control", "raised", "lowered")
HEAD_SHIFT_M = {"control": 0.0, "raised": +0.5, "lowered": -0.5}

# Two rainless stretches, placed in the last two scored years so that every
# store has settled before they start, and long enough that a draining
# catchment is unambiguously the only thing happening. Day indices from the
# record start; the criterion reads them back from the `_regime` column.
DROUGHTS = ((SPINUP_DAYS + 8 * 365 + 150, 90), (SPINUP_DAYS + 9 * 365 + 160, 75))

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

    # Precipitation: seasonal occurrence, gamma depths on wet days.
    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(N_STEPS) < p_wet
    depth = rng.gamma(shape=0.7, scale=13.0, size=N_STEPS)
    pr = np.where(wet, depth, 0.0)
    for start, length in DROUGHTS:
        pr[start:start + length] = 0.0

    # Temperature: annual cycle plus a persistent (AR1) weather anomaly.
    seasonal = 9.0 - 12.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(N_STEPS)
    innovation = rng.normal(0.0, 2.6, N_STEPS)
    for t in range(1, N_STEPS):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    tas = seasonal + noise

    # Potential evapotranspiration: temperature driven with a daylength factor.
    daylength = 1.0 + 0.35 * np.cos(2 * np.pi * (doy - 172) / 365)
    pet = np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength

    # The external head the exchange is judged against. Metres; a fixed datum
    # of 10 m, a 0.15 m annual cycle and a 0.35 m cycle of 20 days, so that the
    # head crosses its local mean about every ten days inside the windows and a
    # model following it reverses there. Prescribed, not hidden: this column
    # has no underscore and every model receives it.
    gwh = (10.0 + 0.15 * np.sin(2 * np.pi * (doy - 100) / 365)
           + 0.35 * np.sin(2 * np.pi * day / 20.0)
           + np.where(day >= SPINUP_DAYS, HEAD_SHIFT_M[variant], 0.0))

    # The criterion is scored on the rainless windows, so they are labelled
    # rather than rediscovered from the rain. Columns beginning with an
    # underscore are the probe's own annotation: the harness strips them before
    # the model sees the forcing, so no model is told which stretch is scored.
    regime = np.full(N_STEPS, "wet", dtype=object)
    for start, length in DROUGHTS:
        regime[start:start + length] = "dry"

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
            "gwh": np.round(gwh, 4),
            "_regime": regime,
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    frame, static = generate(20260912)
    up, _ = generate(20260912, "raised")
    same = all((frame[c] == up[c]).all() for c in ("time", "pr", "tas", "pet", "_regime"))
    spin_same = (frame["gwh"][:SPINUP_DAYS] == up["gwh"][:SPINUP_DAYS]).all()
    print(f"variants differ only in gwh: {same}; spinup gwh identical: {spin_same}; "
          f"scored shift {float((up['gwh'] - frame['gwh'])[SPINUP_DAYS:].mean()):+.2f} m")
    annual_p = frame["pr"].sum() / (N_STEPS / 365)
    annual_pet = frame["pet"].sum() / (N_STEPS / 365)
    dry = int((frame["pr"] == 0).sum())
    scored = int((frame["_regime"] == "dry").sum())
    print(f"{len(frame)} steps, {annual_p:.0f} mm/yr precipitation, "
          f"{annual_pet:.0f} mm/yr potential ET, {dry} rainless days "
          f"({scored} of them in the two scored windows), "
          f"{int((frame['tas'] < 0).sum())} sub-zero days")

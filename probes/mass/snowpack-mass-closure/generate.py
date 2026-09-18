"""Forcing generator for mass/snowpack-mass-closure.

Deterministic given a seed, and never committed as data. Generating the case
at run time is what stops a model from memorising it, and shipping the
generator instead of a file is what lets a reviewer see exactly what the
model will be given.

The case contains two accumulation-storage-melt snow cycles. Temperatures
stay well below the snow threshold during accumulation and storage and well
above it during melt, so the probe tests snowpack bookkeeping rather than a
particular rain-snow partition threshold.

Requirements the harness enforces:
  - generate(seed) returns (DataFrame, dict)
  - the frame has a 'time' column
  - it has exactly PERIOD_DAYS + SPINUP_DAYS rows
  - the same seed produces byte-identical output

Run this file directly to sanity check what it produces.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

STAGE_DAYS = 60
N_CYCLES = 2
PERIOD_DAYS = 3 * STAGE_DAYS * N_CYCLES
SPINUP_DAYS = 0
N_STEPS = PERIOD_DAYS + SPINUP_DAYS

# Catchment attributes handed to every model, and the place criteria look up
# named bounds such as soil_capacity_mm.
STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 44.0,
}


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    day = np.arange(N_STEPS)
    doy = day % 365

    # Two identical stage sequences are labelled separately in time. The
    # leading underscore marks `_regime` as a criterion annotation rather
    # than a forcing variable supplied to the model.
    regime = np.tile(
        np.repeat(["accumulation", "storage", "melt"], STAGE_DAYS),
        N_CYCLES,
    )
    accumulation = regime == "accumulation"
    melt = regime == "melt"

    # Precipitation occurs only during accumulation. Event occurrence and
    # depths vary with the seed, while storage and melt stages remain dry so
    # changes in snow storage can be attributed unambiguously.
    wet = (rng.random(N_STEPS) < 0.30) & accumulation
    depth = rng.gamma(shape=0.7, scale=13.0, size=N_STEPS)
    pr = np.where(wet, depth, 0.0)

    # With no spinup, the first model output is also the harness reference
    # state. Keep the first day dry so an empty snowpack is unchanged before
    # the first scored snowfall.
    pr[0] = 0.0

    # Temperature uses the same simple persistent weather variability as the
    # reference generators, but is kept away from the phase threshold:
    # cold stages are at most -3 degC and melt stages at least +4 degC.
    noise = np.zeros(N_STEPS)
    innovation = rng.normal(0.0, 0.8, N_STEPS)
    for t in range(1, N_STEPS):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]

    tas = np.where(melt, 7.0, -6.0) + noise
    tas = np.where(melt, np.maximum(tas, 4.0), np.minimum(tas, -3.0))

    # Potential evapotranspiration follows the temperature-driven form used
    # by the existing synthetic mass probes.
    daylength = 1.0 + 0.35 * np.cos(2 * np.pi * (doy - 172) / 365)
    pet = np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
            "_regime": regime,
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    frame, _ = generate(20260903)
    cold = frame["_regime"].isin(["accumulation", "storage"])
    melt = frame["_regime"] == "melt"
    print(
        f"{len(frame)} steps, {frame['pr'].sum():.0f} mm precipitation, "
        f"cold tas <= {frame.loc[cold, 'tas'].max():.1f} degC, "
        f"melt tas >= {frame.loc[melt, 'tas'].min():.1f} degC"
    )

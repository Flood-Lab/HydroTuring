"""Seeded synthetic forcing for the snowmelt energy-water coupling probe.

Deterministic given a seed, and never committed as data. One snow year in
three stages, of which only the last is scored.

    spinup        30 d  cold and dry, the pack starts empty
    accumulation 150 d  snowfall well below any rain-snow threshold
    melt          60 d  no precipitation, cold at first and then well above
                        freezing                                   <- scored

The scored block carries no precipitation. That is the one property of the
case the identity still rests on: the contract has neither a melt flux nor a
snowfall flux, so ice arriving during the block could not be told from ice
melting, and the criterion refuses a block that is not dry.

The block deliberately **opens on a cold pack** and warms through zero inside
the scored stretch. An earlier design ripened the pack in an unscored stage so
that the surface residual would be fusion alone, which made ripeness a
precondition no criterion could observe and no case could guarantee for a
model whose snow physics differs from the reference's. The criterion now
carries the cold content as a term of the balance instead, so the case no
longer has to engineer it away -- and warming the pack inside the scored block
is what gives that term something to measure.

The stage labels travel in a `_regime` column. The harness strips every
underscore-prefixed column before staging, so the model is given the weather
and not the answer.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_DAYS = 210
SPINUP_DAYS = 30
N_STEPS = PERIOD_DAYS + SPINUP_DAYS

# (label, days). The order is the physics: a pack has to be built before it can
# go, and the scored block holds both halves of what happens to it -- the cold
# opening where it warms, and the end where it melts.
STAGES = (
    ("spinup", 30),
    ("accumulation", 150),
    ("melt", 60),
)

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 47.0,
}


def _stage_index() -> np.ndarray:
    labels = []
    for name, days in STAGES:
        labels.extend([name] * days)
    return np.array(labels, dtype=object)


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    regime = _stage_index()

    # Every draw happens before any stage branch, so that two seeds differ in
    # the weather and never in the shape of the experiment.
    cold_base = rng.uniform(-20.0, -16.0)      # degC, the accumulation stage
    warm_base = rng.uniform(8.0, 11.0)         # degC, the end of the melt stage
    snowfall_mean = rng.uniform(5.0, 6.5)      # mm/day over the accumulation
    rn_winter = rng.uniform(8.0, 18.0)         # W/m2, a low sun and a bright pack
    rn_melt_peak = rng.uniform(70.0, 100.0)    # W/m2, by the end of the melt
    anomaly = np.zeros(N_STEPS)
    innovation = rng.normal(0.0, 1.6, N_STEPS)
    for t in range(1, N_STEPS):
        anomaly[t] = 0.70 * anomaly[t - 1] + innovation[t]
    wet = rng.random(N_STEPS)
    depth = rng.gamma(shape=2.2, scale=1.0, size=N_STEPS)

    is_spinup = regime == "spinup"
    is_accum = regime == "accumulation"
    is_melt = regime == "melt"

    # --- temperature ------------------------------------------------------
    # The accumulation stage sits clearly below any plausible rain-snow
    # threshold, so the phase split never depends on a model's own PXTEMP. The
    # scored block starts there too and climbs through zero, so a pack that
    # carries cold content has to spend it where the criterion is watching.
    tas = np.empty(N_STEPS)
    tas[is_spinup] = cold_base - 2.0
    tas[is_accum] = cold_base
    tas[is_melt] = np.linspace(cold_base, warm_base, int(is_melt.sum()))
    tas = tas + anomaly
    tas[is_accum | is_spinup] = np.clip(tas[is_accum | is_spinup], -24.0, -12.0)
    # The block has to open on a pack that is actually cold, or the seed gives
    # the cold-content term nothing to measure. The AR(1) anomaly alone can
    # lift the first scored day above zero, so the opening stretch is held
    # below freezing on every seed.
    scored = np.flatnonzero(is_melt)
    tas[scored[:20]] = np.minimum(tas[scored[:20]], -14.0)
    # And it has to end warm enough to take the pack off, or the fusion term
    # has nothing to measure either. Both ends are held rather than trusted to
    # the draw, so no seed produces a block that tests only half the balance.
    tas[scored[-20:]] = np.maximum(tas[scored[-20:]], 4.0)

    # --- precipitation ----------------------------------------------------
    # Snowfall on about four days in five of the accumulation stage, and not a
    # millimetre anywhere else. The scored block is dry by construction.
    pr = np.zeros(N_STEPS)
    snowing = is_accum & (wet < 0.8)
    pr[snowing] = depth[snowing] * snowfall_mean / (0.8 * 2.2)

    # --- net radiation ----------------------------------------------------
    # Low over a bright winter pack, climbing through the melt season as the
    # sun rises and the albedo falls.
    rn = np.empty(N_STEPS)
    rn[is_spinup | is_accum] = rn_winter
    rn[is_melt] = np.linspace(rn_winter, rn_melt_peak, int(is_melt.sum()))
    rn = rn + 0.8 * anomaly

    # --- potential evapotranspiration -------------------------------------
    pet = np.maximum(0.0, 0.13 * (tas + 5.0))

    time = pd.date_range("2000-10-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
            "rn": np.round(rn, 6),
            # Stripped before the model sees it. The criteria read it.
            "_regime": regime,
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    import sys

    LAMBDA_F = 3.337e5
    for seed in (20260913, 11111111, 99999999):
        frame, static = generate(seed)
        accum = frame[frame["_regime"] == "accumulation"]
        melt = frame[frame["_regime"] == "melt"]
        snowfall = accum["pr"].sum()
        print(
            f"seed {seed}: {len(frame)} steps | snowfall {snowfall:.0f} mm | "
            f"scored block {len(melt)} d, T {melt['tas'].min():.1f}..{melt['tas'].max():.1f} C, "
            f"rn {melt['rn'].min():.0f}..{melt['rn'].max():.0f} W/m2 | melting all of that "
            f"pack demands {LAMBDA_F * snowfall / (len(melt) * 86400.0):.1f} W/m2"
        )
        if melt["pr"].sum() > 0:
            sys.exit("the scored block is not dry; ice arriving cannot be told from ice melting")
        if melt["tas"].iloc[0] >= 0:
            sys.exit("the scored block does not open on a cold pack; nothing exercises cold content")
        if melt["tas"].iloc[-1] <= 3:
            sys.exit("the scored block does not warm enough to melt anything")

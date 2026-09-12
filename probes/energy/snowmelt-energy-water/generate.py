"""Seeded synthetic forcing for the snowmelt energy-water coupling probe.

Deterministic given a seed, and never committed as data. One snow year in four
stages, of which only the last is scored.

    spinup        30 d  cold and dry, the pack starts empty
    accumulation 110 d  snowfall well below any rain-snow threshold
    ripening      40 d  no precipitation, temperature rising through zero
    melt          60 d  no precipitation, well above freezing        <- scored

Two properties of the scored stage are what make the probe possible, and both
are built here rather than assumed of the model.

No precipitation falls during it, so snowfall is zero whatever threshold a
model uses internally, and the pack's loss is melt plus sublimation and
nothing else. The contract carries neither a melt flux nor a snowfall flux, so
this is what makes melt observable at all.

The pack is ripe before it opens. Energy spent warming a sub-freezing pack
towards zero does no melting and appears in no contract variable, so the
ripening stage spends that cold content where it is not scored. What is left
over the melt stage is fusion and nothing else.

The stage labels travel in a `_regime` column. The harness strips every
underscore-prefixed column before staging, so the model is given the weather
and not the answer: it has to tell the stages apart from the forcing itself.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_DAYS = 210
SPINUP_DAYS = 30
N_STEPS = PERIOD_DAYS + SPINUP_DAYS

# (label, days). The order is the physics: a pack has to be built before it can
# be ripened, and ripened before what it does next is only melt.
STAGES = (
    ("spinup", 30),
    ("accumulation", 110),
    ("ripening", 40),
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
    cold_base = rng.uniform(-7.0, -5.0)        # degC, the accumulation stage
    warm_base = rng.uniform(4.5, 7.0)          # degC, the melt stage
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
    is_ripen = regime == "ripening"
    is_melt = regime == "melt"

    # --- temperature ------------------------------------------------------
    # The cold stages sit clearly below, and the melt stage clearly above, any
    # plausible rain-snow threshold, so the phase split never depends on a
    # model's own PXTEMP. The anomaly is damped in the melt stage so that no
    # seed drops a scored day back below freezing and reopens the cold content
    # this design exists to close.
    tas = np.empty(N_STEPS)
    tas[is_spinup] = cold_base - 2.0
    tas[is_accum] = cold_base
    tas[is_ripen] = np.linspace(cold_base, 1.5, int(is_ripen.sum()))
    tas[is_melt] = warm_base + np.linspace(0.0, 2.5, int(is_melt.sum()))
    tas = tas + np.where(is_melt, 0.35, 1.0) * anomaly
    tas[is_melt] = np.maximum(tas[is_melt], 3.0)
    # Bounded on both sides. Below -2 C so no seed puts a scored-stage
    # threshold in doubt, and above -10 C because the cold content of the pack
    # scales with how cold it is and the ripening stage has to be able to
    # spend it. See README, "The pack is ripe before the block opens".
    tas[is_accum | is_spinup] = np.clip(tas[is_accum | is_spinup], -10.0, -2.0)

    # --- precipitation ----------------------------------------------------
    # Snowfall on about four days in five of the accumulation stage, and not a
    # millimetre anywhere else. The scored stage is dry by construction.
    pr = np.zeros(N_STEPS)
    snowing = is_accum & (wet < 0.8)
    pr[snowing] = depth[snowing] * snowfall_mean / (0.8 * 2.2)

    # --- net radiation ----------------------------------------------------
    # Low over a bright winter pack, climbing through the melt season as the
    # sun rises and the albedo falls.
    rn = np.empty(N_STEPS)
    rn[is_spinup | is_accum] = rn_winter
    rn[is_ripen] = np.linspace(rn_winter, 55.0, int(is_ripen.sum()))
    rn[is_melt] = np.linspace(55.0, rn_melt_peak, int(is_melt.sum()))
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
    C_ICE = 2100.0              # J kg-1 K-1, for the cold content of the pack
    LAMBDA_S = 2.501e6 + LAMBDA_F   # sublimation takes vaporisation plus fusion
    K_H = 3.0                   # W m-2 K-1, the bulk exchange over a pack
    SUBL_SHARE = 0.35           # share of demand a pack meets by sublimating
    GROUND_SHARE = 0.10         # share of net radiation conducted into the ground
    MIN_MARGIN = 2.0

    def ripening_margin(frame):
        """How much more energy the ripening stage delivers than the pack needs.

        A lower bound on both sides, deliberately. The cold content is taken at
        the *coldest* accumulation day rather than the mean, because a pack
        integrates its season and the coldest day is the conservative end; and
        the energy delivered has the pack's own sublimation subtracted, which
        leaves the surface as latent heat and so cannot also ripen the pack.
        Computing either the generous way inflates the margin, and the margin
        is the whole guarantee that a model carrying real cold content is not
        failed for arriving at the scored block still cold.
        """
        accum = frame[frame["_regime"] == "accumulation"]
        ripen = frame[frame["_regime"] == "ripening"]
        cold_content = C_ICE * accum["pr"].sum() * abs(accum["tas"].min())

        rn = ripen["rn"].to_numpy()
        tas = ripen["tas"].to_numpy()
        sublimation = SUBL_SHARE * np.maximum(0.0, 0.13 * (tas + 5.0))
        latent = LAMBDA_S * sublimation / 86400.0
        sensible = K_H * (0.0 - tas)
        delivered = float((((1.0 - GROUND_SHARE) * rn - sensible - latent).sum()) * 86400.0)
        return delivered / cold_content

    worst = (float("inf"), None)
    for seed in range(0, 8000, 23):
        frame, static = generate(seed)
        if len(frame) != N_STEPS:
            sys.exit(f"seed {seed}: {len(frame)} rows, expected {N_STEPS}")
        melt = frame[frame["_regime"] == "melt"]
        if melt["pr"].sum() > 0:
            sys.exit(f"seed {seed}: the scored stage is not dry; melt is not observable")
        if melt["tas"].min() <= 0:
            sys.exit(f"seed {seed}: a scored day is below freezing; cold content reopens")
        margin = ripening_margin(frame)
        if margin < MIN_MARGIN:
            sys.exit(
                f"seed {seed}: the ripening stage delivers only {margin:.2f}x the "
                "pack's cold content; a model that carries one would meet the "
                "scored block still cold and fail for being right"
            )
        worst = min(worst, (margin, seed))

    for seed in (20260913, 11111111, 99999999):
        frame, static = generate(seed)
        accum = frame[frame["_regime"] == "accumulation"]
        melt = frame[frame["_regime"] == "melt"]
        snowfall = accum["pr"].sum()
        demand = LAMBDA_F * snowfall / (len(melt) * 86400.0)
        print(
            f"seed {seed}: {len(frame)} steps | snowfall {snowfall:.0f} mm | "
            f"melt stage {len(melt)} d, T {melt['tas'].min():.1f}..{melt['tas'].max():.1f} C, "
            f"rn {melt['rn'].mean():.0f} W/m2 | a pack that size demands "
            f"{demand:.1f} W/m2 | ripening delivers {ripening_margin(frame):.2f}x "
            "its cold content"
        )
    print(f"348 seeds checked; worst ripening margin {worst[0]:.2f}x at seed {worst[1]}")

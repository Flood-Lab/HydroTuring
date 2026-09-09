"""Seeded synthetic forcing for the latent-heat / ET coherence probe.

Deterministic given a seed, and never committed as data. Generating the case
at run time is what stops a model from memorising it, and shipping the
generator instead of a NetCDF file is what lets a reviewer see exactly what
the model will be given.

The water side is the record of `mass/catchment-closure`: the same occurrence
and depth process for precipitation, the same annual cycle and AR(1) anomaly
for temperature, so that a model already adapted to that probe sees nothing
unfamiliar. Two things are added.

Net radiation is generated from a clear-sky cycle at the catchment's latitude,
attenuated by a seeded cloudiness process, with an albedo that switches when
the air is below freezing and an outgoing longwave term that follows air
temperature. It averages under 20 W/m2 in the depth of winter and goes below
zero on a handful of overcast days, which is deliberate: it is what forces the
energy criteria to carry an absolute floor rather than a bare percentage.

The shortwave cycle peaks at the summer solstice. The first release had the
sign of its seasonal term flipped, so the year's strongest radiation and
demand fell in the wet winter months and the soil never filled; the exact
model's runoff was then almost pure baseflow and its correlation with recent
rain was a coin toss against the `non_degenerate` threshold. See the README.

Potential evapotranspiration is then derived from that net radiation by the
Priestley-Taylor equation with alpha = 1.26, rather than being an independent
temperature function. A reviewer will ask why the demand and the available
energy in a coherence probe are consistent with each other, and this is the
answer: they are not two forcings, they are one.

The consequence that matters for the probe is that PET stays small but
non-zero on cold, clear days. Without it no snowpack would ever sublimate,
and the phase-change half of the criterion would have nothing to score.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = 10
SPINUP_DAYS = 365
N_STEPS = PERIOD_YEARS * 365 + SPINUP_DAYS

# Priestley-Taylor, and the psychrometric constant at sea level.
PT_ALPHA = 1.26
GAMMA_KPA_PER_C = 0.0665
GROUND_SHARE = 0.10

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 40.0,
}


def _saturation_slope(tas: np.ndarray) -> np.ndarray:
    """d(e_s)/dT in kPa per degC, the Tetens form used everywhere."""
    es = 0.6108 * np.exp(17.27 * tas / (tas + 237.3))
    return 4098.0 * es / np.square(tas + 237.3)


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    day = np.arange(N_STEPS)
    doy = day % 365

    # Precipitation: seasonal occurrence, gamma depths on wet days.
    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(N_STEPS) < p_wet
    depth = rng.gamma(shape=0.7, scale=13.0, size=N_STEPS)
    pr = np.where(wet, depth, 0.0)

    # Temperature: annual cycle plus a persistent (AR1) weather anomaly.
    seasonal = 9.0 - 12.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(N_STEPS)
    innovation = rng.normal(0.0, 2.6, N_STEPS)
    for t in range(1, N_STEPS):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    tas = seasonal + noise

    # Cloudiness: its own AR(1), so that a bright winter day and a dull summer
    # one both occur and net radiation is not a pure function of the calendar.
    cloud = np.zeros(N_STEPS)
    shock = rng.normal(0.0, 0.45, N_STEPS)
    for t in range(1, N_STEPS):
        cloud[t] = 0.55 * cloud[t - 1] + shock[t]
    clearness = np.clip(0.72 + 0.18 * cloud, 0.35, 1.0)

    # Daily-mean incoming shortwave at 40 degN, and the surface it lands on.
    rsds = (245.0 + 135.0 * np.cos(2 * np.pi * (doy - 172) / 365)) * clearness
    albedo = np.where(tas < STATIC["snow_threshold_degC"], 0.60, 0.23)
    longwave_out = np.clip(32.0 + 1.2 * tas, 15.0, 95.0) * (0.55 + 0.45 * clearness)
    rn = (1.0 - albedo) * rsds - longwave_out

    # Priestley-Taylor demand on the available energy.
    slope = _saturation_slope(tas)
    lam = 2.501e6 - 2361.0 * tas
    available = np.maximum(rn * (1.0 - GROUND_SHARE), 0.0)
    pet = PT_ALPHA * (slope / (slope + GAMMA_KPA_PER_C)) * available * 86400.0 / lam

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
            "rn": np.round(rn, 6),
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    frame, static = generate(20260903)
    years = N_STEPS / 365
    dry_frozen = ((frame["pr"] == 0) & (frame["tas"] < 0)).sum()
    print(
        f"{len(frame)} steps, {frame['pr'].sum() / years:.0f} mm/yr precipitation, "
        f"{frame['pet'].sum() / years:.0f} mm/yr potential ET, "
        f"rn mean {frame['rn'].mean():.1f} W/m2 "
        f"(min {frame['rn'].min():.1f}, max {frame['rn'].max():.1f}), "
        f"{int((frame['rn'] < 0).sum())} days of negative rn, "
        f"{int(dry_frozen)} dry sub-zero days"
    )

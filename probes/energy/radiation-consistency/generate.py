"""Seeded hourly forcing for a warm, snow-free gray surface.

Timestamps use local solar time. tas and rlds are instantaneous, as are the
requested ts and rlus. Water/energy inputs pr, pet and rn are means over
[time, time + 1 h). Net radiation is prescribed independently of longwave."""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_DAYS = 30
SPINUP_DAYS = 2
STEPS_PER_DAY = 24
N_STEPS = (PERIOD_DAYS + SPINUP_DAYS) * STEPS_PER_DAY

# Stefan-Boltzmann constant, W m-2 K-4 (CODATA 2018).
STEFAN_BOLTZMANN = 5.670374419e-8
# Shared bounds for shortwave transmission and derived cloud fraction.
TRANSMISSION = (0.45, 1.0)

STATIC = {
    "area_km2": 250.0,
    "soil_capacity_mm": 180.0,
    "canopy_capacity_mm": 0.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 0.0,
}


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    time = pd.date_range("2000-03-20T06:00:00", periods=N_STEPS, freq="h")
    hour = time.hour.to_numpy()
    daytime = (hour >= 6) & (hour < 18)
    solar_day = np.arange(N_STEPS) // STEPS_PER_DAY
    n_days = PERIOD_DAYS + SPINUP_DAYS

    # Exact hourly mean of the positive sine between sunrise at 06:00 and
    # sunset at 18:00, as in energy/surface-energy-closure. Cloud attenuation
    # is drawn at every hour, night included, because it also sets how much
    # the sky radiates downward.
    omega = np.pi / 12.0
    daylight = np.where(
        daytime,
        (np.cos(omega * (hour - 6)) - np.cos(omega * (hour - 5))) / omega,
        0.0,
    )
    peak_shortwave = rng.uniform(420.0, 540.0, n_days)[solar_day]
    transmission = rng.uniform(*TRANSMISSION, N_STEPS)
    longwave_cooling = rng.uniform(25.0, 40.0, n_days)[solar_day]
    rn = peak_shortwave * transmission * daylight - longwave_cooling

    # Sample warm air at the radiation timestamp, not as an hourly mean.
    weather = rng.uniform(-2.0, 2.0, n_days)[solar_day]
    tas = 23.0 + weather + 5.0 * np.cos(2.0 * np.pi * (hour - 14.0) / 24.0)

    # The same cloud draw dims the sun and raises sky emissivity. Overcast
    # conditions close three quarters of the gap from clear sky to a black body.
    cloud = (TRANSMISSION[1] - transmission) / (TRANSMISSION[1] - TRANSMISSION[0])
    eps_clear = rng.uniform(0.74, 0.82, n_days)[solar_day]
    eps_sky = eps_clear + (1.0 - eps_clear) * 0.75 * cloud
    rlds = eps_sky * STEFAN_BOLTZMANN * (tas + 273.15) ** 4

    # Water-side inputs for the reference model, in mm/day at every hour.
    pet = 0.12 + 3.2 * daylight * transmission
    wet = rng.random(N_STEPS) < 0.025
    pr = np.where(wet, rng.uniform(12.0, 48.0, N_STEPS), 0.0)

    # Supply one fixed surface emissivity for the entire case.
    static = dict(STATIC)
    static["eps"] = round(float(rng.uniform(0.95, 0.99)), 4)

    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
            "rn": np.round(rn, 6),
            "rlds": np.round(rlds, 6),
        }
    )
    return forcing, static


if __name__ == "__main__":
    frame, static = generate(20260911)
    print(f"{len(frame)} hourly rows: {SPINUP_DAYS} spinup + {PERIOD_DAYS} scored days")
    print(f"Air temperature: {frame['tas'].min():.1f} to {frame['tas'].max():.1f} degC")
    print(f"Downward longwave: {frame['rlds'].min():.0f} to {frame['rlds'].max():.0f} W/m2")
    print(f"Net radiation: {frame['rn'].min():.0f} to {frame['rn'].max():.0f} W/m2")
    print(f"Surface emissivity: {static['eps']}")

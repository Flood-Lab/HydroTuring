"""A warm, dry soil layer under shared radiative and meteorological forcing."""

from __future__ import annotations

import numpy as np
import pandas as pd

SPINUP_HOURS = 24
HEATING_HOURS = 12
RECOVERY_HOURS = 60
N_STEPS = SPINUP_HOURS + HEATING_HOURS + RECOVERY_HOURS


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    depth = float(rng.uniform(0.1, 0.35))
    solid_capacity = float(rng.uniform(1.8e6, 2.2e6))
    heating = float(rng.uniform(250.0, 400.0))
    temperature = float(rng.uniform(288.15, 298.15))
    porosity, water_content = 0.464, 0.005
    # Prescribed mixture properties, before any model runs. The pore-air term
    # uses the effective coefficient in the documented dry-soil test material.
    water_capacity, pore_air_capacity = 4.188e6, 1004.64
    volumetric_capacity = (
        (1.0 - porosity) * solid_capacity + water_content * water_capacity
        + (porosity - water_content) * pore_air_capacity
    )

    time = pd.date_range("2000-07-01T06:00:00", periods=N_STEPS, freq="h")
    shortwave = np.zeros(N_STEPS)
    shortwave[SPINUP_HOURS:SPINUP_HOURS + HEATING_HOURS] = heating
    forcing = pd.DataFrame({
        "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "pr": 0.0,
        "tas": temperature - 273.15,
        "pet": 0.0,
        # The prescribed forcing is constant within each one-hour interval.
        "rsds": shortwave,
        "rlds": 5.670374419e-8 * temperature**4,
        "sfcWind": 2.0,
        "huss": 0.0,
        "ps": 101325.0,
        "_phase": (["spinup"] * SPINUP_HOURS + ["heating"] * HEATING_HOURS
                   + ["recovery"] * RECOVERY_HOURS),
    })
    static = {
        "area_km2": 1.0,
        "soil_layer_depth_m": depth,
        "soil_heat_capacity_areal": volumetric_capacity * depth,
        "soil_temperature_initial": temperature,
        "soil_solid_heat_capacity": solid_capacity,
        "soil_porosity": porosity,
        "soil_water_content_initial": water_content,
        "soil_water_heat_capacity": water_capacity,
        "soil_pore_air_heat_capacity": pore_air_capacity,
    }
    return forcing, static

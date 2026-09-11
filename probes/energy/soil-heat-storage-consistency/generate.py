"""Hourly heating and recovery of a homogeneous layer, without phase change."""

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
    volumetric_capacity = float(rng.uniform(1.6e6, 2.4e6))
    sensible_exchange = float(rng.uniform(10.0, 15.0))
    top_conductance = float(rng.uniform(4.0, 8.0))
    bottom_conductance = float(rng.uniform(0.5, 1.5))
    heating = float(rng.uniform(120.0, 240.0))
    temperature = float(rng.uniform(288.15, 298.15))

    time = pd.date_range("2000-07-01T06:00:00", periods=N_STEPS, freq="h")
    rn = np.zeros(N_STEPS)
    rn[SPINUP_HOURS:SPINUP_HOURS + HEATING_HOURS] = heating
    forcing = pd.DataFrame({
        "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "pr": 0.0,
        "tas": temperature - 273.15,
        "pet": 0.0,
        # The prescribed forcing is constant within each one-hour interval.
        "rn": rn,
        "_phase": (["spinup"] * SPINUP_HOURS + ["heating"] * HEATING_HOURS
                   + ["recovery"] * RECOVERY_HOURS),
    })
    static = {
        "area_km2": 1.0,
        "soil_layer_depth_m": depth,
        "soil_heat_capacity_areal": volumetric_capacity * depth,
        "soil_temperature_initial": temperature,
        "soil_deep_temperature": temperature,
        "soil_sensible_exchange": sensible_exchange,
        "soil_top_conductance": top_conductance,
        "soil_bottom_conductance": bottom_conductance,
    }
    return forcing, static

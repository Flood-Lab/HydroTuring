"""Forcing generator for {id}.

Deterministic given a seed, and never committed as data. Generating the case
at run time is what stops a model from memorising it, and shipping the
generator instead of a file is what lets a reviewer see exactly what the
model will be given.

Requirements the harness enforces:
  - generate(seed) returns (DataFrame, dict)
  - the frame has a 'time' column
  - it has exactly PERIOD_YEARS * 365 + SPINUP_DAYS rows
  - the same seed produces byte-identical output

Run this file directly to sanity check what it produces.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = {period_years}
SPINUP_DAYS = {spinup_days}
N_STEPS = int(PERIOD_YEARS * 365) + SPINUP_DAYS

# Catchment attributes handed to every model, and the place criteria look up
# named bounds such as soil_capacity_mm.
STATIC = {{
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    # TODO: add whatever your probe's criteria and models need.
}}


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    day = np.arange(N_STEPS)
    doy = day % 365

    # TODO: build your forcing. Everything below is a placeholder.
    # Keep it physically plausible: a model should not be able to tell it is
    # being tested from the statistics of what it is handed.
    pr = np.zeros(N_STEPS)
    tas = np.zeros(N_STEPS)
    pet = np.zeros(N_STEPS)

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {{
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
        }}
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    frame, static = generate(20260903)
    years = N_STEPS / 365
    print(f"{{len(frame)}} steps")
    print(f"{{frame['pr'].sum() / years:.0f}} mm/yr precipitation")
    print(f"{{frame['pet'].sum() / years:.0f}} mm/yr potential ET")

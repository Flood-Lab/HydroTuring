"""Forcing generator for {id}: ordinary years, then a year outside them.

Deterministic given a seed. The record is one continuous series with a
`_regime` column naming each step. The leading underscore is not decoration:
the harness strips underscore-prefixed columns before staging, so the model
receives the weather and not the label. It has to notice the anomaly from the
forcing, which is the thing being tested.

Requirements the harness enforces:
  - generate(seed) returns (DataFrame, dict)
  - the frame has a 'time' column
  - it has exactly PERIOD_YEARS * 365 + SPINUP_DAYS rows
  - the same seed produces byte-identical output

Run this file directly to see what the two regimes look like.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = {period_years}
SPINUP_DAYS = {spinup_days}
N_STEPS = int(PERIOD_YEARS * 365) + SPINUP_DAYS

# Where the ordinary record ends and the anomaly begins, counted from the
# start of the scored window. Everything before it, spinup included, is
# labelled `ordinary`, so the model meets the anomaly having seen nothing
# like it.
ANOMALY_START = SPINUP_DAYS + int((PERIOD_YEARS - 1) * 365)

STATIC = {{
    "area_km2": 250.0,
    "soil_capacity_mm": 320.0,
    "canopy_capacity_mm": 2.0,
    "degree_day_factor_mm_per_C_day": 3.2,
    "baseflow_coefficient": 0.006,
    "snow_threshold_degC": 0.0,
    "latitude_deg": 40.0,
}}


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    day = np.arange(N_STEPS)
    doy = day % 365

    # --- the ordinary years -------------------------------------------------
    # TODO: your own climate. Keep it plausible; a model should not be able to
    # tell it is being tested from the statistics of the ordinary stretch.
    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(N_STEPS) < p_wet
    pr = np.where(wet, rng.gamma(shape=0.7, scale=13.0, size=N_STEPS), 0.0)

    seasonal = 9.0 - 12.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(N_STEPS)
    innovation = rng.normal(0.0, 2.6, N_STEPS)
    for t in range(1, N_STEPS):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    tas = seasonal + noise

    daylength = 1.0 + 0.35 * np.cos(2 * np.pi * (doy - 172) / 365)
    pet = np.maximum(0.0, 0.13 * (tas + 5.0)) * daylength

    # --- the anomaly --------------------------------------------------------
    # TODO: this is the probe. Make the last stretch leave the range of
    # everything before it, in a way that is physically coherent rather than
    # merely large. A record storm needs the synoptic conditions that carry
    # it; a multi-year drought needs the temperature and the demand that go
    # with it. An anomaly that is only a scaled version of an ordinary year
    # tests scaling, not extrapolation.
    #
    # The placeholder below is a wet extreme: state what yours is in README.md
    # and say, in one sentence, why nothing in the ordinary stretch prepares a
    # model for it.
    anomaly = slice(ANOMALY_START, N_STEPS)
    ordinary_max = float(pr[:ANOMALY_START].max())
    pr[anomaly] = pr[anomaly] * 2.5
    storm = ANOMALY_START + int(rng.integers(60, 300))
    pr[storm] = 3.0 * ordinary_max

    regime = np.where(day < ANOMALY_START, "ordinary", "anomaly")

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {{
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
            # Stripped before the model sees it. The criteria read it.
            "_regime": regime,
        }}
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    frame, static = generate(20260903)
    ordinary = frame[frame["_regime"] == "ordinary"]
    anomaly = frame[frame["_regime"] == "anomaly"]
    print(f"{{len(frame)}} steps: {{len(ordinary)}} ordinary, {{len(anomaly)}} anomaly")
    print(f"ordinary: {{ordinary['pr'].sum() / (len(ordinary) / 365):.0f}} mm/yr, "
          f"wettest day {{ordinary['pr'].max():.1f}} mm")
    print(f"anomaly:  {{anomaly['pr'].sum() / (len(anomaly) / 365):.0f}} mm/yr, "
          f"wettest day {{anomaly['pr'].max():.1f}} mm")
    if anomaly["pr"].max() <= ordinary["pr"].max():
        print("WARNING: the anomaly does not leave the range of the ordinary "
              "record. This probe is not testing extrapolation yet.")

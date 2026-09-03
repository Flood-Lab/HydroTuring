"""Forcing generator for {id}: a different catchment on every seed.

Deterministic given a seed. The spatial dimension of this probe is the seed:
each one draws a catchment from attribute space and the weather that belongs
with it. Because the harness requires every seed to pass, a model that closes
its budget in the middle of the distribution and leaks in the corners fails on
the seed that draws a corner.

The whole design is one decision: which region of attribute space is inside
the hull models are normally fitted to, and which is outside. Write that down
in README.md before you tune the sampler.

Requirements the harness enforces:
  - generate(seed) returns (DataFrame, dict)
  - the frame has a 'time' column
  - it has exactly PERIOD_YEARS * 365 + SPINUP_DAYS rows
  - the same seed produces byte-identical output

Run this file directly to see the spread of catchments the sampler produces.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_YEARS = {period_years}
SPINUP_DAYS = {spinup_days}
N_STEPS = int(PERIOD_YEARS * 365) + SPINUP_DAYS

# The share of seeds drawn from outside the hull. Too low and the corners go
# unsampled; too high and the probe stops being about extrapolation and
# becomes a test of the corners alone.
OUT_OF_HULL_FRACTION = 0.35


# The attributes that define a catchment here, each with the range models are
# normally fitted over. TODO: these are the probe. Replace them with ranges you
# can defend, and say in README.md where the boundary comes from — a published
# catchment sample, a training set, an argument. A reviewer will ask why this
# line and not another one.
IN_HULL = {{
    "area_km2": (80.0, 1200.0),
    "soil_capacity_mm": (200.0, 450.0),
    "canopy_capacity_mm": (1.0, 3.0),
    "degree_day_factor_mm_per_C_day": (2.0, 4.5),
    "baseflow_coefficient": (0.003, 0.012),
    "latitude_deg": (35.0, 52.0),
    "aridity_index": (0.6, 1.4),
    "mean_annual_pr_mm": (700.0, 1400.0),
}}

# Where each attribute goes when it is the one pushed out of range.
OUT_OF_HULL = {{
    "soil_capacity_mm": (60.0, 140.0),          # thin, or fractured
    "baseflow_coefficient": (0.03, 0.09),       # flashy
    "latitude_deg": (58.0, 68.0),               # snow dominated
    "aridity_index": (1.8, 2.8),                # strongly water limited
    "mean_annual_pr_mm": (250.0, 450.0),        # dry
}}

# The share of seeds drawn from outside the hull. Too low and the corners go
# unsampled; too high and the probe stops being about extrapolation and becomes
# a test of the corners alone.
OUT_OF_HULL_FRACTION = 0.35


def draw_catchment(rng: np.random.Generator) -> tuple[dict, str | None]:
    """One catchment, and the name of the attribute pushed out of range.

    Exactly one attribute leaves the hull at a time, deliberately. Pushing
    several at once compounds into combinations no real catchment occupies — an
    arid climate on a soil too thin to hold anything sits nowhere near the
    Budyko curve — and a probe built on those is not testing extrapolation, it
    is testing whether a model can be right about a place that does not exist.
    Which single attribute is out of range is also what you report when the
    probe catches something.
    """
    static = {{
        name: float(rng.uniform(*bounds)) for name, bounds in IN_HULL.items()
    }}
    static["snow_threshold_degC"] = 0.0

    outside = None
    if rng.random() < OUT_OF_HULL_FRACTION:
        outside = str(rng.choice(sorted(OUT_OF_HULL)))
        static[outside] = float(rng.uniform(*OUT_OF_HULL[outside]))
    return static, outside


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)

    # Draw the catchment first, then weather consistent with it. Generating
    # weather that does not belong with the attributes produces a case no model
    # can be right about, which reads as a model failure and is a probe failure.
    static, outside = draw_catchment(rng)

    day = np.arange(N_STEPS)
    doy = day % 365

    # TODO: make the climate follow from the attributes. Here aridity is the
    # ratio of annual potential ET to annual precipitation, which is what puts
    # the catchment somewhere defensible on the Budyko curve; latitude sets the
    # temperature range and how much of the year falls below freezing.
    annual_pr = static["mean_annual_pr_mm"]
    annual_pet = static["aridity_index"] * annual_pr

    p_wet = 0.25 + 0.12 * np.cos(2 * np.pi * (doy - 30) / 365)
    wet = rng.random(N_STEPS) < p_wet
    shape = 0.7
    scale = annual_pr / 365.0 / max(float(p_wet.mean()), 1e-6) / shape
    pr = np.where(wet, rng.gamma(shape=shape, scale=scale, size=N_STEPS), 0.0)

    warmth = 26.0 - 0.42 * static["latitude_deg"]
    seasonal = warmth - 12.0 * np.cos(2 * np.pi * (doy - 15) / 365)
    noise = np.zeros(N_STEPS)
    innovation = rng.normal(0.0, 2.6, N_STEPS)
    for t in range(1, N_STEPS):
        noise[t] = 0.72 * noise[t - 1] + innovation[t]
    tas = seasonal + noise

    daylength = 1.0 + 0.35 * np.cos(2 * np.pi * (doy - 172) / 365)
    shape_pet = np.maximum(0.0, tas + 5.0) * daylength
    pet = shape_pet * annual_pet / max(float(shape_pet.sum()) / (N_STEPS / 365.0), 1e-6)

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {{
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
        }}
    )
    return forcing, static


if __name__ == "__main__":
    # The sanity check that matters: every draw must be a catchment an exact
    # model could be right about. If the out-of-hull rows sit somewhere no real
    # catchment does, the probe fails honest models and the fault is here.
    print(f"{{'seed':>10}}  {{'out of hull':<22}}  {{'P mm/yr':>8}}  {{'PET mm/yr':>9}}  "
          f"{{'aridity':>7}}  {{'soil mm':>7}}")
    for s in range(20260901, 20260913):
        frame, static = generate(s)
        rng = np.random.default_rng(s)
        _, outside = draw_catchment(rng)
        years = N_STEPS / 365
        print(f"{{s:>10}}  {{outside or '-':<22}}  "
              f"{{frame['pr'].sum() / years:>8.0f}}  "
              f"{{frame['pet'].sum() / years:>9.0f}}  "
              f"{{static['aridity_index']:>7.2f}}  "
              f"{{static['soil_capacity_mm']:>7.0f}}")

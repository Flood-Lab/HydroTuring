"""Seeded synthetic forcing for mass/routing-network-closure.

one pr/tas/pet series is supplied to the model. The proposed extension is on
the network side: static metadata declares a fixed Y network and compatible
models report reach-indexed inflow, outflow, and channel storage.

The two headwater sub-catchments have different contributing areas so the same
meteorological forcing can produce distinct volumetric headwater hydrographs.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SPINUP_DAYS = 365
PERIOD_DAYS = 90
N_STEPS = SPINUP_DAYS + PERIOD_DAYS


def _triangle(n: int, start: float, duration: float, peak: float) -> np.ndarray:
    """Daily triangular precipitation event in mm/day."""
    x = np.arange(n, dtype=float) + 0.5
    half = duration / 2.0
    out = np.zeros(n, dtype=float)

    rising = (x >= start) & (x < start + half)
    falling = (x >= start + half) & (x < start + duration)

    out[rising] = peak * (x[rising] - start) / half
    out[falling] = peak * (start + duration - x[falling]) / half
    return np.maximum(out, 0.0)


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    """Return 365 spin-up + 90 scored daily forcing and fixed network metadata."""
    rng = np.random.default_rng(seed)

    pr = np.zeros(N_STEPS, dtype=float)

    # Intermittent spin-up storms establish non-empty model states.
    for base in range(0, SPINUP_DAYS, 30):
        duration = float(rng.uniform(2.0, 5.0))
        start = min(
            float(base + rng.uniform(0.0, 8.0)),
            SPINUP_DAYS - duration,
        )
        peak = float(rng.uniform(10.0, 35.0))
        pr[:SPINUP_DAYS] += _triangle(
            SPINUP_DAYS, start, duration, peak
        )

    # One scored storm occurs in the first 20 days; >=70 dry days follow.
    duration = float(rng.uniform(3.0, 6.0))
    start = float(rng.uniform(0.0, 10.0))
    peak = float(rng.uniform(15.0, 40.0))
    pr[SPINUP_DAYS:] = _triangle(PERIOD_DAYS, start, duration, peak)

    # Warm, rain-only forcing; PET remains moderate.
    tas = np.full(N_STEPS, float(rng.uniform(8.0, 15.0)))  # degC
    pet = np.full(N_STEPS, float(rng.uniform(1.5, 3.5)))   # mm/day

    forcing = pd.DataFrame(
        {
            "time": pd.date_range(
                "2000-01-01", periods=N_STEPS, freq="D"
            ).strftime("%Y-%m-%d"),
            "pr": pr,
            "tas": tas,
            "pet": pet,
        }
    )

    # Proposed minimal network extension to static.json.
    # A and B are headwater sub-catchments; C receives A+B through junction J.
    static = {
        "routing_network": {
            "reaches": ["A", "B", "C"],
            "junctions": [
                {"id": "J", "upstream": ["A", "B"], "downstream": "C"}
            ],
            "headwaters": ["A", "B"],
            "outlet": "C",
            "contributing_area_km2": {
                "A": 60.0,
                "B": 90.0,
                "C": 150.0,
            },
            "groundwater_channel_exchange": False,
            "additional_lateral_channel_inflow": False,
        }
    }

    return forcing, static


if __name__ == "__main__":
    frame, static = generate(20260918)
    scored = frame.iloc[SPINUP_DAYS:]
    print(
        f"{len(frame)} daily steps; {PERIOD_DAYS} scored days; "
        f"scored precipitation={scored.pr.sum():.1f} mm; "
        f"areas={static['routing_network']['contributing_area_km2']}"
    )

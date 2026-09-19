"""Forcing generator for mass/routing-network-closure.

Deterministic given a seed, and never committed as data. The case isolates a
reach-resolved routing component: prescribed headwater discharge is its only
water input, and precipitation and potential evaporation are both zero so the
declared control-volume equation has no omitted source or sink.

Requirements expected by the proposed case:
  - generate(seed) returns (DataFrame, dict)
  - the frame has a 'time' column
  - it has exactly PERIOD_DAYS + SPINUP_DAYS rows
  - the same seed produces byte-identical output

Run this file directly to sanity check what it produces.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_DAYS = 90
SPINUP_DAYS = 365
N_STEPS = PERIOD_DAYS + SPINUP_DAYS
SECONDS_PER_DAY = 86_400.0

AREA_A_KM2 = 60.0
AREA_B_KM2 = 90.0
AREA_C_KM2 = AREA_A_KM2 + AREA_B_KM2
MIN_DRY_TAIL_DAYS = 70

# The capacity bounds are deliberately loose relative to generated event
# volumes. They catch invented stores without prescribing a routing equation.
STATIC = {
    "area_km2": AREA_C_KM2,
    "routing_network": {
        "reaches": ["A", "B", "C"],
        "junctions": [
            {"id": "J", "upstream": ["A", "B"], "downstream": "C"}
        ],
        "headwaters": ["A", "B"],
        "outlet": "C",
        "contributing_area_km2": {
            "A": AREA_A_KM2,
            "B": AREA_B_KM2,
            "C": AREA_C_KM2,
        },
        "channel_capacity_m3": {
            "A": 4.0e6,
            "B": 6.0e6,
            "C": 1.0e7,
        },
        "channel_evaporation": False,
        "groundwater_channel_exchange": False,
        "additional_lateral_channel_inflow": False,
        "withdrawal": False,
    },
}


def _triangle(
    n: int,
    start: int,
    duration: int,
    total_mm: float,
) -> np.ndarray:
    """Return a non-negative pulse whose discrete sum is total_mm."""
    pulse = np.zeros(n, dtype=float)
    weights = np.bartlett(duration + 2)[1:-1]
    weights /= weights.sum()
    pulse[start : start + duration] = total_mm * weights
    return pulse


def _add_spinup_events(
    rng: np.random.Generator,
    target: np.ndarray,
) -> None:
    """Add intermittent inflow events during the one-year spin-up."""
    for block_start in range(10, SPINUP_DAYS - 20, 35):
        start = block_start + int(rng.integers(0, 8))
        duration = int(rng.integers(3, 7))
        total_mm = float(rng.uniform(12.0, 45.0))
        target += _triangle(N_STEPS, start, duration, total_mm)


def _depth_to_discharge(depth_mm_day: np.ndarray, area_km2: float) -> np.ndarray:
    """Convert mm/day over area_km2 to interval-mean m3/s."""
    return depth_mm_day * 1.0e-3 * area_km2 * 1.0e6 / SECONDS_PER_DAY


def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)

    inflow_a_mm = np.zeros(N_STEPS, dtype=float)
    inflow_b_mm = np.zeros(N_STEPS, dtype=float)
    _add_spinup_events(rng, inflow_a_mm)
    _add_spinup_events(rng, inflow_b_mm)

    # Independent scored events occur early enough to leave at least 70 days
    # without prescribed inflow. The recession bound uses that dry tail.
    start_a = SPINUP_DAYS + int(rng.integers(1, 7))
    duration_a = int(rng.integers(3, 7))
    total_a_mm = float(rng.uniform(30.0, 65.0))
    inflow_a_mm += _triangle(N_STEPS, start_a, duration_a, total_a_mm)

    start_b = SPINUP_DAYS + int(rng.integers(8, 14))
    duration_b = int(rng.integers(3, 7))
    total_b_mm = float(rng.uniform(30.0, 65.0))
    inflow_b_mm += _triangle(N_STEPS, start_b, duration_b, total_b_mm)

    last_event_end = max(
        start_a + duration_a - SPINUP_DAYS,
        start_b + duration_b - SPINUP_DAYS,
    )
    assert PERIOD_DAYS - last_event_end >= MIN_DRY_TAIL_DAYS

    q_in_a = _depth_to_discharge(inflow_a_mm, AREA_A_KM2)
    q_in_b = _depth_to_discharge(inflow_b_mm, AREA_B_KM2)

    # forcing.csv keeps its standard columns, but both atmospheric water terms
    # are zero: prescribed channel inflow is the only water source and there is
    # no unreported open-water evaporation sink in the routing control volume.
    pr = np.zeros(N_STEPS, dtype=float)
    pet = np.zeros(N_STEPS, dtype=float)
    tas = np.full(N_STEPS, float(rng.uniform(8.0, 15.0)))

    time = pd.date_range("2000-01-01", periods=N_STEPS, freq="D")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d"),
            "pr": np.round(pr, 6),
            "tas": np.round(tas, 6),
            "pet": np.round(pet, 6),
            "q_in_A": np.round(q_in_a, 9),
            "q_in_B": np.round(q_in_b, 9),
        }
    )
    return forcing, dict(STATIC)


if __name__ == "__main__":
    frame, static = generate(20260903)
    repeated, repeated_static = generate(20260903)
    assert frame.equals(repeated) and static == repeated_static

    scored = frame.iloc[SPINUP_DAYS:]
    volume_a = scored["q_in_A"].sum() * SECONDS_PER_DAY
    volume_b = scored["q_in_B"].sum() * SECONDS_PER_DAY
    depth_a = volume_a / (AREA_A_KM2 * 1.0e3)
    depth_b = volume_b / (AREA_B_KM2 * 1.0e3)

    print(f"{len(frame)} steps")
    print(f"{depth_a:.1f} mm scored inflow to reach A")
    print(f"{depth_b:.1f} mm scored inflow to reach B")
    print(f"{int((scored[['q_in_A', 'q_in_B']] > 0).any(axis=1).sum())} inflow days")
    print(f"network: {static['routing_network']['reaches']}")

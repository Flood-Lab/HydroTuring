#!/usr/bin/env python3
"""HydroTuring adapter for reference_rating_drift: the routed bucket whose
gauge reads a slowly decaying memory of the last peak.

The routing and the discharge are honest. The fault is the stage: it is
solved from a peak-weighted signal rather than from the flow at the step,

    h(t) = f( q(t) + alpha * peak(t) ),   peak(t) = max(q(t), beta * peak(t-1))

with `beta` close to one, so the memory of a flood lingers for weeks. After a
large event the gauge stays high through everything that follows, and a small
event at the end of the season reads *above* a larger event at the start of
it. That is a rating that has decoupled from its flow, and it breaks
monotonicity against the running maximum without inverting the loop — the two
faults this probe separates.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

COLUMNS = [
    'time', 'pr', 'evspsbl', 'mrro', 'dis',
    'mrso', 'snw', 'canopy', 'channel', 'stage',
]

MODEL = {"name": "reference_rating_drift", "version": "1.0.0"}

EVAP_SHAPE = 0.5
SECONDS_PER_DAY = 86400.0

# How strongly the gauge is weighted toward the remembered peak, and how long
# the memory lasts. Both are large on purpose: this model has to fail.
#
# The decay is per *day*, and the memory has to outlast the gap between the
# probe's floods for the fault to show. The generator spaces its events 78 to
# 180 days apart, so a memory that halves in a month is gone by the time the
# next event arrives and the gauge would track the flow after all. At 0.997 a
# day the memory keeps a 70% hold across a 120-day dry spell, which is what
# makes a small late flood read above a large early one.
PEAK_WEIGHT = 2.5
PEAK_DECAY = 0.997

TIMESTEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}


def _stage(apparent_q_m3s, static):
    """Manning normal depth for a flow, in metres.

    The gauge is solved from a flow, not from a store, for the same reason
    the honest reference model is: a stage is a length read in a channel
    cross-section, and the flow through that section is what sets it. The
    fault this model carries is in *which* flow it is solved from, not in how
    the flow becomes a length.
    """
    area_km2 = static["area_km2"]
    width_m = static.get("width_m", 18.0)
    slope = static.get("slope", 0.0015)
    manning_n = static.get("manning_n", 0.035)
    q = max(float(apparent_q_m3s), 0.0)
    if q <= 0.0 or width_m <= 0.0 or slope <= 0.0:
        return 0.0
    return float((q * manning_n / (width_m * slope ** 0.5)) ** 0.6)


def simulate(forcing, static, dt_days=1.0):
    soil_cap = static["soil_capacity_mm"]
    canopy_cap = static["canopy_capacity_mm"]
    ddf = static["degree_day_factor_mm_per_C_day"]
    k_base = static["baseflow_coefficient"]
    t_snow = static["snow_threshold_degC"]
    area_km2 = static["area_km2"]

    soil = 0.5 * soil_cap
    swe = 0.0
    canopy = 0.0
    store = 0.0
    peak = 0.0
    rows = []

    for step in forcing:
        pr_rate, tas, pet_rate = step["pr"], step["tas"], step["pet"]
        pr = pr_rate * dt_days
        pet = pet_rate * dt_days

        snowfall = pr if tas < t_snow else 0.0
        rain = 0.0 if tas < t_snow else pr

        swe += snowfall
        melt = min(swe, ddf * max(tas - t_snow, 0.0) * dt_days)
        swe -= melt

        water_in = rain + melt
        intercepted = min(canopy_cap - canopy, water_in)
        canopy += intercepted
        throughfall = water_in - intercepted

        canopy_evap = min(canopy, pet)
        canopy -= canopy_evap
        pet_left = pet - canopy_evap

        soil += throughfall
        surface = max(0.0, soil - soil_cap)
        soil -= surface
        baseflow = k_base * soil * dt_days
        soil -= baseflow
        soil_evap = min(soil, pet_left * min(1.0, soil / (EVAP_SHAPE * soil_cap)))
        soil -= soil_evap

        store += surface + baseflow
        routed = store * (1.0 - math.exp(-dt_days / 2.0))
        store -= routed

        q_m3s = routed / dt_days * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY

        # The fault: the gauge is driven mostly by a lingering peak, so it
        # keeps reading high long after the flow has fallen away.
        peak = max(q_m3s, PEAK_DECAY ** dt_days * peak)
        apparent_q = q_m3s + PEAK_WEIGHT * peak
        stage = _stage(apparent_q, static)

        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (canopy_evap + soil_evap) / dt_days,
            "mrro": routed / dt_days,
            "dis": q_m3s,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
            "channel": store,
            "stage": stage,
        })
    return rows


def read_request(path: Path) -> tuple[dict, Path]:
    request = json.loads(path.read_text())
    return request, path.parent


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet"):
            if key in row:
                row[key] = float(row[key])
    return rows


def write_result(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    request, io_dir = read_request(request_path)
    forcing = read_forcing(io_dir / request["input"]["forcing"])
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    timestep = request.get("timestep", "PT1D")
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    rows = simulate(forcing, static, TIMESTEP_DAYS[timestep])

    write_result(io_dir / request["output"]["table"], rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows)}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

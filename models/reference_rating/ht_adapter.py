#!/usr/bin/env python3
"""HydroTuring adapter for reference_rating: the honest bucket with a gauge.

The catchment is the exact bucket. Its yield enters a reach with two parts:
a deep main channel, which is the slow store and the thing the staff gauge
stands in, and a shallow floodplain, which responds almost immediately to
whatever reaches it.

The two time constants are the entire point, and they are why a single store
cannot pass this probe. If the discharge and the stage are both functions of
one state variable, then at equal state they are equal, and at equal
discharge they are equal too: the two limbs of the rating coincide to machine
precision and there is no loop at any discharge. Giving the floodplain a
short residence time and the channel a long one separates them. On the rise
the floodplain is carrying water the channel has not taken up yet, so the
discharge climbs while the gauge is still low; on the fall the floodplain has
emptied, so the same discharge is carried with the channel still full and the
gauge reads higher. That is the direction the probe asks for.
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

MODEL = {"name": "reference_rating", "version": "1.0.0"}

EVAP_SHAPE = 0.5  # soil moisture at which evaporation reaches its potential rate

SECONDS_PER_DAY = 86400.0

# Residence times in days. The floodplain is the fast path and the channel
# the slow one; their separation is what makes the two limbs miss each other.
#
# Both are chosen so the loop a daily step can resolve clears the 2%-of-span
# floor the criterion applies, on every gate seed, with room to spare. An
# earlier version gave the floodplain a 0.35-day residence, so at a daily step
# it drained within a single row: the fast path was invisible, the gauge saw
# only the channel, and the rating it drew was single-valued at every seed.
# A loop is only readable at the step the model actually writes at, so the
# fast constant has to be slow enough to hold part of the flood in transit
# while the channel is still filling — and the channel slow enough that the
# two limbs separate by more than the floor. At these values the smallest loop
# over the three gate seeds is 3.5% of the rating, against a 2% floor.
FLOODPLAIN_RESIDENCE_D = 8.0
CHANNEL_RESIDENCE_D = 60.0

# The share of the yield that goes down the fast path. Most of a storm runs
# over the floodplain rather than into the channel, and that imbalance is what
# makes the channel lag the wave: the gauge has to fall well behind the flow
# for the two limbs to separate by more than the rating's own noise floor.
FLOODPLAIN_SHARE = 0.9

# How much of the floodplain's release passes the gauge as well. The gauge
# stands in the channel, so its reading is the channel's own release first;
# this term is the share of the bank flow that re-enters the reach and passes
# the same section. It is small, and it is there to keep the rating monotone
# in discharge: with none of it the gauge reads the channel alone, which lags
# the flow by enough that the binned rating dips a third of its span below its
# own running maximum on a few seeds. A little coupling to the flow removes
# that without flattening the loop, which comes from the store.
STAGE_FAST_SHARE = 0.08

TIMESTEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}


def _q_m3s(rate_mm_per_day: float, static: dict) -> float:
    """A depth rate over the catchment, in m3/s."""
    area_km2 = static["area_km2"]
    return max(rate_mm_per_day, 0.0) * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY


def manning_depth(q_m3s: float, static: dict) -> float:
    """The depth a steady flow makes in the reach's cross-section, in metres.

    The gauge stands in the main channel, so what it reads is the depth *that
    section* needs to carry the flow through it, through Manning's
    normal-depth relation:

        h = ( Q * n / (w * sqrt(S)) )^(3/5)

    The flow in the section is the channel's own release, `q_slow`, and not
    the reach's total discharge: the floodplain is the fast path running over
    the bank, and its water is not in the bed the gauge stands in. That
    distinction is the whole loop. On the rise the floodplain is carrying
    water the channel has not taken up yet, so the reach discharges more than
    the channel is passing and the gauge reads low; on the fall the floodplain
    has emptied and the same reach discharge is carried by a full channel, so
    the gauge reads high.

    An earlier version read the gauge off the *volume* the channel store holds,
    spread over the reach bed. That is not a rating a survey would recognise:
    the store is a depth over a 250 km2 catchment, and spreading it over 4.5 km
    of bed turned a few tens of millimetres of storage into a forty-metre
    stage. Manning's relation is the honest bridge from a flow to a length, and
    it puts the rating in the reach's real range.
    """
    width_m = float(static.get("width_m", 18.0))
    slope = float(static.get("slope", 0.0015))
    manning_n = float(static.get("manning_n", 0.035))
    q = max(float(q_m3s), 0.0)
    if q <= 0.0 or width_m <= 0.0 or slope <= 0.0:
        return 0.0
    return float((q * manning_n / (width_m * slope ** 0.5)) ** 0.6)


def simulate(forcing, static, dt_days=1.0):
    """The exact bucket, drained through a channel and a floodplain."""
    soil_cap = static["soil_capacity_mm"]
    canopy_cap = static["canopy_capacity_mm"]
    ddf = static["degree_day_factor_mm_per_C_day"]
    k_base = static["baseflow_coefficient"]
    t_snow = static["snow_threshold_degC"]
    area_km2 = static["area_km2"]

    k_fast = 1.0 - math.exp(-dt_days / FLOODPLAIN_RESIDENCE_D)
    k_slow = 1.0 - math.exp(-dt_days / CHANNEL_RESIDENCE_D)

    soil = 0.5 * soil_cap
    swe = 0.0
    canopy = 0.0
    fast = 0.0
    slow = 0.0
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

        # The yield splits between the two paths. Most of it goes over the
        # floodplain, which is shallow and passes it quickly; the rest enters
        # the channel, which is deep and holds it.
        yield_mm = surface + baseflow
        to_fast = FLOODPLAIN_SHARE * yield_mm
        to_slow = yield_mm - to_fast
        fast += to_fast
        slow += to_slow

        q_fast = fast * k_fast
        q_slow = slow * k_slow
        fast -= q_fast
        slow -= q_slow

        q_total = q_fast + q_slow
        dis_m3s = q_total / dt_days * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY
        # The reach's reported store is everything it is holding; the gauge
        # stands in the channel, so it reads the channel's own release, plus
        # the small share of the bank flow that re-enters the reach. That is
        # what makes the rating hysteretic rather than algebraic: the store the
        # probe pairs on carries the floodplain too, so at equal store the
        # channel is emptier on the rise than on the fall, and the gauge reads
        # lower there.
        gauge_q = q_slow + STAGE_FAST_SHARE * q_fast
        stage = manning_depth(_q_m3s(gauge_q / dt_days, static), static)

        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (canopy_evap + soil_evap) / dt_days,
            "mrro": q_total / dt_days,
            "dis": dis_m3s,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
            "channel": fast + slow,
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

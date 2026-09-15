#!/usr/bin/env python3
"""HydroTuring adapter for reference_leaky_router: the exact bucket with a
routing kernel that loses a fixed share of the water it carries.

The kernel is a three-day triangle that sums to 0.9 rather than to 1. A tenth
of every day's generated runoff enters the channel store and is never
released. The loss is invisible to a whole-catchment budget, because the water
is subtracted from the reach and not from the soil: the bucket still generates
the runoff it always did, and only the release is short. What it does show is
a channel store that grows across every step of a recession, when nothing is
entering the reach and the store can only drain.

It is distinct from `reference_stuck_router`, which shares the same kernel but
is pinned by a different criterion: that one is caught by the store *bound*
(`routing_conservation`) once enough water has accumulated. This one is
pinned by the store's *direction* on recession steps (`recession_drainage`),
which shows up long before any bound is reached.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLUMNS = ["time", "pr", "evspsbl", "mrro", "mrso", "snw", "canopy", "channel"]

MODEL = {"name": "reference_leaky_router", "version": "1.0.0"}


EVAP_SHAPE = 0.5  # soil moisture at which evaporation reaches its potential rate

# Steps the contract can name, as a fraction of a day.
TIMESTEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}


def simulate(forcing, static, dt_days=1.0):
    """The bucket, with a routing kernel that does not sum to one.

    Interception, degree-day snow, saturation-excess runoff and linear
    baseflow generate the runoff; the release is a three-day triangle scaled
    by 0.9, so a tenth of the flow is retained in the channel at every step.
    `channel` is reported as generated minus released, which is what the
    contract asks, and it therefore rises wherever the reach is not being fed.
    """
    soil_cap = static["soil_capacity_mm"]
    canopy_cap = static["canopy_capacity_mm"]
    ddf = static["degree_day_factor_mm_per_C_day"]
    k_base = static["baseflow_coefficient"]
    t_snow = static["snow_threshold_degC"]

    soil = 0.5 * soil_cap
    swe = 0.0
    canopy = 0.0
    rows = []
    generated = []
    cum_gen = cum_out = 0.0
    # A three-day kernel that sums to 0.9: the missing tenth never leaves.
    W = [0.45, 0.315, 0.135]

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

        generated.append(surface + baseflow)
        n = len(generated)
        routed = sum(W[k] * generated[n - 1 - k] for k in range(min(len(W), n)))
        cum_gen += surface + baseflow
        cum_out += routed
        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (canopy_evap + soil_evap) / dt_days,
            "mrro": routed / dt_days,
            "channel": cum_gen - cum_out,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
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

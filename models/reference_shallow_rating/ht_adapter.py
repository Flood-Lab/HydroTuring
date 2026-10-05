#!/usr/bin/env python3
"""HydroTuring adapter: the exact bucket with a rating from another reach.

This is a deliberately broken model, in the sense `docs/writing-a-probe.md`
asks for when no existing model fails a new probe: every quantity it reports
is the exact bucket's, computed the same way and conserving water to the same
floating point, with one thing changed. The gauge reads a depth, it rises
with the flow, and the depth it rises to is the one a reach five times wider
would need:

    stage = (Q n / (5 w sqrt(S)))**0.6

Everything else is honest. Rain becomes snow, interception, soil, baseflow
and evaporation exactly as `reference_bucket` computes them; `mrro` is the
bucket's runoff; `dis` is that runoff over the catchment area; `channel` is
identically zero because the bucket does not route.

Why this is the interesting lie, and why it is *this* model rather than a
gauge pinned at a constant. A stage and a discharge are two readings of one
cross-section. Given the discharge, the width of the section and gravity,
there is exactly one depth at which that section can carry that flow without
the water outrunning the wave it is made of — the subcritical branch of

    Fr = Q / (w * d**1.5 * sqrt(g)) <= 1

A rating drawn for a wider section reports a shallower depth, and depth is
what the velocity is divided by, so the Froude number the pair implies is
`(5)**0.9` times the honest one: the pair says the declared reach delivers
the flow faster than gravity allows, on most of the record and through every
flood.

It is a defect no budget and no rating test can see:

* the water is exact, so no budget moves — `mass/*` passes;
* the stage varies (`cv` 0.70) and rises monotonically with the flow, so
  `non_degenerate` and `rating_monotonic` pass;
* the bucket keeps one state and reads its gauge off it, so its rating is
  single-valued and `rating_loop` reads no loop — which is the honest
  statement for a model that holds no water in transit, and passes.

Two criteria read the two numbers against each other. `uniform_flow_friction`
does it where the flow is steady, and fails this gauge on
`momentum/uniform-flow-friction-consistency`: a borrowed width is wrong at
every flow, steady ones included. `froude_subcritical` does it on a record that
is never steady, and asks which regime the pair implies. A gauge that does not
move at all is caught by `non_degenerate` on
`momentum/stage-discharge-monotonic` for a different reason, so this model
deliberately does not use that shape.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "gwex", "mrso", "snw", "canopy", "channel", "stage"]

MODEL = {"name": "reference_shallow_rating", "version": "1.0.0"}

EVAP_SHAPE = 0.5  # soil moisture at which evaporation reaches its potential rate

SECONDS_PER_DAY = 86400.0

# The width the rating is drawn for, as a multiple of the one declared. Five
# is a mis-specification of the kind a survey makes — a rating curve reused
# from a wider, shallower reach — and it puts the implied Froude number
# `5**0.9` (~4.3) times the honest one, so the pair goes supercritical on
# most of the record rather than at one marginal step.
RATING_WIDTH_FACTOR = 5.0

# `reference_bucket`'s own fallbacks, so that the width factor is the only
# difference between the two models on every case, including the ones that
# declare no reach geometry. A different default here would make this model
# degenerate (a zero width gives zero depth at every step) rather than subtly
# wrong, which is the shape the docstring says it deliberately avoids.
DEFAULT_WIDTH_M = 18.0
DEFAULT_SLOPE = 0.0015
DEFAULT_MANNING_N = 0.035

# Steps the contract can name, as a fraction of a day. Forcing and reported
# fluxes are rates in mm per day at every step; the depth moved in one step
# is the rate times this.
TIMESTEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}


def stage_of(runoff_rate_mm_day: float, static: dict) -> float:
    """The level the gauge reports, in metres: normal depth for a wider reach.

    Manning's law, evaluated with the declared slope, roughness and bed, and
    with a width five times the one the case declares. The law is the right
    one and the geometry is not, which is exactly the failure a stage and a
    discharge are jointly able to expose and neither can expose alone: read
    against the flow it is a plausible rating, and read against the reach it
    is a depth that section cannot carry.

    The rating is a single-valued function of the flow, so the two limbs of a
    flood coincide and there is no loop to read. That is honest for a bucket
    that holds no water in transit, and it is the shape this model has to
    have: the defect has to survive every budget and every rating test.

    `stage` is an elevation on the case's fixed vertical datum, so the depth the
    borrowed rating makes is reported above the declared bed rather than as the
    level itself.
    """
    return float(static.get("bed_elevation_m", 0.0)) + _manning(
        max(runoff_rate_mm_day, 0.0), static
    )


def _manning(flow_rate_mm_day: float, static: dict) -> float:
    area_km2 = float(static.get("area_km2", 0.0))
    width_m = float(static.get("width_m", DEFAULT_WIDTH_M)) * RATING_WIDTH_FACTOR
    slope = float(static.get("slope", DEFAULT_SLOPE))
    manning_n = float(static.get("manning_n", DEFAULT_MANNING_N))
    # The same refusal the bucket makes, for the same reason: this is the
    # wide-rectangular relation, so a case declaring another section is one this
    # gauge cannot draw. The borrowed width is the only difference between the
    # two models, and a section this model silently read as rectangular would be
    # a second one.
    shape = str(static.get("cross_section_shape", "rectangular")).strip().lower()
    if shape != "rectangular":
        raise ValueError("reference_shallow_rating stage requires a rectangular section")
    if area_km2 <= 0.0 or width_m <= 0.0 or slope <= 0.0:
        return 0.0
    q_m3s = max(flow_rate_mm_day, 0.0) * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY
    if q_m3s <= 0.0:
        return 0.0
    return (q_m3s * manning_n / (width_m * slope ** 0.5)) ** 0.6


def discharge_m3s(runoff_rate_mm_day: float, static: dict) -> float:
    """The flow through the reach, in m3/s, from the catchment's runoff.

    Honest, and deliberately so: the defect is confined to the gauge. If the
    discharge were wrong too, the probe would be catching a mass error that
    the rest of the suite already catches.
    """
    area_km2 = float(static.get("area_km2", 0.0))
    return max(runoff_rate_mm_day, 0.0) * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY


def simulate(forcing, static, dt_days=1.0):
    """The exact bucket, step for step, reporting a gauge from another reach.

    Identical to `reference_bucket`'s `simulate` except for the `stage` cell:
    the water is conserved to the same floating point, so a failure here is a
    failure of the gauge and not of the bookkeeping.
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

    for step in forcing:
        pr_rate, tas, pet_rate = step["pr"], step["tas"], step["pet"]
        abstr_rate = step.get("abstr", 0.0)
        pr = pr_rate * dt_days
        pet = pet_rate * dt_days

        # The human term, first call on the store.
        want = max(abstr_rate, 0.0) * dt_days
        removed = min(soil, want)
        soil -= removed

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

        # The human term, second call: the day's outflow.
        rest = want - removed
        divert = min(surface, rest)
        surface -= divert
        removed += divert
        rest -= divert
        divert = min(baseflow, rest)
        baseflow -= divert
        removed += divert

        runoff = (surface + baseflow) / dt_days
        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (canopy_evap + soil_evap) / dt_days,
            "mrro": runoff,
            "dis": discharge_m3s(runoff, static),
            "gwex": -removed / dt_days,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
            # The bucket does not route, so nothing is in transit. Reported,
            # not omitted: it is a statement about the model.
            "channel": 0.0,
            # The one thing this model gets wrong.
            "stage": stage_of(runoff, static),
        })
    return rows


def read_request(path: Path) -> tuple[dict, Path]:
    request = json.loads(path.read_text())
    return request, path.parent


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet", "abstr"):
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
        json.dumps({
            "status": "ok",
            "model": MODEL,
            "n_steps": len(rows),
            "notes": {
                "stage": (
                    f"Manning normal depth drawn for a section "
                    f"{RATING_WIDTH_FACTOR:g}x the declared width"
                ),
                "channel": "identically zero; the bucket does not route",
            },
        }, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

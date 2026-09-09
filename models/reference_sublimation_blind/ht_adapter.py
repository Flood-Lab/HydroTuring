#!/usr/bin/env python3
"""HydroTuring adapter. Standard library only, to show that the /io contract
needs no scientific Python stack and could be written in any language.

Coherent for liquid water, blind to the phase change.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

MODE = "sublimation_blind"
MODEL = {"name": "reference_sublimation_blind", "version": "1.0.0"}

COLUMNS = [
    "time", "pr", "evspsbl", "mrro", "sbl", "hfls", "hfss", "hfg",
    "mrso", "snw", "canopy", "channel",
]

EVAP_SHAPE = 0.5   # soil moisture at which evaporation reaches its potential rate
SUBL_SHARE = 0.35  # share of the remaining demand a snowpack can meet
GROUND_SHARE = 0.10  # share of net radiation conducted into the ground

# Latent heat of vaporisation as A + B*T (T in degC) and of fusion, J kg-1.
LAMBDA_A, LAMBDA_B = 2.501e6, -2361.0
LAMBDA_F = 3.337e5
LAMBDA_CONST = 2.45e6  # the constant a careless model would use
CLIMATOLOGICAL_EF = 0.65  # evaporative fraction of a model whose heads never meet
ENERGY_LEAK = 0.15
SECONDS_PER_DAY = 86400.0

TIMESTEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}


def partition_energy(rn, tas, liquid_mm, sublimated_mm, dt_days):
    """Turn the water that left the surface into the energy that carried it.

    `reference_coupled` is exact: every kilogram is converted at the latent
    heat of the phase change it actually underwent, and the sensible flux is
    what net radiation has left once the latent and ground fluxes are taken.
    Solving for H that way is the physical statement that the surface has no
    other place to put the energy, and it is why this model closes the energy
    budget by construction, exactly as `reference_bucket` closes the water
    budget by construction. The discrimination lives in the other four.
    """
    seconds = dt_days * SECONDS_PER_DAY
    lam_v = LAMBDA_A + LAMBDA_B * tas
    lam_s = LAMBDA_A + LAMBDA_F
    ground = GROUND_SHARE * rn

    if MODE == "constant_lambda":
        latent = LAMBDA_CONST * (liquid_mm + sublimated_mm) / seconds
    elif MODE == "sublimation_blind":
        latent = lam_v * (liquid_mm + sublimated_mm) / seconds
    elif MODE == "two_head":
        # An energy head that never reads the water head: the partition is a
        # climatological evaporative fraction of the available energy, and it
        # has no idea how much water actually left.
        latent = CLIMATOLOGICAL_EF * (rn - ground)
    else:
        latent = (lam_v * liquid_mm + lam_s * sublimated_mm) / seconds

    sensible = rn - ground - latent
    if MODE == "energy_leak":
        sensible -= ENERGY_LEAK * rn
    return latent, sensible, ground


def simulate(forcing, static, dt_days=1.0):
    """The reference bucket, with snow sublimation and a surface energy budget.

    The water side is `reference_bucket` plus one term: a snowpack meets part
    of the evaporative demand by sublimating, which is removed from the pack
    like every other flux, so the water budget still closes to floating point.
    That term exists because without it no step in the record has a phase
    change to be coherent about.
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
        rn = step.get("rn", 0.0)
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

        sublimation = min(swe, SUBL_SHARE * pet_left) if swe > 0.0 else 0.0
        swe -= sublimation
        pet_left -= sublimation

        soil += throughfall
        surface = max(0.0, soil - soil_cap)
        soil -= surface
        baseflow = k_base * soil * dt_days
        soil -= baseflow
        soil_evap = min(soil, pet_left * min(1.0, soil / (EVAP_SHAPE * soil_cap)))
        soil -= soil_evap

        liquid = canopy_evap + soil_evap
        latent, sensible, ground = partition_energy(rn, tas, liquid, sublimation, dt_days)

        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (liquid + sublimation) / dt_days,
            # The sublimating share of the evaporation above, so a
            # criterion never has to infer which kilograms left as ice.
            "sbl": sublimation / dt_days,
            "mrro": (surface + baseflow) / dt_days,
            "hfls": latent,
            "hfss": sensible,
            "hfg": ground,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
            # No routing: runoff leaves the stores and the catchment in the
            # same step, so the water in transit is identically zero.
            # Reported rather than omitted, because it is a statement.
            "channel": 0.0,
        })
    return rows


def read_request(path: Path) -> tuple[dict, Path]:
    request = json.loads(path.read_text())
    return request, path.parent


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet", "rn"):
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

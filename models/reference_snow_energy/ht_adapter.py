#!/usr/bin/env python3
"""HydroTuring adapter. Standard library only, to show that the /io contract
needs no scientific Python stack and could be written in any language.

The exact model for the snowmelt coupling probe: melt is bought with energy.

`reference_coupled` melts a degree-day depth and then closes its energy budget
around whatever evaporated, so the fusion its own water budget reports is never
charged to its own energy budget. This model inverts that. The pack's surface
budget decides how much ice can melt, and the water budget reports exactly that
much, so the two ledgers agree at the phase change by construction.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

MODEL = {"name": "reference_snow_energy", "version": "1.0.0"}

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
SECONDS_PER_DAY = 86400.0

# Bulk turbulent exchange over the pack, W m-2 K-1, and the temperature of a
# melting snow surface. The sensible flux is signed the way the contract signs
# every sink: positive leaves the surface. Air warmer than the pack therefore
# gives a negative flux, which is energy arriving and available to melt ice.
K_H = 3.0
T_SURFACE_SNOW = 0.0

TIMESTEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}


def melt_from_energy(rn, tas, sublimation_mm, swe_mm, dt_days):
    """How much ice this step's surface energy budget can actually melt.

    Net radiation, less what is conducted into the ground, less the latent
    heat carried off by the pack's own sublimation, less what the turbulent
    flux takes or brings, is what remains to break the lattice. Divided by the
    latent heat of fusion it is a mass, and it is capped by the ice that is
    there: energy cannot melt snow that has already gone.
    """
    seconds = dt_days * SECONDS_PER_DAY
    if swe_mm <= 0.0:
        return 0.0

    ground = GROUND_SHARE * rn
    latent_pack = (LAMBDA_A + LAMBDA_F) * sublimation_mm / seconds
    sensible = K_H * (T_SURFACE_SNOW - tas)

    available = rn - ground - latent_pack - sensible
    if available <= 0.0:
        return 0.0
    # Capped by the ice that is there. On the step a pack melts out the cap
    # binds and the energy left over warms the ground that has just been
    # uncovered, which is where the reported sensible flux puts it.
    return min(swe_mm, available * seconds / LAMBDA_F)


def partition_energy(rn, tas, liquid_mm, sublimated_mm, melt_mm, dt_days):
    """Turn the water that changed phase into the energy that carried it.

    Every kilogram is converted at the latent heat of the phase change it
    actually underwent, fusion included, and the sensible flux is what net
    radiation has left once the latent, ground and melt terms are taken.
    Solving for H that way is the physical statement that the surface has no
    other place to put the energy, and it is why this model satisfies the
    melt-fusion identity by construction, exactly as `reference_bucket` closes
    the water budget by construction.
    """
    seconds = dt_days * SECONDS_PER_DAY
    lam_v = LAMBDA_A + LAMBDA_B * tas
    lam_s = LAMBDA_A + LAMBDA_F

    ground = GROUND_SHARE * rn
    latent = (lam_v * liquid_mm + lam_s * sublimated_mm) / seconds
    fusion = LAMBDA_F * melt_mm / seconds
    sensible = rn - ground - latent - fusion
    return latent, sensible, ground


def simulate(forcing, static, dt_days=1.0):
    """The reference bucket, with sublimation and an energy-balance snowpack.

    The water side is `reference_coupled` with one term moved: melt is no
    longer a function of air temperature. The pack sublimates first, because
    that is what fixes the latent flux leaving the snow surface, and what the
    surface budget has left after it is what melts ice.
    """
    soil_cap = static["soil_capacity_mm"]
    canopy_cap = static["canopy_capacity_mm"]
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

        # The pack's own fluxes, in the order the energy budget needs them:
        # sublimation fixes the latent heat leaving the snow surface, and only
        # then is it known what remains to melt.
        sublimation = min(swe, SUBL_SHARE * pet) if swe > 0.0 else 0.0
        swe -= sublimation

        melt = melt_from_energy(rn, tas, sublimation, swe, dt_days)
        swe -= melt

        # While a pack survives the step it is the surface: the soil and the
        # canopy are under it and neither evaporates through it. That is not a
        # simplification, it is what keeps the model coherent. The sensible
        # flux that decided the melt above is then the sensible flux reported
        # below, because the latent flux is the pack's sublimation and nothing
        # else. Let soil evaporation run under a metre of snow and the model
        # would melt ice on one energy budget and report another.
        snow_covered = swe > 0.0

        water_in = rain + melt
        if snow_covered:
            canopy_evap = 0.0
            throughfall = water_in
        else:
            intercepted = min(canopy_cap - canopy, water_in)
            canopy += intercepted
            throughfall = water_in - intercepted
            canopy_evap = min(canopy, pet - sublimation)
            canopy -= canopy_evap
        pet_left = pet - sublimation - canopy_evap

        soil += throughfall
        surface = max(0.0, soil - soil_cap)
        soil -= surface
        baseflow = k_base * soil * dt_days
        soil -= baseflow
        soil_evap = (
            0.0 if snow_covered
            else min(soil, pet_left * min(1.0, soil / (EVAP_SHAPE * soil_cap)))
        )
        soil -= soil_evap

        liquid = canopy_evap + soil_evap
        latent, sensible, ground = partition_energy(
            rn, tas, liquid, sublimation, melt, dt_days
        )

        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (liquid + sublimation) / dt_days,
            # The sublimating share of the evaporation above, so a criterion
            # never has to infer which kilograms left as ice.
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

    rows = simulate(forcing, static, dt_days=TIMESTEP_DAYS[timestep])
    write_result(io_dir / request["output"]["table"], rows)

    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows)}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

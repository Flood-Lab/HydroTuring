#!/usr/bin/env python3
"""Standard-library adapter with a skin satisfying the radiation identity."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

MODE = "radiative"
MODEL = {"name": "reference_radiative", "version": "1.0.0"}

COLUMNS = [
    "time", "pr", "evspsbl", "mrro", "sbl", "hfls", "hfss", "hfg", "rlus",
    "mrso", "snw", "canopy", "channel", "ts",
]

EVAP_SHAPE = 0.5   # soil moisture at which evaporation reaches its potential rate
SUBL_SHARE = 0.35  # share of the remaining demand a snowpack can meet
GROUND_SHARE = 0.10  # share of net radiation conducted into the ground

# Latent heat of vaporisation as A + B*T (T in degC) and of fusion, J kg-1.
LAMBDA_A, LAMBDA_B = 2.501e6, -2361.0
LAMBDA_F = 3.337e5
SECONDS_PER_DAY = 86400.0

# Stefan-Boltzmann constant (CODATA 2018) and bulk heat-transfer coefficient.
STEFAN_BOLTZMANN = 5.670374419e-8
CONDUCTANCE = 20.0  # W m-2 K-1
KELVIN = 273.15

TIMESTEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}


def partition_energy(rn, tas, liquid_mm, sublimated_mm, dt_days):
    """Keep reference_coupled's phase-specific latent heat and closed surface budget."""
    seconds = dt_days * SECONDS_PER_DAY
    lam_v = LAMBDA_A + LAMBDA_B * tas
    lam_s = LAMBDA_A + LAMBDA_F
    ground = GROUND_SHARE * rn
    latent = (lam_v * liquid_mm + lam_s * sublimated_mm) / seconds
    return latent, rn - ground - latent, ground


def radiate(tas, sensible, rlds, eps):
    """Diagnose Ts from H = g_H * (Ts - Ta), then total upward longwave.

    The interval-mean sensible flux is treated as holding at the sampling
    instant. The sign of H sets whether skin is warmer or cooler than air.
    Net radiation remains prescribed independently of longwave.
    MODE changes only the emitting temperature or the reflected-sky term."""
    air_k = tas + KELVIN
    skin_k = air_k + sensible / CONDUCTANCE
    emitter_k = air_k if MODE == "air_emitter" else skin_k
    upward = eps * STEFAN_BOLTZMANN * emitter_k ** 4
    if MODE != "no_reflection":
        upward += (1.0 - eps) * rlds
    return skin_k, upward


def simulate(forcing, static, dt_days=1.0):
    """Run reference_coupled's water and energy core, adding Ts and rlus.

    Snow sublimation remains part of evaporation and closes the water budget."""
    soil_cap = static["soil_capacity_mm"]
    canopy_cap = static["canopy_capacity_mm"]
    ddf = static["degree_day_factor_mm_per_C_day"]
    k_base = static["baseflow_coefficient"]
    t_snow = static["snow_threshold_degC"]
    # Use the supplied emissivity without fitting it to the outputs.
    eps = static["eps"]

    soil = 0.5 * soil_cap
    swe = 0.0
    canopy = 0.0
    rows = []

    for step in forcing:
        pr_rate, tas, pet_rate = step["pr"], step["tas"], step["pet"]
        rn = step.get("rn", 0.0)
        rlds = step["rlds"]
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
        skin_k, upward = radiate(tas, sensible, rlds, eps)

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
            "rlus": upward,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
            # Runoff leaves in the same step; no water remains in routing.
            "channel": 0.0,
            "ts": skin_k,
        })
    return rows


def read_request(path: Path) -> tuple[dict, Path]:
    request = json.loads(path.read_text())
    return request, path.parent


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet", "rn", "rlds"):
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

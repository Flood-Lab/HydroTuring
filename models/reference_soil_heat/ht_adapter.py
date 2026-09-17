#!/usr/bin/env python3
"""One-layer heat reference, using only the Python standard library."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

MODEL = {"name": "reference_soil_heat", "version": "1.0.0"}
COLUMNS = ["time", "pr", "hfls", "hfss", "hfg", "hfg_bottom", "tsoil_layer", "rn"]
TIMESTEP_DAYS = {"PT1H": 1.0 / 24.0}
STEFAN_BOLTZMANN = 5.670374419e-8
ALBEDO = 0.2
EMISSIVITY = 0.96
SENSIBLE_EXCHANGE = 12.0
TOP_CONDUCTANCE = 6.0
BOTTOM_CONDUCTANCE = 1.0


def _thermal_layer(static):
    """Configure the lumped control volume and its fixed material properties."""
    depth = float(static["soil_layer_depth_m"])
    capacity = float(static["soil_heat_capacity_areal"])
    initial = float(static["soil_temperature_initial"])
    if not all(math.isfinite(value) and value > 0 for value in (depth, capacity, initial)):
        raise ValueError("soil layer depth, areal heat capacity and initial temperature must be finite and positive")
    # The modeled layer extends from the surface to this supplied bottom.
    # C_A already integrates Cv over that depth; do not multiply it again.
    return {
        "top_depth_m": 0.0,
        "bottom_depth_m": depth,
        "heat_capacity_areal_j_m2_k": capacity,
        "heat_capacity_volumetric_j_m3_k": capacity / depth,
        "initial_temperature_k": initial,
    }


def _radiative_fluxes(layer, air, shortwave, longwave, deep):
    """Solve the massless skin and return its native flux equations."""
    absorbed = (1.0 - ALBEDO) * shortwave + EMISSIVITY * longwave
    skin = layer
    for _ in range(20):
        residual = (EMISSIVITY * STEFAN_BOLTZMANN * skin**4
                    + SENSIBLE_EXCHANGE * (skin - air)
                    + TOP_CONDUCTANCE * (skin - layer) - absorbed)
        derivative = (4.0 * EMISSIVITY * STEFAN_BOLTZMANN * skin**3
                      + SENSIBLE_EXCHANGE + TOP_CONDUCTANCE)
        change = residual / derivative
        skin -= change
        if abs(change) < 1e-11:
            break
    else:
        raise ValueError("radiative reference skin temperature did not converge")
    radiation = absorbed - EMISSIVITY * STEFAN_BOLTZMANN * skin**4
    sensible = SENSIBLE_EXCHANGE * (skin - air)
    top = TOP_CONDUCTANCE * (skin - layer)
    bottom = BOTTOM_CONDUCTANCE * (layer - deep)
    return radiation, sensible, top, bottom


def _simulate_radiation(forcing, static, seconds, temperature_scale, substep_seconds):
    """Integrate temperature and physical fluxes with matching RK4 quadrature.

    The fixed exchange coefficients and lower-temperature reservoir are this
    reference model's structure, not conditions imposed on other models.
    Wind and humidity are available forcing; this dry reference represents
    exchange with a fixed coefficient and has no evaporation.
    """
    material = _thermal_layer(static)
    capacity = material["heat_capacity_areal_j_m2_k"]
    initial = material["initial_temperature_k"]
    deep = initial
    count = max(1, math.ceil(seconds / substep_seconds))
    step_seconds = seconds / count
    layer = initial
    rows = []
    for step in forcing:
        air = float(step["tas"]) + 273.15
        shortwave, longwave = float(step["rsds"]), float(step["rlds"])

        def rates(temperature):
            fluxes = _radiative_fluxes(temperature, air, shortwave, longwave, deep)
            return (fluxes[2] - fluxes[3]) / capacity, fluxes

        integrated = [0.0] * 4
        for _ in range(count):
            k1, f1 = rates(layer)
            k2, f2 = rates(layer + 0.5 * step_seconds * k1)
            k3, f3 = rates(layer + 0.5 * step_seconds * k2)
            k4, f4 = rates(layer + step_seconds * k3)
            layer += step_seconds * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            for n in range(4):
                integrated[n] += step_seconds * (f1[n] + 2.0 * f2[n] + 2.0 * f3[n] + f4[n]) / 6.0
        radiation, sensible, top, bottom = (value / seconds for value in integrated)
        rows.append({
            "time": step["time"], "pr": step["pr"], "hfls": 0.0,
            "hfss": sensible, "hfg": top, "hfg_bottom": bottom,
            "tsoil_layer": initial + temperature_scale * (layer - initial),
            # Optional native output for inspecting the surface budget. This
            # is not a forcing column or a newly required HydroTuring variable.
            "rn": radiation,
        })
    return rows


def simulate(forcing, static, dt_days=1.0 / 24.0, temperature_scale=1.0,
             *, substep_seconds=60.0):
    """Use incoming shortwave/longwave and the explicitly configured layer."""
    return _simulate_radiation(forcing, static, dt_days * 86400.0,
                               temperature_scale, substep_seconds)


def main(model=MODEL, simulate_model=simulate) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    io_dir = request_path.parent
    request = json.loads(request_path.read_text())
    with (io_dir / request["input"]["forcing"]).open(newline="") as fh:
        forcing = list(csv.DictReader(fh))
    for step in forcing:
        for key in ("pr", "tas", "rsds", "rlds"):
            if key in step:
                step[key] = float(step[key])
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    timestep = request["timestep"]
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    rows = simulate_model(forcing, static, TIMESTEP_DAYS[timestep])

    output_path = io_dir / request["output"]["table"]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok", "model": model, "n_steps": len(rows),
        "thermal_scope": "fixed homogeneous layer; no water transport or phase change",
        "radiation_input": "incoming rsds and rlds",
        "thermal_layer": _thermal_layer(static),
        "time_convention": {
            "time": "interval start",
            "hfg": "interval-mean downward flux at the soil surface",
            "hfg_bottom": "interval-mean downward flux at thermal_layer.bottom_depth_m",
            "tsoil_layer": "interval-end mean temperature over the configured layer",
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

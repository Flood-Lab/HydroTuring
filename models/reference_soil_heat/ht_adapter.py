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
    first_air = float(forcing[0]["tas"]) + 273.15
    capacity = float(static.get("soil_heat_capacity_areal", 400000.0))
    initial = float(static.get("soil_temperature_initial", first_air))
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
    """Use incoming radiation, or the existing prescribed-Rn case format."""
    if "rn" in forcing[0]:
        return _simulate_net_radiation(forcing, static, dt_days, temperature_scale)
    if "rsds" not in forcing[0] or "rlds" not in forcing[0]:
        raise ValueError("soil heat reference requires rn, or both rsds and rlds")
    return _simulate_radiation(forcing, static, dt_days * 86400.0,
                               temperature_scale, substep_seconds)


def _simulate_net_radiation(forcing, static, dt_days, temperature_scale):
    """Integrate native flux laws and the layer temperature analytically.

    With u=T_layer-T_deep, a the sensible exchange coefficient, Ks the top
    conductance and Kb the bottom conductance, the massless skin obeys
    Rn=a*(Ts-Ta)+Ks*(Ts-T_layer). Substitution gives
    C_A*du/dt=eta*(Rn+a*(Ta-T_deep))-(a*eta+Kb)*u, eta=Ks/(a+Ks).

    Rn and air temperature are constant in each supplied interval. The exact
    mean u determines the mean native fluxes; no flux uses a storage residual.
    temperature_scale changes only reported temperature for the two controls.
    """
    # Other surface-energy probes do not prescribe a thermal layer. Use the
    # reference's fixed 0.2 m layer (Cv=2 MJ m-3 K-1), initialized using only
    # the first forcing row. Explicit case settings always take precedence.
    first_air = float(forcing[0]["tas"]) + 273.15
    capacity = float(static.get("soil_heat_capacity_areal", 400000.0))
    initial = float(static.get("soil_temperature_initial", first_air))
    deep = float(static.get("soil_deep_temperature", first_air))
    a = float(static.get("soil_sensible_exchange", 12.0))
    ks = float(static.get("soil_top_conductance", 6.0))
    kb = float(static.get("soil_bottom_conductance", 1.0))
    seconds = dt_days * 86400.0
    eta = ks / (a + ks)
    loss = a * eta + kb
    tau = capacity / loss
    fraction = -math.expm1(-seconds / tau)
    mean_fraction = tau * fraction / seconds
    u = initial - deep
    rows = []
    for step in forcing:
        rn = float(step["rn"])
        air = float(step["tas"]) + 273.15 - deep
        equilibrium = eta * (rn + a * air) / loss
        end = u + (equilibrium - u) * fraction
        mean = equilibrium + (u - equilibrium) * mean_fraction
        skin_mean = (rn + a * air + ks * mean) / (a + ks)
        rows.append({
            "time": step["time"],
            "pr": step["pr"],
            "hfls": 0.0,
            "hfss": a * (skin_mean - air),
            "hfg": ks * (skin_mean - mean),
            "hfg_bottom": kb * mean,
            "tsoil_layer": initial + temperature_scale * (deep + end - initial),
            "rn": rn,
        })
        u = end
    return rows


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
        for key in ("pr", "tas", "rn", "rsds", "rlds"):
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
        "radiation_input": "prescribed rn" if "rn" in forcing[0] else "incoming rsds and rlds",
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

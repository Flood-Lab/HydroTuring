#!/usr/bin/env python3
"""Analytic one-layer heat reference, using only the Python standard library."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

MODEL = {"name": "reference_soil_heat", "version": "1.0.0"}
COLUMNS = ["time", "pr", "hfls", "hfss", "hfg", "hfg_bottom", "tsoil_layer"]
TIMESTEP_DAYS = {"PT1H": 1.0 / 24.0}


def simulate(forcing, static, dt_days=1.0 / 24.0, temperature_scale=1.0):
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
        for key in ("pr", "tas", "rn"):
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
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

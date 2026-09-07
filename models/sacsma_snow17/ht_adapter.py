#!/usr/bin/env python3
"""HydroTuring adapter for SAC-SMA + Snow-17 + unit hydrograph.

The physics is in sacsma_snow17.py, a port of the NWS Fortran. This file
is the units layer and the bookkeeping: it feeds the probe's forcing to
Snow-17 and SAC-SMA at the case's step, routes the channel inflow through
the gamma unit hydrograph, and reports every store the model has.

Choices that make it a conservative reference, all stated:

* Snow-17's SCF, the gauge-catch multiplier on snowfall, is 1.0. At its
  calibrated values above one the model manufactures snow to correct for
  gauge undercatch; a physical reference is fed what fell.
* SAC-SMA's SIDE, the ratio of deep to channel baseflow, sends water out
  of the catchment for good. It is kept at a small value and the loss is
  declared as a negative `gwex`, so the budget can close over what the
  model says it did with the water.
* RIVA, riparian evapotranspiration drawn from channel inflow, is part of
  the reported evaporation, as the Fortran counts it.
* The catchment's soil capacity sets UZTWM + UZFWM + LZTWM, scaled from
  the default set keeping their ratios; the lower-zone free water
  capacities keep their defaults. Canopy: SAC-SMA has no interception
  store, so `canopy` is identically zero.

Stores reported (catchment-average, weighting the pervious-area stores by
PAREA and the ADIMP store by ADIMP, as the runoff components are):
`mrso` = UZTWC + UZFWC + LZTWC (+ ADIMC on its area), `gw` = LZFSC + LZFPC,
`snw` = Snow-17's total water equivalent (WE + LIQW + lagged excess +
storage), `channel` = channel inflow generated but not yet released by
the unit hydrograph.

Steps: SAC-SMA takes the step in days and Snow-17 in whole hours, as the
Fortran does; PT1D and PT1H are the steps the model is defined at.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sys
from pathlib import Path

from sacsma_snow17 import SacState, SnowState, gamma_uh, sac1, sac_storage, snow17

COLUMNS = ["time", "pr", "evspsbl", "mrro", "gwex", "mrso", "snw", "canopy", "gw", "channel"]
MODEL = {"name": "sacsma_snow17", "version": "1.0.0"}
STEP_HOURS = {"PT1D": 24, "PT1H": 1}

# A mid-range parameter set from the NWS calibration guidance (Anderson
# 2002) and CAMELS-scale calibrations; the catchment's capacity rescales the
# tension stores. SIDE is small so the declared deep loss is exercised.
SAC = {
    "uztwm": 60.0, "uzfwm": 30.0, "uzk": 0.3, "pctim": 0.01, "adimp": 0.02, "riva": 0.02,
    "zperc": 80.0, "rexp": 2.0, "lztwm": 130.0, "lzfsm": 30.0, "lzfpm": 120.0,
    "lzsk": 0.06, "lzpk": 0.008, "pfree": 0.15, "side": 0.02, "rserv": 0.3,
}
SNOW = {
    "scf": 1.0, "mfmax": 1.2, "mfmin": 0.4, "uadj": 0.06, "si": 500.0, "nmf": 0.15,
    "tipm": 0.1, "mbase": 0.0, "pxtemp": 0.0, "plwhc": 0.04, "daygm": 0.1,
}
ADC = [0.05, 0.15, 0.26, 0.45, 0.5, 0.56, 0.61, 0.65, 0.69, 0.82, 1.0]  # areal depletion curve
UH = {"shape": 2.5, "scale_days": 1.0}


def surface_pressure_hpa(elev_m: float) -> float:
    """The NCAR driver's fit of surface pressure to elevation."""
    e = elev_m / 100.0
    return 33.86 * (29.9 - 0.335 * e + 0.00022 * e ** 2.4)


def parameters(static: dict) -> tuple[dict, dict]:
    sac = dict(SAC)
    if "soil_capacity_mm" in static:
        tension = sac["uztwm"] + sac["uzfwm"] + sac["lztwm"]
        factor = float(static["soil_capacity_mm"]) / tension
        for key in ("uztwm", "uzfwm", "lztwm"):
            sac[key] *= factor
    snow = dict(SNOW)
    if "snow_threshold_degC" in static:
        snow["pxtemp"] = float(static["snow_threshold_degC"])
        snow["mbase"] = float(static["snow_threshold_degC"])
    return sac, snow


def simulate(forcing: list[dict], static: dict, timestep: str) -> tuple[list[dict], dict]:
    hours = STEP_HOURS[timestep]
    dt_days = hours / 24.0
    sac, snow = parameters(static)
    lat = float(static.get("latitude_deg", 40.0))
    elev = float(static.get("elevation_m", 500.0))
    pa = surface_pressure_hpa(elev)
    weights = gamma_uh(UH["shape"], UH["scale_days"], dt_days)

    sac_state = SacState(uztwc=0.5 * sac["uztwm"], uzfwc=0.2 * sac["uzfwm"], lztwc=0.5 * sac["lztwm"],
                         lzfsc=0.2 * sac["lzfsm"], lzfpc=0.5 * sac["lzfpm"], adimc=0.5 * (sac["uztwm"] + sac["lztwm"]))
    snow_state = SnowState(nexlag=5 // hours + 2)

    generated: list[float] = []
    rows = []
    for step in forcing:
        when = dt.datetime.fromisoformat(str(step["time"]))
        pr_rate, tas, pet_rate = step["pr"], step["tas"], step["pet"]
        pcp = pr_rate * dt_days
        pet = pet_rate * dt_days

        raim, _snowfall = snow17(hours, when.year, when.month, when.day, pcp, tas, lat, snow, pa, ADC, snow_state)
        fluxes = sac1(dt_days, raim, pet, sac, sac_state)

        generated.append(fluxes["tci"])
        n = len(generated)
        routed = sum(weights[k] * generated[n - 1 - k] for k in range(min(len(weights), n)))
        soil, lower = sac_storage(sac, sac_state)
        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": fluxes["tet"] / dt_days,
            "mrro": routed / dt_days,
            "gwex": -fluxes["bfncc"] / dt_days,  # deep baseflow leaves the catchment
            "mrso": soil,
            "snw": snow_state.total(),
            "canopy": 0.0,
            "gw": lower,
            "channel": 0.0,  # filled below from the hydrograph's contents
        })
    cum_in = cum_out = 0.0
    for row, g in zip(rows, generated):
        cum_in += g
        cum_out += row["mrro"] * dt_days
        row["channel"] = cum_in - cum_out
    notes = {
        "port": "sacsma_snow17.py, checked against the f2py build of the Fortran",
        "parameters": {"sac": sac, "snow17": snow, "unit_hydrograph": UH},
        "canopy": "identically zero; SAC-SMA has no interception store",
        "gwex": "minus the non-channel baseflow (SIDE); deep groundwater leaving the catchment",
    }
    return rows, notes


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet"):
            if key in row:
                row[key] = float(row[key])
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent
    forcing = read_forcing(io_dir / request["input"]["forcing"])
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    timestep = request.get("timestep", "PT1D")
    if timestep not in STEP_HOURS:
        raise SystemExit(f"unsupported timestep {timestep!r}: Snow-17 is defined on whole hours")
    rows, notes = simulate(forcing, static, timestep)
    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows), "notes": notes}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

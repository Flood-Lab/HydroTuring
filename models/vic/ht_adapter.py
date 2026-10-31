#!/usr/bin/env python3
"""HydroTuring adapter for VIC 5.1.0, classic driver (UW-Hydro/VIC, MIT).

Runs upstream's vic_classic.exe, compiled unmodified from tag 5.1.0, on one grid cell.
It writes the global parameter, soil, vegetation and forcing files VIC reads, runs it
once, reads its ASCII output and maps it onto the /io contract. README.md lists every
parameter and what every choice moves; this docstring says only what the code does.

Configuration
-------------
FULL_ENERGY = TRUE, QUICK_FLUX = FALSE, FROZEN_SOIL = FALSE, ten thermal nodes (the
settings of upstream's tests/science/global_param.classic.5.0_snotel.txt). Soil, the
vegetation library and the tiles are upstream's Stehekin sample (cell 86340, pinned data
commit), with capacities from the case: the layer depths scale so the column holds
`soil_capacity_mm`, and VEG_LAI_WATER_FACTOR and VEG_LAI_SNOW_MULTIPLIER scale so the
cell's canopy holds `canopy_capacity_mm`. VIC_HT_FULL_ENERGY=0 gives the water-balance
mode, for the README's diagnostic.

The step
--------
Rows go in one to one: MODEL_STEPS_PER_DAY = SNOW_STEPS_PER_DAY = FORCE_STEPS_PER_DAY =
the row's step, and ATMOS_STEPS_PER_DAY follows it below the hour, because the driver's
default of 24 must divide evenly by the model step (get_global_param.c:839).
VIC 5.1.0 has no daily forcing and the classic driver no daily model step. A snow step
under 4 a day stops at get_global_param.c:738, a forcing step other than the snow step at
:979, and MODEL_STEPS_PER_DAY = 1 with sub-daily forcing at read_atmos_data.c:60, though
GlobalParam.md still offers that setting. So a daily row needs the `substep` the manifest
declares: its values are held over one record per substep, precipitation split evenly, and
VIC aggregates the steps back to the day: a flux is the day's mean, the stores and the
outgoing longwave and radiative temperature are the last step's (AGGFREQ, OUTVAR). Nothing
is read from a later row. VIC_HT_DAILY_RECORDS=<n> overrides the count, for the README's
sensitivity table. It is not set in the evaluated image.

Forcing VIC needs that the probes do not generate, MOCKED from each row alone, with no
calendar and nothing from another row (the same mock as models/summa):
  SWDOWN   the row's pet through Priestley-Taylor at a fixed 20 C to a net radiation,
           through an albedo of 0.23
  LWDOWN   the probe's `rlds` where it supplies one; otherwise the emission of a
           reference surface at air temperature, so no longwave deficit
  VP       RH_REF (0.70) of saturation at the row's air temperature
  PRESSURE a standard atmosphere at the case's elevation
  WIND     2 m/s at the sample's WIND_H
VIC then computes its own potential evaporation from those: it never reads `pet`.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

MODEL = {"name": "vic", "version": "5.1.0-14a371a.4"}
FULL_ENERGY = os.environ.get("VIC_HT_FULL_ENERGY", "1") == "1"
COLUMNS = ["time", "pr", "evspsbl", "mrro"] + (["hfls", "hfss", "hfg", "rlus", "ts"] if FULL_ENERGY else []) + \
          ["mrso", "snw", "canopy"]
DAILY_RECORDS = int(os.environ.get("VIC_HT_DAILY_RECORDS", "0"))

NOT_REPORTED = {
    "gw": "VIC has no store below the third soil layer",
    "channel": "the classic driver has no routing",
    "dis": "the classic driver has no routing",
    "gwex": "VIC has no exchange with the outside",
    "sbl": "VIC's vapour fluxes are net, negative on deposition",
    "snm": "VIC has no output for liquid leaving the snow module with rain passed through",
    "OUT_SURFSTOR": "no lakes, so identically zero; it stays in the budget in notes.water_budget",
}

VIC_EXE = os.environ.get("VIC_EXE", "/opt/vic/vic_classic.exe")
WORK = Path(os.environ.get("VIC_WORKDIR", "/tmp/vic"))
SAMPLE = Path(os.environ.get("VIC_SAMPLE", "/model/stehekin"))
DIAG = os.environ.get("VIC_HT_DIAG", "0") == "1"

STEPS_PER_DAY = {"PT1D": 1, "PT1H": 24, "PT15M": 96, "PT5M": 288, "PT1M": 1440}

RH_REF = float(os.environ.get("VIC_HT_RH", "0.70"))
WIND_M_S = float(os.environ.get("VIC_HT_WIND", "2.0"))
ALBEDO_REF = float(os.environ.get("VIC_HT_ALBEDO", "0.23"))
PARAMETER_SET = os.environ.get("VIC_HT_PARAMS", "stehekin")
VEG_SET = os.environ.get("VIC_HT_VEG", "cell")
PT_ALPHA = 1.26
PT_TREF_C = 20.0
SIGMA = 5.670374419e-8
LW_EMISSIVITY_REF = 1.0

CELL_ID = "86340"  # the sample's first cell
# VIC's own defaults (initialize_parameters.c), which the ramp's half-width keeps.
SNOW_RAMP_HALF_WIDTH_C = 0.5

# Soil values that differ in the alternative set: the runoff and baseflow parameters of
# the worked example in upstream's docs/Documentation/Drivers/Classic/SoilParam.md.
ALTERNATIVE_SOIL = {"infilt": 0.01, "Ds": 1.0e-4, "Dsmax": 3.05, "Ws": 0.93}

OUTPUT_VARS = ["OUT_PREC", "OUT_EVAP", "OUT_RUNOFF", "OUT_BASEFLOW", "OUT_SWE", "OUT_SNOW_CANOPY",
               "OUT_WDEW", "OUT_SURFSTOR", "OUT_WATER_ERROR", "OUT_PET", "OUT_SNOW_MELT",
               "OUT_SUB_SNOW", "OUT_SUB_CANOP", "OUT_RAINF", "OUT_SNOWF",
               "OUT_LATENT", "OUT_LATENT_SUB", "OUT_SENSIBLE", "OUT_GRND_FLUX", "OUT_DELTAH", "OUT_FUSION",
               "OUT_ADV_SENS", "OUT_DELTACC", "OUT_RFRZ_ENERGY", "OUT_SWNET", "OUT_LWNET", "OUT_LWDOWN",
               "OUT_SURF_TEMP", "OUT_RAD_TEMP", "OUT_ENERGY_ERROR"]
# One instant's values, on a row of several steps, as the contract's rlus and ts are.
END_OF_ROW = {"OUT_LWNET", "OUT_LWDOWN", "OUT_RAD_TEMP"}
NLAYER = 3


def parse_time(text: str) -> datetime:
    return datetime.fromisoformat(text.strip().replace(" ", "T"))


def standard_pressure_pa(elevation_m: float) -> float:
    return 101325.0 * (1.0 - 2.25577e-5 * elevation_m) ** 5.25588


def pt_w_per_mm_day(tc: float, pressure_pa: float) -> float:
    """Net radiation (W m-2) that Priestley-Taylor turns into 1 mm/day at tc."""
    es = 0.6108 * math.exp(17.27 * tc / (tc + 237.3))
    slope = 4098.0 * es / (tc + 237.3) ** 2
    gamma = 0.000665 * pressure_pa / 1000.0
    lam = 2.501e6 - 2361.0 * tc
    return lam * (slope + gamma) / (PT_ALPHA * slope) / 86400.0


def mock_atmosphere(tas: float, pet: float, rlds, pt_factor: float, pressure_pa: float):
    """SWDOWN, LWDOWN (W m-2), VP (kPa) from one row. See the module docstring."""
    rn = pet * pt_factor
    sw = max(rn, 0.0) / (1.0 - ALBEDO_REF)
    lw = rlds if rlds is not None else LW_EMISSIVITY_REF * SIGMA * (tas + 273.15) ** 4 + min(rn, 0.0)
    es_pa = 611.2 * math.exp(17.67 * tas / (tas + 243.5))
    return sw, max(lw, 0.0), RH_REF * es_pa / 1000.0


def sample_soil_row() -> list[str]:
    for line in (SAMPLE / "Stehekin_soil.txt").read_text().splitlines():
        fields = line.split()
        if len(fields) > 2 and fields[1] == CELL_ID:
            return fields
    raise RuntimeError(f"cell {CELL_ID} is not in the Stehekin soil file")


def sample_tiles() -> list[tuple[int, str, list[float]]]:
    """(class, the class's Cv-and-roots line, its monthly LAI) for the sample cell."""
    lines = (SAMPLE / "Stehekin_vegparam.txt").read_text().splitlines()
    for i, line in enumerate(lines):
        fields = line.split()
        if len(fields) == 2 and fields[0] == CELL_ID:
            count = int(fields[1])
            tiles = []
            for j in range(count):
                head = lines[i + 1 + 2 * j]
                lai = [float(v) for v in lines[i + 2 + 2 * j].split()]
                tiles.append((int(head.split()[0]), head, lai))
            return tiles
    raise RuntimeError(f"cell {CELL_ID} is not in the Stehekin vegetation file")


def write_soil(path: Path, static: dict, latitude: float, longitude: float, tas: list[float],
               rows_per_year: int) -> dict:
    fields = sample_soil_row()
    # 12*Nlayer+17 columns plus July_Tavg in the sample; see SoilParam.md for the indices.
    bulk = [float(fields[i]) for i in range(7 * NLAYER + 12, 8 * NLAYER + 12)]
    dens = [float(fields[i]) for i in range(8 * NLAYER + 12, 9 * NLAYER + 12)]
    depth_index = range(4 * NLAYER + 10, 5 * NLAYER + 10)
    depths = [float(fields[i]) for i in depth_index]
    porosity = [1.0 - b / d for b, d in zip(bulk, dens)]
    sample_capacity = sum(z * p * 1000.0 for z, p in zip(depths, porosity))
    capacity = float(static.get("soil_capacity_mm", sample_capacity))
    scale = capacity / sample_capacity

    fields[2], fields[3] = "%.4f" % latitude, "%.4f" % longitude
    fields[4 * NLAYER + 9] = "%.3f" % float(static.get("elevation_m", 0.0))
    for i in depth_index:
        fields[i] = "%.6f" % (float(fields[i]) * scale)
    for i in range(3 * NLAYER + 9, 4 * NLAYER + 9):
        fields[i] = "%.6f" % (float(fields[i]) * scale)
    if PARAMETER_SET == "alternative":
        for name, index in (("infilt", 4), ("Ds", 5), ("Dsmax", 6), ("Ws", 7)):
            fields[index] = "%.6g" % ALTERNATIVE_SOIL[name]
    elif PARAMETER_SET != "stehekin":
        raise SystemExit(f"VIC_HT_PARAMS must be stehekin or alternative, not {PARAMETER_SET!r}")
    if FULL_ENERGY:
        # A column's annual mean temperature, taken from the first year of the record.
        head = tas[:rows_per_year]
        fields[5 * NLAYER + 10] = "%.4f" % (sum(head) / len(head))
    fields[9 * NLAYER + 12] = "0"  # off_gmt: the record's clock is the local one
    path.write_text(" ".join(fields) + "\n")
    return {"soil_capacity_mm": capacity, "depth_scale": scale, "sample_capacity_mm": sample_capacity,
            "layer_depths_m": [round(z * scale, 6) for z in depths]}


def write_vegetation(path: Path, static: dict) -> dict:
    capacity = float(static.get("canopy_capacity_mm", -1.0))
    tiles = sample_tiles()
    if VEG_SET == "grass":
        tiles = [(10, tiles[3][1].replace(tiles[3][1].split()[1], "1.000000", 1), tiles[3][2])]
    elif VEG_SET != "cell":
        raise SystemExit(f"VIC_HT_VEG must be cell or grass, not {VEG_SET!r}")
    if capacity == 0.0:
        path.write_text(f"{CELL_ID} 0\n")
        return {"vegetation": "bare soil: the case has no canopy", "params": {}}
    with open(path, "w") as fh:
        fh.write(f"{CELL_ID} {len(tiles)}\n")
        for _, head, lai in tiles:
            fh.write(head + "\n" + " ".join("%.3f" % v for v in lai) + "\n")
    params = {}
    if capacity > 0.0:
        cover = [float(head.split()[1]) for _, head, _ in tiles]
        cell_lai = max(sum(c * lai[m] for c, (_, _, lai) in zip(cover, tiles)) for m in range(12))
        # Wdmax = factor x LAI, in mm; the snow store's ceiling is 4 x multiplier x LAI, in m.
        params = {"VEG_LAI_WATER_FACTOR": capacity / cell_lai,
                  "VEG_LAI_SNOW_MULTIPLIER": capacity / 1000.0 / (4.0 * cell_lai)}
    return {"vegetation": "the sample cell's tiles and monthly LAI" if VEG_SET == "cell"
                          else "one grass tile (class 10) from the sample cell", "params": params}


def write_forcing(path: Path, pr, tas, pet, rlds, dt_days: float, per_row: int, pressure_pa: float,
                  pt_factor: float) -> None:
    with open(path, "w") as fh:
        for i in range(len(pr)):
            sw, lw, vp = mock_atmosphere(tas[i], pet[i], rlds[i] if rlds else None, pt_factor, pressure_pa)
            # VIC reads precipitation as a depth per forcing record.
            depth = pr[i] * dt_days / per_row
            line = "%.9g %.9g %.9g %.9g %.9g %.9g %.9g\n" % (
                depth, tas[i], sw, lw, pressure_pa / 1000.0, vp, WIND_M_S)
            fh.write(line * per_row)


def write_global(path: Path, start: datetime, n: int, model_steps: int, snow_steps: int, runoff_steps: int,
                 atmos_steps: int, daily_output: bool) -> None:
    lines = [
        "NLAYER %d" % NLAYER,
        "NODES %d" % (10 if FULL_ENERGY else 3),
        "MODEL_STEPS_PER_DAY %d" % model_steps,
        "SNOW_STEPS_PER_DAY %d" % snow_steps,
        "RUNOFF_STEPS_PER_DAY %d" % runoff_steps,
        "ATMOS_STEPS_PER_DAY %d" % atmos_steps,
        "STARTYEAR %d" % start.year, "STARTMONTH %d" % start.month, "STARTDAY %d" % start.day,
        "STARTSEC %d" % (start.hour * 3600 + start.minute * 60 + start.second),
        "NRECS %d" % n,
        "FULL_ENERGY %s" % ("TRUE" if FULL_ENERGY else "FALSE"),
        "FROZEN_SOIL FALSE",
        *(["QUICK_FLUX FALSE"] if FULL_ENERGY else []),
        "FORCING1 %s/forcing_" % WORK,
        "FORCE_FORMAT ASCII",
        "FORCE_TYPE PREC", "FORCE_TYPE AIR_TEMP", "FORCE_TYPE SWDOWN", "FORCE_TYPE LWDOWN",
        "FORCE_TYPE PRESSURE", "FORCE_TYPE VP", "FORCE_TYPE WIND",
        "FORCE_STEPS_PER_DAY %d" % snow_steps,
        "FORCEYEAR %d" % start.year, "FORCEMONTH %d" % start.month, "FORCEDAY %d" % start.day,
        "FORCESEC %d" % (start.hour * 3600 + start.minute * 60 + start.second),
        "GRID_DECIMAL 4",
        "WIND_H 10.0",
        "SOIL %s/soil.txt" % WORK,
        "BASEFLOW ARNO",
        "JULY_TAVG_SUPPLIED FALSE",
        "ORGANIC_FRACT FALSE",
        "VEGLIB %s/Stehekin_veglib.txt" % SAMPLE,
        "VEGPARAM %s/vegparam.txt" % WORK,
        "ROOT_ZONES 3",
        "VEGPARAM_LAI TRUE",
        "LAI_SRC FROM_VEGPARAM",
        "SNOW_BAND 1",
        "CONSTANTS %s/constants.txt" % WORK,
        "RESULT_DIR %s/out" % WORK,
        "OUTFILE ht",
        "AGGFREQ NSTEPS %d" % (model_steps if daily_output else 1),
        "COMPRESS FALSE",
        "OUT_FORMAT ASCII",
        "OUTVAR OUT_SOIL_MOIST %.10g * *",
    ] + ["OUTVAR %s %%.10g * * %s" % (name, "AGG_TYPE_END" if name in END_OF_ROW else "*") for name in OUTPUT_VARS]
    path.write_text("\n".join(lines) + "\n")


def read_output(path: Path, at_least: int, at_most: int) -> dict[str, list[float]]:
    names = [f"OUT_SOIL_MOIST_{k}" for k in range(NLAYER)] + OUTPUT_VARS
    columns = {name: [] for name in names}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line[0].isalpha():
            continue
        parts = line.split()
        if len(parts) < len(names) + 3 or not parts[0].isdigit():
            continue
        for name, value in zip(names, parts[-len(names):]):
            columns[name].append(float(value))
    if not at_least <= len(columns["OUT_EVAP"]) <= at_most:
        raise RuntimeError(f"VIC wrote {len(columns['OUT_EVAP'])} rows for {at_least} to {at_most}")
    return columns


def simulate(rows: list[dict], static: dict, timestep: str, substep: str | None = None):
    started = time.monotonic()
    steps_per_day = STEPS_PER_DAY[timestep]
    dt_days = 1.0 / steps_per_day
    n = len(rows)
    pr = [float(r["pr"]) for r in rows]
    tas = [float(r["tas"]) for r in rows]
    pet = [float(r["pet"]) for r in rows]
    rlds = [float(r["rlds"]) for r in rows] if "rlds" in rows[0] else None
    start = parse_time(rows[0]["time"])
    latitude = float(static.get("latitude_deg", 48.3125))
    longitude = 0.0
    elevation = float(static.get("elevation_m", 0.0))
    pressure_pa = standard_pressure_pa(elevation)
    pt_factor = pt_w_per_mm_day(PT_TREF_C, pressure_pa)

    if WORK.exists():
        shutil.rmtree(WORK)
    (WORK / "out").mkdir(parents=True)

    soil = write_soil(WORK / "soil.txt", static, latitude, longitude, tas, 365 * steps_per_day)
    veg = write_vegetation(WORK / "vegparam.txt", static)
    threshold = float(static.get("snow_threshold_degC", 0.0))
    constants = {"SNOW_MAX_SNOW_TEMP": threshold + SNOW_RAMP_HALF_WIDTH_C,
                 "SNOW_MIN_RAIN_TEMP": threshold - SNOW_RAMP_HALF_WIDTH_C, **veg["params"]}
    if "eps" in static:
        constants["EMISS_GRND"] = float(static["eps"])
        constants["EMISS_VEG"] = float(static["eps"])
        constants["EMISS_SNOW"] = float(static["eps"])
    (WORK / "constants.txt").write_text("".join("%s %.10g\n" % kv for kv in constants.items()))

    if timestep == "PT1D":
        records = DAILY_RECORDS or (STEPS_PER_DAY[substep] if substep else 0)
        if records < 4:
            raise SystemExit("VIC 5.1.0 has no daily forcing and needs a daily row held over at least "
                             "4 records; the manifest's substep says how many")
        model_steps, snow_steps, per_row = records, records, records
    else:
        model_steps, snow_steps, per_row = steps_per_day, steps_per_day, 1
    # Upstream's Stehekin setup solves runoff hourly; a finer model step takes it along.
    runoff_steps = max(24, model_steps)
    atmos_steps = max(24, model_steps)
    tag = "%.4f_%.4f" % (latitude, longitude)
    # make_dmy.c stops any run that does not end on a day boundary. A case that does not
    # gets rainless copies of its last row after the end, and they are dropped from the
    # output, so they can reach no row that is reported.
    start_seconds = start.hour * 3600 + start.minute * 60 + start.second
    # In records, so a daily row that does not start at midnight is padded in hours.
    record_s = 86400 // model_steps
    leftover = (start_seconds + n * per_row * record_s) % 86400
    padding = (86400 - leftover) // record_s if leftover else 0
    write_forcing(WORK / f"forcing_{tag}", pr, tas, pet, rlds, dt_days, per_row, pressure_pa, pt_factor)
    if padding:
        with open(WORK / f"forcing_{tag}", "a") as fh:
            sw, lw, vp = mock_atmosphere(tas[-1], pet[-1], rlds[-1] if rlds else None, pt_factor, pressure_pa)
            fh.write("0 %.9g %.9g %.9g %.9g %.9g %.9g\n" % (tas[-1], sw, lw, pressure_pa / 1000.0, vp, WIND_M_S) * padding)
    nrecs = n * per_row + padding
    write_global(WORK / "global.txt", start, nrecs, model_steps, snow_steps, runoff_steps, atmos_steps,
                 daily_output=timestep == "PT1D" and model_steps > 1)
    staged = time.monotonic()

    proc = subprocess.run([VIC_EXE, "-g", str(WORK / "global.txt")], cwd=WORK, capture_output=True, text=True)
    log = proc.stdout + proc.stderr
    (WORK / "vic.log").write_text(log)
    outputs = sorted((WORK / "out").glob("ht*"))
    if proc.returncode != 0 or not outputs:
        tail = "\n".join(log.strip().splitlines()[-25:])
        raise RuntimeError(f"VIC failed (exit {proc.returncode}):\n{tail}")
    ran = time.monotonic()
    # The padding is dropped: at a record per row it is its own rows, and on a held row a
    # part-filled row VIC may or may not write.
    out = read_output(outputs[0], n, n + (padding if per_row == 1 else padding // per_row + 1))
    out = {name: values[:n] for name, values in out.items()}

    per_day = float(steps_per_day)
    moist = [sum(out[f"OUT_SOIL_MOIST_{k}"][i] for k in range(NLAYER)) for i in range(n)]
    result = {
        "evspsbl": [v * per_day for v in out["OUT_EVAP"]],
        "mrro": [(a + b) * per_day for a, b in zip(out["OUT_RUNOFF"], out["OUT_BASEFLOW"])],
        "mrso": moist,
        "snw": list(out["OUT_SWE"]),
        "canopy": [a + b for a, b in zip(out["OUT_WDEW"], out["OUT_SNOW_CANOPY"])],
    }
    if FULL_ENERGY:
        result["hfls"] = [a + b for a, b in zip(out["OUT_LATENT"], out["OUT_LATENT_SUB"])]
        result["hfss"] = list(out["OUT_SENSIBLE"])
        # Upward fluxes in VIC's own sign. VIC's ground flux leaves its thin surface layer
        # at the first node; that layer's storage and phase change are above it, and
        # AGENTS.md puts them in the flux at the actual surface.
        result["hfg"] = [a + b + c for a, b, c in zip(out["OUT_GRND_FLUX"], out["OUT_DELTAH"], out["OUT_FUSION"])]
        # VIC reports net longwave at the surface, so upward is what net leaves of the sky.
        result["rlus"] = [a - b for a, b in zip(out["OUT_LWDOWN"], out["OUT_LWNET"])]
        result["ts"] = list(out["OUT_RAD_TEMP"])  # already in kelvin, unlike OUT_SURF_TEMP

    storage = [moist[i] + result["snw"][i] + result["canopy"][i] + out["OUT_SURFSTOR"][i] for i in range(n)]
    residuals = [
        (pr[i] - result["evspsbl"][i] - result["mrro"][i]) * dt_days - (storage[i] - storage[i - 1])
        for i in range(1, n)
    ]
    pet_total = sum(pet) * dt_days
    notes = {
        "vic": {"release": "5.1.0", "commit": "14a371a8dbee7ea834152f8ac65e2d6d54a9c30a",
                "driver": "classic", "full_energy": FULL_ENERGY, "frozen_soil": False,
                "baseflow": "ARNO", "parameter_set": PARAMETER_SET,
                "soil_and_vegetation": f"VIC_sample_data@3f2dcb4 classic/Stehekin cell {CELL_ID}"},
        "timestep": timestep,
        "rows_added_after_the_end_and_dropped": padding,
        "vic_steps": {"model": model_steps, "snow": snow_steps, "forcing_records_per_row": per_row,
                      "runoff": runoff_steps, "atmos": atmos_steps},
        "steps_per_row": per_row,
        "catchment": {**soil, **veg, "constants": constants},
        "mock_inputs": {
            "SWDOWN": f"pet x Priestley-Taylor at a fixed {PT_TREF_C} C / (1 - {ALBEDO_REF})",
            "LWDOWN": "the case's rlds" if rlds else "emission of a reference surface at air temperature, plus negative net radiation",
            "VP": f"relative humidity {RH_REF} at the row's air temperature",
            "PRESSURE": "standard atmosphere at the case's elevation",
            "WIND": f"{WIND_M_S} m/s at WIND_H = 10 m",
        },
        "pet_not_read": "VIC computes its own potential evaporation; the row's pet reaches it only through SWDOWN",
        "states": {"mrso": "OUT_SOIL_MOIST summed over the three layers (liquid and ice)",
                   "snw": "OUT_SWE (ground pack)",
                   "canopy": "OUT_WDEW + OUT_SNOW_CANOPY"},
        "fluxes": {"evspsbl": "OUT_EVAP (includes sublimation)", "mrro": "OUT_RUNOFF + OUT_BASEFLOW",
                   "hfls": "OUT_LATENT + OUT_LATENT_SUB", "hfss": "OUT_SENSIBLE",
                   "hfg": "OUT_GRND_FLUX + OUT_DELTAH + OUT_FUSION",
                   "rlus": "OUT_LWDOWN - OUT_LWNET", "ts": "OUT_RAD_TEMP"},
        "water_budget": {
            "max_abs_step_residual_mm": max((abs(v) for v in residuals), default=0.0),
            "cumulative_residual_mm": sum(residuals),
            "max_abs_vic_water_error_mm": max(abs(v) for v in out["OUT_WATER_ERROR"]),
            "input_mm": sum(pr) * dt_days,
            "max_surfstor_mm": max(out["OUT_SURFSTOR"]),
            "evspsbl_over_row_pet": sum(result["evspsbl"]) * dt_days / pet_total if pet_total > 0 else None,
            "vic_potential_over_row_pet": sum(out["OUT_PET"]) / pet_total if pet_total > 0 else None,
            "max_abs_vic_energy_error_w_m2": max((abs(v) for v in out["OUT_ENERGY_ERROR"] if abs(v) < 1e10), default=None),
        },
        "seconds": {"stage": round(staged - started, 3), "vic": round(ran - staged, 3)},
    }
    return result, notes, out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    started = time.monotonic()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent
    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        rows = list(csv.DictReader(fh))
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    timestep = str(request.get("timestep", "PT1D"))
    if timestep not in STEPS_PER_DAY:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    substep = request.get("substep")

    result, notes, raw = simulate(rows, static, timestep, substep)

    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(COLUMNS)
        for i, row in enumerate(rows):
            writer.writerow([row["time"], row["pr"]] + ["%.10g" % result[c][i] for c in COLUMNS[2:]])

    if DIAG:
        names = sorted(raw)
        with open(io_dir / "diag.csv", "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["time"] + names)
            for i, row in enumerate(rows):
                writer.writerow([row["time"]] + [raw[k][i] for k in names])
        shutil.copy(WORK / "vic.log", io_dir / "vic.log")
        shutil.copy(WORK / "global.txt", io_dir / "global.txt")

    notes["seed"] = int(request.get("seed", 0))
    notes["wall_seconds"] = round(time.monotonic() - started, 3)
    (io_dir / request["output"]["run"]).write_text(json.dumps(
        {"status": "ok", "model": MODEL, "n_steps": len(rows),
         # No vegetation tiles at all on a case with no canopy, so that store is a structural zero.
         "zero_stores": {"canopy": "the case has no canopy, so VIC gets no vegetation tiles"}
                        if float(static.get("canopy_capacity_mm", 1.0)) == 0.0 else {},
         "not_reported": NOT_REPORTED, "notes": notes}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

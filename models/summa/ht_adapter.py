#!/usr/bin/env python3
"""HydroTuring adapter for SUMMA v4.0.0 (CH-Earth/summa, GPL-3.0-or-later).

SUMMA (Clark et al. 2015, WRR) solves the coupled water and energy
conservation equations of a vegetation canopy, a layered snowpack, a layered
soil column and an aquifer with one implicit solver. This adapter runs the
compiled release as one lumped HRU, reads SUMMA's own output file and maps it
onto the /io contract. It does not change a line of SUMMA.

What SUMMA is given
-------------------
Decisions, default parameter tables (localParamInfo.txt, basinParamInfo.txt)
and Noah-MP lookup tables are the files of the one complete setup the
release ships, test_ngen/gauge_01073000 (CAMELS basin 01073000), copied into
the image unedited. Its decisions: the homegrown backward-Euler solver with
analytic derivatives, mixed-form Richards' equation with free drainage into a
big-bucket aquifer, Ball-Berry stomata, Beer's-law canopy radiation, monthly
LAI and SAI from the vegetation table, a time-delay histogram for within-GRU
routing, enthalpy-form energy equations. Its catchment attributes are kept
wherever the probe says nothing (USGS mixed forest, ROSETTA loam, 10 m
measurement height, 16 m canopy top, flat for radiation). From static.json,
only what maps onto a SUMMA quantity directly:

* latitude and area;
* `snow_threshold_degC` -> `tempCritRain`;
* `soil_capacity_mm` -> the depth of the soil column, capacity / theta_sat,
  so the most water the column can hold is the capacity; the shipped layer
  thicknesses are kept down to that depth (a remainder under half the layer
  above is merged into it) and the rooting depth is kept inside it;
* `canopy_capacity_mm` -> `refInterceptCapRain` and `refInterceptCapSnow`,
  SUMMA's capacities per unit of vegetation area, divided by the largest
  monthly LAI + SAI the table gives the class, so neither phase can be held
  above the capacity by table design; a capacity of zero is a bare surface
  (USGS class 19, which SUMMA gives no leaves or stems).

The cold state is the shipped one (283.16 K soil at 0.3 volumetric water, no
snow, no canopy water) except the aquifer, which starts empty: under the
shipped aquiferBaseflowRate the shipped 0.4 m drains within the first hours
and would put 400 mm of runoff into the spinup.

Forcing SUMMA needs that the probes do not generate is mocked from each row
alone, with no calendar and nothing from another row (README.md has what
every choice moves):

* the net radiation a row carries: the probe's `rn` where it supplies one,
  otherwise what Priestley-Taylor (alpha 1.26) needs to produce the row's
  `pet` at the row's temperature;
* shortwave and longwave: split so that a reference surface at air
  temperature, with albedo 0.23 and SUMMA's soil emissivity 0.96, would have
  exactly that net radiation and no longwave deficit. Downwelling longwave is
  the reference surface's own emission, plus the net radiation where that is
  negative; shortwave carries a positive net radiation through the albedo.
  The split is linear in the net radiation, so a day of hourly rows delivers
  the same radiant energy as the same day as one row. SUMMA then computes its
  own net radiation with its own albedo, emissivity and surface temperature,
  which is not the probe's `rn`; run.json records how far apart they are;
* humidity: 70 percent relative humidity at air temperature; pressure: a
  standard atmosphere at sea level; wind: 2 m/s at the measurement height.

Time
----
The rows go in at the step they came at, as SUMMA's data step; nothing is
resampled. SUMMA's forcing stamps are period-ending. An hour is stamped at
its end. A day is stamped at 23:00 rather than 24:00: SUMMA's clear-sky
geometry (CLRSKY_RAD, called from derivforce) starts a 24-hour window at the
stamp hour minus 24, so a midnight stamp puts the whole window in the day
before, finds no daylight in it and returns a zero cosine of the zenith
angle, and vegSWavRad then discards every watt of shortwave. Stamped at 23:00
the window holds the day's daylight. SUMMA itself reads the clock (solar
geometry; the phenology's monthly LAI and SAI); the adapter's mocks do not.

What is reported
----------------
Fluxes are step means, as rates in mm/day or W m-2; states are end of step.

* `pr`       the forcing, echoed as the text it arrived as
* `evspsbl`  total evapotranspiration including snow and canopy sublimation
             (SUMMA's scalarTotalET leaves sublimation out), positive upward;
             net deposition makes it negative
* `sbl`      snow plus canopy sublimation, positive upward; frost deposition
             is negative
* `mrro`     routed runoff (averageRoutedRunoff): surface runoff plus aquifer
             baseflow after the time-delay histogram
* `channel`  cumulative instantaneous minus routed runoff: water inside the
             histogram, which SUMMA normalises to sum to one
* `hfls`     -scalarLatHeatTotal and `hfss` -scalarSenHeatTotal (SUMMA's
             energy fluxes are positive downward)
* `hfg`      scalarGroundNetNrgFlux: the net energy flux into the snow-soil
             column through its top. With no snow that is the soil surface;
             with snow it is the snow surface, so snowpack heat storage and
             melt are inside it, as the correction rule for storage above the
             soil surface requires. The canopy's net energy flux is not in it
             (the canopy is not ground); run.json reports its size
* `mrso`     scalarTotalSoilWat (liquid plus ice) plus the water SUMMA stores
             by compressing the soil matrix, the cumulative scalarSoilCompress
             since the cold state. SUMMA's soil balance carries that term
             outside the volumetric water content; left out, it shows as a
             step residual of up to 3 mm when a dry column rewets and a drift
             of -4 mm over two rainless years
* `snw`      scalarSWE, ice plus liquid in the pack, including snow without a
             layer. scalarSfcMeltPond is not added: it records melt SUMMA has
             already passed to the soil in the same step, and adding it
             counts that water twice
* `canopy`   scalarCanopyLiq + scalarCanopyIce
* `gw`       scalarAquiferStorage, m to mm

Developer switches: SUMMA_HT_DIAG=1 writes SUMMA's raw series next to the
request (only possible where that directory is writable, never under the
harness); SUMMA_HT_EXTRA_VARS adds variables to SUMMA's output file;
SUMMA_HT_ALBEDO, SUMMA_HT_RH, SUMMA_HT_WIND, SUMMA_HT_SPLIT (neutral|clear),
SUMMA_HT_VEG and SUMMA_HT_DAILY_END_HOUR reproduce the README's sensitivity
table. None is set in the evaluated image.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np

MODEL = {"name": "summa", "version": "4.0.0-f787fa5.1"}
COLUMNS = ["time", "pr", "evspsbl", "mrro", "sbl", "hfls", "hfss", "hfg",
           "mrso", "snw", "canopy", "gw", "channel"]

SUMMA_EXE = os.environ.get("SUMMA_EXE", "/opt/summa/bin/summa.exe")
SHIPPED = Path(os.environ.get("SUMMA_SETTINGS", "/model/summa_settings"))
WORK = Path(os.environ.get("SUMMA_WORKDIR", "/tmp/summa"))
DIAG = os.environ.get("SUMMA_HT_DIAG", "0") == "1"

STEP_SECONDS = {"PT1D": 86400, "PT1H": 3600, "PT15M": 900, "PT5M": 300, "PT1M": 60}
SECONDS_PER_DAY = 86400.0
MISSING = -9000.0  # SUMMA writes -9999 where a diagnostic does not apply

# The shipped test case's catchment (attributes_tiled_by_hru.nc,
# trialParams_default_tiled_by_hru.nc, coldstate_tiled_by_hru.nc).
VEG_TYPE = int(os.environ.get("SUMMA_HT_VEG", "15"))        # USGS mixed forest
BARE_VEG_TYPE = 19                                            # USGS barren
SOIL_TYPE = 3                                                 # ROSETTA loam
MEASUREMENT_HEIGHT_M = 10.0
TAN_SLOPE = 0.1
CONTOUR_LENGTH_M = 100.0
CANOPY_TOP_M = 16.0
SHIPPED_LAYERS_M = [0.025, 0.075, 0.15, 0.25, 0.5, 0.5, 1.0, 1.5]
COLD_TEMPERATURE_K = 283.16
COLD_THETA = 0.3
COLD_MATRIC_HEAD_M = -1.0
DEFAULT_ROOTING_DEPTH_M = 2.0  # localParamInfo.txt

# Mock atmosphere, from each forcing row alone.
ALBEDO_REF = float(os.environ.get("SUMMA_HT_ALBEDO", "0.23"))
RH_REF = float(os.environ.get("SUMMA_HT_RH", "0.70"))
WIND_M_S = float(os.environ.get("SUMMA_HT_WIND", "2.0"))
SPLIT = os.environ.get("SUMMA_HT_SPLIT", "neutral")
DAILY_END_HOUR = int(os.environ.get("SUMMA_HT_DAILY_END_HOUR", "23"))
PT_ALPHA = 1.26
SIGMA = 5.670374419e-8
SOIL_EMISSIVITY = 0.96  # vegNrgFlux.f90

OUTPUT_VARS = [
    "scalarTotalET", "scalarSnowSublimation", "scalarCanopySublimation",
    "scalarSurfaceRunoff", "scalarAquiferBaseflow", "scalarRainfall", "scalarSnowfall",
    "scalarLatHeatTotal", "scalarSenHeatTotal", "scalarGroundNetNrgFlux",
    "scalarCanopyNetNrgFlux", "scalarNetRadiation", "scalarGroundAbsorbedSolar",
    "scalarCanopyAbsorbedSolar", "scalarLWNetGround", "scalarLWNetCanopy",
    "scalarGroundAdvectiveHeatFlux", "scalarCanopyAdvectiveHeatFlux",
    "scalarSWE", "scalarSfcMeltPond", "scalarCanopyLiq", "scalarCanopyIce",
    "scalarTotalSoilWat", "scalarSoilCompress", "scalarAquiferStorage", "scalarSurfaceTemp",
    "scalarCosZenith", "scalarLAI", "scalarSAI",
    "balanceCasNrg", "balanceVegNrg", "balanceSnowNrg", "balanceSoilNrg",
    "balanceVegMass", "balanceSnowMass", "balanceSoilMass", "balanceAqMass",
    "averageInstantRunoff", "averageRoutedRunoff",
]


# --- shipped tables ------------------------------------------------------------


def rosetta_soil(soil_type: int) -> dict[str, float]:
    text = (SHIPPED / "TBL_SOILPARM.TBL").read_text().splitlines()
    start = next(i for i, line in enumerate(text) if line.strip() == "ROSETTA")
    header = re.findall(r"[A-Za-z_]+", text[start + 1].split("'")[1])
    for line in text[start + 2:]:
        fields = line.split()
        if fields and fields[0] == str(soil_type):
            return dict(zip(header, (float(v) for v in fields[1:1 + len(header)])))
    raise ValueError(f"soil type {soil_type} not in the ROSETTA table")


def monthly_vai(veg_type: int) -> np.ndarray:
    """LAI + SAI by month for a USGS vegetation class (TBL_MPTABLE.TBL)."""
    block = (SHIPPED / "TBL_MPTABLE.TBL").read_text().split("MODIFIED_IGBP_MODIS_NOAH")[0]

    def table(name: str) -> np.ndarray:
        match = re.search(rf"\b{name}\s*=\s*(.*?)(?:\n\s*\n)", block, re.S)
        return np.array([float(v) for v in re.findall(r"-?\d+\.?\d*", match.group(1))]).reshape(12, -1)

    return table("LAIM")[:, veg_type - 1] + table("SAIM")[:, veg_type - 1]


# --- mock atmosphere -----------------------------------------------------------


def mock_atmosphere(tas: np.ndarray, pet: np.ndarray, rn: np.ndarray | None,
                    elevation_m: float) -> dict[str, np.ndarray]:
    """SUMMA's missing forcing, from each row alone. See the module docstring."""
    tk = tas + 273.15
    pressure = np.full_like(tas, 101325.0 * (1.0 - 2.25577e-5 * elevation_m) ** 5.25588)
    ea = RH_REF * 611.2 * np.exp(17.67 * tas / (tas + 243.5))
    spechum = 0.622 * ea / (pressure - 0.378 * ea)
    if rn is None:
        es_kpa = 0.6108 * np.exp(17.27 * tas / (tas + 237.3))
        slope = 4098.0 * es_kpa / (tas + 237.3) ** 2
        gamma = 0.000665 * pressure / 1000.0
        lam = 2.501e6 - 2361.0 * tas
        target = pet / SECONDS_PER_DAY * lam * (slope + gamma) / (PT_ALPHA * slope)
    else:
        target = rn.copy()
    emitted = SIGMA * tk ** 4
    if SPLIT == "clear":
        # Sensitivity only: Brutsaert (1975) clear-sky longwave where the
        # target is positive. Not linear in the target, so not step-consistent.
        lw_clear = 1.24 * (ea / 100.0 / tk) ** (1.0 / 7.0) * emitted
        positive = target > 0.0
        sw = np.where(positive, (target - SOIL_EMISSIVITY * (lw_clear - emitted)) / (1.0 - ALBEDO_REF), 0.0)
        lw = np.where(positive, lw_clear, emitted + target / SOIL_EMISSIVITY)
    else:
        sw = np.maximum(target, 0.0) / (1.0 - ALBEDO_REF)
        lw = emitted + np.minimum(target, 0.0) / SOIL_EMISSIVITY
    return {"SWRadAtm": np.maximum(sw, 0.0), "LWRadAtm": np.maximum(lw, 0.0), "airpres": pressure,
            "spechum": spechum, "windspd": np.full_like(tas, WIND_M_S), "target": target}


# --- SUMMA files -----------------------------------------------------------------


def parse_time(text: str) -> datetime:
    return datetime.fromisoformat(text.strip().replace(" ", "T"))


def soil_layers(depth_m: float) -> list[float]:
    """The shipped layer thicknesses down to `depth_m`; a thin remainder is merged up."""
    layers, top = [], 0.0
    for d in SHIPPED_LAYERS_M:
        if top + d >= depth_m - 1e-9:
            layers.append(depth_m - top)
            break
        layers.append(d)
        top += d
    else:
        layers[-1] += depth_m - top
    if len(layers) > 1 and layers[-1] < 0.5 * layers[-2]:
        remainder = layers.pop()
        layers[-1] += remainder
    return layers


def netcdf_vars(nc, spec):
    for name, dtype, dims, value in spec:
        var = nc.createVariable(name, dtype, dims)
        var[:] = value


def write_forcing(path: Path, stamps, reference, step_s, pr, tas, atmos) -> None:
    with netCDF4.Dataset(path, "w") as nc:
        nc.createDimension("time", len(stamps))
        nc.createDimension("hru", 1)
        t = nc.createVariable("time", "f8", ("time",))
        t.units = f"seconds since {reference:%Y-%m-%d %H:%M}"
        t.calendar = "standard"
        t[:] = [(s - reference).total_seconds() for s in stamps]
        nc.createVariable("data_step", "f8").assignValue(float(step_s))
        series = {"pptrate": pr / SECONDS_PER_DAY, "airtemp": tas + 273.15,
                  **{k: atmos[k] for k in ("SWRadAtm", "LWRadAtm", "windspd", "airpres", "spechum")}}
        netcdf_vars(nc, [("hruId", "i4", ("hru",), 1)] +
                    [(k, "f8", ("time", "hru"), v[:, None]) for k, v in series.items()])


def write_attributes(path: Path, static: dict, veg_type: int) -> None:
    with netCDF4.Dataset(path, "w") as nc:
        nc.createDimension("hru", 1)
        nc.createDimension("gru", 1)
        netcdf_vars(nc, [
            ("hruId", "i4", ("hru",), 1), ("gruId", "i4", ("gru",), 1),
            ("hru2gruId", "i4", ("hru",), 1), ("downHRUindex", "i4", ("hru",), 0),
            ("longitude", "f8", ("hru",), 0.0),
            ("latitude", "f8", ("hru",), float(static.get("latitude_deg", 40.0))),
            ("elevation", "f8", ("hru",), float(static.get("elevation_m", 0.0))),
            ("HRUarea", "f8", ("hru",), float(static["area_km2"]) * 1.0e6),
            ("tan_slope", "f8", ("hru",), TAN_SLOPE),
            ("contourLength", "f8", ("hru",), CONTOUR_LENGTH_M),
            ("slopeTypeIndex", "i4", ("hru",), 1), ("soilTypeIndex", "i4", ("hru",), SOIL_TYPE),
            ("vegTypeIndex", "i4", ("hru",), veg_type),
            ("mHeight", "f8", ("hru",), MEASUREMENT_HEIGHT_M),
        ])


def write_trial_params(path: Path, params: dict[str, float]) -> None:
    with netCDF4.Dataset(path, "w") as nc:
        nc.createDimension("hru", 1)
        netcdf_vars(nc, [("hruId", "i4", ("hru",), 1)] +
                    [(k, "f8", ("hru",), v) for k, v in params.items()])


def write_cold_state(path: Path, layers: list[float], step_s: int) -> None:
    n = len(layers)
    with netCDF4.Dataset(path, "w") as nc:
        for dim, size in (("hru", 1), ("scalarv", 1), ("midSoil", n), ("midToto", n), ("ifcToto", n + 1)):
            nc.createDimension(dim, size)
        scalar = ("scalarv", "hru")
        netcdf_vars(nc, [
            ("hruId", "i4", ("hru",), 1),
            ("dt_init", "f8", scalar, float(min(3600, step_s))),
            ("nSoil", "i4", scalar, n), ("nSnow", "i4", scalar, 0),
            *[(k, "f8", scalar, 0.0) for k in ("scalarCanopyIce", "scalarCanopyLiq", "scalarSnowDepth",
                                               "scalarSWE", "scalarSfcMeltPond", "scalarAquiferStorage",
                                               "scalarSnowAlbedo")],
            ("scalarCanairTemp", "f8", scalar, COLD_TEMPERATURE_K),
            ("scalarCanopyTemp", "f8", scalar, COLD_TEMPERATURE_K),
            ("mLayerTemp", "f8", ("midToto", "hru"), np.full((n, 1), COLD_TEMPERATURE_K)),
            ("mLayerVolFracIce", "f8", ("midToto", "hru"), np.zeros((n, 1))),
            ("mLayerVolFracLiq", "f8", ("midToto", "hru"), np.full((n, 1), COLD_THETA)),
            ("mLayerMatricHead", "f8", ("midSoil", "hru"), np.full((n, 1), COLD_MATRIC_HEAD_M)),
            ("iLayerHeight", "f8", ("ifcToto", "hru"), np.concatenate(([0.0], np.cumsum(layers)))[:, None]),
            ("mLayerDepth", "f8", ("midToto", "hru"), np.array(layers)[:, None]),
        ])


def write_file_manager(path: Path, start: datetime, end: datetime) -> None:
    entries = [
        ("controlVersion", "SUMMA_FILE_MANAGER_V3.0.0"),
        ("simStartTime", f"{start:%Y-%m-%d %H:%M}"), ("simEndTime", f"{end:%Y-%m-%d %H:%M}"),
        ("tmZoneInfo", "localTime"), ("outFilePrefix", "ht"),
        ("settingsPath", f"{WORK}/settings/"), ("forcingPath", f"{WORK}/forcing/"),
        ("outputPath", f"{WORK}/output/"),
        ("initConditionFile", "coldState.nc"), ("attributeFile", "attributes.nc"),
        ("trialParamFile", "trialParams.nc"), ("forcingListFile", "forcingFileList.txt"),
        ("decisionsFile", "modelDecisions.txt"), ("outputControlFile", "outputControl.txt"),
        ("globalHruParamFile", "localParamInfo.txt"), ("globalGruParamFile", "basinParamInfo.txt"),
        ("vegTableFile", "TBL_VEGPARM.TBL"), ("soilTableFile", "TBL_SOILPARM.TBL"),
        ("generalTableFile", "TBL_GENPARM.TBL"), ("noahmpTableFile", "TBL_MPTABLE.TBL"),
    ]
    path.write_text("".join(f"{k:<20} '{v}'\n" for k, v in entries))


# --- run -------------------------------------------------------------------------


def simulate(rows: list[dict], columns: list[str], static: dict, timestep: str):
    started = time.monotonic()
    step_s = STEP_SECONDS[timestep]
    dt_days = step_s / SECONDS_PER_DAY
    pr = np.array([float(r["pr"]) for r in rows])
    tas = np.array([float(r["tas"]) for r in rows])
    pet = np.array([float(r["pet"]) for r in rows])
    rn = np.array([float(r["rn"]) for r in rows]) if "rn" in columns else None

    starts = [parse_time(r["time"]) for r in rows]
    shift = timedelta(hours=DAILY_END_HOUR) if step_s == 86400 else timedelta(seconds=step_s)
    stamps = [s + shift for s in starts]
    # The reference is the first day of the record, so the numbers on the time
    # axis do not depend on the calendar year.
    reference = datetime(starts[0].year, starts[0].month, starts[0].day)

    soil = rosetta_soil(SOIL_TYPE)
    depth_m = float(static["soil_capacity_mm"]) / 1000.0 / soil["theta_sat"]
    layers = soil_layers(depth_m)
    canopy_capacity = float(static.get("canopy_capacity_mm", 0.0))
    veg_type = VEG_TYPE if canopy_capacity > 0.0 else BARE_VEG_TYPE
    vai_max = float(monthly_vai(veg_type).max())
    params = {
        "tempCritRain": float(static.get("snow_threshold_degC", 0.0)) + 273.15,
        "heightCanopyTop": CANOPY_TOP_M,
        "rootingDepth": min(DEFAULT_ROOTING_DEPTH_M, depth_m),
    }
    if canopy_capacity > 0.0 and vai_max > 0.0:
        params["refInterceptCapRain"] = canopy_capacity / vai_max
        params["refInterceptCapSnow"] = canopy_capacity / vai_max
    atmos = mock_atmosphere(tas, pet, rn, float(static.get("elevation_m", 0.0)))

    if WORK.exists():
        shutil.rmtree(WORK)
    for sub in ("settings", "forcing", "output"):
        (WORK / sub).mkdir(parents=True)
    for name in ("TBL_VEGPARM.TBL", "TBL_SOILPARM.TBL", "TBL_GENPARM.TBL", "TBL_MPTABLE.TBL",
                 "localParamInfo.txt", "basinParamInfo.txt", "modelDecisions.txt"):
        shutil.copy(SHIPPED / name, WORK / "settings" / name)
    extra = [v for v in os.environ.get("SUMMA_HT_EXTRA_VARS", "").split(",") if v]
    (WORK / "settings" / "outputControl.txt").write_text(
        "outputPrecision | double\n" + "".join(f"{v} | 1\n" for v in OUTPUT_VARS + extra))
    (WORK / "settings" / "forcingFileList.txt").write_text("'forcing.nc'\n")
    write_forcing(WORK / "forcing" / "forcing.nc", stamps, reference, step_s, pr, tas, atmos)
    write_attributes(WORK / "settings" / "attributes.nc", static, veg_type)
    write_trial_params(WORK / "settings" / "trialParams.nc", params)
    write_cold_state(WORK / "settings" / "coldState.nc", layers, step_s)
    write_file_manager(WORK / "settings" / "fileManager.txt", stamps[0], stamps[-1])
    staged = time.monotonic()

    proc = subprocess.run(
        [SUMMA_EXE, "-m", str(WORK / "settings" / "fileManager.txt"), "-p", "never", "-r", "never"],
        cwd=WORK, capture_output=True, text=True)
    log = proc.stdout + proc.stderr
    (WORK / "summa.log").write_text(log)
    outputs = sorted((WORK / "output").glob("*timestep*.nc"))
    if proc.returncode != 0 or not outputs or "FATAL" in log.upper():
        tail = "\n".join(log.strip().splitlines()[-25:])
        raise RuntimeError(f"SUMMA failed (exit {proc.returncode}):\n{tail}")
    ran = time.monotonic()

    out: dict[str, np.ndarray] = {}
    with netCDF4.Dataset(outputs[0]) as nc:
        for name in OUTPUT_VARS:
            data = np.ma.filled(nc.variables[name][:].astype(float), np.nan)
            out[name] = data.reshape(data.shape[0], -1)[:, 0]
    if len(out["scalarTotalET"]) != len(rows):
        raise RuntimeError(f"SUMMA wrote {len(out['scalarTotalET'])} steps for {len(rows)} rows")

    per_day = SECONDS_PER_DAY
    sublimation = -(out["scalarSnowSublimation"] + out["scalarCanopySublimation"]) * per_day
    instant = out["averageInstantRunoff"] * 1000.0 * per_day
    routed = out["averageRoutedRunoff"] * 1000.0 * per_day
    result = {
        "evspsbl": -out["scalarTotalET"] * per_day + sublimation,
        "mrro": routed,
        "sbl": sublimation,
        "hfls": -out["scalarLatHeatTotal"],
        "hfss": -out["scalarSenHeatTotal"],
        "hfg": out["scalarGroundNetNrgFlux"],
        # Liquid and ice, plus the elastic storage SUMMA's soil balance keeps
        # outside them (kg m-2 s-1 over the step, accumulated from the cold state).
        "mrso": out["scalarTotalSoilWat"] + np.cumsum(out["scalarSoilCompress"] * step_s),
        "snw": out["scalarSWE"],
        "canopy": out["scalarCanopyLiq"] + out["scalarCanopyIce"],
        "gw": out["scalarAquiferStorage"] * 1000.0,
        "channel": np.cumsum((instant - routed) * dt_days),
    }
    for name, values in result.items():
        if not np.isfinite(values).all():
            raise RuntimeError(f"SUMMA produced non-finite {name} on {int((~np.isfinite(values)).sum())} steps")

    notes = describe(result, out, atmos, pr, dt_days, timestep, veg_type, layers, depth_m, params,
                     rn is not None)
    notes["seconds"] = {"stage": round(staged - started, 2), "summa": round(ran - staged, 2)}
    diag = {**out, "rn_target": atmos["target"], "SWRadAtm": atmos["SWRadAtm"],
            "LWRadAtm": atmos["LWRadAtm"], "spechum": atmos["spechum"]}
    return result, notes, diag


def describe(result, out, atmos, pr, dt_days, timestep, veg_type, layers, depth_m, params, has_rn):
    """What run.json says about a case: the mapping, and the budgets as SUMMA kept them."""
    storage = sum(result[k] for k in ("mrso", "snw", "canopy", "gw", "channel"))
    fluxes = (pr - result["evspsbl"] - result["mrro"]) * dt_days
    step_residual = fluxes[1:] - np.diff(storage)

    surface = result["hfls"] + result["hfss"] + result["hfg"]
    target, summa_rn = atmos["target"], out["scalarNetRadiation"]
    advective = out["scalarGroundAdvectiveHeatFlux"] + out["scalarCanopyAdvectiveHeatFlux"]
    balance = {}
    for name in ("balanceCasNrg", "balanceVegNrg", "balanceSnowNrg", "balanceSoilNrg",
                 "balanceVegMass", "balanceSnowMass", "balanceSoilMass", "balanceAqMass"):
        values = out[name][np.isfinite(out[name]) & (np.abs(out[name]) < -MISSING)]
        balance[name] = float(np.abs(values).max()) if values.size else 0.0

    return {
        "summa": {"release": "v4.0.0", "commit": "f787fa5e63d3c67030721a85f2f244d22a8d691e",
                  "solver": "homegrown backward Euler (num_method homegrown, fDerivMeth analytic)",
                  "setup": "test_ngen/gauge_01073000/settings/SUMMA: decisions, localParamInfo, "
                           "basinParamInfo and Noah-MP tables unedited"},
        "timestep": timestep,
        "forcing_stamps": ("period-ending, a day stamped at 23:00 (a 24:00 stamp makes SUMMA's "
                           "CLRSKY_RAD return a zero zenith cosine and discard all shortwave)"
                           if dt_days == 1.0 else "period-ending"),
        "catchment": {
            "vegetation_type_usgs": veg_type, "soil_type_rosetta": SOIL_TYPE,
            "soil_depth_m": round(depth_m, 4), "soil_layers_m": [round(x, 4) for x in layers],
            "trial_parameters": {k: round(v, 6) for k, v in params.items()},
            "measurement_height_m": MEASUREMENT_HEIGHT_M,
        },
        "mock_inputs": {
            "net_radiation": "the probe's rn" if has_rn else
                             f"Priestley-Taylor (alpha {PT_ALPHA}) inverted from pet",
            "radiation_split": (f"{SPLIT}: reference surface at air temperature, albedo {ALBEDO_REF}, "
                                f"emissivity {SOIL_EMISSIVITY}"),
            "SWRadAtm": "positive net radiation / (1 - albedo)" if SPLIT != "clear" else
                        "positive net radiation minus the clear-sky longwave deficit, / (1 - albedo)",
            "LWRadAtm": "reference emission plus negative net radiation / emissivity" if SPLIT != "clear" else
                        "Brutsaert clear sky where net radiation is positive",
            "spechum": f"relative humidity {RH_REF}",
            "airpres": "standard atmosphere at sea level",
            "windspd": f"{WIND_M_S} m/s at the measurement height",
        },
        "states": {
            "mrso": "scalarTotalSoilWat (liquid + ice) + cumulative scalarSoilCompress (elastic storage "
                    "since the cold state)",
            "snw": "scalarSWE",
            "canopy": "scalarCanopyLiq + scalarCanopyIce",
            "gw": "scalarAquiferStorage (big bucket), m to mm",
            "channel": "cumulative averageInstantRunoff - averageRoutedRunoff (time-delay histogram)",
        },
        "fluxes": {
            "evspsbl": "-(scalarTotalET + scalarSnowSublimation + scalarCanopySublimation)",
            "sbl": "-(scalarSnowSublimation + scalarCanopySublimation); negative is deposition",
            "mrro": "averageRoutedRunoff", "hfls": "-scalarLatHeatTotal", "hfss": "-scalarSenHeatTotal",
            "hfg": "scalarGroundNetNrgFlux (top of the snow-soil column; canopy energy not included)",
        },
        "diagnostics_whole_record": {
            "water_max_abs_step_residual_mm": float(np.abs(step_residual).max()) if step_residual.size else 0.0,
            "water_cumulative_residual_mm": float(step_residual.sum()),
            "water_input_mm": float((pr * dt_days).sum()),
            "soil_elastic_storage_change_mm": float(np.sum(out["scalarSoilCompress"]) * dt_days * SECONDS_PER_DAY),
            "mean_target_net_radiation_w_m2": float(target.mean()),
            "mean_summa_net_radiation_w_m2": float(summa_rn.mean()),
            "mean_abs_target_minus_summa_net_radiation_w_m2": float(np.abs(target - summa_rn).mean()),
            "mean_hfls_hfss_hfg_w_m2": float(surface.mean()),
            "mean_summa_net_radiation_minus_hfls_hfss_hfg_w_m2": float((summa_rn - surface).mean()),
            "mean_advective_heat_w_m2": float(advective.mean()),
            "mean_canopy_net_energy_w_m2": float(out["scalarCanopyNetNrgFlux"].mean()),
            "summa_balance_max_abs": balance,
            "deposition_steps": int((result["sbl"] < 0.0).sum()),
            "negative_evspsbl_steps": int((result["evspsbl"] < 0.0).sum()),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    started = time.monotonic()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent
    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        reader = csv.DictReader(fh)
        rows, columns = list(reader), list(reader.fieldnames or [])
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    timestep = str(request.get("timestep", "PT1D"))
    if timestep not in STEP_SECONDS:
        raise SystemExit(f"unsupported timestep {timestep!r}")

    result, notes, diag = simulate(rows, columns, static, timestep)

    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(COLUMNS)
        for i, row in enumerate(rows):
            writer.writerow([row["time"], row["pr"]] + ["%.12g" % result[c][i] for c in COLUMNS[2:]])

    if DIAG:
        names = sorted(diag)
        with open(io_dir / "diag.csv", "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["time"] + names)
            for i, row in enumerate(rows):
                writer.writerow([row["time"]] + ["%.12g" % diag[n][i] for n in names])
        shutil.copy(WORK / "summa.log", io_dir / "summa.log")
        for nc_file in (WORK / "output").glob("*timestep*.nc"):
            shutil.copy(nc_file, io_dir / "summa_timestep.nc")

    notes["seed"] = int(request.get("seed", 0))
    notes["wall_seconds"] = round(time.monotonic() - started, 2)
    (io_dir / request["output"]["run"]).write_text(json.dumps(
        {"status": "ok", "model": MODEL, "n_steps": len(rows), "notes": notes}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

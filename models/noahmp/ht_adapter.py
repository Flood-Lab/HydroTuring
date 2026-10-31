#!/usr/bin/env python3
"""HydroTuring adapter for Noah-MP v5.2.1 (NCAR/noahmp, UCAR licence).

Noah-MP (Niu et al. 2011, JGR) is the land scheme of the National Water Model,
WRF and NASA LIS. This adapter runs upstream's unmodified Fortran as one column
through ht_noahmp_driver, which does what HRLDAS does for a land point with the
netCDF I/O removed; verify_driver.py compares the two. One Noah-MP step per row,
except on a row the request gives a `substep` for: a daily row holds its values
flat for 24 hourly Noah-MP steps, because one 86400 s step stops upstream's own
checks (README.md). Fluxes are then the mean over the row's steps; the stores,
`tsoil_layer`, `rlus` and `ts` are the last step's. Nothing is read from a later row.

What Noah-MP is given
---------------------
Physics options: upstream's namelist.hrldas_example, unedited except for the
step, the layers and the file paths. Parameters: upstream's NoahmpTable.TBL,
unedited, with the MODIS land classes. Land point: upstream's own single-point
example (Bondville: MODIS class 12 croplands, STAS class 8 silty clay loam,
96 percent maximum vegetation fraction, its initial soil moisture).

From static.json, only what maps onto a Noah-MP quantity:

* `soil_capacity_mm`: the column depth, capacity / porosity, so the most water
  it can hold is the capacity. The first layer keeps upstream's 0.1 m and the
  rest of the depth is shared in upstream's ratio 0.3 : 0.6 : 1.0;
* `canopy_capacity_mm`: CH2OP, the liquid capacity per unit LAI + SAI, divided
  by the class's largest monthly LAI + SAI. Noah-MP's snow interception capacity
  is a formula with no parameter and is not touched. Zero is a bare surface:
  MODIS class 16, no vegetation fraction;
* `eps`: EG, the soil-surface emissivity;
* `soil_layer_depth_m`, `soil_heat_capacity_areal`, `soil_temperature_initial`,
  `soil_water_content_initial`: the first layer's thickness, CSOIL (the solid
  heat capacity that gives this areal capacity at this water content, through
  Noah-MP's own mixture rule), and the initial column;
* `latitude_deg`, for the hemisphere of the phenology; `elevation_m` for the
  pressure of the standard atmosphere.

Forcing the probes do not generate is mocked from each row alone, with no
calendar and nothing from another row. `rsds`, `rlds`, `sfcWind`, `huss` and `ps`
go in directly wherever a probe supplies them. Otherwise:

* the net radiation a row carries: the probe's `rn` where it supplies one,
  otherwise `pet` times the Priestley-Taylor (alpha 1.26) conversion at a fixed
  20 C, as models/summa does;
* shortwave and longwave: split so that a reference surface at air temperature,
  albedo 0.23 and Noah-MP's soil emissivity, would have exactly that net
  radiation. Noah-MP then computes its own net radiation from its own albedo,
  emissivity and surface temperature. That is NOT the probe's `rn`, so `rn` is
  not declared in model.yaml: the adapter uses it to build a sky and nothing
  holds Noah-MP's budget to it. run.json reports both and how far apart they are;
* humidity 70 percent, wind 2 m/s, standard atmosphere;
* the cosine of the zenith angle, which HRLDAS computes from the clock and
  Noah-MP uses for albedo and to discard shortwave at night: 0.5 where the row
  has shortwave, 0 where it has none, so the row's radiation reaches the column
  whatever the hour and the case's calendar;
* the deep soil temperature at 8 m, and the initial soil temperature: the mean
  air temperature of the first year of the record.

Time
----
A row labelled t forces the step ending at t + dt, as HRLDAS's forcing stamps
do. The calendar reaches Noah-MP only through the phenology's monthly LAI and SAI
for that end time; that is the model reading the clock, not a mock.

What is reported
----------------
Fluxes are step means, as rates in mm/day or W m-2; states are end of step.

* `pr`          the forcing, echoed as the text it arrived as
* `evspsbl`     EvapCanopyNet + Transpiration + EvapGroundNet, net of dew and frost
* `mrro`        RunoffSurface + RunoffSubsurface (mm per soil step)
* `sbl`         SublimSnowSfcIce + SublimCanopyIce, Noah-MP's gross sublimation
* `snm`         SnowBotOutflow, the thin-snow pond terms, and rain with no snow layer
* `hfls`        HeatLatentCanopy + HeatLatentGrd + HeatLatentTransp
* `hfss`        HeatSensibleSfc
* `hfg`         HeatGroundTot: into the top of the snow-soil column, so under
                snow it is the snow surface. Noah-MP's own budget also carries
                HeatPrecipAdvSfc and HeatCanStorageChg, which have no column in
                the contract; run.json reports them
* `rlus`        RadLwEmitSfc, the downward longwave plus the net longwave
* `ts`          TemperatureRadSfc, which Noah-MP inverts from RadLwEmitSfc with
                its own emissivity and the reflected sky
* `hfg_bottom`  the Fourier flux through the base of the first soil layer, from
                the conductivity and the end-of-step temperatures of Noah-MP's
                implicit solve (SoilSnowThermalDiffusion)
* `tsoil_layer` TemperatureSoilSnow(1), the first soil layer, end of step
* `mrso`        the column's SoilMoisture (liquid and ice) over layer thickness
* `snw`         SnowWaterEquiv
* `canopy`      CanopyLiqWater + CanopyIce

Under the free-drainage runoff option upstream's example selects, the aquifer
WaterStorageAquifer is a constant 4900 mm that nothing moves: it is not reported
(`not_reported`). There is no routing in Noah-MP, so there is no `channel` or `dis`.

Developer switches: NOAHMP_HT_DIAG=1 writes every driver column next to the
request (only possible where that directory is writable, never under the
harness); NOAHMP_HT_OPTIONS="NAME=VALUE,..." overrides namelist options;
NOAHMP_HT_VEG, NOAHMP_HT_SOIL, NOAHMP_HT_ALBEDO, NOAHMP_HT_RH, NOAHMP_HT_WIND,
NOAHMP_HT_COSZ reproduce the README's sensitivity table. None is set in the
evaluated image.
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
import tempfile
import time
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

MODEL = {"name": "noahmp", "version": "5.2.1-17751dc.5"}
COLUMNS = ["time", "pr", "evspsbl", "mrro", "sbl", "snm", "hfls", "hfss", "hfg", "rlus", "hfg_bottom",
           "mrso", "snw", "canopy", "ts", "tsoil_layer"]

DRIVER = os.environ.get("NOAHMP_DRIVER", "/opt/noahmp/ht_noahmp_driver")
SHIPPED = Path(os.environ.get("NOAHMP_SETTINGS", "/opt/noahmp"))
DIAG = os.environ.get("NOAHMP_HT_DIAG", "0") == "1"

# In the order read_column writes them.
OUTPUT_NAMES = [
    "RadSwDownRefHeight", "RadLwDownRefHeight", "RadSwAbsSfc", "RadSwReflSfc", "RadLwNetSfc", "RadLwEmitSfc",
    "TemperatureRadSfc", "EmissivitySfc", "AlbedoSfc", "TemperatureGrd", "TemperatureCanopy", "HeatSensibleSfc",
    "HeatLatentCanopy", "HeatLatentGrd", "HeatLatentTransp", "HeatLatentIrriEvap", "HeatGroundTot",
    "HeatPrecipAdvSfc", "HeatCanStorageChg", "EnergyBalanceError", "RadSwPenetrateGrd1", "EvapCanopyNet",
    "Transpiration", "EvapGroundNet", "RunoffSurface", "RunoffSubsurface", "WaterBalanceError", "SnowBotOutflow",
    "RainfallGround", "SnowfallGround", "SublimSnowSfcIce", "FrostSnowSfcIce", "SublimCanopyIce",
    "PondSfcThinSnwMelt", "PondSfcThinSnwComb", "PondSfcThinSnwTrans", "SoilWater", "SnowWaterEquiv", "SnowDepth",
    "CanopyLiqWater", "CanopyIce", "WaterStorageAquifer", "NumSnowLayerNeg", "VegFrac", "VegAreaIndEff",
    "TemperatureSoil1", "HeatSoil1Bottom", "HeatCapacSoil1", "RainWithoutPack",
]

# Eight significant digits, about what single precision resolves, and single precision is
# all Noah-MP computes. Twelve wrote 9.85 MB and nine 8.12 MB on the 40 days of minute rows
# of mass/resolution-invariance, whose limit is 8 MB. The rounding is about 1e-5 mm on mrso
# against 7e-4 mm moved per minute row at 1 mm/day, and under 1e-3 W/m2 on the soil heat
# budget against a 1 W/m2 tolerance. A shorter window would fit too, but would score less
# of the month that probe integrates over.
OUTPUT_FORMAT = "%.8g"

STEP_SECONDS = {"PT1D": 86400, "PT1H": 3600, "PT15M": 900, "PT5M": 300, "PT1M": 60}
SECONDS_PER_DAY = 86400.0
SIGMA = 5.670374419e-8

# Noah-MP has a snowpack and a canopy, so no store is an invented zero.
ZERO_STORES: dict[str, str] = {}

NOT_REPORTED = {
    "gw": "WaterStorageAquifer is a constant 4900 mm that free drainage never moves",
    "channel": "Noah-MP has no routing",
    "dis": "Noah-MP has no routing",
    "gwex": "nothing crosses the boundary except drainage, which is in mrro",
}

# Upstream's single-point example, hrldas/HRLDAS_forcing/run/examples/single_point/bondville.dat.
VEG_CROPLAND = int(os.environ.get("NOAHMP_HT_VEG", "12"))  # MODIS croplands
VEG_BARREN = 16                                             # MODIS barren
SOIL_CLASS = int(os.environ.get("NOAHMP_HT_SOIL", "8"))     # STAS silty clay loam
MAX_VEG_PCT = 96.0
MIN_VEG_PCT = 1.0
LAI_INITIAL = 2.0
SOIL_MOISTURE_INITIAL = [0.298, 0.294, 0.271, 0.307]
LATITUDE_DEFAULT = 40.01
LONGITUDE = -88.37

# Upstream's namelist.hrldas_example.
LAYERS_M = [0.10, 0.30, 0.60, 1.00]

# A saturated column holds porosity x depth in single precision, printed to 8 significant
# digits: at 320 mm that is 320.00001, one printed digit past the capacity and ten times the
# bound's 1e-6 mm. So the depth that sets the capacity sits this fraction under it, 3e-4 mm
# at 320 mm. The reported water is never clipped.
CAPACITY_GUARD = 1.0e-6

# Mock atmosphere, from each forcing row alone.
ALBEDO_REF = float(os.environ.get("NOAHMP_HT_ALBEDO", "0.23"))
RH_REF = float(os.environ.get("NOAHMP_HT_RH", "0.70"))
WIND_M_S = float(os.environ.get("NOAHMP_HT_WIND", "2.0"))
COSZ_DAY = float(os.environ.get("NOAHMP_HT_COSZ", "0.5"))
PT_ALPHA = 1.26
PT_TREF_C = 20.0

# ConstantDefineMod.F90, for inverting Noah-MP's mixture rule (SoilThermalPropertyMod.F90).
HEAT_CAPACITY_WATER = 4.188e6
HEAT_CAPACITY_AIR = 1004.64


# --- shipped tables ------------------------------------------------------------


def table_block(name: str) -> str:
    text = (SHIPPED / "NoahmpTable.TBL").read_text()
    start = text.index(f"&{name}")
    return text[start:text.index("\n/", start)]


def table_values(block: str, key: str) -> list[float]:
    match = re.search(rf"^\s*{key}\s*=\s*(.*?)(?:!.*)?$", block, re.M)
    if match is None:
        raise SystemExit(f"NoahmpTable.TBL has no {key}")
    return [float(v) for v in re.findall(r"-?\d+\.?\d*(?:[Ee][+-]?\d+)?", match.group(1))]


def monthly_vai_max(veg: int) -> float:
    """The largest monthly LAI + SAI upstream's table gives a MODIS class."""
    block = table_block("noahmp_modis_parameters")
    months = ("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
    return max(table_values(block, f"LAI_{m}")[veg - 1] + table_values(block, f"SAI_{m}")[veg - 1] for m in months)


def table_scalar_eg() -> float:
    return table_values(table_block("noahmp_rad_parameters"), "EG")[0]


def porosity(soil: int) -> float:
    return table_values(table_block("noahmp_soil_stas_parameters"), "MAXSMC")[soil - 1]


# --- mock atmosphere -----------------------------------------------------------


def standard_pressure(elevation_m: float) -> float:
    return 101325.0 * (1.0 - 2.25577e-5 * elevation_m) ** 5.25588


def mock_spechum(tas: np.ndarray, pressure: np.ndarray) -> np.ndarray:
    ea = RH_REF * 611.2 * np.exp(17.67 * tas / (tas + 243.5))
    return 0.622 * ea / (pressure - 0.378 * ea)


def priestley_taylor_w_per_mm_day(tc: float, pressure: np.ndarray) -> np.ndarray:
    """Net radiation (W m-2) that Priestley-Taylor turns into 1 mm/day of evaporation."""
    es_kpa = 0.6108 * np.exp(17.27 * tc / (tc + 237.3))
    slope = 4098.0 * es_kpa / (tc + 237.3) ** 2
    gamma = 0.000665 * pressure / 1000.0
    lam = 2.501e6 - 2361.0 * tc
    return lam * (slope + gamma) / (PT_ALPHA * slope) / SECONDS_PER_DAY


def atmosphere(rows: list[dict], columns: list[str], elevation_m: float, emissivity: float) -> dict:
    """Noah-MP's forcing, one array per input, and what each came from. See the module docstring."""

    def column(name: str) -> np.ndarray | None:
        return np.array([float(r[name]) for r in rows]) if name in columns else None

    tas = column("tas")
    tk = tas + 273.15
    supplied_ps = column("ps")
    pressure = supplied_ps if supplied_ps is not None else np.full_like(tas, standard_pressure(elevation_m))
    supplied_q = column("huss")
    spechum = supplied_q if supplied_q is not None else mock_spechum(tas, pressure)
    supplied_wind = column("sfcWind")
    wind = supplied_wind if supplied_wind is not None else np.full_like(tas, WIND_M_S)
    rn, pet, rsds, rlds = column("rn"), column("pet"), column("rsds"), column("rlds")

    target = rn
    if target is None and pet is not None and (rsds is None or rlds is None):
        target = pet * priestley_taylor_w_per_mm_day(PT_TREF_C, pressure)
    emitted = SIGMA * tk ** 4
    if rsds is not None and rlds is not None:
        sw, lw, split = rsds, rlds, "the case's rsds and rlds, directly"
    elif rlds is not None:
        lw = rlds
        sw = np.maximum(target + emissivity * (emitted - rlds), 0.0) / (1.0 - ALBEDO_REF)
        split = (f"the case's rlds; shortwave closes a reference surface at air temperature, albedo {ALBEDO_REF}, "
                 f"emissivity {emissivity}, on the target net radiation")
    elif rsds is not None:
        sw = rsds
        lw = np.maximum(emitted + (target - (1.0 - ALBEDO_REF) * rsds) / emissivity, 0.0)
        split = (f"the case's rsds; longwave closes a reference surface at air temperature, albedo {ALBEDO_REF}, "
                 f"emissivity {emissivity}, on the target net radiation")
    else:
        sw = np.maximum(target, 0.0) / (1.0 - ALBEDO_REF)
        lw = np.maximum(emitted + np.minimum(target, 0.0) / emissivity, 0.0)
        split = (f"neutral: a reference surface at air temperature, albedo {ALBEDO_REF}, emissivity {emissivity} "
                 "would have exactly the target net radiation and no longwave deficit")
    if rn is not None and rsds is not None and rlds is not None:
        target_label = "the probe's rn, for comparison only: the case's rsds and rlds drive the column"
    elif rn is not None:
        target_label = "the probe's rn, used only to build the sky; rn is not declared and Noah-MP computes its own"
    elif rsds is not None and rlds is not None:
        target_label = "none: the case supplies rsds and rlds"
    else:
        target_label = f"pet times the Priestley-Taylor (alpha {PT_ALPHA}) conversion at a fixed {PT_TREF_C} C"
    labels = {
        "net_radiation_target": target_label,
        "radiation_split": split,
        "spechum": "the case's huss" if supplied_q is not None else f"relative humidity {RH_REF} at the row's air temperature",
        "pressure": "the case's ps" if supplied_ps is not None else "standard atmosphere at the case's elevation",
        "wind": "the case's sfcWind" if supplied_wind is not None else f"{WIND_M_S} m/s at the 10 m forcing height",
        "cosz": f"{COSZ_DAY} where the row has shortwave, 0 where it has none",
    }
    return {"tk": tk, "sw": sw, "lw": lw, "pressure": pressure, "spechum": spechum, "wind": wind,
            "target": target, "prescribed_rn": rn, "labels": labels}


# --- Noah-MP files -------------------------------------------------------------


def parse_time(text: str) -> datetime:
    return datetime.fromisoformat(text.strip().replace(" ", "T"))


def julian(when: datetime) -> float:
    """CALC_DECLIN's JULIAN: whole days since 1 January, plus the hour as a fraction."""
    days = (when.date() - datetime(when.year, 1, 1).date()).days
    return days + (when.hour + when.minute / 60.0 + when.second / 3600.0) / 24.0


def column_layers(static: dict, theta_sat: float) -> list[float]:
    """Layer thicknesses (module docstring)."""
    if "soil_layer_depth_m" in static:
        return [float(static["soil_layer_depth_m"])] + LAYERS_M[1:]
    if "soil_capacity_mm" not in static:
        return list(LAYERS_M)
    depth = float(static["soil_capacity_mm"]) * (1.0 - CAPACITY_GUARD) / 1000.0 / theta_sat
    top = min(LAYERS_M[0], 0.5 * depth)
    weights = np.array(LAYERS_M[1:])
    return [top] + [float(v) for v in (depth - top) * weights / weights.sum()]


def write_namelist(path: Path, step_s: int, layers: list[float]) -> None:
    """Upstream's example with the Noah-MP step, which is the row's unless it has a substep."""
    text = (SHIPPED / "namelist.hrldas_example").read_text()
    settings = {
        "HRLDAS_SETUP_FILE": '"none"', "INDIR": '"none"', "OUTDIR": '"none"',
        "KDAY": "1", "FORCING_TIMESTEP": str(step_s), "NOAH_TIMESTEP": str(step_s),
        "SOIL_TIMESTEP": "0.0", "OUTPUT_TIMESTEP": str(step_s), "RESTART_FREQUENCY_HOURS": "0",
    }
    settings.update({f"soil_thick_input({i + 1})": repr(d) for i, d in enumerate(layers)})
    for item in filter(None, os.environ.get("NOAHMP_HT_OPTIONS", "").split(",")):
        key, value = item.split("=")
        settings[key] = value
    for key, value in settings.items():
        text, count = re.subn(rf"(?m)^\s*{re.escape(key)}\s*=.*$", f" {key} = {value}", text)
        if count != 1:
            raise SystemExit(f"upstream's namelist has no single {key}")
    path.write_text(text)


def write_control(path: Path, nrows: int, steps_per_row: int, veg: int, tmn: float, skin: float, soil_t: list[float],
                  soil_m: list[float], latitude: float, csoil: float, eg: float, ch2op: float) -> None:
    gvf = (0.0, 0.0) if veg == VEG_BARREN else (MAX_VEG_PCT, MIN_VEG_PCT)
    path.write_text(
        "&control\n"
        f" nrows={nrows}, steps_per_row={steps_per_row}, vegtype={veg}, soiltype={SOIL_CLASS},\n"
        f" latitude={latitude!r}, longitude={LONGITUDE!r}, tmn={tmn!r}, skin={skin!r},\n"
        f" gvfmax={gvf[0]!r}, gvfmin={gvf[1]!r}, lai_initial={0.0 if veg == VEG_BARREN else LAI_INITIAL!r},\n"
        f" soil_t={','.join(repr(v) for v in soil_t)},\n"
        f" soil_m={','.join(repr(v) for v in soil_m)},\n"
        f" csoil_override={csoil!r}, eg_override={eg!r}, ch2op_override={ch2op!r}\n"
        "/\n")


# --- run -----------------------------------------------------------------------


def simulate(rows: list[dict], columns: list[str], static: dict, timestep: str, substep: str | None = None):
    started = time.monotonic()
    step_s = STEP_SECONDS[timestep]
    dt_days = step_s / SECONDS_PER_DAY
    noah_s = STEP_SECONDS[substep] if substep else step_s
    steps_per_row = step_s // noah_s
    noah_days = noah_s / SECONDS_PER_DAY
    pr = np.array([float(r["pr"]) for r in rows])
    tas = np.array([float(r["tas"]) for r in rows])

    theta_sat = porosity(SOIL_CLASS)
    if "soil_porosity" in static and abs(float(static["soil_porosity"]) - theta_sat) > 1.0e-9:
        raise SystemExit(f"the case's porosity {static['soil_porosity']} is not the soil class's {theta_sat}")
    layers = column_layers(static, theta_sat)

    canopy_capacity = float(static.get("canopy_capacity_mm", -1.0))
    veg = VEG_BARREN if canopy_capacity == 0.0 else VEG_CROPLAND
    vai_max = monthly_vai_max(veg)
    ch2op = canopy_capacity / vai_max if canopy_capacity > 0.0 and vai_max > 0.0 else -1.0

    eg = float(static["eps"]) if "eps" in static else -1.0
    emissivity = eg if eg >= 0.0 else table_scalar_eg()
    latitude = float(static.get("latitude_deg", LATITUDE_DEFAULT))

    # First-year air temperature, and nothing later (AGENTS.md, Nothing from the future).
    first_year = tas[: int(round(365.0 * SECONDS_PER_DAY / step_s))]
    if "soil_temperature_initial" in static:
        tmn = float(static["soil_temperature_initial"])
    else:
        tmn = float(first_year.mean()) + 273.15
    skin = tmn if "soil_temperature_initial" in static else float(tas[0]) + 273.15
    theta0 = static.get("soil_water_content_initial")
    soil_m = [float(theta0)] * 4 if theta0 is not None else list(SOIL_MOISTURE_INITIAL)

    csoil = -1.0
    if "soil_heat_capacity_areal" in static:
        # Noah-MP's mixture rule, solved for the solid at the case's water content.
        theta = soil_m[0]
        volumetric = float(static["soil_heat_capacity_areal"]) / layers[0]
        csoil = (volumetric - theta * HEAT_CAPACITY_WATER - (theta_sat - theta) * HEAT_CAPACITY_AIR) / (1.0 - theta_sat)

    atmos = atmosphere(rows, columns, float(static.get("elevation_m", 0.0)), emissivity)
    cosz = np.where(atmos["sw"] > 0.0, COSZ_DAY, 0.0)

    # One record per Noah-MP step: the row's values held, the clock moving on.
    starts = [parse_time(r["time"]) for r in rows]
    ends = [s + timedelta(seconds=noah_s * (k + 1)) for s in starts for k in range(steps_per_row)]
    held = np.column_stack([
        atmos["sw"], atmos["lw"], atmos["tk"], atmos["wind"], atmos["spechum"], atmos["pressure"],
        pr / SECONDS_PER_DAY, cosz])
    table = np.column_stack([np.repeat(held, steps_per_row, axis=0), [julian(e) for e in ends], [e.year for e in ends]])

    work = Path(tempfile.mkdtemp(prefix="noahmp"))
    try:
        shutil.copy(SHIPPED / "NoahmpTable.TBL", work / "NoahmpTable.TBL")
        write_namelist(work / "namelist.hrldas", noah_s, layers)
        write_control(work / "control.nml", len(rows), steps_per_row, veg, tmn, skin, [tmn] * 4, soil_m, latitude,
                      csoil, eg, ch2op)
        table.astype(np.float64).tofile(work / "forcing.bin")
        staged = time.monotonic()
        proc = subprocess.run([DRIVER], cwd=work, capture_output=True, text=True)
        log = proc.stdout + proc.stderr
        if proc.returncode != 0 or "error" in log.lower():
            tail = "\n".join(log.strip().splitlines()[-25:])
            raise RuntimeError(f"Noah-MP failed (exit {proc.returncode}):\n{tail}")
        raw = np.fromfile(work / "output.bin", dtype=np.float64)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    ran = time.monotonic()

    if raw.size != len(rows) * len(OUTPUT_NAMES):
        raise RuntimeError(f"the driver wrote {raw.size} values for {len(rows)} rows")
    out = dict(zip(OUTPUT_NAMES, raw.reshape(len(rows), len(OUTPUT_NAMES)).T))

    per_day = SECONDS_PER_DAY
    ponds = out["PondSfcThinSnwMelt"] + out["PondSfcThinSnwComb"] + out["PondSfcThinSnwTrans"]
    result = {
        "evspsbl": (out["EvapCanopyNet"] + out["Transpiration"] + out["EvapGroundNet"]) * per_day,
        "mrro": (out["RunoffSurface"] + out["RunoffSubsurface"]) / noah_days,
        "sbl": (out["SublimSnowSfcIce"] + out["SublimCanopyIce"]) * per_day,
        "snm": (out["SnowBotOutflow"] + out["RainWithoutPack"]) * per_day + ponds / noah_days,
        "hfls": out["HeatLatentCanopy"] + out["HeatLatentGrd"] + out["HeatLatentTransp"],
        "hfss": out["HeatSensibleSfc"],
        "hfg": out["HeatGroundTot"],
        "rlus": out["RadLwEmitSfc"],
        "hfg_bottom": out["HeatSoil1Bottom"],
        "mrso": out["SoilWater"],
        "snw": out["SnowWaterEquiv"],
        "canopy": out["CanopyLiqWater"] + out["CanopyIce"],
        "ts": out["TemperatureRadSfc"],
        "tsoil_layer": out["TemperatureSoil1"],
    }
    for name, values in result.items():
        if not np.isfinite(values).all():
            raise RuntimeError(f"Noah-MP produced non-finite {name} on {int((~np.isfinite(values)).sum())} steps")

    catchment = {
        "modis_class": veg, "stas_class": SOIL_CLASS, "porosity": theta_sat,
        "soil_layers_m": [round(x, 5) for x in layers],
        "max_vegetation_pct": 0.0 if veg == VEG_BARREN else MAX_VEG_PCT,
        "deep_soil_temperature_K": round(tmn, 4),
        "soil_moisture_initial": soil_m,
        "overrides": {"CH2OP": ch2op if ch2op >= 0.0 else "table", "EG(1)": eg if eg >= 0.0 else "table",
                      "CSOIL": csoil if csoil >= 0.0 else "table"},
    }
    notes = describe(result, out, atmos, pr, dt_days, timestep, steps_per_row, catchment, soil_m, layers)
    notes["seconds"] = {"stage": round(staged - started, 2), "noahmp": round(ran - staged, 2)}
    diag = {**out, "rn_target": atmos["target"] if atmos["target"] is not None else np.full(len(rows), np.nan),
            "rn_prescribed": atmos["prescribed_rn"] if atmos["prescribed_rn"] is not None else np.full(len(rows), np.nan),
            "rn_noahmp": out["RadSwAbsSfc"] - out["RadLwNetSfc"]}
    return result, notes, diag


def describe(result, out, atmos, pr, dt_days, timestep, steps_per_row, catchment, soil_m, layers):
    """What run.json says about a case: the mapping, and the budgets as Noah-MP kept them."""
    storage = result["mrso"] + result["snw"] + result["canopy"]
    initial = sum(m * d * 1000.0 for m, d in zip(soil_m, layers))
    fluxes = (pr - result["evspsbl"] - result["mrro"]) * dt_days
    step_residual = fluxes - np.diff(np.concatenate(([initial], storage)))
    noahmp_rn = out["RadSwAbsSfc"] - out["RadLwNetSfc"]
    surface = result["hfls"] + result["hfss"] + result["hfg"]
    target = atmos["target"]
    prescribed = atmos["prescribed_rn"]
    sunlit =out["RadSwDownRefHeight"] > 0.0  # AlbedoSfc is Noah-MP's undefined value in the dark
    diagnostics = {
        "water_max_abs_step_residual_mm": float(np.abs(step_residual).max()),
        "water_cumulative_residual_mm": float(step_residual.sum()),
        "water_input_mm": float((pr * dt_days).sum()),
        "noahmp_water_balance_error_max_abs_mm": float(np.abs(out["WaterBalanceError"]).max()),
        "noahmp_energy_balance_error_max_abs_w_m2": float(np.abs(out["EnergyBalanceError"]).max()),
        "aquifer_mm_constant": [float(out["WaterStorageAquifer"].min()), float(out["WaterStorageAquifer"].max())],
        # The prescribed rn beside Noah-MP's own, as issue #59 asked: this tells an
        # input mismatch from an energy-budget problem. `target` is the net radiation
        # the sky was built on, which is the probe's rn only where one is prescribed.
        "mean_target_net_radiation_w_m2": float(target.mean()) if target is not None else None,
        "mean_prescribed_rn_w_m2": float(prescribed.mean()) if prescribed is not None else None,
        "mean_noahmp_net_radiation_w_m2": float(noahmp_rn.mean()),
        "mean_abs_target_minus_noahmp_net_radiation_w_m2":
            float(np.abs(target - noahmp_rn).mean()) if target is not None else None,
        "mean_abs_prescribed_rn_minus_noahmp_net_radiation_w_m2":
            float(np.abs(prescribed - noahmp_rn).mean()) if prescribed is not None else None,
        "mean_hfls_hfss_hfg_w_m2": float(surface.mean()),
        "mean_noahmp_net_radiation_minus_hfls_hfss_hfg_w_m2": float((noahmp_rn - surface).mean()),
        "mean_abs_noahmp_net_radiation_minus_hfls_hfss_hfg_w_m2": float(np.abs(noahmp_rn - surface).mean()),
        "mean_precipitation_advected_heat_w_m2": float(out["HeatPrecipAdvSfc"].mean()),
        "mean_canopy_heat_storage_change_w_m2": float(out["HeatCanStorageChg"].mean()),
        "mean_albedo_when_sunlit": float(out["AlbedoSfc"][sunlit].mean()) if sunlit.any() else None,
        "snow_layers_max": int(-out["NumSnowLayerNeg"].min()),
        "negative_evspsbl_steps": int((result["evspsbl"] < 0.0).sum()),
    }
    return {
        "noahmp": {"release": "v5.2.1", "commit": "17751dcd7a2442a3c57138f530f523daf08dc1c5",
                   "hrldas_commit": "cd96df470220f7d7133cdbccd5f9c5355cf173e2",
                   "precision": "upstream's default single precision",
                   "setup": "namelist.hrldas_example options and NoahmpTable.TBL, unedited; land point from the "
                            "single_point example (Bondville)"},
        "timestep": timestep,
        "steps_per_row": steps_per_row,
        "catchment": catchment,
        "mock_inputs": atmos["labels"],
        "prescribed_rn": "not held: rn only builds the sky; Noah-MP computes its own net radiation",
        "states": {"mrso": "SoilMoisture over layer thickness", "snw": "SnowWaterEquiv",
                   "canopy": "CanopyLiqWater + CanopyIce"},
        "fluxes": {"hfg": "HeatGroundTot, top of the snow-soil column",
                   "hfg_bottom": "2 k1 (T1 - T2) / -z2 from SoilSnowThermalDiffusion, end-of-step temperatures",
                   "rlus": "RadLwEmitSfc", "ts": "TemperatureRadSfc"},
        "diagnostics_whole_record": diagnostics,
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
    timestep = str(request.get("timestep", "PT1H"))
    if timestep not in STEP_SECONDS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    substep = request.get("substep")

    result, notes, diag = simulate(rows, columns, static, timestep, substep)

    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(COLUMNS)
        for i, row in enumerate(rows):
            writer.writerow([row["time"], row["pr"]] + [OUTPUT_FORMAT % result[c][i] for c in COLUMNS[2:]])

    if DIAG:
        names = sorted(diag)
        with open(io_dir / "diag.csv", "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["time"] + names)
            for i, row in enumerate(rows):
                writer.writerow([row["time"]] + ["%.12g" % diag[n][i] for n in names])

    notes["seed"] = int(request.get("seed", 0))
    notes["wall_seconds"] = round(time.monotonic() - started, 2)
    (io_dir / request["output"]["run"]).write_text(json.dumps(
        {"status": "ok", "model": MODEL, "n_steps": len(rows), "zero_stores": ZERO_STORES,
         "not_reported": NOT_REPORTED, "notes": notes}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

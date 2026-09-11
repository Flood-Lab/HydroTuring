#!/usr/bin/env python3
"""HydroTuring adapter for CWatM 1.11, the Community Water Model (IIASA).

CWatM is a grid model driven by one settings file and netCDF maps. This
adapter builds the smallest domain CWatM accepts for a lumped catchment, one
active grid cell whose area is the catchment's, writes the forcing and every
parameter map that cell needs into a scratch directory, and runs the model
through its own entry points (parse_configuration, checkifDate, CWATModel,
ModelFrame), stepping the frame one row at a time so the stores can be read
after every step. Nothing in the model's code is patched.

What is reported
----------------
Fluxes as rates in mm per day, states as absolute storages in mm, all taken
from CWatM's own variables after each step (CWatM works in metres of water
over the cell; land-cover variables are its own fraction-weighted sums):

* `pr`       the forcing, echoed
* `evspsbl`  totalET: transpiration, bare-soil evaporation, open-water
             evaporation, interception evaporation and snow evaporation
* `mrro`     runoff: surface runoff, interflow and baseflow after CWatM's
             runoff-concentration lag, i.e. what leaves the cell
* `dis`      the same over the catchment area, m3/s
* `mrso`     sum_soil = sum_w1 + sum_w2 + sum_w3 (+ sum_topwater, zero
             without paddy fields)
* `snw`      SnowCover (the degree-day snow store; it holds no liquid water)
* `canopy`   sum_interceptStor
* `gw`       storGroundwater, the linear groundwater reservoir
* `channel`  gridcell_storage, runoff generated but still inside the
             runoff-concentration lag

CWatM has no exchange with the outside in this configuration (no MODFLOW,
no abstraction, no inflow), so `gwex` is not emitted. The budget
P - ET - Q = d(tws) that CWatM itself defines (tws = storGroundwater +
totalSto + gridcell_storage) is checked every step and its largest residual
is recorded in run.json.

Configuration
-------------
Switched off because no probe prescribes them: irrigation and water demand
(so abstraction, return flows and paddy fields do not exist), lakes and
reservoirs, MODFLOW groundwater, inflow hydrographs, environmental flow,
water quality, glaciers, pySnowClim, and kinematic-wave routing (a single
cell has no river network to route along; the within-cell lag CWatM does
have, runoff concentration, stays on). CWatM still initialises its routing
module, so the channel maps and lakeEvaFactor are written with placeholder
values that no flux reads. Potential evaporation is not computed:
calc_evaporation = False makes CWatM read reference evapotranspiration
(ETMaps) and open-water evaporation (E0Maps) directly, and both are fed the
probe's pet. Open water and sealed surfaces have zero area, so E0Maps is read
but used by no flux. Capillary rise and runoff concentration are on, as in
every 30-arcminute template of CWatM and of CWatM-Earth-30min. Preferential
flow is on as in the CWatM-Earth-30min template
(settings_CWatM_template_30min.ini at e9dfd99), the repository the parameter
maps come from; the pinned model repository's own 30-arcminute templates and
its global 30-arcminute test setup ship it off. That switch is a packaging
choice that decides two verdicts: with preferentialFlow = False,
pet-consistency and response-nonnegativity both pass on the gate seeds
(README.md, Sensitivity).

The cell is forest and grassland only. Land cover, soil hydraulics, rooting,
crop coefficients, groundwater recession, relative elevation, slope and
orographic spread are not in static.json; they take the median over the land
cells of CWatM's own 30-arcminute input maps (iiasa/CWatM-Earth-30min at
e9dfd99), listed in GLOBAL_MEDIAN below with the spread in README.md. The
forest share is the area-weighted forest share of forest plus grassland in
the same maps. Calibration factors take the neutral values the templates'
[CALIBRATION] comments document, not the example calibration; SnowMeltCoef,
for which no neutral value is documented, is a choice (see CALIBRATION_DEFAULTS).

From static.json: the cell area (area_km2), the snow threshold
(snow_threshold_degC -> TempSnow), the interception capacity
(canopy_capacity_mm -> interceptCap of forest and grassland) and the soil
capacity (soil_capacity_mm). CWatM's soil store is bounded by saturation, so
the soil column is sized to hold exactly soil_capacity_mm at saturation: the
two storage depths keep the ratio of the global medians and are scaled
together until the saturated storage of the three layers, weighted over both
land covers, equals the capacity. Crop coefficients and interception
capacities are held constant through the year: a global median of a seasonal
cycle mixes hemispheres and describes no catchment.

The timestep
------------
CWatM 1.11 has no sub-daily step. miscInitial sets DtSec = 86400 as a
literal, and the soil, groundwater, snow, interception and runoff-concentration
rates are per day with nothing that rescales them. The adapter therefore
hands CWatM one row per step at every timestep, as the depth that row carries
(rate x step length), and CWatM's calendar advances one day per row. Its snow
scheme reads that calendar: the melt coefficient follows a northern-hemisphere
sine of the day of year, and extra ice melt acts on snow between days 166 and
259, so at PT1H both run on a calendar 24 times too fast. At PT1H
that is the model applying a day of drainage, percolation and evaporation to
every hour of forcing; the resolution probe measures exactly that, and no
parameter is rescaled to hide it.

The model is deterministic; the request seed is recorded and otherwise
unused.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import datetime as dt
import io
import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

MODEL = {"name": "cwatm", "version": "1.11-5baaadd.1"}
COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "mrso", "snw", "canopy", "gw", "channel"]
TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0, "PT15M": 1.0 / 96.0, "PT5M": 1.0 / 288.0, "PT1M": 1.0 / 1440.0}

CELL_DEG = 0.5  # size of the one-cell domain in degrees; its area is set by CellArea
DZREL_NAMES = ["dzRel0001", "dzRel0005", "dzRel0010", "dzRel0020", "dzRel0030", "dzRel0040",
               "dzRel0050", "dzRel0060", "dzRel0070", "dzRel0080", "dzRel0090", "dzRel0100"]
SOIL_NAMES = (("KSat", "ksat"), ("alpha", "alpha"), ("lambda", "lambda"), ("thetas", "thetas"), ("thetar", "thetar"))

# Median over land cells of CWatM's 30-arcminute input maps (CWatM-Earth-30min
# at e9dfd99; README.md has the script and the spread). Units as CWatM reads
# them: KSat in cm/day, depths in m, recession coefficient per day, crop
# coefficient and interception capacity as the annual mean of the 10-day maps.
GLOBAL_MEDIAN = {
    # soil hydraulics (van Genuchten), layers 1-3; forest has its own maps for layers 1-2
    "ksat1": 3.05714, "ksat2": 3.05714, "ksat3": 2.91572,
    "forest_ksat1": 3.436, "forest_ksat2": 3.436,
    "alpha1": 0.03767, "alpha2": 0.03767, "alpha3": 0.04208,
    "forest_alpha1": 0.03959, "forest_alpha2": 0.03959,
    "lambda1": 0.1577, "lambda2": 0.1577, "lambda3": 0.15418,
    "forest_lambda1": 0.15503, "forest_lambda2": 0.15503,
    "thetas1": 0.44685, "thetas2": 0.44685, "thetas3": 0.42345,
    "forest_thetas1": 0.48756, "forest_thetas2": 0.48756,
    "thetar1": 0.09659, "thetar2": 0.09659, "thetar3": 0.09921,
    "forest_thetar1": 0.09659, "forest_thetar2": 0.09659,
    "percolationImp": 0.16628,
    "storageDepth1": 0.2669, "storageDepth2": 1.0678,
    "cropgrp": 3.0,
    # land cover
    "forest_cropCoefficient": 0.8569, "forest_interceptCap": 0.00042,
    "forest_maxRootDepth": 1.0, "forest_rootFraction1": 0.64026,
    "grassland_cropCoefficient": 0.51557, "grassland_interceptCap": 0.00017,
    "grassland_maxRootDepth": 0.5, "grassland_rootFraction1": 0.67863,
    "irrPaddy_maxRootDepth": 0.5, "irrPaddy_rootFraction1": 0.97896,
    "irrNonPaddy_maxRootDepth": 1.0, "irrNonPaddy_rootFraction1": 0.88463,
    "forest_share": 0.5076,
    # groundwater and topography
    "recessionCoeff": 0.0069, "specificYield": 0.05,
    "elvstd": 57.026, "tanslope": 0.012183,
    "dzRel": {
        "dzRel0001": -41.881, "dzRel0005": -38.392, "dzRel0010": -33.142, "dzRel0020": -25.202,
        "dzRel0030": -17.23, "dzRel0040": -9.602, "dzRel0050": -2.152, "dzRel0060": 0.0,
        "dzRel0070": 7.628, "dzRel0080": 20.105, "dzRel0090": 38.729, "dzRel0100": 115.163,
    },
}

# preferentialFlowConstant, arnoBeta_add, factor_interflow, recessionCoeff_factor
# and runoffConc_factor take the neutral values the [CALIBRATION] comments of the
# 30-arcminute templates document; crop_correct and soildepth_factor are 1 (no
# scaling). SnowMeltCoef has no documented neutral value: CWatM-Earth-30min's
# template sets 0.0027 and the model repository's 30-arcminute tests 0.0034.
# 0.004 m/degC/day is a choice: the value those files carry commented out under
# [SNOW] ("default: 4.0" mm) and the one the model repository's Bhima 1 km test
# sets (pytest/settings/1km/Bhima/settings_Bhima.ini:171). manningsN and
# lakeEvaFactor are read by the routing module's initialisation only; routing
# is off.
CALIBRATION_DEFAULTS = {
    "SnowMeltCoef": 0.004,
    "crop_correct": 1.0,
    "soildepth_factor": 1.0,
    "preferentialFlowConstant": 4.0,
    "arnoBeta_add": 0.1,
    "factor_interflow": 1.0,
    "recessionCoeff_factor": 1.0,
    "runoffConc_factor": 1.0,
    "manningsN": 1.0,
    "lakeEvaFactor": 1.0,
}

# Constants of the shipped settings template that have no counterpart in
# static.json, taken as they are.
TEMPLATE = {
    "GlacierTransportZone": 3,
    "TemperatureLapseRate": 0.0065,
    "SnowFactor": 1.0,
    "SnowSeasonAdj": 0.001,
    "TempMelt": 1.0,
    "TempSnow": 1.0,
    "IceMeltCoef": 0.007,
    "SnowWaterEquivalent": 0.45,
    "Afrost": 0.97,
    "Kfrost": 0.57,
    "FrostIndexThreshold": 56,
    "maxGWCapRise": 5.0,
    "minCropKC": 0.2,
    "minTopWaterLayer": 0.0,
    "forest_arnoBeta": 0.2,
    "grassland_arnoBeta": 0.0,
    "irrPaddy_arnoBeta": 0.2,
    "irrNonPaddy_arnoBeta": 0.2,
    "minInterceptCap": 0.001,
    "forest_runoff_peaktime": 1.0,
    "grassland_runoff_peaktime": 0.5,
    "irrPaddy_runoff_peaktime": 0.5,
    "irrNonPaddy_runoff_peaktime": 0.5,
    "sealed_runoff_peaktime": 0.15,
    "water_runoff_peaktime": 0.01,
    "interflow_runoff_peaktime": 1.0,
    "baseflow_runoff_peaktime": 2.0,
    "NoRoutingSteps": 10,
    "chanBeta": 0.6,
    "chanGradMin": 0.0001,
}

# One elevation zone. CWatM's template splits a cell into seven zones along
# its relative-elevation distribution with a lapse rate; a lumped case has one
# elevation.
NUMBER_SNOW_LAYERS = 1


# --- the domain ---------------------------------------------------------------


def write_grid(path: Path, lat: float, fields: dict, times: np.ndarray | None = None,
               time_units: str | None = None) -> None:
    """A 2 x 2 netCDF grid (CWatM reads the cell size from two coordinates);
    the adapter's mask selects the upper-left cell. Every field is uniform."""
    with Dataset(path, "w", format="NETCDF4") as nc:
        if times is not None:
            nc.createDimension("time", len(times))
            tv = nc.createVariable("time", "f8", ("time",))
            tv.units = time_units
            tv.calendar = "standard"
            tv[:] = times
        nc.createDimension("lat", 2)
        nc.createDimension("lon", 2)
        la = nc.createVariable("lat", "f8", ("lat",))
        la[:] = [lat, lat - CELL_DEG]
        lo = nc.createVariable("lon", "f8", ("lon",))
        lo[:] = [CELL_DEG / 2.0, 1.5 * CELL_DEG]
        for name, value in fields.items():
            if times is None:
                var = nc.createVariable(name, "f8", ("lat", "lon"))
                var[:] = np.full((2, 2), float(value))
            else:
                var = nc.createVariable(name, "f8", ("time", "lat", "lon"))
                series = np.asarray(value, dtype=float)
                var[:] = np.repeat(np.repeat(series[:, None, None], 2, axis=1), 2, axis=2)


def saturated_storage(depth1: float, depth2: float, forest: float, med: dict) -> float:
    """CWatM's soil store at saturation (m) for StorDepth1/2, with the root
    depths landcoverType builds for forest and grassland, fraction-weighted."""
    s1 = max(0.05, depth1 - 0.05)
    s2 = max(0.05, depth2)
    s12 = s1 + s2
    total = 0.0
    for cover, frac in (("forest", forest), ("grassland", 1.0 - forest)):
        if cover == "forest":
            h1 = max(s1, med["forest_maxRootDepth"] - 0.05)
            r1 = min(s12 - 0.05, h1)
            r2 = max(0.05, s12 - r1)
            ts = (med["forest_thetas1"], med["forest_thetas2"], med["thetas3"])  # layer 3 shares the map
        else:
            r1, r2 = s1, s2
            ts = (med["thetas1"], med["thetas2"], med["thetas3"])
        total += frac * (ts[0] * 0.05 + ts[1] * r1 + ts[2] * r2)
    return total


def soil_depths(capacity_m: float, forest: float, med: dict) -> tuple[float, float]:
    """Scale the median storage depths together until saturation holds the capacity."""
    d1, d2 = med["storageDepth1"], med["storageDepth2"]
    lo, hi = 0.0, 50.0
    if saturated_storage(0.0, 0.0, forest, med) > capacity_m:
        raise ValueError(f"soil capacity {capacity_m * 1000:.1f} mm is below CWatM's minimum soil column")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if saturated_storage(d1 * mid, d2 * mid, forest, med) < capacity_m:
            lo = mid
        else:
            hi = mid
    k = 0.5 * (lo + hi)
    return d1 * k, d2 * k


def parameter_maps(med: dict, forest: float, depth1: float, depth2: float) -> dict[str, float]:
    """Every binding CWatM's own setups give as a map, as one-cell maps."""
    values: dict[str, float] = {}
    for layer in (1, 2, 3):
        for cwatm_name, key in SOIL_NAMES:
            values[f"{cwatm_name}{layer}"] = med[f"{key}{layer}"]
            values[f"forest_{cwatm_name}{layer}"] = med[f"forest_{key}{layer}" if layer < 3 else f"{key}{layer}"]
    values.update({
        "ElevationStD": med["elvstd"],
        "cropgroupnumber": med["cropgrp"],
        "tanslope": med["tanslope"],
        "percolationImp": med["percolationImp"],
        "StorDepth1": depth1,
        "StorDepth2": depth2,
        "forest_fracVegCover": forest,
        "irrPaddy_fracVegCover": 0.0,
        "irrNonPaddy_fracVegCover": 0.0,
        "sealed_fracVegCover": 0.0,
        "water_fracVegCover": 0.0,
        "recessionCoeff": med["recessionCoeff"],
        "specificYield": med["specificYield"],
        # routing is off; its initialisation still reads the channel geometry
        "chanGrad": 1.0, "chanMan": 1.0, "chanLength": 1.0, "chanWidth": 1.0, "chanDepth": 1.0,
    })
    for cover in ("forest", "grassland", "irrPaddy", "irrNonPaddy"):
        values[f"{cover}_rootFraction1"] = med[f"{cover}_rootFraction1"]
        values[f"{cover}_maxRootDepth"] = med[f"{cover}_maxRootDepth"]
    return values


def settings_text(work: Path, start: dt.date, n_steps: int, lat: float, static: dict, forest: float) -> str:
    p = {**TEMPLATE, **CALIBRATION_DEFAULTS}
    if "snow_threshold_degC" in static:
        p["TempSnow"] = float(static["snow_threshold_degC"])
    min_icap = p["minInterceptCap"]
    if "canopy_capacity_mm" in static:
        min_icap = min(min_icap, float(static["canopy_capacity_mm"]) / 1000.0)
    ul_lat = lat + CELL_DEG / 2.0
    root = str(work)

    def m(key: str) -> str:
        return f"{root}/maps/{key}.nc"

    soil_block = "\n".join(
        f"{prefix}{name}{layer} = {m(prefix + name + str(layer))}"
        for layer in (1, 2, 3) for name, _ in SOIL_NAMES for prefix in ("", "forest_")
    )
    return f"""[OPTIONS]
TemperatureInKelvin = False
gridSizeUserDefined = True
calc_evaporation = False
includeIrrigation = False
includeWaterDemand = False
usingAllocSegments = False
limitAbstraction = True
calc_environflow = False
preferentialFlow = True
CapillarRise = True
includeRunoffConcentration = True
includeWaterBodies = False
includeRouting = False
inflow = False
waterquality = False
modflow_coupling = False
includeCrops = False
includeGlaciers = False
usepySnowClim = False
staticLandCoverMaps = True
writeNetcdfStack = False
reportMap = False
reportTss = False
calcWaterBalance = False
sumWaterBalance = False

[FILE_PATHS]
PathRoot = {root}
PathOut = {root}/out
PathMaps = {root}
PathMeteo = {root}

[NETCDF_ATTRIBUTES]
institution = HydroTuring
title = CWatM lumped probe case
metaNetcdfFile = {root}/metaNetcdf.xml

[MASK_OUTLET]
MaskMap = 1 1 {CELL_DEG} 0.0 {ul_lat!r}
Gauges = {CELL_DEG / 2.0} {lat!r}
GaugesLocal = True

[TIME-RELATED_CONSTANTS]
StepStart = {start.day:02d}/{start.month:02d}/{start.year:04d}
SpinUp = None
StepEnd = {n_steps}

[INITITIAL CONDITIONS]
load_initial = False
initLoad = {root}/init.nc
save_initial = False
initSave = {root}/init
StepInit = 31/12/2100

[CALIBRATION]
SnowMeltCoef = {p['SnowMeltCoef']!r}
crop_correct = {p['crop_correct']!r}
soildepth_factor = {p['soildepth_factor']!r}
preferentialFlowConstant = {p['preferentialFlowConstant']!r}
arnoBeta_add = {p['arnoBeta_add']!r}
factor_interflow = {p['factor_interflow']!r}
recessionCoeff_factor = {p['recessionCoeff_factor']!r}
runoffConc_factor = {p['runoffConc_factor']!r}
manningsN = {p['manningsN']!r}
lakeEvaFactor = {p['lakeEvaFactor']!r}

[TOPOP]
Ldd = {root}/ldd.nc
ElevationStD = {m('ElevationStD')}
CellArea = {root}/cellarea.nc

[METEO]
PrecipitationMaps = {root}/pr.nc
TavgMaps = {root}/tavg.nc
E0Maps = {root}/ewref.nc
ETMaps = {root}/etref.nc
precipitation_coversion = 0.001
evaporation_coversion = 0.001

[SNOW]
NumberSnowLayers = {NUMBER_SNOW_LAYERS}
GlacierTransportZone = {p['GlacierTransportZone']}
TemperatureLapseRate = {p['TemperatureLapseRate']!r}
SnowFactor = {p['SnowFactor']!r}
SnowSeasonAdj = {p['SnowSeasonAdj']!r}
TempMelt = {p['TempMelt']!r}
TempSnow = {p['TempSnow']!r}
IceMeltCoef = {p['IceMeltCoef']!r}

[FROST]
SnowWaterEquivalent = {p['SnowWaterEquivalent']!r}
Afrost = {p['Afrost']!r}
Kfrost = {p['Kfrost']!r}
FrostIndexThreshold = {p['FrostIndexThreshold']!r}

[VEGETATION]
cropgroupnumber = {m('cropgroupnumber')}

[SOIL]
tanslope = {m('tanslope')}
relativeElevation = {root}/dzrel.nc
{soil_block}
percolationImp = {m('percolationImp')}
maxGWCapRise = {p['maxGWCapRise']!r}
minCropKC = {p['minCropKC']!r}
minTopWaterLayer = {p['minTopWaterLayer']!r}
StorDepth1 = {m('StorDepth1')}
StorDepth2 = {m('StorDepth2')}

[LANDCOVER]
coverTypes = forest, grassland, irrPaddy, irrNonPaddy, sealed, water
coverTypesShort = f, g, i, n, s, w
dynamicLandcover = False
fixLandcoverYear = 2000

[__forest]
forest_arnoBeta = {p['forest_arnoBeta']!r}
forest_minInterceptCap = {min_icap!r}
forest_cropDeplFactor = 0.0
forest_fracVegCover = {m('forest_fracVegCover')}
forest_rootFraction1 = {m('forest_rootFraction1')}
forest_maxRootDepth = {m('forest_maxRootDepth')}
forest_cropCoefficientNC = {root}/kc_forest.nc
forest_interceptCapNC = {root}/icap_forest.nc

[__grassland]
grassland_arnoBeta = {p['grassland_arnoBeta']!r}
grassland_minInterceptCap = {min_icap!r}
grassland_cropDeplFactor = 0.0
grassland_rootFraction1 = {m('grassland_rootFraction1')}
grassland_maxRootDepth = {m('grassland_maxRootDepth')}
grassland_cropCoefficientNC = {root}/kc_grassland.nc
grassland_interceptCapNC = {root}/icap_grassland.nc

[__irrPaddy]
irrPaddy_arnoBeta = {p['irrPaddy_arnoBeta']!r}
irrPaddy_minInterceptCap = {p['minInterceptCap']!r}
irrPaddy_cropDeplFactor = 0.0
irrPaddy_fracVegCover = {m('irrPaddy_fracVegCover')}
irrPaddy_rootFraction1 = {m('irrPaddy_rootFraction1')}
irrPaddy_maxRootDepth = {m('irrPaddy_maxRootDepth')}
irrPaddy_cropCoefficientNC = {root}/kc_grassland.nc
irrPaddy_maxtopwater = 0.05

[__irrNonPaddy]
irrNonPaddy_arnoBeta = {p['irrNonPaddy_arnoBeta']!r}
irrNonPaddy_minInterceptCap = {p['minInterceptCap']!r}
irrNonPaddy_cropDeplFactor = 0.0
irrNonPaddy_fracVegCover = {m('irrNonPaddy_fracVegCover')}
irrNonPaddy_rootFraction1 = {m('irrNonPaddy_rootFraction1')}
irrNonPaddy_maxRootDepth = {m('irrNonPaddy_maxRootDepth')}
irrNonPaddy_cropCoefficientNC = {root}/kc_grassland.nc

[__sealed]
sealed_minInterceptCap = {p['minInterceptCap']!r}
sealed_fracVegCover = {m('sealed_fracVegCover')}

[__open_water]
water_minInterceptCap = 0.0
water_fracVegCover = {m('water_fracVegCover')}

[GROUNDWATER]
recessionCoeff = {m('recessionCoeff')}
specificYield = {m('specificYield')}

[RUNOFF_CONCENTRATION]
forest_runoff_peaktime = {p['forest_runoff_peaktime']!r}
grassland_runoff_peaktime = {p['grassland_runoff_peaktime']!r}
irrPaddy_runoff_peaktime = {p['irrPaddy_runoff_peaktime']!r}
irrNonPaddy_runoff_peaktime = {p['irrNonPaddy_runoff_peaktime']!r}
sealed_runoff_peaktime = {p['sealed_runoff_peaktime']!r}
water_runoff_peaktime = {p['water_runoff_peaktime']!r}
interflow_runoff_peaktime = {p['interflow_runoff_peaktime']!r}
baseflow_runoff_peaktime = {p['baseflow_runoff_peaktime']!r}

[ROUTING]
NoRoutingSteps = {p['NoRoutingSteps']}
chanBeta = {p['chanBeta']!r}
chanGrad = {m('chanGrad')}
chanGradMin = {p['chanGradMin']!r}
chanMan = {m('chanMan')}
chanLength = {m('chanLength')}
chanWidth = {m('chanWidth')}
chanDepth = {m('chanDepth')}

[OUTPUT]
OUT_Dir = {root}/out
"""


def build_case(work: Path, forcing: list[dict], static: dict, dt_days: float,
               med: dict | None = None) -> dict:
    """Write every file CWatM reads for this case. Returns what was chosen."""
    med = GLOBAL_MEDIAN if med is None else med
    lat = min(max(float(static.get("latitude_deg", 45.0)), -89.0), 89.0)
    area_m2 = float(static["area_km2"]) * 1.0e6
    forest = float(med["forest_share"])

    if "soil_capacity_mm" in static:
        depth1, depth2 = soil_depths(float(static["soil_capacity_mm"]) / 1000.0, forest, med)
    else:
        depth1, depth2 = med["storageDepth1"], med["storageDepth2"]
    icap = float(static["canopy_capacity_mm"]) / 1000.0 if "canopy_capacity_mm" in static else None

    start = dt.date.fromisoformat(str(forcing[0]["time"])[:10])
    n = len(forcing)
    days = np.arange(n, dtype=float)
    units = f"days since {start.isoformat()} 00:00:00"
    pr = np.array([r["pr"] for r in forcing]) * dt_days          # mm per row
    tas = np.array([r["tas"] for r in forcing])
    pet = np.array([r["pet"] for r in forcing]) * dt_days        # mm per row
    write_grid(work / "pr.nc", lat, {"pr": pr}, days, units)
    write_grid(work / "tavg.nc", lat, {"tavg": tas}, days, units)
    write_grid(work / "etref.nc", lat, {"ETRef": pet}, days, units)
    write_grid(work / "ewref.nc", lat, {"EWRef": pet}, days, units)

    write_grid(work / "ldd.nc", lat, {"ldd": 5.0})
    write_grid(work / "cellarea.nc", lat, {"cellarea": area_m2})
    write_grid(work / "dzrel.nc", lat, {name: med["dzRel"][name] for name in DZREL_NAMES})
    (work / "maps").mkdir()
    for key, value in parameter_maps(med, forest, depth1, depth2).items():
        write_grid(work / "maps" / f"{key}.nc", lat, {"value": value})
    ten_day = np.arange(37, dtype=float)
    for cover in ("forest", "grassland"):
        write_grid(work / f"kc_{cover}.nc", lat, {"kc": np.full(37, med[f"{cover}_cropCoefficient"])},
                   ten_day, "days since 2000-01-01 00:00:00")
        cap = icap if icap is not None else med[f"{cover}_interceptCap"]
        write_grid(work / f"icap_{cover}.nc", lat, {"icap": np.full(37, cap)},
                   ten_day, "days since 2000-01-01 00:00:00")
    (work / "out").mkdir(exist_ok=True)
    import cwatm
    shutil.copy(Path(cwatm.__file__).parent / "metaNetcdf.xml", work / "metaNetcdf.xml")
    (work / "settings.ini").write_text(settings_text(work, start, n, lat, static, forest))
    return {
        "soil_depths_m": {"StorDepth1": depth1, "StorDepth2": depth2},
        "soil_saturated_storage_mm": 1000.0 * saturated_storage(depth1, depth2, forest, med),
        "forest_fraction": forest,
        "interception_capacity_mm": 1000.0 * icap if icap is not None else None,
        "calendar_start": start.isoformat(),
    }


# --- the model ----------------------------------------------------------------


def scalar(x) -> float:
    return float(np.asarray(x, dtype=float).reshape(-1)[0])


def run_cwatm(settings: Path, n_steps: int) -> tuple[dict[str, np.ndarray], dict]:
    """CWATMexe from cwatm/run_cwatm.py, with the frame stepped here so every
    step's stores can be read."""
    from cwatm.management_modules.globals import Flags, dateVar, globalFlags, settingsfile
    from cwatm.management_modules.configuration import parse_configuration, read_metanetcdf
    from cwatm.management_modules.data_handling import cbinding
    from cwatm.management_modules.timestep import checkifDate
    from cwatm.management_modules.dynamicModel import ModelFrame
    from cwatm.cwatm_model import CWATModel
    from cwatm.run_cwatm import headerinfo

    globalFlags(str(settings), ["-v"], settingsfile, Flags)
    headerinfo()
    parse_configuration(settingsfile[0])
    read_metanetcdf("metaNetcdf.xml")
    checkifDate("StepStart", "StepEnd", "SpinUp", cbinding("PrecipitationMaps"))
    model = CWATModel()
    frame = ModelFrame(model, firstTimestep=dateVar["intStart"], lastTimeStep=dateVar["intEnd"])
    frame.initialize_run()
    v = model.var

    names = ("P", "ET", "Q", "soil", "snow", "canopy", "gw", "channel", "tws")
    out = {k: np.zeros(n_steps) for k in names}
    previous = scalar(v.totalSto) + scalar(v.storGroundwater) + scalar(v.gridcell_storage)
    worst = 0.0
    i = 0
    while frame.currentStep <= model.lastStep:
        frame.step()
        out["P"][i] = scalar(v.Precipitation)
        out["ET"][i] = scalar(v.totalET)
        out["Q"][i] = scalar(v.runoff)
        out["soil"][i] = scalar(v.sum_soil)
        out["snow"][i] = scalar(v.SnowCover)
        out["canopy"][i] = scalar(v.sum_interceptStor)
        out["gw"][i] = scalar(v.storGroundwater)
        out["channel"][i] = scalar(v.gridcell_storage)
        out["tws"][i] = scalar(v.tws)
        now = out["soil"][i] + out["snow"][i] + out["canopy"][i] + out["gw"][i] + out["channel"][i]
        worst = max(worst, abs(now - previous - (out["P"][i] - out["ET"][i] - out["Q"][i])))
        previous = now
        i += 1
    if i != n_steps:
        raise RuntimeError(f"CWatM ran {i} steps for {n_steps} forcing rows")
    reported = out["soil"] + out["snow"] + out["canopy"] + out["gw"] + out["channel"]
    check = {
        "largest_step_residual_mm": 1000.0 * worst,
        "largest_tws_minus_reported_stores_mm": 1000.0 * float(np.max(np.abs(out["tws"] - reported))),
    }
    return out, check


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet"):
            row[key] = float(row[key])
    return rows


def simulate(forcing: list[dict], static: dict, timestep: str, med: dict | None = None) -> tuple[list[dict], dict]:
    dt_days = TIMESTEP_DAYS[timestep]
    work = Path(tempfile.mkdtemp(prefix="cwatm-"))
    try:
        chosen = build_case(work, forcing, static, dt_days, med)
        with contextlib.redirect_stdout(io.StringIO()):
            out, check = run_cwatm(work / "settings.ini", len(forcing))
    finally:
        shutil.rmtree(work, ignore_errors=True)

    area_m2 = float(static["area_km2"]) * 1.0e6
    rows = []
    for i, step in enumerate(forcing):
        mrro = out["Q"][i] * 1000.0 / dt_days
        rows.append({
            "time": step["time"],
            "pr": step["pr"],
            "evspsbl": out["ET"][i] * 1000.0 / dt_days,
            "mrro": mrro,
            "dis": mrro / 1000.0 * area_m2 / 86400.0,
            "mrso": out["soil"][i] * 1000.0,
            "snw": out["snow"][i] * 1000.0,
            "canopy": out["canopy"][i] * 1000.0,
            "gw": out["gw"][i] * 1000.0,
            "channel": out["channel"][i] * 1000.0,
        })
    notes = {
        **chosen,
        "timestep": timestep,
        "cwatm_step": "one CWatM step (a day in its calendar) per forcing row; DtSec is fixed at 86400 in CWatM 1.11",
        "forcing_as_given_to_cwatm": "pr and pet as depth per row (rate x step length); ETMaps and E0Maps both the probe's pet",
        "water_balance_check": check,
        "states": {
            "mrso": "sum_w1 + sum_w2 + sum_w3 + sum_topwater, fraction-weighted over forest and grassland",
            "snw": "SnowCover, one elevation zone; the degree-day scheme holds no liquid water",
            "canopy": "sum_interceptStor",
            "gw": "storGroundwater",
            "channel": "gridcell_storage: runoff inside CWatM's runoff-concentration lag",
        },
        "not_emitted": {"gwex": "no exchange with the outside in this configuration"},
    }
    return rows, notes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    started = time.monotonic()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent
    forcing = read_forcing(io_dir / request["input"]["forcing"])
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    timestep = str(request.get("timestep", "PT1D"))
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")

    rows, notes = simulate(forcing, static, timestep)
    notes["seed"] = int(request.get("seed", 0))

    table = io_dir / request["output"]["table"]
    table.parent.mkdir(parents=True, exist_ok=True)
    with open(table, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "n_steps": len(rows),
        "wall_seconds": round(time.monotonic() - started, 2),
        "notes": notes,
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

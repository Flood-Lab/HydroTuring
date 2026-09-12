#!/usr/bin/env python3
"""HydroTuring adapter for LISFLOOD 5.0.0 (EC Joint Research Centre).

LISFLOOD is a distributed model: a run is one XML settings file, a mask, a
local drain direction map, and parameter maps on that grid. A lumped probe
case is given here as one representative cell at the resolution the
parameters come from: the shipped test catchment's 5 km grid (cell length
5000 m, cell area 25 km2, channel length 5000 m, that catchment's median).
The cell's drain direction is a pit and it carries a channel, so everything
it generates leaves through the model's own routing. Every flux and store is
a depth over the cell, which is a depth over the catchment; the catchment's
stated area enters only the discharge, `dis = mrro * area_km2 / 86.4`. The
adapter writes that domain and a settings file into /tmp, builds LISFLOOD's
model object through the package's Python API, and runs its dynamic
framework step by step, reading the stores after every step rather than
asking the model to write maps.

What is reported
----------------
Fluxes, as rates in mm per day, from the per-step terms LISFLOOD's own water
balance module (waterbalance.py) closes its budget with:

* `pr`      the forcing, echoed.
* `evspsbl` transpiration + evaporation of intercepted water + soil
            evaporation (TaWB + TaInterceptionWB + ESActWB).
* `mrro`    the channel outflow at the outlet over the step (ChanQAvg * DtSec),
            as a depth over the cell.
* `dis`     `mrro` over the catchment's stated area, in m3/s.
* `gwex`    what leaves the reported stores to the outside, negative: the loss
            from the lower groundwater zone to deep groundwater (GwLossWB,
            zero at LISFLOOD's default GwLoss = 0) and, when the case
            prescribes a withdrawal, the water LISFLOOD's water-use module
            actually took from the lower groundwater zone and from the
            channel (abstraction_GW_actual_M3, withdrawal_CH_actual_M3).

States, absolute, in mm over the cell, weighted by land-use fraction exactly
as waterbalance.py weights them:

* `snw`     SnowCover, the mean of the three elevation zones' packs. LISFLOOD's
            degree-day snow holds no liquid water.
* `canopy`  interception storage (CumInterception), plus sealed-surface
            depression storage (zero: no sealed fraction).
* `mrso`    the three soil layers, W1a + W1b + W2.
* `gw`      the upper and lower groundwater zones, UZ + LZ.
* `channel` overland flow storage (WaterDepth) plus channel water
            (ChanM3 over the cell): generated runoff not yet past the outlet.

The per-step budget over those terms (precipitation in; evaporation, outlet
flow, deep loss and withdrawals out; change in the stores) is recomputed
after the run and its largest value reported in run.json as
`budget_residual_mm_max_abs`. LISFLOOD's own reporting of the same check
(repMBTs) is off: it writes four timeseries through PCRaster every step, a
quarter of the run time, and switching it off leaves every output
bit-identical.

A prescribed withdrawal (`abstr`)
---------------------------------
When the forcing carries an `abstr` column (mm/day, net of return flow), the
water-use option (`wateruse`) is switched on for the run, whatever the column
holds, so the natural variant of a paired probe (a column of zeros) and the
irrigated one run the same configuration. Runs without the column are
configured exactly as before. The withdrawal is LISFLOOD's own: each step the
adapter sets the industrial demand map to `abstr * DtDay` before the
water-use module runs, with the industrial consumptive-use fraction at 1 (the
prescription is already net, so nothing returns). The module then splits the
demand as it always does: `FractionGroundwaterUsed` of it is taken from the
lower groundwater zone with no availability check (LZ can fall below its
threshold, and below zero), and the rest from the channel, limited to the
channel water above the environmental-flow threshold (`EFlowThreshold * DtSec`);
what the channel cannot supply is a shortage LISFLOOD records and does not
take elsewhere. `gwex` carries what was actually taken, never the
prescription, and run.json carries the prescribed, withdrawn and short totals.

How the probe's forcing is fed
------------------------------
The meteorological reader (readmeteo.dynamic) is replaced by one that sets
the same five variables it would set from map stacks, in the same units:
precipitation and the three potential evaporations as depths per step (rate
times DtDay), temperature in degC. LISFLOOD wants three potential
evaporations: ET0 (reference crop, for transpiration), ES0 (bare soil) and
EW0 (open water, used here for evaporation of intercepted water). The probe
gives one potential evaporation, and all three are set to it. LISFLOOD
documents `Tavg` as the daily mean even at sub-daily steps; at PT1H the probe's
hourly temperature is fed as given, because rows are never resampled. DtSec
is the case's step.

The catchment
-------------
Parameters come from three places, in this order of preference, and run.json
records which is which:

1. static.json, where the definition is unambiguous: the latitude (LISFLOOD's
   snowmelt season changes sign with hemisphere), the rain-snow threshold
   (TempSnow), the degree-day factor (SnowMeltCoef), the soil capacity (the
   saturated water content of the three layers: the 50 mm top layer is kept
   and the two lower layers are scaled together) and the canopy capacity (LAI
   chosen so that LISFLOOD's interception capacity SMax = 0.935 + 0.498 LAI -
   0.00575 LAI^2 equals it; the LAI is held constant through the year). The
   area only scales `dis`. `baseflow_coefficient` has no LISFLOOD counterpart
   (LISFLOOD's recession is two linear zones with their own time constants)
   and is not used. The latitude, the area and the canopy capacity are
   required: the adapter stops before it stages anything if one is missing.
   The rain-snow threshold, degree-day factor and soil capacity fall back to
   the sources in 2 and 3 when absent, and run.json records which was used.
2. LISFLOOD's documented defaults in src/lisfloodSettings_reference.xml, for the
   calibration parameters (groundwater time constants, percolation, GwLoss,
   LZThreshold, b_Xinanjiang, PowerPrefFlow, CalChanMan), fixed constants,
   the channel routing sub-step (DtSecChannel 3600 s) and most water-use
   constants.
3. The shipped test catchment (tests/data/LF_ETRS89_UseCase), for maps with no
   lumped counterpart and no documented default: soil hydraulic properties and
   depths, crop coefficient and group, overland Manning's n, hillslope
   gradient and elevation spread (catchment means), channel Manning's n
   (median), the channel's bottom width, bankfull depth and gradient and the
   environmental-flow threshold (medians over the catchment's 961 headwater
   cells, whose upstream area is one cell: these grow with upstream area, and
   the representative cell is a pit that drains only its own 25 km2), and the
   groundwater share of water use (mean).

LISFLOOD's snow adds terms keyed to the day of year on top of the degree-day
factor: a seasonal melt coefficient of +-0.5 mm/degC/day at SnowSeasonAdj = 1,
a summer ice-melt term on any pack left between days 165 and 257, and melt
enhanced by the rain depth of the step. They are the model's own physics,
driven by the forcing's dates, not adapter inputs.

The whole cell is the rainfed "other" land-use fraction. Lakes, reservoirs,
polders, transmission loss, open-water evaporation, split and MCT routing,
variable water fraction, rice and drained irrigation, indicators and land-use
change are switched off: no probe prescribes them and each needs inputs a
synthetic catchment does not have. The soil starts at field capacity, the
groundwater zones, snow, interception and overland flow empty, the channel at
half bankfull (LISFLOOD's default); the probe's spinup year does the rest.

Speed
-----
NUMBA_DISABLE_JIT=1: LISFLOOD's soil and interception loops are numba-jitted
with parallel loops over pixels; on one cell they gain nothing, and compiling
them costs about 40 seconds per run because a read-only container cannot keep
the cache. As Python they give bit-identical output. numexpr runs on one
thread and repMBTs is off, both bit-identical. The store sums are taken on
the numpy arrays under LISFLOOD's vegetation-fraction wrappers.

The model is deterministic. The request seed is recorded and otherwise unused.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import os
import sys
import tempfile
import time
import uuid
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
# numexpr evaluates LISFLOOD's vegetation-fraction expressions on one-cell
# arrays about sixty times a step; a thread pool only adds dispatch cost.
os.environ.setdefault("NUMEXPR_MAX_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np  # noqa: E402

MODEL = {"name": "lisflood", "version": "5.0.0-onecell.5"}
COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "gwex", "mrso", "snw", "canopy", "gw", "channel"]
TIMESTEP_SECONDS = {"PT1D": 86400, "PT1H": 3600, "PT15M": 900, "PT5M": 300, "PT1M": 60}

# The representative cell: the shipped test catchment's grid, from which every
# map-derived parameter below comes (pixleng.nc, pixarea.nc, chanlength.nc).
CELL_LENGTH_M = 5000.0
CELL_AREA_M2 = 25.0e6
CHANNEL_LENGTH_M = 5000.0  # chanlength median; the mean is 5073 m

# First day of each ten-day LAI interval, as leafarea.py numbers its map stack.
LAI_INTERVAL_DAYS = [1, 11, 21, 32, 42, 52, 60, 70, 80, 91, 101, 111, 121, 131, 141, 152, 162, 172,
                     182, 192, 202, 213, 223, 233, 244, 254, 264, 274, 284, 294, 305, 315, 325, 335,
                     345, 355]

# LISFLOOD's documented defaults (src/lisfloodSettings_reference.xml at v5.0.0).
REFERENCE_DEFAULTS = {
    "UpperZoneTimeConstant": 10.0,     # days
    "LowerZoneTimeConstant": 100.0,    # days
    "GwPercValue": 0.5,                # mm/day, UZ -> LZ
    "GwLoss": 0.0,                     # mm/day, LZ -> deep groundwater ("closed lower boundary")
    "LZThreshold": 10.0,               # mm
    "b_Xinanjiang": 0.7,
    "PowerPrefFlow": 3.5,
    "CalChanMan": 2.0,
    "SnowMeltCoef": 4.0,               # mm/degC/day, replaced by static.json when given
    "TempSnow": 1.0,                   # degC, replaced by static.json when given
    "TempMelt": 1.0,                   # degC
    "SnowSeasonAdj": 1.0,
    "SnowFactor": 1.0,
    "TemperatureLapseRate": 0.0065,
    "LeafDrainageTimeConstant": 1.0,   # days
    "kdf": 0.72,
    "AvWaterRateThreshold": 5.0,
    "SMaxSealed": 1.0,
    "Afrost": 0.97,
    "Kfrost": 0.57,
    "SnowWaterEquivalent": 0.45,
    "FrostIndexThreshold": 56.0,
    "beta": 0.6,
    "OFDepRef": 5.0,
    "GradMin": 0.001,
    "ChanGradMin": 0.0001,
    "CourantCrit": 0.4,
    "BankFullPerc": 0.5,
    "DtSecChannel": 3600.0,            # s, channel routing sub-step
    "PrScaling": 1.0,
    "CalEvaporation": 1.0,
    "ChanBottomWMult": 1.0,
    "ChanDepthTMult": 1.0,
    "ChanSMult": 1.0,
}

# tests/data/LF_ETRS89_UseCase at v5.0.0, over its 2847-cell mask: means for
# hillslope and soil maps of the "other" land use; medians for channel geometry,
# over the 961 headwater cells (upstream area one cell) where it grows with
# upstream area.
TEST_CATCHMENT = {
    "ElevationStD": 159.9,             # m, elvstd
    "Grad": 0.2361,                    # gradient
    "SoilDepth1": 50.0,                # mm, soildepth1_o (50 everywhere)
    "SoilDepth2": 748.1,               # mm, soildepth2_o
    "SoilDepth3": 1709.0,              # mm, soildepth3_o
    "MapThetaSat1": 0.4516, "MapThetaSat2": 0.4434, "MapThetaSat3": 0.4236,
    "MapThetaRes1": 0.09877, "MapThetaRes2": 0.1028, "MapThetaRes3": 0.1051,
    "MapLambda1": 0.1586, "MapLambda2": 0.1545, "MapLambda3": 0.1450,
    "MapGenuAlpha1": 0.03203, "MapGenuAlpha2": 0.03455, "MapGenuAlpha3": 0.03793,
    "MapKSat1": 2.815, "MapKSat2": 3.004, "MapKSat3": 2.923,
    "MapCropCoef": 0.9994, "MapCropGroupNumber": 2.692, "MapN": 0.09327,
    "ChanMan": 0.04676,                # ec_chanman, median
    "ChanBottomWidth": 4.344,          # m, ec_chanbw, headwater median (all cells 8.147)
    "ChanDepthThreshold": 0.1646,      # m, ec_chanbnkf, headwater median (all cells 0.331)
    "ChanSdXdY": 1.0,                  # chans (1 everywhere)
    "ChanGrad": 0.06579,               # changrad, headwater median (all cells 0.009729)
}

# The water-use option's inputs, used only when the case prescribes `abstr`.
# (value, source)
WATER_USE = {
    "IndustryConsumptiveUseFraction": (1.0, "packaging: abstr is net of return flow, so nothing returns (reference 0.15)"),
    "IndustrialDemandMaps": (0.0, "set every step to abstr * DtDay before the water-use module runs"),
    "DomesticDemandMaps": (0.0, "no demand but the prescribed one"),
    "LivestockDemandMaps": (0.0, "no demand but the prescribed one"),
    "EnergyDemandMaps": (0.0, "no demand but the prescribed one"),
    "FractionGroundwaterUsed": (0.168, "shipped test catchment mean of fracgwusedNew, the map the reference settings name"),
    "FractionNonConventionalWaterUsed": (0.0, "shipped test catchment (fracncused is 0 everywhere)"),
    "EFlowThreshold": (0.0700, "shipped test catchment: median of ad_dis_nat_10 over its 961 headwater "
                               "cells (upstream area one cell), m3/s; all cells 0.2604, per unit upstream "
                               "area on 25 km2 0.0813"),
    "GroundwaterBodies": (1.0, "shipped test catchment median of ad_gwbodies"),
    "FractionLakeReservoirWaterUsed": (0.25, "LISFLOOD reference default (lakes and reservoirs are off)"),
    "WUsePercRemain": (0.5, "LISFLOOD reference default"),
    "maxNoWateruse": (5.0, "LISFLOOD reference default"),
    "DomesticConsumptiveUseFraction": (0.20, "LISFLOOD reference default"),
    "LivestockConsumptiveUseFraction": (0.15, "LISFLOOD reference default"),
    "EnergyConsumptiveUseFraction": (0.06378, "shipped test catchment mean of energyconsumptiveuse"),
    "LeakageFraction": (0.2, "LISFLOOD reference default"),
    "LeakageWaterLoss": (0.75, "LISFLOOD reference default"),
    "LeakageReductionFraction": (0.0, "LISFLOOD reference default"),
    "WaterSavingFraction": (0.0, "LISFLOOD reference default"),
    "IrrigationEfficiency": (0.75, "LISFLOOD reference default"),
    "ConveyanceEfficiency": (0.80, "LISFLOOD reference default"),
    "IrrigationType": (1.0, "LISFLOOD reference default"),
    "IrrigationMult": (1.20, "LISFLOOD reference default"),
    "IrrigationWaterReUseM3": (0.0, "LISFLOOD reference default"),
    "IrrigationWaterReUseNumDays": (143.0, "LISFLOOD reference default"),
}

OPTIONS_OFF = [
    "wateruse", "TransientWaterDemandChange", "useWaterDemandAveYear", "wateruseRegion",
    "groundwaterSmooth", "drainedIrrigation", "riceIrrigation", "openwaterevapo", "varfractionwater",
    "simulateLakes", "simulateReservoirs", "simulatePolders", "TransLoss", "SplitRouting",
    "MCTRouting", "dynamicWave", "inflow", "indicator", "InitLisflood", "ColdStart",
    "TransientLandUseChange", "simulatePF", "cropsEPIC", "simulateWaterLevels",
    "readNetcdfStack", "writeNetcdfStack", "writeNetcdf",
    "repDischargeTs", "repStateMaps", "repEndMaps", "repDischargeMaps", "repStateUpsGauges",
    "repRateUpsGauges", "repMeteoUpsGauges", "repsimulateLakes", "repsimulateReservoirs",
    "repE2O1", "repE2O2", "repMBTs",
]
OPTIONS_ON = ["gridSizeUserDefined"]


# --- the catchment ------------------------------------------------------------


def lai_for_canopy_capacity(capacity_mm: float) -> float:
    """The LAI at which LISFLOOD's interception capacity equals the catchment's.

    soilloop.py: SMax = 0 for LAI <= 0.1, else 0.935 + 0.498 LAI - 0.00575 LAI^2.
    A capacity below what LAI = 0.1 gives cannot be represented and becomes no
    interception at all.
    """
    if capacity_mm < 0.935 + 0.498 * 0.1 - 0.00575 * 0.01:
        return 0.1
    c = capacity_mm - 0.935
    return (0.498 - math.sqrt(0.498 ** 2 - 4.0 * 0.00575 * c)) / (2.0 * 0.00575)


# static.json keys the adapter cannot run without and will not invent: the
# latitude sets the hemisphere of LISFLOOD's seasonal snowmelt coefficient, the
# area scales the discharge and the canopy capacity sets the LAI. Every probe
# supplies them. The snow and soil keys fall back to LISFLOOD's reference and
# the test catchment's values when absent, and run.json names the source used.
REQUIRED_STATIC = ("latitude_deg", "area_km2", "canopy_capacity_mm")


def require_static(static: dict) -> None:
    missing = [key for key in REQUIRED_STATIC if key not in static]
    if missing:
        raise SystemExit(f"static.json has no {', '.join(missing)}; LISFLOOD needs each of "
                         f"{', '.join(REQUIRED_STATIC)} and the adapter does not invent them")


def catchment_parameters(static: dict) -> tuple[dict, dict]:
    p = dict(REFERENCE_DEFAULTS)
    p.update(TEST_CATCHMENT)
    source = {k: "LISFLOOD reference default" for k in REFERENCE_DEFAULTS}
    source.update({k: "shipped test catchment" for k in TEST_CATCHMENT})

    if "snow_threshold_degC" in static:
        p["TempSnow"] = float(static["snow_threshold_degC"])
        source["TempSnow"] = "static.json snow_threshold_degC"
    if "degree_day_factor_mm_per_C_day" in static:
        p["SnowMeltCoef"] = float(static["degree_day_factor_mm_per_C_day"])
        source["SnowMeltCoef"] = "static.json degree_day_factor_mm_per_C_day"

    ws_top = p["MapThetaSat1"] * p["SoilDepth1"]
    ws_lower = p["MapThetaSat2"] * p["SoilDepth2"] + p["MapThetaSat3"] * p["SoilDepth3"]
    if "soil_capacity_mm" in static:
        scale = max(float(static["soil_capacity_mm"]) - ws_top, 1.0) / ws_lower
        p["SoilDepth2"] *= scale
        p["SoilDepth3"] *= scale
        source["SoilDepth2"] = source["SoilDepth3"] = (
            f"shipped test catchment x {scale:.4f}, so saturated storage equals static.json soil_capacity_mm"
        )

    canopy = float(static["canopy_capacity_mm"])
    p["LAI"] = lai_for_canopy_capacity(canopy)
    source["LAI"] = "static.json canopy_capacity_mm through LISFLOOD's SMax(LAI), constant in time"
    return p, source


def parse_time(text: str) -> dt.datetime:
    value = dt.datetime.fromisoformat(text.strip().replace("Z", "+00:00"))
    return value.replace(tzinfo=None)


def write_domain(work: Path, static: dict, params: dict, forcing: list[dict], timestep: str,
                 water_use: bool) -> Path:
    """One representative cell, a pit with a channel, and the settings file that points at them."""
    import netCDF4
    import pcraster as pcr
    from lisflood.global_modules.add1 import generateName

    latitude = float(static["latitude_deg"])
    maps = work / "maps"
    out = work / "out"
    lai_dir = maps / "lai"
    for d in (maps, out, lai_dir):
        d.mkdir(parents=True, exist_ok=True)

    # The clone is a lat/lon cell centred on the catchment's latitude, which is
    # all LISFLOOD reads from its coordinates (the hemisphere of the snowmelt
    # season). Cell length and area are given explicitly (gridSizeUserDefined).
    cell_deg = 0.1
    pcr.setclone(1, 1, cell_deg, -cell_deg / 2.0, latitude + cell_deg / 2.0)
    one = np.ones((1, 1))

    def pcr_map(name: str, kind, value: float) -> str:
        path = maps / f"{name}.map"
        pcr.report(pcr.numpy2pcr(kind, one * value, -9999), str(path))
        return str(path)

    with netCDF4.Dataset(maps / "template.nc", "w") as nc:
        nc.createDimension("lat", 1)
        nc.createDimension("lon", 1)
        nc.createVariable("lat", "f8", ("lat",))[:] = [latitude]
        nc.createVariable("lon", "f8", ("lon",))[:] = [0.0]
        nc.createVariable("mask", "f4", ("lat", "lon"))[:] = one

    for prefix in ("laio", "laif", "laii"):
        for day in LAI_INTERVAL_DAYS:
            pcr.report(pcr.numpy2pcr(pcr.Scalar, one * params["LAI"], -9999),
                       generateName(str(lai_dir / prefix), day))

    bindings = {
        "MapsCaching": "False", "numCPUs_parallelNumba": "1", "OutputMapsChunks": "1",
        "OutputMapsDataType": "float64", "NetCDFTimeChunks": "auto",
        "CalendarConvention": "proleptic_gregorian",
        "CalendarDayStart": parse_time(forcing[0]["time"]).strftime("%d/%m/%Y %H:%M"),
        "StepStart": "1", "StepEnd": str(len(forcing)), "timestepInit": "1",
        "DtSec": str(TIMESTEP_SECONDS[timestep]), "NumDaysSpinUp": "0",
        "MaskMap": pcr_map("mask", pcr.Boolean, 1),
        "Ldd": pcr_map("ldd", pcr.Ldd, 5),
        "Channels": pcr_map("chan", pcr.Boolean, 1),
        "PixelLengthUser": pcr_map("pixleng", pcr.Scalar, CELL_LENGTH_M),
        "PixelAreaUser": pcr_map("pixarea", pcr.Scalar, CELL_AREA_M2),
        "ChanLength": pcr_map("chanlength", pcr.Scalar, CHANNEL_LENGTH_M),
        "netCDFtemplate": str(maps / "template.nc"),
        "LAIOtherMaps": str(lai_dir / "laio"), "LAIForestMaps": str(lai_dir / "laif"),
        "LAIIrrigationMaps": str(lai_dir / "laii"),
        "OtherFraction": pcr_map("fracother", pcr.Scalar, 1.0),
        "ForestFraction": pcr_map("fracforest", pcr.Scalar, 0.0),
        "IrrigationFraction": pcr_map("fracirrigated", pcr.Scalar, 0.0),
        "RiceFraction": pcr_map("fracrice", pcr.Scalar, 0.0),
        "DirectRunoffFraction": pcr_map("fracsealed", pcr.Scalar, 0.0),
        "WaterFraction": pcr_map("fracwater", pcr.Scalar, 0.0),
        "DrainedFraction": "0",
        # Forest and irrigated soils: required bindings, zero area; given the "other" values.
        "SoilDepth1Forest": repr(params["SoilDepth1"]), "SoilDepth2Forest": repr(params["SoilDepth2"]),
        "SoilDepth3Forest": repr(params["SoilDepth3"]),
        "MapThetaSat1Forest": repr(params["MapThetaSat1"]), "MapThetaSat2Forest": repr(params["MapThetaSat2"]),
        "MapThetaRes1Forest": repr(params["MapThetaRes1"]), "MapThetaRes2Forest": repr(params["MapThetaRes2"]),
        "MapLambda1Forest": repr(params["MapLambda1"]), "MapLambda2Forest": repr(params["MapLambda2"]),
        "MapGenuAlpha1Forest": repr(params["MapGenuAlpha1"]), "MapGenuAlpha2Forest": repr(params["MapGenuAlpha2"]),
        "MapKSat1Forest": repr(params["MapKSat1"]), "MapKSat2Forest": repr(params["MapKSat2"]),
        "MapForestCropCoef": repr(params["MapCropCoef"]), "MapForestCropGroupNumber": repr(params["MapCropGroupNumber"]),
        "MapForestN": repr(params["MapN"]),
        "MapIrrigationCropCoef": repr(params["MapCropCoef"]),
        "MapIrrigationCropGroupNumber": repr(params["MapCropGroupNumber"]),
        # Initial state: soil at field capacity (-9999), everything else empty,
        # channel at half bankfull (-9999, LISFLOOD's default).
        "OFDirectInitValue": "0", "OFOtherInitValue": "0", "OFForestInitValue": "0",
        "SnowCoverAInitValue": "0", "SnowCoverBInitValue": "0", "SnowCoverCInitValue": "0",
        "FrostIndexInitValue": "0",
        "CumIntInitValue": "0", "CumIntForestInitValue": "0", "CumIntIrrigationInitValue": "0",
        "CumIntSealedInitValue": "0",
        "UZInitValue": "0", "UZForestInitValue": "0", "UZIrrigationInitValue": "0",
        "DSLRInitValue": "1", "DSLRForestInitValue": "1", "DSLRIrrigationInitValue": "1",
        "LZInitValue": "0", "LZAvInflowMap": "0", "TotalCrossSectionAreaInitValue": "-9999",
        "ThetaInit1Value": "-9999", "ThetaInit2Value": "-9999", "ThetaInit3Value": "-9999",
        "ThetaForestInit1Value": "-9999", "ThetaForestInit2Value": "-9999", "ThetaForestInit3Value": "-9999",
        "ThetaIrrigationInit1Value": "-9999", "ThetaIrrigationInit2Value": "-9999",
        "ThetaIrrigationInit3Value": "-9999",
        "PrevDischarge": "-9999", "PrevDischargeAvg": "-9999", "CumQInit": "0",
        # Read only on code paths switched off; present so no branch fails on a missing key.
        "CrossSection2AreaInitValue": "-9999", "PrevSideflowInitValue": "-9999", "AvgDis": "-9999",
        "QSplitMult": "2.0", "CalChanMan2": "3.0", "WaterDepthInitValue": "0", "DefineEndofYear": "304",
        "LZInflowCUMInit": "0", "TimeSinceStartPrerunChunkInit": "0",
        "SeepTopToSubBAverageOtherMap": "-9999", "SeepTopToSubBAverageForestMap": "-9999",
        "SeepTopToSubBAverageIrrigationMap": "-9999", "cumSeepTopToSubBOtherInit": "0",
        "cumSeepTopToSubBForestInit": "0", "cumSeepTopToSubBIrrigationInit": "0",
        # Never read: the meteorological reader is replaced (see simulate).
        "PrecipitationMaps": str(maps / "pr"), "TavgMaps": str(maps / "ta"), "ET0Maps": str(maps / "et"),
        "ES0Maps": str(maps / "es"), "E0Maps": str(maps / "e0"),
    }
    for key, value in params.items():
        if key != "LAI":
            bindings.setdefault(key, repr(float(value)))
    options_off = list(OPTIONS_OFF)
    options_on = list(OPTIONS_ON)
    if water_use:
        for key, (value, _) in WATER_USE.items():
            bindings[key] = repr(float(value))
        # loadmap("WUseRegion").astype(int) needs a map, not a constant: one region.
        bindings["WUseRegion"] = pcr_map("wregion", pcr.Nominal, 1)
        options_off.remove("wateruse")
        options_on.append("wateruse")

    lines = ['<?xml version="1.0" encoding="UTF-8"?>', "<lfsettings>", "<lfoptions>"]
    lines += [f'  <setoption choice="0" name="{name}"/>' for name in options_off]
    lines += [f'  <setoption choice="1" name="{name}"/>' for name in options_on]
    lines += ["</lfoptions>", "<lfuser>",
              f'  <textvar name="PathOut" value="{out}"/>',
              '  <textvar name="ReportSteps" value="1..9999"/>',
              '  <textvar name="FilterSteps" value="0"/>',
              '  <textvar name="EnsMembers" value="1"/>',
              '  <textvar name="nrCores" value="1"/>',
              "</lfuser>", "<lfbinding>"]
    lines += [f'  <textvar name="{k}" value="{v}"/>' for k, v in bindings.items()]
    lines += ["</lfbinding>", "</lfsettings>"]
    path = work / "settings.xml"
    path.write_text("\n".join(lines) + "\n")
    return path


# --- the model ----------------------------------------------------------------


def first(value) -> float:
    """The one cell's value of a LISFLOOD pixel array (or a scalar)."""
    return float(np.asarray(value).reshape(-1)[0])


def plain(value) -> np.ndarray:
    """The numpy array under LISFLOOD's vegetation-fraction wrapper, or the array itself."""
    return np.asarray(getattr(value, "values", value))


def storage_terms(m, veg_axis: int) -> tuple[float, float, float, float, float]:
    """Soil, snow, canopy, groundwater, channel: waterbalance.py's stores, in mm over the cell,
    weighted by land-use fraction as that module weights them, on plain numpy arrays."""
    sf = plain(m.SoilFraction)
    soil = np.sum(sf * (plain(m.W1a) + plain(m.W1b) + plain(m.W2)), axis=veg_axis)
    canopy = np.sum(sf * plain(m.CumInterception), axis=veg_axis) + plain(m.DirectRunoffFraction) * plain(m.CumInterSealed)
    groundwater = np.sum(sf * plain(m.UZ), axis=veg_axis) + plain(m.LZ)
    channel = plain(m.ChanM3) * plain(m.M3toMM) + plain(m.WaterDepth)
    return first(soil), first(m.SnowCover), first(canopy), first(groundwater), first(channel)


def simulate(forcing: list[dict], static: dict, timestep: str) -> tuple[list[dict], dict]:
    if timestep not in TIMESTEP_SECONDS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    require_static(static)
    started = time.monotonic()
    params, sources = catchment_parameters(static)
    water_use = "abstr" in forcing[0]
    work = Path(tempfile.mkdtemp(prefix="lisflood-"))
    settings_path = write_domain(work, static, params, forcing, timestep, water_use)

    from lisflood.global_modules.settings import CDFFlags, LisSettings, MaskInfo
    from lisflood.global_modules.zusatz import DynamicFramework
    from lisflood.main import LisfloodModel

    settings = LisSettings(str(settings_path), ["-v"])
    CDFFlags(uuid.uuid4())
    dt_day = TIMESTEP_SECONDS[timestep] / 86400.0
    pr = np.array([float(r["pr"]) for r in forcing])
    tas = np.array([float(r["tas"]) for r in forcing])
    pet = np.array([float(r["pet"]) for r in forcing])
    abstr = np.array([float(r["abstr"]) for r in forcing]) if water_use else np.zeros(len(forcing))
    records: list[tuple[float, ...]] = []
    veg_axis = [0]

    class SteppedLisflood(LisfloodModel):
        def dynamic(self):
            super().dynamic()
            evaporation = self.TaWB + self.TaInterceptionWB + self.ESActWB
            outflow_m3s = np.where(self.AtLastPointC, self.ChanQAvg, 0.0)
            m3_to_mm = first(self.M3toMM)
            if water_use:
                from_groundwater = first(self.abstraction_GW_actual_M3) * m3_to_mm
                from_channel = (first(self.withdrawal_CH_actual_M3) + first(self.LakeAbstractionM3)
                                + first(self.ReservoirAbstractionM3)
                                - first(self.returnflow_GwAbs2Channel_M3_routStep) * self.NoRoutSteps) * m3_to_mm
                short = first(self.areatotal_shortage_SW_M3) * m3_to_mm
            else:
                from_groundwater = from_channel = short = 0.0
            records.append((
                first(evaporation), first(outflow_m3s) * self.DtSec * m3_to_mm, first(self.GwLossWB),
                *storage_terms(self, veg_axis[0]), first(self.TotalPrecipitationWB),
                from_groundwater, from_channel, short, first(self.LZ),
            ))

    model = SteppedLisflood()
    veg_axis[0] = model.SoilFraction.dims.index("vegetation")
    initial_storage = sum(storage_terms(model, veg_axis[0]))
    mask = MaskInfo.instance()

    def read_forcing_rows():
        # readmeteo.py's variables, in its units: depths per step, degC.
        i = model.currentTimeStep() - model.firstTimeStep()
        model.Precipitation = mask.in_zero() + pr[i] * dt_day * model.PrScaling
        model.Tavg = mask.in_zero() + tas[i]
        demand = mask.in_zero() + pet[i] * dt_day * model.CalEvaporation
        model.ETRef = demand
        model.ESRef = demand.copy()
        model.EWRef = demand.copy()
        if water_use:
            # The water-use module's own demand input, in its units (mm per step).
            model.IndustrialDemandMM = mask.in_zero() + abstr[i] * dt_day

    model.readmeteo_module.dynamic = read_forcing_rows
    initialised = time.monotonic()
    framework = DynamicFramework(model, firstTimestep=settings.model_steps[0],
                                 lastTimeStep=settings.model_steps[1])
    framework.rquiet = True
    framework.rtrace = False
    framework.run()
    if len(records) != len(forcing):
        raise RuntimeError(f"LISFLOOD ran {len(records)} steps for {len(forcing)} forcing rows")

    # waterbalance.py's budget, per step, from the same terms: precipitation in,
    # evaporation, outlet flow, deep loss and withdrawals out, change in the stores.
    storage = [initial_storage] + [sum(r[3:8]) for r in records]
    residual_max = max(
        abs(r[8] - r[0] - r[1] - r[2] - r[9] - r[10] - (storage[i + 1] - storage[i]))
        for i, r in enumerate(records)
    )

    area_km2 = float(static["area_km2"])
    rows = []
    for step, rec in zip(forcing, records):
        evaporation, outflow_mm, loss, soil, snow, canopy, groundwater, channel, _, gw_take, ch_take, _, _ = rec
        mrro = outflow_mm / dt_day
        rows.append({
            "time": step["time"],
            "pr": step["pr"],
            "evspsbl": evaporation / dt_day,
            "mrro": mrro,
            "dis": mrro * area_km2 / 86.4,
            "gwex": (0.0 - loss - gw_take - ch_take) / dt_day,
            "mrso": soil,
            "snw": snow,
            "canopy": canopy,
            "gw": groundwater,
            "channel": channel,
        })

    notes = {
        "timestep": timestep,
        "dt_seconds": TIMESTEP_SECONDS[timestep],
        "domain": {
            "cells": 1,
            "cell_length_m": CELL_LENGTH_M,
            "cell_area_km2": CELL_AREA_M2 / 1.0e6,
            "channel_length_m": CHANNEL_LENGTH_M,
            "source": "the shipped test catchment's 5 km grid; chanlength median",
            "ldd": "pit with a channel",
            "land_use": "rainfed 'other' fraction 1.0",
            "catchment_area_km2": area_km2,
            "catchment_area_enters": "only dis = mrro * area_km2 / 86.4",
            "latitude_deg": float(static["latitude_deg"]),
            "latitude_source": "static.json latitude_deg; the clone cell is centred on it, and it sets "
                               "the hemisphere of LISFLOOD's seasonal snowmelt coefficient",
        },
        "static_keys_unused": sorted(k for k in ("baseflow_coefficient",) if k in static),
        "options_off": [o for o in OPTIONS_OFF if not (water_use and o == "wateruse")],
        "options_on": OPTIONS_ON + (["wateruse"] if water_use else []),
        "parameters": {k: {"value": v, "source": sources[k]} for k, v in params.items()},
        "forcing": "readmeteo replaced: Precipitation = pr*DtDay, Tavg = tas (hourly at PT1H), "
                   "ET0 = ES0 = EW0 = pet*DtDay",
        "states": {
            "snw": "SnowCover (mean of three elevation zones; no liquid water)",
            "canopy": "CumInterception weighted by fraction, plus sealed depression storage (zero)",
            "mrso": "W1a + W1b + W2 weighted by fraction",
            "gw": "UZ weighted by fraction + LZ",
            "channel": "overland flow storage (WaterDepth) + channel water (ChanM3) over the cell",
        },
        "fluxes": {
            "evspsbl": "TaWB + TaInterceptionWB + ESActWB",
            "mrro": "ChanQAvg * DtSec at the outlet, over the cell",
            "dis": "mrro * area_km2 / 86.4",
            "gwex": "-(GwLossWB + abstraction_GW_actual_M3 + withdrawal_CH_actual_M3 over the cell)",
        },
        "budget_residual_mm_max_abs": residual_max,
        "numba_disable_jit": os.environ.get("NUMBA_DISABLE_JIT"),
        "numexpr_threads": os.environ.get("NUMEXPR_NUM_THREADS"),
        "initialise_seconds": round(initialised - started, 2),
        "run_seconds": round(time.monotonic() - initialised, 2),
    }
    if water_use:
        prescribed = float(np.sum(abstr) * dt_day)
        from_groundwater = float(sum(r[9] for r in records))
        from_channel = float(sum(r[10] for r in records))
        notes["water_use"] = {
            "enabled_because": "the forcing carries an abstr column",
            "demand": "industrial demand map set to abstr * DtDay each step; consumptive fraction 1",
            "sources": "FractionGroundwaterUsed of the demand from the lower groundwater zone (no "
                       "availability check); the rest from channel water above EFlowThreshold * DtSec",
            "parameters": {k: {"value": v, "source": s} for k, (v, s) in WATER_USE.items()},
            "prescribed_mm": prescribed,
            "withdrawn_from_groundwater_mm": from_groundwater,
            "withdrawn_from_channel_mm": from_channel,
            "channel_shortage_mm": float(sum(r[11] for r in records)),
            "share_withdrawn": (from_groundwater + from_channel) / prescribed if prescribed > 0 else None,
            "lz_min_mm": float(min(r[12] for r in records)),
            "lz_min_note": "LISFLOOD's lower groundwater zone LZ at the end of each step, over the whole record",
        }
    return rows, notes


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet", "abstr"):
            if key in row:
                row[key] = float(row[key])
    return rows


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

    rows, notes = simulate(forcing, static, str(request["timestep"]))
    notes["seed"] = request.get("seed")

    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
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

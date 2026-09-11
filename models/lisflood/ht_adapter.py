#!/usr/bin/env python3
"""HydroTuring adapter for LISFLOOD 5.0.0 (EC Joint Research Centre).

LISFLOOD is a distributed model: a run is one XML settings file, a mask, a
local drain direction map, and parameter maps on that grid. A lumped probe
case is given here as the smallest valid domain LISFLOOD accepts: one square
cell whose area is the catchment's, whose drain direction is a pit, and
which carries a channel, so everything the cell generates leaves through the
model's own routing. The adapter writes that domain and a settings file into
/tmp, builds LISFLOOD's model object through the package's Python API, and
runs its dynamic framework step by step, reading the stores after every step
rather than asking the model to write maps.

What is reported
----------------
Fluxes, as rates in mm per day, from the per-step terms LISFLOOD's own water
balance module (waterbalance.py) closes its budget with:

* `pr`      the forcing, echoed.
* `evspsbl` transpiration + evaporation of intercepted water + soil
            evaporation (TaWB + TaInterceptionWB + ESActWB).
* `mrro`    the channel outflow at the outlet over the step (ChanQAvg * DtSec),
            as a depth over the cell.
* `dis`     the same outflow in m3/s.
* `gwex`    minus the loss from the lower groundwater zone to deep
            groundwater (GwLossWB), an exchange with the outside. Identically
            zero at LISFLOOD's documented default GwLoss = 0.

States, absolute, in mm over the cell, weighted by land-use fraction exactly
as waterbalance.py weights them:

* `snw`     SnowCover, the mean of the three elevation zones' packs. LISFLOOD's
            degree-day snow holds no liquid water.
* `canopy`  interception storage (CumInterception), plus sealed-surface
            depression storage (zero: no sealed fraction).
* `mrso`    the three soil layers, W1a + W1b + W2.
* `gw`      the upper and lower groundwater zones, UZ + LZ.
* `channel` overland flow storage (WaterDepth) plus channel water
            (ChanM3 over the cell area): generated runoff not yet past the
            outlet.

The per-step budget over those same terms (LISFLOOD's own precipitation in,
evaporation, outlet flow and deep loss out, change in the stores) is
recomputed after the run and its largest value reported in run.json as
`budget_residual_mm_max_abs`. LISFLOOD's own reporting of the same check
(repMBTs) is off: it writes four timeseries through PCRaster every step, a
quarter of the run time, and switching it off leaves every output bit-identical.

How the probe's forcing is fed
------------------------------
The meteorological reader (readmeteo.dynamic) is replaced by one that sets
the same five variables it would set from map stacks, in the same units:
precipitation and the three potential evaporations as depths per step (rate
times DtDay), temperature in degC. LISFLOOD wants three potential
evaporations: ET0 (reference crop, for transpiration), ES0 (bare soil) and
EW0 (open water, used here for evaporation of intercepted water). The probe
gives one potential evaporation, and all three are set to it. Rows are fed at
the step the case runs at: DtSec is the case's step and nothing is resampled.

The catchment
-------------
Parameters come from three places, in this order of preference, and run.json
records which is which:

1. static.json, where the definition is unambiguous: the area (cell area and
   length), the latitude (LISFLOOD's snowmelt season changes sign with
   hemisphere), the rain-snow threshold (TempSnow), the degree-day factor
   (SnowMeltCoef), the soil capacity (the saturated water content of the three
   layers: the 50 mm top layer is kept and the two lower layers are scaled
   together) and the canopy capacity (LAI chosen so that LISFLOOD's interception
   capacity SMax = 0.935 + 0.498 LAI - 0.00575 LAI^2 equals it; the LAI is held
   constant through the year).
2. LISFLOOD's documented defaults in src/lisfloodSettings_reference.xml, for the
   calibration parameters (groundwater time constants, percolation, GwLoss,
   LZThreshold, b_Xinanjiang, PowerPrefFlow, CalChanMan) and fixed constants.
3. The shipped test catchment (tests/data/LF_ETRS89_UseCase), for maps with no
   lumped counterpart and no documented default: soil hydraulic properties and
   depths, crop coefficient and group, overland Manning's n, hillslope
   gradient and elevation spread (catchment means), and channel geometry
   (medians, because channel dimensions grow with upstream area and the mean
   is set by the few large-river cells).

The whole cell is the rainfed "other" land-use fraction. Water use, lakes,
reservoirs, polders, transmission loss, open-water evaporation, split and MCT
routing, variable water fraction, rice irrigation and land-use change are
switched off: no probe prescribes them and each needs inputs a synthetic
catchment does not have. The soil starts at field capacity, the groundwater
zones, snow, interception and overland flow empty, the channel at half
bankfull (LISFLOOD's default); the probe's spinup year does the rest.

Speed
-----
Several probes score a ten-year daily record within 60 seconds per run, and
the image is amd64 (conda-forge has no linux-aarch64 PCRaster), so it runs
emulated on an Apple-silicon host. Four settings serve that budget:

* NUMBA_DISABLE_JIT=1. LISFLOOD's soil and interception loops are
  numba-jitted with parallel loops over pixels; on one cell they gain
  nothing, and compiling them costs about 40 seconds per run because a
  read-only container cannot keep the cache. As Python, and with a warm
  cache, they run at the same speed and every output is bit-identical.
* numexpr on one thread, and repMBTs off (above): outputs bit-identical.
* DtSecChannel = 21600 s instead of LISFLOOD's reference 3600 s: channel
  routing four sub-steps a day instead of twenty-four at the daily step, and
  once an hour at the hourly step either way. This one changes the routed
  numbers; README.md compares every probe at both settings.

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

MODEL = {"name": "lisflood", "version": "5.0.0-onecell.2"}

# Channel routing sub-step. LISFLOOD's reference settings use 3600 s; routing
# at 24 sub-steps a day puts a ten-year daily record at about 100 s under amd64
# emulation, beyond the 60 s budget of the probes scored on it. 21600 s (four
# sub-steps a day, LISFLOOD's reference model step) fits. At an hourly step
# the model step is shorter and routing runs once a step either way.
CHANNEL_SUBSTEP_SECONDS = 21600.0
COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "gwex", "mrso", "snw", "canopy", "gw", "channel"]
TIMESTEP_SECONDS = {"PT1D": 86400, "PT1H": 3600, "PT15M": 900, "PT5M": 300, "PT1M": 60}

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
    "DtSecChannel": 3600.0,
    "PrScaling": 1.0,
    "CalEvaporation": 1.0,
    "ChanBottomWMult": 1.0,
    "ChanDepthTMult": 1.0,
    "ChanSMult": 1.0,
}

# tests/data/LF_ETRS89_UseCase at v5.0.0, over its 2847-cell mask: means for
# hillslope and soil maps of the "other" land use, medians for channel geometry.
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
    "ChanBottomWidth": 8.147,          # m, ec_chanbw, median
    "ChanDepthThreshold": 0.331,       # m, ec_chanbnkf, median
    "ChanSdXdY": 1.0,                  # chans (1 everywhere)
    "ChanGrad": 0.009729,              # changrad, median
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


def catchment_parameters(static: dict) -> tuple[dict, dict]:
    p = dict(REFERENCE_DEFAULTS)
    p.update(TEST_CATCHMENT)
    source = {k: "LISFLOOD reference default" for k in REFERENCE_DEFAULTS}
    source.update({k: "shipped test catchment" for k in TEST_CATCHMENT})
    p["DtSecChannel"] = CHANNEL_SUBSTEP_SECONDS
    source["DtSecChannel"] = (
        "packaging choice: 21600 s instead of the reference 3600 s, so a ten-year daily "
        "record fits a 60 s probe budget under amd64 emulation (README.md has the sensitivity)"
    )

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

    canopy = float(static.get("canopy_capacity_mm", 0.935 + 0.498 * 2.84 - 0.00575 * 2.84 ** 2))
    p["LAI"] = lai_for_canopy_capacity(canopy)
    source["LAI"] = "static.json canopy_capacity_mm through LISFLOOD's SMax(LAI), constant in time"
    return p, source


def first(value) -> float:
    """The one cell's value of a LISFLOOD pixel array (or a scalar)."""
    return float(np.asarray(value).reshape(-1)[0])


def storage_terms(m) -> tuple[float, float, float, float, float]:
    """Soil, snow, canopy, groundwater, channel: waterbalance.py's stores, in mm over the cell,
    weighted by land-use fraction as that module weights them."""
    ax = m.SoilFraction.dims.index("vegetation")
    sf = m.SoilFraction
    soil = np.sum(sf * (m.W1a + m.W1b + m.W2), ax)
    canopy = np.sum(sf * m.CumInterception, ax) + m.DirectRunoffFraction * m.CumInterSealed
    groundwater = np.sum(sf * m.UZ, ax) + m.LZ
    channel = m.ChanM3 * m.M3toMM + m.WaterDepth
    return first(soil), first(m.SnowCover), first(canopy), first(groundwater), first(channel)


def parse_time(text: str) -> dt.datetime:
    value = dt.datetime.fromisoformat(text.strip().replace("Z", "+00:00"))
    return value.replace(tzinfo=None)


def write_domain(work: Path, static: dict, params: dict, forcing: list[dict], timestep: str) -> Path:
    """One cell, a pit with a channel, and the settings file that points at them."""
    import netCDF4
    import pcraster as pcr
    from lisflood.global_modules.add1 import generateName

    area_m2 = float(static["area_km2"]) * 1.0e6
    length_m = math.sqrt(area_m2)
    latitude = float(static.get("latitude_deg", 45.0))
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
        "PixelLengthUser": pcr_map("pixleng", pcr.Scalar, length_m),
        "PixelAreaUser": pcr_map("pixarea", pcr.Scalar, area_m2),
        "ChanLength": pcr_map("chanlength", pcr.Scalar, length_m),
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

    lines = ['<?xml version="1.0" encoding="UTF-8"?>', "<lfsettings>", "<lfoptions>"]
    lines += [f'  <setoption choice="0" name="{name}"/>' for name in OPTIONS_OFF]
    lines += [f'  <setoption choice="1" name="{name}"/>' for name in OPTIONS_ON]
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


def simulate(forcing: list[dict], static: dict, timestep: str) -> tuple[list[dict], dict]:
    if timestep not in TIMESTEP_SECONDS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    started = time.monotonic()
    params, sources = catchment_parameters(static)
    work = Path(tempfile.mkdtemp(prefix="lisflood-"))
    settings_path = write_domain(work, static, params, forcing, timestep)

    from lisflood.global_modules.settings import CDFFlags, LisSettings, MaskInfo
    from lisflood.global_modules.zusatz import DynamicFramework
    from lisflood.main import LisfloodModel

    settings = LisSettings(str(settings_path), ["-v"])
    CDFFlags(uuid.uuid4())
    dt_day = TIMESTEP_SECONDS[timestep] / 86400.0
    pr = np.array([float(r["pr"]) for r in forcing])
    tas = np.array([float(r["tas"]) for r in forcing])
    pet = np.array([float(r["pet"]) for r in forcing])
    records: list[tuple[float, ...]] = []

    class SteppedLisflood(LisfloodModel):
        def dynamic(self):
            super().dynamic()
            evaporation = self.TaWB + self.TaInterceptionWB + self.ESActWB
            outflow_m3s = np.where(self.AtLastPointC, self.ChanQAvg, 0.0)
            records.append((
                first(evaporation), first(outflow_m3s * self.DtSec * self.M3toMM), first(outflow_m3s),
                first(self.GwLossWB), *storage_terms(self), first(self.TotalPrecipitationWB),
            ))

    model = SteppedLisflood()
    initial_storage = sum(storage_terms(model))
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
    # evaporation, outlet flow and deep loss out, change in the stores.
    storage = [initial_storage] + [sum(r[4:9]) for r in records]
    residual_max = max(
        abs(r[9] - r[0] - r[1] - r[3] - (storage[i + 1] - storage[i])) for i, r in enumerate(records)
    )

    area_km2 = float(static["area_km2"])
    rows = []
    for step, rec in zip(forcing, records):
        evaporation, outflow_mm, outflow_m3s, loss, soil, snow, canopy, groundwater, channel, _ = rec
        rows.append({
            "time": step["time"],
            "pr": step["pr"],
            "evspsbl": evaporation / dt_day,
            "mrro": outflow_mm / dt_day,
            "dis": outflow_m3s,
            "gwex": (0.0 - loss) / dt_day,  # a loss, so negative; 0.0 - 0.0 avoids writing -0.0
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
            "cell_area_km2": area_km2,
            "cell_length_m": math.sqrt(area_km2 * 1.0e6),
            "ldd": "pit with a channel",
            "land_use": "rainfed 'other' fraction 1.0",
        },
        "options_off": OPTIONS_OFF,
        "options_on": OPTIONS_ON,
        "parameters": {k: {"value": v, "source": sources.get(k, "")} for k, v in params.items()},
        "forcing": "readmeteo replaced: Precipitation = pr*DtDay, Tavg = tas, ET0 = ES0 = EW0 = pet*DtDay",
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
            "gwex": "-GwLossWB (loss from the lower zone to deep groundwater)",
        },
        "budget_residual_mm_max_abs": residual_max,
        "channel_substep_seconds": params["DtSecChannel"],
        "numba_disable_jit": os.environ.get("NUMBA_DISABLE_JIT"),
        "numexpr_threads": os.environ.get("NUMEXPR_NUM_THREADS"),
        "initialise_seconds": round(initialised - started, 2),
        "run_seconds": round(time.monotonic() - initialised, 2),
    }
    return rows, notes


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet"):
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

    rows, notes = simulate(forcing, static, str(request.get("timestep", "PT1D")))
    notes["seed"] = int(request.get("seed", 0))

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

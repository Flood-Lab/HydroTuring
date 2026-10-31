#!/usr/bin/env python3
"""HydroTuring adapter for PRMS 6.0.0 (USGS, CC0).

Runs the official Linux executable from prms_6.0.0_linux.zip as one HRU draining
to one stream segment. It writes the control, parameter and climate-by-HRU files
PRMS reads, runs it once, reads its statvar file and maps it onto the /io contract.
No PRMS code is changed or reimplemented. README.md lists every parameter and what
every choice moves; this docstring says only what the code does.

Inputs PRMS gets
----------------
Through its climate_hru module, the forcing rows directly:
  precip   the row's depth in mm (precip_units = mm)
  potet    the row's pet as a depth in inches, so PRMS's own PET modules do not run
  tmax, tmin  MOCKED: both equal to tas (no diurnal range), because the probes give
              one daily temperature. HT_PRMS_DTR=<degC> sets a symmetric range.
  swrad    MOCKED from the row's pet alone: net radiation by inverting Priestley-Taylor
           at a fixed 20 C (the same mock as models/summa), carried through an albedo
           of 0.23, in Langleys per row. Only snowcomp reads it. HT_PRMS_SOLRAD=ddsolrad
           uses PRMS's own degree-day module instead (calendar and tmax, so not a mock
           and not the evaluated configuration).

At a step finer than a day (PRMS has none, see below)
-----------------------------------------------------
PRMS stops on any data step other than 24 h. The resolution probe still runs every
model at its next finer step, so hourly rows are fed one per PRMS "day": the depth
per row goes in, and the per-day rate constants are put in the units of the row
(fractions 1-(1-k)^dt; amounts and first-order rates x dt; K_coef / dt). The dates
PRMS sees then advance one day per row, and snowcomp's compiled half-day energy
terms stay per day. The manifest does not declare PT1H for that reason.
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
from datetime import date, datetime, timedelta
from pathlib import Path

MODEL = {"name": "prms", "version": "6.0.0-ht.3"}

NOT_REPORTED = {
    "basin_imperv_stor": "no impervious area (hru_percent_imperv = 0), so identically zero",
    "depression storage": "dprst_flag = 0",
}
COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "gwex", "snm", "sbl",
           "mrso", "snw", "canopy", "gw", "channel"]

PRMS_EXE = os.environ.get("PRMS_EXE", "/opt/prms/bin/prms")
WORK = Path(os.environ.get("PRMS_WORKDIR", "/tmp/prms"))
DIAG = os.environ.get("HT_PRMS_DIAG", "0") == "1"

MM_PER_INCH = 25.4
ACRES_PER_KM2 = 247.10538146717
STEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0, "PT15M": 1.0 / 96.0, "PT5M": 1.0 / 288.0, "PT1M": 1.0 / 1440.0}

# Sagehen Creek (upstream's projects/sagehen/input/sagehen.params): capillary share
# of soil_moist_max + sat_threshold, medians over its 128 HRUs, 3.6759 and 4.5706 in.
CAPILLARY_SHARE = 3.6759 / (3.6759 + 4.5706)
TRUE_AREA = os.environ.get("HT_PRMS_TRUE_AREA", "0") == "1"
REFERENCE_AREA_KM2 = 1.0
PRINT_QUANTUM_IN = 1.0e-6  # statvar writes doubles with C's %lf
GWSTOR_INIT_IN = 2.0  # gwstor_init, PRMS's declared default
COVDEN_DEFAULT = 0.5  # covden_sum and covden_win, PRMS's declared defaults
ALBEDO_REF = 0.23
PT_ALPHA = 1.26
PT_TREF_C = 20.0
LANGLEY_PER_WM2_DAY = 86400.0 / 41840.0

# The alternative parameter set: every value is upstream's Sagehen example (median
# over HRUs where it varies). Used only for the README's sensitivity table.
SAGEHEN = {
    "soil_rechr_max_frac": 0.568, "ssr2gw_rate": 0.0499, "ssr2gw_exp": 0.75,
    "fastcoef_lin": 0.4, "pref_flow_den": 0.1, "smidx_coef": 0.0004,
    "gwflow_coef": 0.0037, "potet_sublim": 0.75, "snowinfil_max": 2.75,
    "rad_trncf": 0.2331, "melt_force": 90,
}
PARAMETER_SET = os.environ.get("HT_PRMS_PARAMS", "default")
DTR_C = float(os.environ.get("HT_PRMS_DTR", "0"))
SOLRAD = os.environ.get("HT_PRMS_SOLRAD", "mock")
K_COEF_H = float(os.environ.get("HT_PRMS_K", "1.0"))  # PRMS's declared default
AET_FLAG = int(os.environ.get("HT_PRMS_AET_FLAG", "1"))  # upstream's sagehen controls set 1

# (name, statVar_element)
STAT_VARS = [
    ("basin_ppt", 1), ("basin_rain", 1), ("basin_snow", 1), ("basin_actet", 1),
    ("basin_stflow_in", 1), ("basin_stflow_out", 1), ("basin_cfs", 1),
    ("basin_soil_moist", 1), ("basin_ssstor", 1), ("basin_gwstor", 1), ("basin_gwsink", 1),
    ("basin_intcp_stor", 1), ("basin_imperv_stor", 1), ("basin_pweqv", 1),
    ("basin_segment_storage", 1), ("basin_snowevap", 1), ("basin_intcp_evap", 1),
    ("basin_sroff", 1), ("basin_ssflow", 1), ("basin_gwflow", 1), ("basin_potet", 1), ("basin_gwin", 1),
    ("snowmelt", 1), ("net_rain", 1), ("net_snow", 1), ("pk_precip", 1),
    ("pptmix_nopack", 1), ("intcp_form", 1), ("hru_intcpevap", 1), ("snow_evap", 1),
]


def pt_w_per_mm_day(tc: float, pressure_pa: float = 101325.0) -> float:
    """Net radiation (W m-2) that Priestley-Taylor turns into 1 mm/day at tc."""
    es = 0.6108 * math.exp(17.27 * tc / (tc + 237.3))
    slope = 4098.0 * es / (tc + 237.3) ** 2
    gamma = 0.000665 * pressure_pa / 1000.0
    lam = 2.501e6 - 2361.0 * tc
    return lam * (slope + gamma) / (PT_ALPHA * slope) / 86400.0


def parse_day(text: str) -> date:
    return datetime.fromisoformat(text.strip().replace(" ", "T")).date()


def prms_area_km2(static: dict) -> float:
    """The hru_area PRMS is given: a fixed reference, unless HT_PRMS_TRUE_AREA=1.

    One HRU on one segment has no process that depends on its area. PRMS works in depths
    and multiplies by area only to report volumes, in single precision. Giving it the
    case's area changes nothing but which way a value rounds at statvar's sixth decimal,
    so `dis` alone takes the case's area. README.md has the measurement.
    """
    if TRUE_AREA:
        return float(static.get("area_km2", 250.0))
    return REFERENCE_AREA_KM2


def rate_fraction(k: float, dt: float) -> float:
    return 1.0 - (1.0 - k) ** dt


def build_parameters(static: dict, dt: float, chosen: dict) -> list[tuple[str, str, int, list]]:
    """Every parameter this adapter sets; anything absent is PRMS's declared default."""
    area_acres = prms_area_km2(static) * ACRES_PER_KM2
    total_in = float(static["soil_capacity_mm"]) / MM_PER_INCH if "soil_capacity_mm" in static else None
    # Parameter names here follow PRMS's dimension table: one, nhru, nssr, ngw,
    # nsegment, nmonths. Type 1 is integer and 2 is real.
    p: list[tuple[str, str, int, list]] = []

    def add(name, dim, kind, *values):
        p.append((name, dim, kind, list(values)))

    add("elev_units", "one", 1, 1)
    add("temp_units", "one", 1, 1)
    add("precip_units", "one", 1, 1)
    add("print_type", "one", 1, 0)
    add("print_freq", "one", 1, 0)
    add("runoff_units", "one", 1, 0)
    add("outlet_sta", "one", 1, 0)
    add("hru_area", "nhru", 2, area_acres)
    add("hru_elev", "nhru", 2, float(static.get("elevation_m", 0.0)))
    add("hru_lat", "nhru", 2, float(static.get("latitude_deg", 40.0)))
    add("hru_type", "nhru", 1, 1)
    add("hru_percent_imperv", "nhru", 2, 0.0)

    threshold = float(static.get("snow_threshold_degC", 0.0))
    p.append(("tmax_allsnow", "nmonths", 2, [threshold] * 12))
    p.append(("tmax_allrain_offset", "nmonths", 2, [1.0] * 12))

    if total_in is not None:
        add("soil_moist_max", "nhru", 2, total_in * CAPILLARY_SHARE)
        add("sat_threshold", "nhru", 2, total_in * (1.0 - CAPILLARY_SHARE))

    if "canopy_capacity_mm" in static:
        capacity_in = float(static["canopy_capacity_mm"]) / MM_PER_INCH
        if capacity_in > 0.0:
            # Storage per unit of canopy, so the basin-mean store tops out at the capacity.
            stor = capacity_in / COVDEN_DEFAULT
            add("cov_type", "nhru", 1, 3)
            add("srain_intcp", "nhru", 2, stor)
            add("wrain_intcp", "nhru", 2, stor)
            add("snow_intcp", "nhru", 2, stor)
        else:
            add("cov_type", "nhru", 1, 0)
            add("covden_sum", "nhru", 2, 0.0)
            add("covden_win", "nhru", 2, 0.0)
            add("srain_intcp", "nhru", 2, 0.0)
            add("wrain_intcp", "nhru", 2, 0.0)
            add("snow_intcp", "nhru", 2, 0.0)
            add("rad_trncf", "nhru", 2, 1.0)  # tm6-b7: 0.9917 exp(-2.7557 covden) + 0.0083

    # Rate constants: the declared defaults, or the alternative set, in the units of the row.
    rates = {"gwflow_coef": 0.015, "ssr2gw_rate": 0.1, "slowcoef_lin": 0.015, "fastcoef_lin": 0.1}
    linear_rates = {"slowcoef_sq": 0.1, "fastcoef_sq": 0.8}
    amounts = {"snowinfil_max": 2.0}
    settle = 0.1
    cecn = 5.0
    if PARAMETER_SET == "sagehen":
        for key in rates:
            rates[key] = SAGEHEN.get(key, rates[key])
        amounts["snowinfil_max"] = SAGEHEN["snowinfil_max"]
        add("soil_rechr_max_frac", "nhru", 2, SAGEHEN["soil_rechr_max_frac"])
        add("ssr2gw_exp", "nssr", 2, SAGEHEN["ssr2gw_exp"])
        add("pref_flow_den", "nhru", 2, SAGEHEN["pref_flow_den"])
        add("smidx_coef", "nhru", 2, SAGEHEN["smidx_coef"])
        add("potet_sublim", "nhru", 2, SAGEHEN["potet_sublim"])
        add("melt_force", "nhru", 1, SAGEHEN["melt_force"])
        if float(static.get("canopy_capacity_mm", 1.0)) > 0.0:
            add("rad_trncf", "nhru", 2, SAGEHEN["rad_trncf"])
    elif PARAMETER_SET != "default":
        raise SystemExit(f"HT_PRMS_PARAMS must be default or sagehen, not {PARAMETER_SET!r}")

    dims = {"gwflow_coef": "ngw", "ssr2gw_rate": "nssr"}
    for name, value in rates.items():
        add(name, dims.get(name, "nhru"), 2, rate_fraction(value, dt))
    for name, value in linear_rates.items():
        add(name, "nhru", 2, value * dt)
    for name, value in amounts.items():
        add(name, "nhru", 2, value * dt)
    add("settle_const", "nhru", 2, rate_fraction(settle, dt))
    p.append(("cecn_coef", "nmonths", 2, [cecn * dt] * 12))

    # One segment, one HRU, no observed-flow replacement.
    add("hru_segment", "nhru", 1, 1)
    add("tosegment", "nsegment", 1, 0)
    add("segment_type", "nsegment", 1, 0)
    add("obsin_segment", "nsegment", 1, 0)
    add("segment_flow_init", "nsegment", 2, 0.0)
    add("obsout_segment", "nsegment", 1, 0)
    add("K_coef", "nsegment", 2, min(K_COEF_H / dt, 24.0))
    add("x_coef", "nsegment", 2, 0.2)
    chosen["K_coef_hours_written"] = min(K_COEF_H / dt, 24.0)
    return p


def write_parameter_file(path: Path, parameters, dims: dict) -> None:
    lines = ["HydroTuring one-HRU parameter file", "Version: 1.7", "** Dimensions **"]
    for name, size in dims.items():
        lines += ["####", name, str(size)]
    lines.append("** Parameters **")
    for name, dim, kind, values in parameters:
        lines += ["####", name, "1", dim, str(len(values)), str(kind)]
        lines += [("%d" % v) if kind == 1 else ("%.12g" % v) for v in values]
    path.write_text("\n".join(lines) + "\n")


def write_climate(path: Path, name: str, stamps: list[date], values: list[float]) -> None:
    with open(path, "w") as fh:
        fh.write(f"Written by the HydroTuring adapter\n{name} 1\n########################################\n")
        for d, v in zip(stamps, values):
            fh.write("%d %d %d 0 0 0 %.9g\n" % (d.year, d.month, d.day, v))


def control_block(name: str, kind: int, *values) -> str:
    return "####\n%s\n%d\n%d\n%s\n" % (name, len(values), kind, "\n".join(str(v) for v in values))


def write_control(path: Path, first: date, last: date, solrad: str, transfers: bool) -> None:
    blocks = [
        control_block("model_mode", 4, "PRMS5"),
        control_block("start_time", 1, first.year, first.month, first.day, 0, 0, 0),
        control_block("end_time", 1, last.year, last.month, last.day, 0, 0, 0),
        control_block("data_file", 4, "./input/ht.data"),
        control_block("param_file", 4, "./input/ht.params"),
        control_block("tmax_day", 4, "./input/tmax.day"),
        control_block("tmin_day", 4, "./input/tmin.day"),
        control_block("precip_day", 4, "./input/precip.day"),
        control_block("potet_day", 4, "./input/potet.day"),
        control_block("model_output_file", 4, "./output/ht.out"),
        control_block("temp_module", 4, "climate_hru"),
        control_block("precip_module", 4, "climate_hru"),
        control_block("et_module", 4, "climate_hru"),
        control_block("solrad_module", 4, solrad),
        control_block("transp_module", 4, "transp_tindex"),
        control_block("srunoff_module", 4, "srunoff_smidx"),
        control_block("strmflow_module", 4, "muskingum"),
        control_block("cascade_flag", 1, 0),
        control_block("cascadegw_flag", 1, 0),
        control_block("dprst_flag", 1, 0),
        control_block("soilzone_aet_flag", 1, AET_FLAG),
        control_block("print_debug", 1, -1),
        control_block("parameter_check_flag", 1, 0),
        control_block("statsON_OFF", 1, 1),
        control_block("stat_var_file", 4, "./output/statvar.dat"),
        control_block("nstatVars", 1, len(STAT_VARS)),
        control_block("statVar_element", 4, *[e for _, e in STAT_VARS]),
        control_block("statVar_names", 4, *[n for n, _ in STAT_VARS]),
    ]
    if solrad == "climate_hru":
        blocks.append(control_block("swrad_day", 4, "./input/swrad.day"))
    if transfers:
        blocks.append(control_block("gwr_transferON_OFF", 1, 1))
        blocks.append(control_block("gwr_transfer_file", 4, "./input/gwr_transfer.wu"))
    path.write_text("HydroTuring control file\n" + "".join(blocks))


def read_statvar(path: Path) -> dict[str, list[float]]:
    with open(path) as fh:
        count = int(fh.readline())
        names = [fh.readline().split()[0] for _ in range(count)]
        columns: dict[str, list[float]] = {n: [] for n in names}
        for line in fh:
            parts = line.split()
            if len(parts) != count + 7:
                continue
            for name, value in zip(names, parts[7:]):
                columns[name].append(float(value))
    return columns


def simulate(rows: list[dict], fieldnames: list[str], static: dict, timestep: str):
    started = time.monotonic()
    dt = STEP_DAYS[timestep]
    area_km2 = float(static.get("area_km2", 250.0))
    pr = [float(r["pr"]) for r in rows]
    tas = [float(r["tas"]) for r in rows]
    pet = [float(r["pet"]) for r in rows]
    n = len(rows)

    first = parse_day(rows[0]["time"])
    if timestep == "PT1D":
        stamps = [parse_day(r["time"]) for r in rows]
        for a, b in zip(stamps, stamps[1:]):
            if (b - a).days != 1:
                raise SystemExit(f"PRMS needs consecutive days; {a} is followed by {b}")
    else:
        # One PRMS "day" per row: see the module docstring.
        stamps = [first + timedelta(days=i) for i in range(n)]
    last = stamps[-1]

    if WORK.exists():
        shutil.rmtree(WORK)
    (WORK / "input").mkdir(parents=True)
    (WORK / "output").mkdir()

    chosen: dict = {}
    parameters = build_parameters(static, dt, chosen)
    dims = {"one": 1, "nhru": 1, "nssr": 1, "ngw": 1, "nsegment": 1, "nmonths": 12,
            "nobs": 1, "nrain": 0, "ntemp": 0}

    abstr = [max(float(r["abstr"]), 0.0) for r in rows] if "abstr" in fieldnames else None
    use_transfers = abstr is not None and any(v > 0.0 for v in abstr)
    if use_transfers:
        dims["nwateruse"] = 1
        dims["nexternal"] = 1
    write_parameter_file(WORK / "input" / "ht.params", parameters, dims)

    half_range = DTR_C / 2.0
    write_climate(WORK / "input" / "tmax.day", "tmax", stamps, [t + half_range for t in tas])
    write_climate(WORK / "input" / "tmin.day", "tmin", stamps, [t - half_range for t in tas])
    write_climate(WORK / "input" / "precip.day", "precip", stamps, [p * dt for p in pr])
    write_climate(WORK / "input" / "potet.day", "potet", stamps, [e * dt / MM_PER_INCH for e in pet])
    pt_factor = pt_w_per_mm_day(PT_TREF_C)
    if SOLRAD == "mock":
        swrad = [max(e, 0.0) * pt_factor / (1.0 - ALBEDO_REF) * LANGLEY_PER_WM2_DAY * dt for e in pet]
        write_climate(WORK / "input" / "swrad.day", "swrad", stamps, swrad)
        solrad_module = "climate_hru"
    elif SOLRAD == "ddsolrad":
        solrad_module = "ddsolrad"
    else:
        raise SystemExit(f"HT_PRMS_SOLRAD must be mock or ddsolrad, not {SOLRAD!r}")

    with open(WORK / "input" / "ht.data", "w") as fh:
        fh.write("HydroTuring data file: no observed values\nrunoff 1\n#####\n")
        for d in stamps:
            fh.write("%d %d %d 0 0 0 -999.0\n" % (d.year, d.month, d.day))

    write_control(WORK / "ht.control", first=stamps[0], last=last, solrad=solrad_module, transfers=False)
    staged = time.monotonic()

    def run_prms():
        proc = subprocess.run([PRMS_EXE, "./ht.control"], cwd=WORK, capture_output=True, text=True)
        log = proc.stdout + proc.stderr
        (WORK / "prms.log").write_text(log)
        if proc.returncode != 0 or not (WORK / "output" / "statvar.dat").exists():
            errors = [line for line in log.splitlines() if "ERROR" in line or "storage" in line]
            tail = "\n".join((errors + log.strip().splitlines()[-3:])[-25:])
            raise RuntimeError(f"PRMS failed (exit {proc.returncode}):\n{tail}")

    run_prms()
    removed = [0.0] * n
    if use_transfers:
        # PRMS's water-use module takes `abstr` from the groundwater reservoir to an
        # external location, as a rate in cfs over the PRMS step, and an event holds
        # until the next one. PRMS stops rather than overdraw the reservoir, and that
        # is the only store it can take from besides the stream. So the natural run's
        # recharge (basin_gwin, which a withdrawal from below cannot change) is replayed
        # through gwflow's own order of operations to find what the reservoir can give
        # each day, and PRMS is run again with that. What it could not supply is neither
        # removed nor declared, and run.json says how much that was.
        recharge = read_statvar(WORK / "output" / "statvar.dat")["basin_gwin"]
        coef = dict((name, values[0]) for name, _, _, values in parameters)["gwflow_coef"]
        stored = GWSTOR_INIT_IN
        taken_in = []
        for i in range(n):
            stored += recharge[i]
            want = abstr[i] * dt / MM_PER_INCH
            take = min(want, max(0.99 * stored - 1.0e-5, 0.0))
            taken_in.append(take)
            stored = (stored - take) * (1.0 - coef)
        cfs_per_inch = prms_area_km2(static) * ACRES_PER_KM2 * 43560.0 / 12.0 / 86400.0  # Cfs_conv x hru_area
        with open(WORK / "input" / "gwr_transfer.wu", "w") as fh:
            fh.write("HydroTuring: groundwater withdrawal to an external location\n####\n")
            previous = 0.0
            for d, take in zip(stamps, taken_in):
                if take != previous:
                    fh.write("%d %d %d 1 4 1 %.9g\n" % (d.year, d.month, d.day, take * cfs_per_inch))
                    previous = take
        removed = [take * MM_PER_INCH / dt for take in taken_in]
        write_control(WORK / "ht.control", first=stamps[0], last=last, solrad=solrad_module, transfers=True)
        run_prms()
    ran = time.monotonic()

    log = (WORK / "prms.log").read_text()
    out = read_statvar(WORK / "output" / "statvar.dat")
    if len(out["basin_actet"]) != n:
        raise RuntimeError(f"PRMS wrote {len(out['basin_actet'])} steps for {n} rows\n" + log[-2000:])

    to_flux = MM_PER_INCH / dt  # inches per row to mm/day
    evspsbl = [v * to_flux for v in out["basin_actet"]]
    mrro = [v * to_flux for v in out["basin_stflow_out"]]
    dis = [m * area_km2 / 86.4 for m in mrro]
    sink = [v * to_flux for v in out["basin_gwsink"]]
    gwex = [-(s + removed[i]) for i, s in enumerate(sink)]
    # PRMS holds actual ET at or below potential by construction (soilzone's avail_potet),
    # but statvar prints inches to six decimals, so a day at potential can read one print
    # quantum above the row's pet, and over decades that sums to a ratio just over one.
    # Only a value within that quantum of the row's pet is taken as equal to it.
    snapped = 0
    for i in range(n):
        if abs(out["basin_actet"][i] - pet[i] * dt / MM_PER_INCH) <= PRINT_QUANTUM_IN:
            evspsbl[i] = pet[i]
            snapped += 1
    snm = []
    sbl = []
    for i in range(n):
        delivered = out["snowmelt"][i]
        # srunoff's compute_infil: rain reaches the ground outside the pack on a
        # rain/snow mix with no pack, or on a rain day with no pack to take it.
        if out["pptmix_nopack"][i] > 0.0:
            delivered += out["net_rain"][i]
        elif delivered <= 0.0 and out["basin_pweqv"][i] <= 0.0 and out["net_snow"][i] <= 0.0 \
                and out["net_rain"][i] > 0.0 and out["pk_precip"][i] <= 0.0:
            delivered += out["net_rain"][i]
        snm.append(delivered * to_flux)
        ice = out["snow_evap"][i]
        if out["intcp_form"][i] > 0.0:
            ice += out["hru_intcpevap"][i]
        sbl.append(ice * to_flux)

    result = {
        "evspsbl": evspsbl, "mrro": mrro, "dis": dis, "gwex": gwex, "snm": snm, "sbl": sbl,
        "mrso": [(a + b) * MM_PER_INCH for a, b in zip(out["basin_soil_moist"], out["basin_ssstor"])],
        "snw": [v * MM_PER_INCH for v in out["basin_pweqv"]],
        "canopy": [v * MM_PER_INCH for v in out["basin_intcp_stor"]],
        "gw": [v * MM_PER_INCH for v in out["basin_gwstor"]],
        "channel": [v * MM_PER_INCH for v in out["basin_segment_storage"]],
    }

    # The same quantum at the top of the interception store: a canopy at capacity can print
    # up to half a quantum above it. Only a value within that of the capacity is moved.
    canopy_snapped = 0
    if "canopy_capacity_mm" in static:
        capacity_mm = float(static["canopy_capacity_mm"])
        for i, value in enumerate(result["canopy"]):
            if capacity_mm < value <= capacity_mm + PRINT_QUANTUM_IN * MM_PER_INCH:
                result["canopy"][i] = capacity_mm
                canopy_snapped += 1

    storage = [sum(result[k][i] for k in ("mrso", "snw", "canopy", "gw", "channel")) for i in range(n)]
    residuals = [
        (pr[i] + gwex[i] - evspsbl[i] - mrro[i]) * dt - (storage[i] - storage[i - 1])
        for i in range(1, n)
    ]
    notes = {
        "prms": {"release": "6.0.0 (prms_6.0.0_linux.zip, official executable, unmodified)",
                 "modules": {"temp": "climate_hru", "precip": "climate_hru", "et": "climate_hru",
                             "solrad": solrad_module, "transp": "transp_tindex", "srunoff": "srunoff_smidx",
                             "soilzone": "soilzone", "strmflow": "muskingum", "snow": "snowcomp",
                             "intcp": "intcp", "gw": "gwflow"},
                 "parameter_set": PARAMETER_SET,
                 "soilzone_aet_flag": AET_FLAG},
        "timestep": timestep,
        "calendar": "the rows' own dates" if timestep == "PT1D"
                    else "one PRMS day per row from the first row's date; PRMS has no finer step",
        "rate_constants_per_row": timestep != "PT1D",
        "mock_inputs": {
            "tmax": f"tas + {half_range} C",
            "tmin": f"tas - {half_range} C",
            "swrad": ("pet x Priestley-Taylor conversion at a fixed 20 C / (1 - 0.23), Langleys per row; "
                      "row's forcing alone, no calendar") if SOLRAD == "mock"
                     else "not mocked: PRMS's ddsolrad (tmax and calendar)",
        },
        "states": {
            "mrso": "basin_soil_moist (capillary) + basin_ssstor (gravity and preferential)",
            "snw": "basin_pweqv", "canopy": "basin_intcp_stor", "gw": "basin_gwstor",
            "channel": "basin_segment_storage (Muskingum: cumulative in - out)",
        },
        "evspsbl_rows_snapped_to_pet_within_print_quantum": snapped,
        "canopy_rows_snapped_to_capacity_within_print_quantum": canopy_snapped,
        "fluxes": {"evspsbl": "basin_actet", "mrro": "basin_stflow_out",
                   "dis": "mrro over the case's area_km2; PRMS itself is given a fixed reference hru_area",
                   "gwex": "- basin_gwsink - the abstr PRMS removed through gwr_transfer",
                   "snm": "snowmelt + net_rain where compute_infil passes it around the pack",
                   "sbl": "snow_evap + hru_intcpevap when intcp_form is snow"},
        "water_budget": {
            "max_abs_step_residual_mm": max(abs(v) for v in residuals) if residuals else 0.0,
            "cumulative_residual_mm": sum(residuals),
            "input_mm": sum(p * dt for p in pr),
            "max_imperv_stor_in": max(out["basin_imperv_stor"]),
            "max_basin_gwsink_in": max(out["basin_gwsink"]),
        },
        "abstr": {"applied": bool(use_transfers),
                  "prescribed_mm": sum(abstr) * dt if abstr else 0.0,
                  "removed_mm": sum(removed) * dt,
                  "source": "groundwater reservoir (gwr_transfer), capped at what it holds"},
        "seconds": {"stage": round(staged - started, 3), "prms": round(ran - staged, 3)},
    }
    notes.update(chosen)
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
        reader = csv.DictReader(fh)
        rows, fieldnames = list(reader), list(reader.fieldnames or [])
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    timestep = str(request.get("timestep", "PT1D"))
    if timestep not in STEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")

    result, notes, raw = simulate(rows, fieldnames, static, timestep)

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
        shutil.copy(WORK / "prms.log", io_dir / "prms.log")

    notes["seed"] = int(request.get("seed", 0))
    notes["wall_seconds"] = round(time.monotonic() - started, 3)
    (io_dir / request["output"]["run"]).write_text(json.dumps(
        {"status": "ok", "model": MODEL, "n_steps": len(rows),
         "zero_stores": {"canopy": "the case has no canopy, so PRMS gets no vegetation cover"}
                        if float(static.get("canopy_capacity_mm", 1.0)) == 0.0 else {},
         "not_reported": NOT_REPORTED, "notes": notes}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""HydroTuring adapter for mHM v5.13.4.

The physics is upstream's, compiled unmodified: ht_mhm_driver.f90 is upstream's
driver with mhm_eval's loop written out and one record written per row. This
file writes the namelists and netCDF forcing mHM always reads and turns those
records into the contract's rows. README.md states every choice and how much it
moves results.

Domain. A sub-basin of upstream's test domain, cut at build time by
derive_domain.py. Forcing is uniform over it, MPR runs from its 500 m grids, the
river network is upstream's, and the stated area enters only `dis`.

Step. mHM runs at `timestep = 1` hour, every shipped configuration's value, with
the rows at the case's step as its forcing file: daily rows take its own
disaggregation to hours, hourly rows its hourly-forcing branch. The other value
upstream's namelist names, 24, does not work (README.md, "timestep = 24").

Budget, in mm over the sub-basin. `mrro` is what passes the outlet according to
mRM, and `channel` is what has been generated and not yet passed it.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np

COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "snm", "mrso", "snw", "canopy", "gw", "channel"]
MODEL = {"name": "mhm", "version": "5.13.4-ht.1", "revision": "5f2de41bf244e9bd0d4b7a27806e0c96ff8186db"}

STEPS_PER_ROW = {"PT1D": 24, "PT1H": 1}
DOMAIN = Path("/model/domain")
UPSTREAM = Path("/model/upstream")

# Upstream's own process selection in test_domain/mhm.nml, with PET as an input.
# Case 5 = 0: PET is input, corrected by aspect. Case 8 = 3: adaptive routing with
# slope-dependent celerity. Every geological unit in the shipped table is
# non-karstic, so karstic loss is off and no water leaves the domain below ground.
PROCESSES = {1: 1, 2: 1, 3: 1, 4: 1, 5: 0, 6: 1, 7: 1, 8: 3, 9: 1}

# Soil. mHM's capacity is theta_s * depth summed over horizons, so the depth of
# the second horizon carries the catchment's stated capacity; upstream's first
# horizon, its tillage depth, stays. These two are that sum's slope over this
# sub-basin's unsealed share, mm of saturated soil per mm of each horizon, read
# from the driver's own record with upstream's parameters (README.md).
TILLAGE_MM = 200.0
SATURATED_PER_MM_H1 = 87.70943384 / 200.0
SATURATED_PER_MM_H2 = 176.80420444 / (711.11111111 - 200.0)
# The sealed retention store counts in `mrso`: imperviousStorageCapacity 0.5 mm over
# the sub-basin's 1.68% sealed share.
SEALED_RETENTION_MM = 0.0084
# Upstream's soil table's usual depth, for a case that states no capacity.
DEFAULT_DEPTH_MM = 1000.0
MIN_SECOND_HORIZON_MM = 10.0

# Upstream's day/night ratios are monthly tables. Each is flat here at its annual
# mean, for the reason the LAI table is, so daily rows are split with no calendar:
# rain evenly, PET at upstream's constant 10% by night, temperature +-1.84 degC.
FNIGHT_PREC = 0.5
FNIGHT_PET = 0.1
FNIGHT_TEMP = -1.84

# For the sensitivity runs README.md reports. A case never sets these.
OPT_PARAM = os.environ.get("HT_MHM_PARAMETER_NML", "")
OPT_PROCESS = os.environ.get("HT_MHM_PROCESS", "")
OPT_DOMAIN = os.environ.get("HT_MHM_DOMAIN", "")
OPT_KEEP = os.environ.get("HT_MHM_KEEP", "")
OPT_BINARY = os.environ.get("HT_MHM_BINARY", "/model/ht_mhm_driver")
OPT_FNIGHT = os.environ.get("HT_MHM_FNIGHT_TEMP", "")


def parse_time(text: str) -> datetime:
    return datetime.fromisoformat(text.strip())


def write_forcing(run: Path, name: str, values: np.ndarray, times: np.ndarray, units: str, meta: dict) -> None:
    folder = run / "meteo" / name
    folder.mkdir(parents=True)
    nx, ny, res = meta["ncols_l1"], meta["nrows_l1"], meta["resolution_m"]
    (folder / "header.txt").write_text(
        f"ncols {nx}\nnrows {ny}\nxllcorner {meta['xll']}\nyllcorner {meta['yll']}\n"
        f"cellsize {res}\nNODATA_value -9999\n"
    )
    with netCDF4.Dataset(folder / f"{name}.nc", "w", format="NETCDF4") as nc:
        nc.createDimension("time", len(times))
        nc.createDimension("xc", nx)
        nc.createDimension("yc", ny)
        var = nc.createVariable("time", "i4", ("time",))
        var.units = units
        var.calendar = "standard"
        var[:] = times
        nc.createVariable("xc", "f8", ("xc",))[:] = meta["xll"] + (np.arange(nx) + 0.5) * res
        nc.createVariable("yc", "f8", ("yc",))[:] = meta["yll"] + (ny - 1 - np.arange(ny) + 0.5) * res
        var = nc.createVariable(name, "f8", ("time", "yc", "xc"), fill_value=-9999.0)
        # Uniform forcing: the same value over every cell of the grid.
        var[:] = np.broadcast_to(values[:, None, None], (len(times), ny, nx))


def namelist(run: Path, domain: Path, start: datetime, end: datetime, depth_mm: float, meta: dict) -> str:
    processes = dict(PROCESSES)
    for item in filter(None, OPT_PROCESS.split(",")):
        key, value = item.split("=")
        processes[int(key)] = int(value)
    process_lines = "\n".join(f"processCase({k}) = {v}" for k, v in processes.items())
    res = meta["resolution_m"]
    twelve = lambda v: ", ".join([str(v)] * 12)
    fnight_temp = float(OPT_FNIGHT) if OPT_FNIGHT else FNIGHT_TEMP
    return f"""
&project_description
project_details="HydroTuring case"
setup_description="sub-basin of the mHM test domain, uniform forcing"
simulation_type="benchmark"
Conventions="XXX"
contact="none"
mHM_details="mHM v5.13.4"
history="none"
/
&mainconfig
iFlag_cordinate_sys = 0
nDomains = 1
resolution_Hydrology(1) = {res}
L0Domain(1) = 1
write_restart = .FALSE.
read_opt_domain_data(1) = 0
/
&mainconfig_mhm_mrm
mhm_file_RestartIn(1) = "{run}/mHM_restart_in.nc"
mrm_file_RestartIn(1) = "{run}/mRM_restart_in.nc"
resolution_Routing(1) = {res}
timestep = 1
read_restart = .FALSE.
optimize = .FALSE.
optimize_restart = .FALSE.
opti_method = 1
opti_function = 10
/
&mainconfig_mrm
/
&directories_general
dirConfigOut = "{run}/"
dirCommonFiles = "{domain}/morph/"
dir_Morpho(1) = "{domain}/morph/"
dir_LCover(1) = "{domain}/luse/"
mhm_file_RestartOut(1) = "{run}/mHM_restart.nc"
mrm_file_RestartOut(1) = "{run}/mRM_restart.nc"
dir_Out(1) = "{run}/out/"
file_LatLon(1) = "{domain}/latlon.nc"
/
&directories_mHM
inputFormat_meteo_forcings = "nc"
dir_Precipitation(1) = "{run}/meteo/pre/"
dir_Temperature(1) = "{run}/meteo/tavg/"
dir_ReferenceET(1) = "{run}/meteo/pet/"
time_step_model_inputs(1) = 0
/
&directories_mRM
dir_Gauges(1) = "{run}/gauge/"
dir_Total_Runoff(1) = "{run}/out/"
dir_Bankfull_Runoff(1) = "{run}/"
/
&processSelection
{process_lines}
/
&LCover
nLCoverScene = 1
LCoverYearStart(1) = 1700
LCoverYearEnd(1) = 2600
LCoverfName(1) = "lc_1991.asc"
/
&time_periods
warming_Days(1) = 0
eval_Per(1)%yStart = {start.year}
eval_Per(1)%mStart = {start.month}
eval_Per(1)%dStart = {start.day}
eval_Per(1)%yEnd = {end.year}
eval_Per(1)%mEnd = {end.month}
eval_Per(1)%dEnd = {end.day}
/
&soildata
iFlag_soilDB = 1
tillageDepth = {TILLAGE_MM}
nSoilHorizons_mHM = 2
soil_Depth(1) = {TILLAGE_MM}
soil_Depth(2) = {depth_mm:.6f}
/
&LAI_data_information
timeStep_LAI_input = 0
inputFormat_gridded_LAI = "nc"
/
&LCover_MPR
fracSealed_cityArea = 0.6
/
&evaluation_gauges
nGaugesTotal = 1
NoGauges_domain(1) = 1
Gauge_id(1,1) = {meta["gauge_id"]}
gauge_filename(1,1) = "gauge.txt"
/
&inflow_gauges
nInflowGaugesTotal = 0
NoInflowGauges_domain(1) = 0
InflowGauge_Headwater(1,1) = .FALSE.
/
&panEvapo
evap_coeff = {twelve(1.0)}
/
&nightDayRatio
read_meteo_weights = .FALSE.
fnight_prec = {twelve(FNIGHT_PREC)}
fnight_pet = {twelve(FNIGHT_PET)}
fnight_temp = {twelve(fnight_temp)}
/
&Optimization
nIterations = 7
seed = 1235876
sce_ngs = 2
/
"""


# Nothing through upstream's netCDF writers: the driver's record replaces them.
OUTPUTS = """
&NLoutputResults
output_deflate_level = 1
output_double_precision = .true.
output_time_reference = 0
timeStep_model_outputs = -1
/
"""
OUTPUTS_MRM = """
&NLoutputResults
output_deflate_level_mrm = 1
output_double_precision_mrm = .true.
output_time_reference_mrm = 0
timeStep_model_outputs_mrm = -1
/
"""


def soil_depth(static: dict) -> tuple[float, dict]:
    """The second horizon's lower depth that gives the stated capacity."""
    capacity = static.get("soil_capacity_mm")
    if capacity is None:
        return DEFAULT_DEPTH_MM, {"soil_capacity_mm": None, "second_horizon_lower_depth_mm": DEFAULT_DEPTH_MM}
    capacity = float(capacity)
    depth = TILLAGE_MM + (capacity - SEALED_RETENTION_MM - SATURATED_PER_MM_H1 * TILLAGE_MM) / SATURATED_PER_MM_H2
    # A second horizon has to exist. A capacity this tillage depth already exceeds
    # cannot be met, and the record below says what it came to.
    depth = max(depth, TILLAGE_MM + MIN_SECOND_HORIZON_MM)
    return depth, {"soil_capacity_mm": capacity, "second_horizon_lower_depth_mm": depth}


def simulate(forcing: list[dict], static: dict, timestep: str) -> tuple[list[dict], dict]:
    domain = Path(OPT_DOMAIN) if OPT_DOMAIN else DOMAIN
    meta = json.loads((domain / "domain.json").read_text())
    steps = STEPS_PER_ROW[timestep]
    days = steps / 24.0
    n = len(forcing)
    start = parse_time(forcing[0]["time"])
    pr = np.array([float(row["pr"]) for row in forcing])
    tas = np.array([float(row["tas"]) for row in forcing])
    pet = np.array([float(row["pet"]) for row in forcing])

    if start.hour != 0 or (n * steps) % 24 != 0:
        raise SystemExit("mHM simulates whole days: a case must start at midnight and hold whole days")
    n_days = n * steps // 24
    end = start + timedelta(days=n_days - 1)

    depth_mm, checks = soil_depth(static)

    with tempfile.TemporaryDirectory() as tmp:
        run = Path(tmp)
        (run / "out").mkdir()
        (run / "gauge").mkdir()
        # The gauge mRM requires, all nodata: nothing is calibrated against it, and
        # its discharge at that location is what the driver records.
        with open(run / "gauge" / "gauge.txt", "w") as fh:
            fh.write("HydroTuring no observation DAILY\nnodata -9999\nn 1 measurements per day [1, 1440]\n")
            fh.write(f"start {start.year} {start.month:02d} {start.day:02d} 00 00\n")
            fh.write(f"end {end.year} {end.month:02d} {end.day:02d} 00 00\n")
            for k in range(n_days):
                d = start + timedelta(days=k)
                fh.write(f"{d.year} {d.month:02d} {d.day:02d} 00 00 -9999.000\n")

        if steps == 24:
            times = np.array([(parse_time(row["time"]) - start).days for row in forcing])
            units = f"days since {start:%Y-%m-%d %H:%M:%S}"
        else:
            times = np.array([round((parse_time(row["time"]) - start).total_seconds() / 3600) for row in forcing])
            units = f"hours since {start:%Y-%m-%d %H:%M:%S}"
        # mHM reads a depth per forcing row: mm/day at daily rows, mm/hour at hourly.
        write_forcing(run, "pre", pr * days, times, units, meta)
        write_forcing(run, "tavg", tas, times, units, meta)
        write_forcing(run, "pet", pet * days, times, units, meta)

        (run / "mhm.nml").write_text(namelist(run, domain, start, end, depth_mm, meta))
        (run / "mhm_outputs.nml").write_text(OUTPUTS)
        (run / "mrm_outputs.nml").write_text(OUTPUTS_MRM)
        parameter = OPT_PARAM or str(UPSTREAM / "test_domain" / "mhm_parameter.nml")
        done = subprocess.run(
            [OPT_BINARY, "--nml", str(run / "mhm.nml"), "--parameter", parameter,
             "--mhm_output", str(run / "mhm_outputs.nml"), "--mrm_output", str(run / "mrm_outputs.nml"), str(run)],
            cwd=run, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            env={**os.environ, "HT_STEPS_PER_ROW": str(steps)},
        )
        if OPT_KEEP:
            import shutil
            (run / "mhm_out.txt").write_text(done.stdout)
            shutil.copytree(run, OPT_KEEP, dirs_exist_ok=True)
        if done.returncode != 0:
            raise SystemExit("mhm failed:\n" + "\n".join(done.stdout.splitlines()[-30:]))

        record = np.fromfile(run / "ht_rows.bin", dtype=np.float64).reshape(-1, 11)
        outlet_m3s = np.fromfile(run / "ht_outlet.bin", dtype=np.float64)
        record_meta = (run / "ht_meta.txt").read_text()

    if record.shape[0] != n or outlet_m3s.shape[0] != n:
        raise SystemExit(f"mHM wrote {record.shape[0]} rows for {n} forcing rows")

    soil, sealed, unsat, sat, snow, canopy = record[:, :6].T
    aet, generated, effective, melt, demand = record[:, 6:].T
    area_mhm_km2 = meta["area_km2"]
    out_mm = outlet_m3s * 86.4 * days / area_mhm_km2
    transit = np.cumsum(generated - out_mm)

    area_km2 = float(static["area_km2"])
    rows = [
        {
            "time": forcing[i]["time"],
            "pr": forcing[i]["pr"],
            "evspsbl": aet[i] / days,
            "mrro": out_mm[i] / days,
            "dis": out_mm[i] / days * area_km2 / 86.4,
            "snm": effective[i] / days,
            "mrso": soil[i] + sealed[i],
            "snw": snow[i],
            "canopy": canopy[i],
            "gw": unsat[i] + sat[i],
            "channel": transit[i],
        }
        for i in range(n)
    ]

    storage = soil + sealed + unsat + sat + snow + canopy
    pr_mm = pr * days
    residual = pr_mm[1:] - aet[1:] - generated[1:] - np.diff(storage)
    fields = {}
    for line in record_meta.splitlines():
        parts = line.split()
        fields["_".join(parts[:-1])] = float(parts[-1])
    checks.update({
        "driver": fields,
        "land_budget_residual_share_of_rain": float(residual.sum() / max(pr_mm[1:].sum(), 1e-9)),
        "max_row_land_budget_residual_mm": float(np.abs(residual).max()),
        "transit_mm_min_max": [float(transit.min()), float(transit.max())],
        "pet_after_aspect_correction_over_forcing": float(demand.sum() / max(pet.sum() * days, 1e-9)),
        "generated_runoff_mm": float(generated.sum()),
        "outlet_runoff_mm": float(out_mm.sum()),
        "melt_mm": float(melt.sum()),
    })
    return rows, checks


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


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
    if timestep not in STEPS_PER_ROW:
        raise SystemExit(f"unsupported timestep {timestep!r}: mHM reads daily or hourly forcing and nothing finer")
    rows, checks = simulate(forcing, static, timestep)

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
        "not_reported": {
            "stage": "mRM has no water level",
            "hfls, hfss, hfg and other energy terms": "mHM has no energy balance",
        },
        "notes": {
            "domain": "a 286 km2 sub-basin of upstream's test domain, at 4 km; the stated area enters only dis",
            "mrso": "soil horizons over the unsealed share, plus the retention store of the sealed share",
            "gw": "the unsaturated and saturated reservoirs below the soil horizons, over the unsealed share",
            "channel": "cumulative generated runoff less cumulative outlet discharge from mRM",
            "mrro": "mRM's routed discharge at the outlet, as a depth over the sub-basin",
            "snm": "preEffect: rain and melt leaving the snow module; that rain has already passed interception",
            "seed": "unused; the model is deterministic",
        },
        "checks": checks,
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

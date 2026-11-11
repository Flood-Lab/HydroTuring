#!/usr/bin/env python3
"""HydroTuring adapter for mizuRoute v3.1.1.

The adapter writes what upstream's stand-alone executable reads, a control file,
a river network and a runoff file, runs it, and reads the history file back. It
computes no routing. README.md explains each choice.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np

MODEL = {"name": "mizuroute", "version": "3.1.1-ht.2"}
EXECUTABLE = Path("/model/route_runoff")
PARAM_NML = Path("/model/param.nml.default")
WORK = Path("/tmp/mizu")

STEP_SECONDS = {"PT1H": 3600, "PT15M": 900, "PT5M": 300, "PT1M": 60}

# Upstream's route_opt digits and the prefix of their history variables.
SCHEMES = {"1": "IRF", "2": "KWT", "3": "KW", "4": "MC", "5": "DW"}
# route/settings/SAMPLE.control sets route_opt to 5. The environment variable
# exists for the scheme comparison in README.md; the harness never sets it.
ROUTE_OPT = os.environ.get("HT_MIZU_ROUTE_OPT", "5")

# The headwater reach does no routing under upstream's default hw_drain_point of 2:
# its REACH_Q is its local runoff. Its length enters only KWT, and it is the same
# in both reaches of a paired experiment.
SOURCE_REACH_M = 100.0
SOURCE_SEG, TESTED_SEG = 1, 2


def write_network(path: Path, static: dict) -> None:
    length = float(static["reach_length_m"])
    width = float(static["width_m"])
    slope = float(static["slope"])
    manning_n = float(static["manning_n"])
    area_m2 = float(static["area_km2"]) * 1.0e6
    with netCDF4.Dataset(path, "w", format="NETCDF4_CLASSIC") as nc:
        nc.createDimension("seg", 2)
        nc.createDimension("hru", 1)
        segs = {
            "segId": ("i4", [SOURCE_SEG, TESTED_SEG]),
            "downSegId": ("i4", [TESTED_SEG, 0]),
            "length": ("f8", [SOURCE_REACH_M, length]),
            "slope": ("f8", [slope, slope]),
            "width": ("f8", [width, width]),
            "man_n": ("f8", [manning_n, manning_n]),
            # Vertical banks: the declared section is rectangular.
            "sideSlope": ("f8", [0.0, 0.0]),
        }
        for name, (kind, values) in segs.items():
            nc.createVariable(name, kind, ("seg",))[:] = values
        nc.createVariable("HRUid", "i4", ("hru",))[:] = [1]
        nc.createVariable("hruSegId", "i4", ("hru",))[:] = [SOURCE_SEG]
        nc.createVariable("area", "f8", ("hru",))[:] = [area_m2]


def write_runoff(path: Path, times: list[datetime], q_in: np.ndarray, area_km2: float) -> str:
    """The prescribed discharge in the only form the executable reads: a depth rate.

    q_in is already a discharge. mizuRoute multiplies a runoff depth rate by the
    HRU area to get one, so the depth rate given to it is q_in over that same area,
    in mm/s, and the area cancels; `area_km2` has no effect on `dis`.
    """
    units = f"seconds since {times[0]:%Y-%m-%d %H:%M:%S}"
    seconds = [(t - times[0]).total_seconds() for t in times]
    with netCDF4.Dataset(path, "w", format="NETCDF4_CLASSIC") as nc:
        nc.createDimension("time", None)
        nc.createDimension("hru", 1)
        var = nc.createVariable("time", "f8", ("time",))
        var.units = units
        var.calendar = "standard"
        var[:] = seconds
        nc.createVariable("hru", "i4", ("hru",))[:] = [1]
        runoff = nc.createVariable("runoff", "f8", ("time", "hru"), fill_value=-9999.0)
        runoff.units = "mm/s"
        runoff[:, 0] = q_in / (area_km2 * 1.0e6) * 1000.0
    return units


def write_control(path: Path, times: list[datetime], dt: int, units: str, scheme: str) -> None:
    lines = [
        ("ancil_dir", f"{WORK}/"),
        ("input_dir", f"{WORK}/"),
        ("output_dir", f"{WORK}/output/"),
        ("case_name", "ht"),
        ("sim_start", f"{times[0]:%Y-%m-%d %H:%M:%S}"),
        ("sim_end", f"{times[-1]:%Y-%m-%d %H:%M:%S}"),
        ("route_opt", ROUTE_OPT),
        # q_in is a river inflow, so it does not pass through a hillslope hydrograph.
        ("doesBasinRoute", "0"),
        ("dt_qsim", str(dt)),
        ("is_lake_sim", "F"),
        ("is_flux_wm", "F"),
        ("fname_ntopOld", "network.nc"),
        ("dname_sseg", "seg"),
        ("dname_nhru", "hru"),
        ("fname_qsim", "runoff.nc"),
        ("vname_qsim", "runoff"),
        ("vname_time", "time"),
        ("dname_time", "time"),
        ("vname_hruid", "hru"),
        ("dname_hruid", "hru"),
        ("units_qsim", "mm/s"),
        ("dt_ro", str(dt)),
        ("ro_calendar", "standard"),
        ("ro_time_units", units),
        ("ro_time_stamp", "start"),
        ("is_remap", "F"),
        # Reach width, roughness and side slope come from the network file.
        ("hydGeometryOption", "0"),
        ("topoNetworkOption", "1"),
        ("computeReachList", "1"),
        ("param_nml", "param.nml"),
        ("varname_area", "area"),
        ("varname_length", "length"),
        ("varname_slope", "slope"),
        ("varname_HRUid", "HRUid"),
        ("varname_hruSegId", "hruSegId"),
        ("varname_segId", "segId"),
        ("varname_downSegId", "downSegId"),
        ("varname_width", "width"),
        ("varname_man_n", "man_n"),
        ("varname_sideSlope", "sideSlope"),
        ("newFileFrequency", "single"),
        ("outputFrequency", "1"),
        ("pio_netcdf_type", "netcdf"),
        ("basRunoff", "F"),
        ("dlayRunoff", "F"),
        ("sumUpstreamRunoff", "F"),
        (f"{scheme}routedRunoff", "T"),
        (f"{scheme}volume", "T"),
        ("restart_write", "never"),
    ]
    path.write_text("".join(f"<{key}> {value} !\n" for key, value in lines))


def read_history(path: Path, scheme: str, n_steps: int) -> tuple[np.ndarray, np.ndarray]:
    """Outflow of the declared reach, and the water held in both reaches."""
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_mask(False)
        reach_ids = list(np.asarray(nc.variables["reachID"][:]).ravel())
        tested = reach_ids.index(TESTED_SEG)

        def table(name: str) -> np.ndarray:
            variable = nc.variables[name]
            data = np.asarray(variable[:], dtype=float)
            if variable.dimensions[0] != "time":
                data = data.T
            return data

        routed = table(f"{scheme}routedRunoff")
        volume = table(f"{scheme}volume")
    if routed.shape[0] != n_steps:
        raise RuntimeError(f"history file holds {routed.shape[0]} steps, expected {n_steps}")
    return routed[:, tested], volume.sum(axis=1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent

    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        forcing = list(csv.DictReader(fh))
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    shape = str(static["cross_section_shape"]).strip().lower()
    if shape != "rectangular":
        raise ValueError(f"the adapter maps a rectangular section only, not {shape!r}")
    timestep = request["timestep"]
    dt = STEP_SECONDS[timestep]
    times = [datetime.fromisoformat(row["time"]) for row in forcing]
    if any(b - a != timedelta(seconds=dt) for a, b in zip(times, times[1:])):
        raise ValueError(f"forcing rows are not {timestep} apart")
    q_in = np.array([float(row["q_in"]) for row in forcing])
    area_km2 = float(static["area_km2"])
    scheme = SCHEMES[ROUTE_OPT]

    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "output").mkdir(exist_ok=True)
    shutil.copy(PARAM_NML, WORK / "param.nml")
    write_network(WORK / "network.nc", static)
    units = write_runoff(WORK / "runoff.nc", times, q_in, area_km2)
    write_control(WORK / "control.txt", times, dt, units, scheme)

    proc = subprocess.run(
        [str(EXECUTABLE), str(WORK / "control.txt")],
        cwd=WORK, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    histories = sorted((WORK / "output").glob("ht.h.*.nc"))
    if proc.returncode != 0 or not histories:
        sys.stderr.write(proc.stdout[-4000:])
        raise RuntimeError(f"route_runoff exited {proc.returncode} with {len(histories)} history files")

    dis, volume = read_history(histories[0], scheme, len(forcing))
    # m3 over the same area the inflow was put over, in mm: the catchment-average
    # form `channel` is defined in, from upstream's own volume and nothing else.
    channel = volume / (area_km2 * 1.0e6) * 1000.0

    columns = ["time", "dis", "channel"]
    output = io_dir / request["output"]["table"]
    with open(output, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        for i, row in enumerate(forcing):
            writer.writerow([row["time"], repr(float(dis[i])), repr(float(channel[i]))])

    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "route_opt": ROUTE_OPT,
        "scheme": scheme,
        "dt_qsim_s": dt,
        "network": {
            "source_reach_m": SOURCE_REACH_M,
            "reach_length_m": float(static["reach_length_m"]),
            "hw_drain_point": 2,
            "doesBasinRoute": 0,
        },
        "not_reported": {
            # Outside the control volume, which starts at the inflow: unknown, not zero.
            "snw": "mizuRoute takes q_in and never sees weather, so snow is outside its control volume",
            "canopy": "outside the control volume, as snw",
            "mrso": "outside the control volume, as snw",
            "gw": "outside the control volume, as snw",
        },
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

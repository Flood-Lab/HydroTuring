#!/usr/bin/env python3
"""HydroTuring adapter for RVIC 1.1.1.

Writes the files `rvic parameters` and `rvic convolution` read, and runs upstream's
own routines. Two things are added to the run: a recorder on `Rvar.get_rof`, which
`convolution_run` already calls every step, because RVIC's history tape has no storage
field and the ring is otherwise visible only in a restart file; and the whole-second
rounding in `WholeSecondStep`. README.md explains each choice.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import logging
import os
import sys
import traceback
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np

MODEL = {"name": "rvic", "version": "1.1.1-ht.1"}
WORK = Path("/tmp/rvic")
DT_S = 3600

# Grid spacing of the two cells. RVIC reads the distance from Flow_Distance and not
# from the coordinates, so this value only has to give two distinct cells.
RESOLUTION_DEG = 1.0
LATS = np.array([40.0])
LONS = np.array([-100.0, -100.0 + RESOLUTION_DEG])
# VIC direction code 3 is east: source cell into outlet cell.
FLOW_EAST = 3

# From samples/configs/rvic.parameters.pnw.cfg, except SUBSET_DAYS. The sample's 10 days
# shortens the ring, and that path passes a float as an array shape, which every NumPy
# with a wheel for this Python refuses. Equal to BASIN_FLOWDAYS, finish_params skips the
# subset: nothing is truncated, and the ring is 1200 hours long.
BASIN_FLOWDAYS = 50
CELL_FLOWDAYS = 2
SUBSET_DAYS = BASIN_FLOWDAYS


def normal_depth(q: float, width: float, slope: float, n: float) -> float:
    """Manning's rectangular normal depth, by bisection."""

    def discharge(h: float) -> float:
        area = width * h
        return area * (area / (width + 2.0 * h)) ** (2.0 / 3.0) * slope ** 0.5 / n

    lo, hi = 0.0, 1.0
    while discharge(hi) < q:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if discharge(mid) < q else (lo, mid)
    return 0.5 * (lo + hi)


def routing_parameters(q0: float, static: dict) -> dict:
    """RVIC's one velocity and diffusivity, from the section at one discharge.

    RVIC takes these as constants of the parameter file. Upstream's samples set them
    by hand. Here they are the kinematic celerity dQ/dA and the low-inertia hydraulic
    diffusivity Q/(2BS) of the declared section at the first row's q_in, the only
    discharge known before the run without reading ahead.
    """
    width = float(static["width_m"])
    slope = float(static["slope"])
    n = float(static["manning_n"])
    h = normal_depth(q0, width, slope, n)
    area, perimeter = width * h, width + 2.0 * h
    velocity = q0 * (5.0 / (3.0 * area) - 4.0 / (3.0 * width * perimeter))
    diffusion = q0 / (2.0 * width * slope)
    return {
        "reference_q_m3s": q0,
        "normal_depth_m": h,
        "velocity_m_s": float(os.environ.get("HT_RVIC_VELOCITY", velocity)),
        "diffusion_m2_s": float(os.environ.get("HT_RVIC_DIFFUSION", diffusion)),
    }


def config_number(value: float) -> str:
    """RVIC's config parser takes "1.0" for an integer and then fails to convert it."""
    return str(int(value)) if value == int(value) else repr(value)


def write_inputs(static: dict, times: list, q_in: np.ndarray, params: dict) -> None:
    area_m2 = float(static["area_km2"]) * 1.0e6
    length = float(static["reach_length_m"])

    with netCDF4.Dataset(WORK / "domain.nc", "w", format="NETCDF3_CLASSIC") as nc:
        nc.createDimension("lat", len(LATS))
        nc.createDimension("lon", len(LONS))
        nc.createVariable("lat", "f8", ("lat",))[:] = LATS
        nc.createVariable("lon", "f8", ("lon",))[:] = LONS
        nc.createVariable("mask", "i4", ("lat", "lon"))[:] = 1
        nc.createVariable("frac", "f8", ("lat", "lon"))[:] = 1.0
        area = nc.createVariable("area", "f8", ("lat", "lon"))
        area.units = "m2"
        area[:] = area_m2

    with netCDF4.Dataset(WORK / "routing.nc", "w", format="NETCDF3_CLASSIC") as nc:
        nc.createDimension("lat", len(LATS))
        nc.createDimension("lon", len(LONS))
        nc.createVariable("lat", "f8", ("lat",))[:] = LATS
        nc.createVariable("lon", "f8", ("lon",))[:] = LONS
        nc.createVariable("Flow_Direction", "i4", ("lat", "lon"))[:] = [[FLOW_EAST, 0]]
        nc.createVariable("Basin_ID", "i4", ("lat", "lon"))[:] = 1
        # make_grid_uh gives a cell uh_box convolved with the river response of the
        # cell it drains to, which at the outlet is the Green's function at the
        # outlet's own Flow_Distance. So the outlet carries the reach length, and the
        # source cell's value enters nothing it delivers to this outlet.
        distance = nc.createVariable("Flow_Distance", "f8", ("lat", "lon"))
        distance.units = "m"
        distance[:] = length

    (WORK / "pour_points.csv").write_text(f"lons,lats\n{LONS[1]},{LATS[0]}\n")
    # An impulse in the first hour: q_in is a river inflow, so no hillslope delay.
    (WORK / "uh_box.csv").write_text(f"time,uh\n0,1.0\n{DT_S},0.0\n")

    # RVIC's clock adds the step to a float number of days since year 1, and it stops
    # if a forcing time differs from that clock by a microsecond. So the forcing times
    # are written as that same float sum, in RVIC's own units; they are the row times.
    from rvic.core.share import SECSPERDAY, TIMEUNITS

    clock = [float(netCDF4.date2num(times[0], TIMEUNITS, calendar="standard"))]
    for _ in times[1:]:
        clock.append(clock[-1] + float(DT_S) / float(SECSPERDAY))
    with netCDF4.Dataset(WORK / "forcing.nc", "w", format="NETCDF3_CLASSIC") as nc:
        nc.createDimension("time", None)
        nc.createDimension("lat", len(LATS))
        nc.createDimension("lon", len(LONS))
        nc.createVariable("lat", "f8", ("lat",))[:] = LATS
        nc.createVariable("lon", "f8", ("lon",))[:] = LONS
        var = nc.createVariable("time", "f8", ("time",))
        var.units = TIMEUNITS
        var.calendar = "standard"
        var[:] = clock
        runoff = nc.createVariable("inflow", "f8", ("time", "lat", "lon"))
        runoff.units = "kg m-2 s-1"
        field = np.zeros((len(times), len(LATS), len(LONS)))
        # q_in is a discharge. RVIC routes a flux per unit area and multiplies back by
        # the same cell area on output, so it is put over that area and the area cancels.
        field[:, 0, 0] = q_in / area_m2 * 1000.0
        runoff[:] = field

    (WORK / "parameters.cfg").write_text(
        "[OPTIONS]\n"
        "LOG_LEVEL:INFO\nVERBOSE:False\nCLEAN:False\n"
        "CASEID:ht\nGRIDID:ht\n"
        f"CASE_DIR:{WORK}/params_case\nTEMP_DIR:{WORK}/params_case/temp\n"
        "REMAP:False\nAGGREGATE:False\nAGG_PAD:0\n"
        "NETCDF_FORMAT:NETCDF3_CLASSIC\nNETCDF_ZLIB:False\nNETCDF_COMPLEVEL:4\nNETCDF_SIGFIGS:None\n"
        f"SUBSET_DAYS:{SUBSET_DAYS}\nCONSTRAIN_FRACTIONS:False\nSEARCH_FOR_CHANNEL:False\n"
        f"[POUR_POINTS]\nFILE_NAME:{WORK}/pour_points.csv\n"
        f"[UH_BOX]\nFILE_NAME:{WORK}/uh_box.csv\nHEADER_LINES=1\n"
        f"[ROUTING]\nFILE_NAME:{WORK}/routing.nc\n"
        "LONGITUDE_VAR:lon\nLATITUDE_VAR:lat\nFLOW_DISTANCE_VAR:Flow_Distance\n"
        "FLOW_DIRECTION_VAR:Flow_Direction\nBASIN_ID_VAR:Basin_ID\n"
        f"VELOCITY:{config_number(params['velocity_m_s'])}\n"
        f"DIFFUSION:{config_number(params['diffusion_m2_s'])}\n"
        f"OUTPUT_INTERVAL:{DT_S}\nBASIN_FLOWDAYS:{BASIN_FLOWDAYS}\nCELL_FLOWDAYS:{CELL_FLOWDAYS}\n"
        f"[DOMAIN]\nFILE_NAME:{WORK}/domain.nc\n"
        "LONGITUDE_VAR:lon\nLATITUDE_VAR:lat\nLAND_MASK_VAR:mask\nFRACTION_VAR:frac\nAREA_VAR:area\n"
    )


def convolution_config(param_file: str, times: list) -> str:
    path = WORK / "convolution.cfg"
    path.write_text(
        "[OPTIONS]\n"
        "LOG_LEVEL:INFO\nVERBOSE:False\n"
        f"CASE_DIR:{WORK}/conv_case/\nCASEID:ht\nCASESTR:hydroturing\nCALENDAR:standard\n"
        f"RUN_TYPE:drystart\nRUN_STARTDATE:{times[0]:%Y-%m-%d-%H}\n"
        f"STOP_OPTION:nsteps\nSTOP_N:{len(times)}\nSTOP_DATE:-999\n"
        # RVIC has no "never": a restart due after the last step is not written.
        f"REST_OPTION:nsteps\nREST_N:{len(times) + 10}\nREST_DATE:-999\nREST_NCFORM:NETCDF3_CLASSIC\n"
        "[HISTORY]\n"
        f"RVICHIST_NTAPES=1\nRVICHIST_MFILT:{len(times) + 1}\nRVICHIST_NDENS:2\nRVICHIST_NHTFRQ:1\n"
        "RVICHIST_AVGFLAG=I\nRVICHIST_OUTTYPE:array\nRVICHIST_NCFORM:NETCDF3_CLASSIC\n"
        "RVICHIST_NETCDF_ZLIB:False\nRVICHIST_NETCDF_COMPLEVEL:4\nRVICHIST_NETCDF_SIGFIGS:None\n"
        "RVICHIST_UNITS:m3/s\n"
        f"[DOMAIN]\nFILE_NAME:{WORK}/domain.nc\n"
        "LONGITUDE_VAR:lon\nLATITUDE_VAR:lat\nAREA_VAR:area\nLAND_MASK_VAR:mask\nFRACTION_VAR:frac\n"
        "[INITIAL_STATE]\nFILE_NAME:None\n"
        f"[PARAM_FILE]\nFILE_NAME:{param_file}\n"
        f"[INPUT_FORCINGS]\nDATL_PATH:{WORK}/\nDATL_FILE:forcing.nc\n"
        "TIME_VAR:time\nLATITUDE_VAR:lat\nDATL_LIQ_FLDS:inflow\nSTART:\nEND:\n"
    )
    return str(path)


def run_rvic(times: list) -> tuple:
    """Upstream's parameter step, then upstream's convolution, recording the ring."""
    import rvic.convolution as upstream
    from rvic.convolution import convolution_final, convolution_init, convolution_run
    from rvic.parameters import parameters

    class WholeSecondStep(upstream.DataModel):
        """Upstream's forcing reader, with its step rounded to the whole second.

        It takes the step from two times in days since year 1, which comes out as
        3599.9999966 s at this date, and then refuses to run because that is not a
        divisor of the impulse response's 3600 s. Nothing else is changed.
        """

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.secs_per_step = float(round(self.secs_per_step))

    upstream.DataModel = WholeSecondStep

    parameters(str(WORK / "parameters.cfg"))
    param_files = glob.glob(str(WORK / "params_case" / "params" / "*.nc"))
    if len(param_files) != 1:
        raise RuntimeError(f"expected one parameter file, found {param_files}")

    hist_tapes, data_model, rout_var, time_handle, directories = convolution_init(
        convolution_config(param_files[0], times)
    )
    outflow, in_transit = [], []
    get_rof = rout_var.get_rof

    def recording_rof():
        # Called by convolution_run right after each convolution. Rvar.get_storage is
        # not used: it sums the ring over outlets and not over time, so with one outlet
        # it returns the first entry, the outflow, and not the water held. Entry 0 leaves
        # during this step, so what is still in transit at the end of it is the rest.
        rof = get_rof()
        outflow.append(float(rof["LIQ"][0]))
        in_transit.append(float(rout_var.ring["LIQ"][1:].sum()))
        return rof

    rout_var.get_rof = recording_rof
    convolution_run(hist_tapes, data_model, rout_var, time_handle, directories)
    convolution_final(time_handle, hist_tapes)
    logging.shutdown()
    return np.array(outflow), np.array(in_transit)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent

    with open(str(io_dir / request["input"]["forcing"]), newline="") as fh:
        forcing = list(csv.DictReader(fh))
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    shape = str(static["cross_section_shape"]).strip().lower()
    if shape != "rectangular":
        raise ValueError("the adapter maps a rectangular section only, not %r" % shape)
    if request["timestep"] != "PT1H":
        raise ValueError("this adapter's impulse response is hourly")
    times = [datetime.strptime(row["time"], "%Y-%m-%dT%H:%M:%S") for row in forcing]
    if any(b - a != timedelta(seconds=DT_S) for a, b in zip(times, times[1:])):
        raise ValueError("forcing rows are not hourly")
    q_in = np.array([float(row["q_in"]) for row in forcing])

    WORK.mkdir(parents=True, exist_ok=True)
    params = routing_parameters(float(q_in[0]), static)
    write_inputs(static, times, q_in, params)
    outflow, in_transit = run_rvic(times)
    if len(outflow) != len(forcing):
        raise RuntimeError("RVIC ran %d steps, expected %d" % (len(outflow), len(forcing)))

    area_m2 = float(static["area_km2"]) * 1.0e6
    # The ring is a flux in kg m-2 s-1 over the outlet cell, whose area is `area_km2`.
    # Per step that is a depth, and kg m-2 is mm.
    dis = outflow * area_m2 / 1000.0
    channel = in_transit * DT_S

    columns = ["time", "dis", "snw", "canopy", "channel"]
    with open(str(io_dir / request["output"]["table"]), "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        for i, row in enumerate(forcing):
            writer.writerow([row["time"], repr(float(dis[i])), "0.0", "0.0", repr(float(channel[i]))])

    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "routing": params,
        "reach_length_m": float(static["reach_length_m"]),
        "impulse_response": {
            "uh_box": "impulse in the first hour",
            "basin_flowdays": BASIN_FLOWDAYS,
            "cell_flowdays": CELL_FLOWDAYS,
            "subset_days": SUBSET_DAYS,
        },
        "zero_stores": ["snw", "canopy"],
        "not_reported": {"mrso": "a router has no soil", "gw": "a router has no aquifer"},
    }))
    return 0


if __name__ == "__main__":
    try:
        code = main()
    except Exception:
        # RVIC's logger replaces sys.stderr, so a failure would otherwise say nothing.
        traceback.print_exc(file=sys.__stderr__)
        code = 1
    raise SystemExit(code)

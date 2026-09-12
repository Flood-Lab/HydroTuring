"""Reproduce the supplemental full HRLDAS/Noah-MP single-point experiment.

Requires Linux, git, make, GNU Fortran, NetCDF C/Fortran and Python netCDF4.
Choose a fresh scratch directory; no model sources or outputs enter this repo.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import numpy as np
from netCDF4 import Dataset

from score_full import evaluate

HRLDAS_COMMIT = "cd96df470220f7d7133cdbccd5f9c5355cf173e2"
NOAHMP_COMMIT = "17751dcd7a2442a3c57138f530f523daf08dc1c5"
HERE = Path(__file__).resolve().parent


def build(source: Path, compiler: str, netcdf: Path, log_path: Path) -> None:
    config = (source / "hrldas/arch/user_build_options.gfortran.serial").read_text()
    updates = {
        "COMPILERF90": compiler,
        "F90FLAGS": "-O1 -g -fconvert=big-endian -fbounds-check -fno-range-check -fallow-argument-mismatch -fbacktrace",
        "CPPFLAGS": "-P -traditional -D_GFORTRAN_",
        "NETCDFMOD": f"-I{netcdf / 'include'}",
        "NETCDFLIB": f"-L{netcdf / 'lib'} -Wl,-rpath,{netcdf / 'lib'} -lnetcdff -lnetcdf",
        "LIBJASPER": "", "INCJASPER": "",
    }
    lines = []
    for line in config.splitlines():
        key = line.split("=", 1)[0].strip()
        lines.append(f" {key} = {updates[key]}" if key in updates else line)
    (source / "hrldas/user_build_options").write_text("\n".join(lines) + "\n")
    # Official executable targets in Makefile order. The unrelated GRIB
    # converter needs Jasper; this case uses the official text converter.
    targets = ["hrldas/Utility_routines", "noahmp/utility", "noahmp/src",
               "noahmp/drivers/hrldas", "urban/wrf", "hrldas/IO_code", "hrldas/run"]
    with log_path.open("w") as log:
        for target in targets:
            subprocess.run(["make", "-j1"], cwd=source / target, stdout=log,
                           stderr=subprocess.STDOUT, check=True)


def prepare_case(source: Path, case: Path, compiler: str, netcdf: Path) -> None:
    forcing = case / "forcing"
    forcing.mkdir(parents=True)
    example = source / "hrldas/HRLDAS_forcing/run/examples/single_point"
    header = (example / "bondville.dat").read_text().splitlines()[:54]
    updates = {
        "latitude": "0.0", "longitude": "0.0", "vegetation_category": "16",
        "soil_category": "8", "deep_soil_temperature": "295.0", "elevation": "0.0",
        "maximum_vegetation_pct": "0.0", "minimum_vegetation_pct": "0.0",
        "snow_depth": "0.0", "snow_water_equivalent": "0.0", "canopy_water": "0.0",
        "skin_temperature": "295.0", "soil_temperature": "295.0, 295.0, 295.0, 295.0",
        "soil_moisture": "0.20, 0.20, 0.20, 0.20", "leaf_area_index": "0.0",
        "have_relative_humidity": ".false.", "temperature_offset": "0.0",
        "temperature_scale": "1.0", "pressure_scale": "1.0", "precipitation_scale": "1.0",
    }
    for n, line in enumerate(header):
        key = line.split("=", 1)[0].strip()
        if key in updates:
            header[n] = f" {key} = {updates[key]}"
    header[50:54] = [
        "! Synthetic controlled forcing, not Bondville observations",
        "! yyyy mm dd hh mi | wind | airT K | specific humidity | pressure Pa | SW | LW | rain mm/s",
        "! 24h spinup, 12h shortwave heating, 60h recovery; no precipitation",
        "! Surface heat fluxes are calculated by the complete model",
    ]
    start = datetime(2020, 6, 20, 6)
    longwave = 5.670374419e-8 * 295.0 ** 4
    rows = []
    # The extra forcing endpoint serves the driver's final interpolation.
    for step in range(1153):
        timestamp = start + timedelta(seconds=step * 300)
        shortwave = 300.0 if 288 < step <= 432 else 0.0
        rows.append(f"{timestamp:%Y %m %d %H %M} 2.0 295.0 0.010 101325.0 "
                    f"{shortwave:.6f} {longwave:.9f} 0.0")
    (forcing / "bondville.dat").write_text("\n".join(header + rows) + "\n")
    converter = case / "create_point_data.exe"
    subprocess.run([compiler, f"-I{netcdf / 'include'}", "-o", str(converter),
                    str(example / "create_point_data.f90"), f"-L{netcdf / 'lib'}",
                    f"-Wl,-rpath,{netcdf / 'lib'}", "-lnetcdff", "-lnetcdf"], check=True)
    with (case / "create_forcing.log").open("w") as log:
        subprocess.run([str(converter)], cwd=forcing, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    for variant in ("unmodified", "diagnostic"):
        run = case / variant
        output = run / "output"
        output.mkdir(parents=True)
        namelist = (source / "hrldas/run/namelist.hrldas_example").read_text()
        settings = {
            "HRLDAS_SETUP_FILE": f'"{forcing / "hrldas_setup_single_point.nc"}"',
            "INDIR": f'"{forcing}"', "OUTDIR": f'"{output}"',
            "START_YEAR": "2020", "START_MONTH": "6", "START_DAY": "20", "START_HOUR": "6",
            "KDAY": "4", "TEMP_TIME_SCHEME_OPTION": "1", "FORCING_TIMESTEP": "300",
            "NOAH_TIMESTEP": "300", "SOIL_TIMESTEP": "300", "OUTPUT_TIMESTEP": "300",
            "SPLIT_OUTPUT_COUNT": "0", "RESTART_FREQUENCY_HOURS": "0", "NOAHMP_OUTPUT": "1",
        }
        for key, value in settings.items():
            namelist, count = re.subn(rf"(?m)^\s*{key}\s*=.*$", f" {key} = {value}", namelist)
            if count != 1:
                raise ValueError(f"Missing or repeated namelist setting: {key}")
        (run / "namelist.hrldas").write_text(namelist)
        shutil.copy2(source / "noahmp/parameters/NoahmpTable.TBL", run / "NoahmpTable.TBL")


def run_model(source: Path, run: Path) -> float:
    executable = run / "hrldas.exe"
    shutil.copy2(source / "hrldas/run/hrldas.exe", executable)
    started = time.monotonic()
    with (run / "run.log").open("w") as log:
        subprocess.run([str(executable)], cwd=run, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    return time.monotonic() - started


def compare_outputs(case: Path) -> dict:
    name = "202006200600.LDASOUT_DOMAIN1"
    checked = 0
    differences = []
    with Dataset(case / "unmodified/output" / name) as a, Dataset(case / "diagnostic/output" / name) as b:
        a.set_auto_mask(False)
        b.set_auto_mask(False)
        if set(a.variables) != set(b.variables):
            raise ValueError("Instrumentation changed the set of native output variables")
        for key in a.variables:
            left, right = a.variables[key][:], b.variables[key][:]
            equal = np.array_equal(left, right, equal_nan=True) if left.dtype.kind in "fc" else np.array_equal(left, right)
            if not equal:
                differences.append(key)
            checked += 1
        records = len(a.dimensions["Time"])
    result = {"variables_compared": checked, "records_including_initial": records,
              "changed_variables": differences, "all_native_output_values_equal": not differences}
    if records != 1153 or differences:
        raise ValueError(f"Full-run output comparison failed: {result}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True, help="fresh scratch directory")
    parser.add_argument("--fc", default="gfortran")
    parser.add_argument("--netcdf-prefix", type=Path, required=True)
    args = parser.parse_args()
    workdir = args.workdir.resolve()
    if workdir.exists() and any(workdir.iterdir()):
        parser.error("choose a fresh --workdir so previous model outputs are preserved")
    workdir.mkdir(parents=True, exist_ok=True)
    compiler = shutil.which(args.fc)
    if compiler is None:
        parser.error("GNU Fortran compiler not found")
    source = workdir / "source"
    subprocess.run(["git", "clone", "https://github.com/NCAR/hrldas", str(source)], check=True)
    subprocess.run(["git", "checkout", HRLDAS_COMMIT], cwd=source, check=True)
    subprocess.run(["git", "submodule", "update", "--init", "--recursive"], cwd=source, check=True)
    noah_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source / "noahmp", text=True).strip()
    if noah_commit != NOAHMP_COMMIT:
        raise ValueError("The checked-out Noah-MP version differs from the documented experiment")
    case = workdir / "case"
    build(source, compiler, args.netcdf_prefix.resolve(), workdir / "build_unmodified.log")
    prepare_case(source, case, compiler, args.netcdf_prefix.resolve())
    unmodified_seconds = run_model(source, case / "unmodified")
    subprocess.run([sys.executable, str(HERE / "instrument_full.py"), "--source", str(source)], check=True)
    build(source, compiler, args.netcdf_prefix.resolve(), workdir / "build_diagnostic.log")
    diagnostic_seconds = run_model(source, case / "diagnostic")
    comparison = compare_outputs(case)
    report = evaluate(case / "diagnostic")
    provenance = {
        "hrldas_commit": HRLDAS_COMMIT, "noahmp_commit": NOAHMP_COMMIT,
        "compiler": subprocess.check_output([compiler, "--version"], text=True).splitlines()[0],
        "precision": "official default single precision", "execution": "serial CPU",
        "forcing": "synthetic warm bare-soil single point; official point-data converter",
        "unmodified_run_seconds": unmodified_seconds, "diagnostic_run_seconds": diagnostic_seconds,
        "instrumentation_comparison": comparison,
        "validation_report": "case/diagnostic/full_model_validation.json",
    }
    (workdir / "reproduction.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2))
    print("Fixed-capacity probe applicable:", report["fixed_capacity_criterion_replay"]["applicable_to_full_run"])


if __name__ == "__main__":
    main()

"""Optional native Noah-MP thermal-component check; see docs/noahmp-soil-heat-validation.md."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import subprocess
from urllib.request import urlopen

import pandas as pd

from hydroturing.criteria import get
from hydroturing.protocol import Case, RunResult

COMMIT = "17751dcd7a2442a3c57138f530f523daf08dc1c5"  # official NCAR/noahmp v5.2.1
SOURCE = f"https://raw.githubusercontent.com/NCAR/noahmp/{COMMIT}"
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MODULES = [
    "Machine", "ConstantDefineMod", "ForcingVarType", "ConfigVarType",
    "EnergyVarType", "WaterVarType", "BiochemVarType", "NoahmpVarType",
    "SoilThermalPropertyMod", "MatrixSolverTriDiagonalMod",
    "SoilSnowThermalDiffusionMod", "SoilSnowTemperatureSolverMod",
]
SOURCE_NOTICE = """The source of this material is the Research Applications Laboratory at the
National Center for Atmospheric Research, a program of the University Corporation
for Atmospheric Research (UCAR) pursuant to a Cooperative Agreement with the
National Science Foundation; ©2007 University Corporation for Atmospheric
Research. All Rights Reserved.
"""


def instrument_solver(source: str) -> str:
    """Expose the implicit lower-face Fourier flux without changing temperature equations."""
    signature = (
        "subroutine SoilSnowTemperatureSolver(noahmp, TimeStep, "
        "MatLeft1, MatLeft2, MatLeft3, MatRight)"
    )
    declaration = "    real(kind=kind_noahmp), intent(in)    :: TimeStep"
    location = "    ! deallocate local arrays to avoid memory leaks"
    for target in (signature, declaration, location):
        if source.count(target) != 1:
            raise ValueError(f"pinned Noah-MP source does not contain expected insertion: {target}")
    source = source.replace(signature, signature[:-1] + ", HeatTopLayerBottom)")
    source = source.replace(
        declaration,
        "    real(kind=kind_noahmp), optional, intent(out) :: HeatTopLayerBottom\n" + declaration,
    )
    return source.replace(location, """    ! HydroTuring diagnostic only: native Fourier relation evaluated after
    ! the implicit solution. No snow; OptSnowSoilTempTime=1 (prescribed flux).
    if (present(HeatTopLayerBottom)) then
       HeatTopLayerBottom = 2.0 * noahmp%energy%state%ThermConductSoilSnow(1) * &
          (TemperatureSoilSnow(1)-TemperatureSoilSnow(2)) / &
          (-noahmp%config%domain%DepthSnowSoilLayer(2))
    endif

""" + location)


def build_and_run(workdir: Path, compiler: str) -> str:
    """Keep downloaded source, build products and output entirely in the scratch directory."""
    for name in MODULES:
        folder = "utility" if name == "Machine" else "src"
        with urlopen(f"{SOURCE}/{folder}/{name}.F90", timeout=30) as response:
            (workdir / f"{name}.F90").write_bytes(response.read())
    with urlopen(f"{SOURCE}/LICENSE.txt", timeout=30) as response:
        (workdir / "LICENSE.txt").write_bytes(response.read())
    (workdir / "NOTICE.txt").write_text(SOURCE_NOTICE, encoding="utf-8")
    original = (workdir / "SoilSnowTemperatureSolverMod.F90").read_text(encoding="utf-8")
    instrumented = workdir / "SoilSnowTemperatureSolverInstrumented.F90"
    instrumented.write_text(instrument_solver(original), encoding="utf-8")
    executable = workdir / ("thermal_component_driver.exe" if os.name == "nt" else "thermal_component_driver")
    files = [str(workdir / f"{name}.F90") for name in MODULES[:-1]]
    files += [str(instrumented), str(HERE / "thermal_component_driver.F90")]
    try:
        version = subprocess.check_output([compiler, "--version"], text=True).splitlines()[0]
    except FileNotFoundError:
        raise SystemExit(
            f"GNU Fortran was not found at '{compiler}'. Provide --fc /path/to/gfortran, "
            "or use --skip-build with existing output CSV files."
        ) from None
    subprocess.run(
        [compiler, "-cpp", "-DDOUBLE_PREC", "-O2", "-ffree-line-length-none", "-fcheck=bounds",
         *files, "-o", str(executable)], cwd=workdir, check=True,
    )
    for step in (300, 3600):
        with (workdir / f"thermal_component_output_{step}s.csv").open("w", encoding="utf-8") as output:
            subprocess.run([str(executable), str(step)], cwd=workdir, stdout=output, check=True)
    (workdir / "compiler.txt").write_text(version + "\n", encoding="utf-8")
    return version


def score_outputs(path: Path) -> dict:
    data = pd.read_csv(path, encoding="utf-8")
    if len(data) != 96:
        raise ValueError(f"{path.name}: expected 24 spinup + 72 scored hourly rows")
    times = pd.date_range("2001-07-01", periods=len(data), freq="h")
    initial = float(data.t_start_k.iloc[0])
    case = Case(
        probe_id="validation/noahmp-soil-heat-component", seed=0,
        forcing=pd.DataFrame({"time": times, "_phase": data.phase}),
        static={"soil_heat_capacity_areal": float(data.heat_capacity_area_j_m2_k.iloc[0]),
                "soil_temperature_initial": initial},
        spinup_steps=24, timestep="PT1H",
    )
    native = pd.DataFrame({
        "time": times, "hfg": data.g_top_mean_w_m2,
        "hfg_bottom": data.g_bottom_mean_w_m2, "tsoil_layer": data.t_end_k,
    })
    results = {}
    for variant in ("native", "frozen_temperature", "half_temperature_change", "wrong_hour_start_flux"):
        table = native.copy()
        if variant == "frozen_temperature":
            table["tsoil_layer"] = initial
        elif variant == "half_temperature_change":
            table["tsoil_layer"] = initial + 0.5 * (table.tsoil_layer - initial)
        elif variant == "wrong_hour_start_flux":
            table["hfg_bottom"] = data.g_bottom_hour_start_w_m2
        run = RunResult(case, table, {"source": "Noah-MP thermal component"}, 0.0)
        result = get("soil_heat_storage")(run, None, {"threshold": 0.05, "floor": 1.0})
        results[variant] = asdict(result)
    return {
        "source_csv": path.name, "rows": len(data),
        "temperature_range_k": [float(data.t_end_k.min()), float(data.t_end_k.max())],
        "heat_capacity_area_j_m2_k": case.static["soil_heat_capacity_areal"],
        "cases": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True,
                        help="scratch directory outside the repository, or within its _tmp/")
    parser.add_argument("--fc", default="gfortran", help="GNU Fortran compiler executable")
    parser.add_argument("--skip-build", action="store_true",
                        help="score existing thermal_component_output_{300,3600}s.csv; no compiler/network")
    args = parser.parse_args()
    workdir = args.workdir.resolve()
    if workdir == ROOT or (workdir.is_relative_to(ROOT) and not workdir.is_relative_to(ROOT / "_tmp")):
        parser.error("--workdir must be outside the repository or within its _tmp/ directory")
    workdir.mkdir(parents=True, exist_ok=True)
    compiler = "not invoked (--skip-build)"
    if not args.skip_build:
        compiler = build_and_run(workdir, args.fc)
    report = {
        "scope": "native thermal-component output evaluated by the soil_heat_storage criterion; not a full-model run or gate",
        "source_version": "NCAR/noahmp v5.2.1", "source_commit": COMMIT,
        "source_url": f"https://github.com/NCAR/noahmp/tree/{COMMIT}", "compiler": compiler,
        "runs": {str(step): score_outputs(workdir / f"thermal_component_output_{step}s.csv")
                 for step in (300, 3600)},
    }
    destination = workdir / "hydroturing_criterion_validation.json"
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    expected = {"native": "pass", "frozen_temperature": "fail", "half_temperature_change": "fail"}
    unexpected = []
    for step, run in report["runs"].items():
        for name, result in run["cases"].items():
            print(f"{step}s {name}: {result['status'].upper()}; {result['message']}")
            if name in expected and result["status"] != expected[name]:
                unexpected.append(f"{step}s {name}")
    print(f"Results: {destination}")
    if unexpected:
        raise SystemExit("Unexpected native/temperature-control result: " + ", ".join(unexpected))


if __name__ == "__main__":
    main()

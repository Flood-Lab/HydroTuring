#!/usr/bin/env python3
"""Run the local MODFLOW 6 groundwater/riv experiment through /io."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
from pathlib import Path
import shutil
import tempfile

import flopy
import numpy as np

MODEL = {"name": "modflow6", "version": "0.1.0"}
COLUMNS = ["time", "gwex", "gw_to_sw", "sw_to_gw", "gw"]
NROW = 10
NCOL = 10
DELR = 100.0
DELC = 100.0
TOP = 20.0
BOTM = 0.0
K = 1.0
RIVER_ROW = 5
RIVER_COL = 5


def resolve_mf6_exe() -> str:
    """Resolve the platform-specific MODFLOW 6.7.0 executable."""
    override = os.environ.get("MODFLOW6_EXE")
    if override:
        return override

    system = platform.system().lower()
    if system == "windows":
        candidates = [Path(r"D:\hydrowang\modflow\bin\mf6.exe")]
        executable_names = ("mf6.exe", "mf6")
    elif system == "linux":
        candidates = [
            Path("/opt/modflow6/bin/mf6"),
            Path("/opt/modflow6/mf6"),
            Path("/usr/local/bin/mf6"),
        ]
        executable_names = ("mf6",)
    else:
        candidates = []
        executable_names = ("mf6", "mf6.exe")

    root = os.environ.get("MODFLOW6_ROOT")
    if root:
        root_path = Path(root)
        for name in executable_names:
            candidates.extend(root_path.rglob(name))

    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    for name in executable_names:
        found = shutil.which(name)
        if found:
            return found

    release = "https://github.com/MODFLOW-ORG/modflow6/releases/tag/6.7.0"
    raise RuntimeError(
        f"MODFLOW 6.7.0 executable not found for {system or 'unknown'}; "
        f"set MODFLOW6_EXE or MODFLOW6_ROOT (release: {release})"
    )


def run_modflow(forcing: list[dict], static: dict) -> list[dict]:
    mf6_exe = resolve_mf6_exe()
    if not Path(mf6_exe).exists():
        raise RuntimeError(f"MODFLOW 6 executable not found: {mf6_exe}")
    area_m2 = float(static["area_km2"]) * 1.0e6
    storage_coefficient = float(static["aquifer_storage_coefficient"])
    specific_yield = float(static["aquifer_specific_yield"])
    river_conductance = float(static["river_conductance_m2_per_day"])
    river_bottom_offset = float(static["river_bottom_offset_m"])
    initial_head = float(static.get("aquifer_initial_head_m", 10.0))

    with tempfile.TemporaryDirectory(prefix="hydroturing-mf6-") as workdir:
        sim = flopy.mf6.MFSimulation(sim_name="ht_gw", sim_ws=workdir, exe_name=mf6_exe)
        nper = len(forcing)
        flopy.mf6.ModflowTdis(
            sim, time_units="days", nper=nper,
            perioddata=[(1.0, 1, 1.0)] * nper,
        )
        flopy.mf6.ModflowIms(
            sim, complexity="SIMPLE", outer_dvclose=1e-8,
            inner_dvclose=1e-8,
        )
        gwf = flopy.mf6.ModflowGwf(sim, modelname="ht_gw", save_flows=True)
        flopy.mf6.ModflowGwfdis(
            gwf, nlay=1, nrow=NROW, ncol=NCOL, delr=DELR, delc=DELC,
            top=TOP, botm=BOTM,
        )
        flopy.mf6.ModflowGwfic(gwf, strt=initial_head)
        flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=K)
        flopy.mf6.ModflowGwfsto(
            gwf, iconvert=1, storagecoefficient=False,
            ss=1.0e-5, sy=specific_yield,
            steady_state={0: False},
            transient={period: True for period in range(nper)},
        )
        flopy.mf6.ModflowGwfrcha(
            gwf,
            recharge={period: float(row["gw_recharge"]) / 1000.0 for period, row in enumerate(forcing)},
        )
        riv = {
            period: [[(0, RIVER_ROW, RIVER_COL), float(row["sw_stage_m"]),
                      river_conductance, float(row["sw_stage_m"]) - river_bottom_offset]]
            for period, row in enumerate(forcing)
        }
        flopy.mf6.ModflowGwfriv(gwf, stress_period_data=riv, save_flows=True)
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord="ht_gw.hds",
            budget_filerecord="ht_gw.cbc",
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )
        sim.write_simulation(silent=True)
        success, output = sim.run_simulation(silent=True)
        if not success:
            raise RuntimeError("MODFLOW 6 failed:\n" + "\n".join(output[-20:]))

        head_file = gwf.output.head()
        budget_file = gwf.output.budget()
        heads = head_file.get_alldata()
        rows = []
        for period, forcing_row in enumerate(forcing):
            q_records = budget_file.get_data(text="RIV", kstpkper=(0, period))[0]
            # MODFLOW RIV q is positive into the aquifer. The probe's signed
            # components are positive GW -> river and negative river -> GW.
            q_m3_day = float(q_records["q"].sum())
            q_mm_day = q_m3_day / area_m2 * 1000.0
            gw_to_sw = max(-q_mm_day, 0.0)
            sw_to_gw = min(-q_mm_day, 0.0)
            gwex = gw_to_sw + sw_to_gw
            mean_head = float(np.mean(heads[period, 0]))
            gw_mm = specific_yield * max(mean_head - BOTM, 0.0) * 1000.0
            rows.append({
                "time": forcing_row["time"],
                "gwex": gwex,
                "gw_to_sw": gw_to_sw,
                "sw_to_gw": sw_to_gw,
                "gw": gw_mm,
            })
        head_file.close()
        budget_file.close()
        return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    io_dir = request_path.parent
    request = json.loads(request_path.read_text())
    with open(io_dir / request["input"]["forcing"], newline="") as handle:
        forcing = list(csv.DictReader(handle))
    for row in forcing:
        row["gw_recharge"] = float(row["gw_recharge"])
        row["sw_stage_m"] = float(row["sw_stage_m"])
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    rows = run_modflow(forcing, static)
    output = io_dir / request["output"]["table"]
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows)}, indent=2)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

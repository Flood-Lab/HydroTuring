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
        candidates = []
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
            ss=storage_coefficient, sy=specific_yield,
            save_flows=True,
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

        budget_file = gwf.output.budget()
        rows = []
        gw_mm = specific_yield * max(initial_head - BOTM, 0.0) * 1000.0
        for period, forcing_row in enumerate(forcing):
            q_records = budget_file.get_data(text="RIV", kstpkper=(0, period))[0]
            # MODFLOW RIV q is positive into the aquifer, i.e. the river
            # losing water to groundwater: sw_to_gw (>= 0) per AGENTS.md's
            # gwex-positive-into-the-aquifer convention.
            q_m3_day = float(q_records["q"].sum())
            q_mm_day = q_m3_day / area_m2 * 1000.0
            sw_to_gw = max(q_mm_day, 0.0)
            gw_to_sw = min(q_mm_day, 0.0)
            gwex = gw_to_sw + sw_to_gw

            # MODFLOW's own STO-SS/STO-SY budget terms, not a head-derived
            # approximation: they are exact even with a non-uniform head
            # field, and use aquifer_storage_coefficient (Ss) as well as
            # specific yield. Positive q is water released FROM storage
            # (storage decreasing), so the storage change is its negative.
            # Unlike RIV/RCHA, MODFLOW 6 writes STO-SS/STO-SY as a dense
            # per-cell array (one value per model cell), not a list of
            # discrete (node, q) records, because storage change applies
            # everywhere rather than at a few boundary cells. So it has no
            # named "q" field: the array itself is the flow.
            sto_ss = budget_file.get_data(text="STO-SS", kstpkper=(0, period))[0]
            sto_sy = budget_file.get_data(text="STO-SY", kstpkper=(0, period))[0]
            storage_release_m3_day = float(sto_ss.sum() + sto_sy.sum())
            gw_mm += -storage_release_m3_day / area_m2 * 1000.0

            rows.append({
                "time": forcing_row["time"],
                "gwex": gwex,
                "gw_to_sw": gw_to_sw,
                "sw_to_gw": sw_to_gw,
                "gw": gw_mm,
            })
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

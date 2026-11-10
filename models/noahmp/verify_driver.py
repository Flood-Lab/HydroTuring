#!/usr/bin/env python3
"""Compare ht_noahmp_driver with upstream's own entry point, hrldas.exe.

Both run one land point on the same synthetic forcing, with the same namelist and
table, at one Noah-MP step per forcing step. hrldas.exe reads the forcing through
upstream's single-point converter and writes LDASOUT; the driver takes the same
numbers, the cosine of the zenith angle HRLDAS computed included. Every series
both can report is compared. Run inside the image:

    python3 /model/verify_driver.py <scratch directory>
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np

SHIPPED = Path("/opt/noahmp")
DRIVER = SHIPPED / "ht_noahmp_driver"
STEP_S = 3600
N_STEPS = 24 * 150
START = datetime(2001, 1, 1, 0)
SKIN = 263.7
SOIL_T = [266.1, 274.0, 276.9, 279.9]
SOIL_M = [0.298, 0.294, 0.271, 0.307]
LAYERS = [0.10, 0.30, 0.60, 1.00]
NODES = [0.05, 0.25, 0.70, 1.50]


def forcing() -> dict[str, np.ndarray]:
    """A winter into spring at Bondville's latitude, with snow, rain and a thaw."""
    rng = np.random.default_rng(59)
    hours = np.arange(N_STEPS)
    hour_of_day = hours % 24
    day = hours / 24.0
    tk = 268.0 + 16.0 * day / 150.0 + 4.0 * np.cos(2.0 * np.pi * (hour_of_day - 14.0) / 24.0)
    sw = np.maximum(0.0, 500.0 * np.sin(np.pi * (hour_of_day - 6.0) / 12.0)) * rng.uniform(0.4, 1.0, N_STEPS)
    lw = 0.8 * 5.67e-8 * tk ** 4
    rain = np.where(rng.random(N_STEPS) < 0.05, rng.uniform(0.5, 6.0, N_STEPS), 0.0) / 3600.0
    ea = 0.7 * 611.2 * np.exp(17.67 * (tk - 273.15) / (tk - 273.15 + 243.5))
    p = np.full(N_STEPS, 99000.0)
    return {"wind": 2.0 + 2.0 * rng.random(N_STEPS), "tk": tk, "q": 0.622 * ea / (p - 0.378 * ea), "p": p,
            "sw": sw, "lw": lw, "rain": rain}


def namelist_set(text: str, settings: dict[str, str]) -> str:
    for key, value in settings.items():
        text, count = re.subn(rf"(?m)^\s*{re.escape(key)}\s*=.*$", f" {key} = {value}", text)
        if count != 1:
            raise ValueError(f"namelist has no single {key}")
    return text


def run_hrldas(work: Path, f: dict[str, np.ndarray], veg: int, gvf: tuple[float, float]) -> netCDF4.Dataset:
    forcing_dir = work / "forcing"
    forcing_dir.mkdir(parents=True)
    header = (SHIPPED / "bondville.dat").read_text().splitlines()[:54]
    updates = {
        "vegetation_category": str(veg), "maximum_vegetation_pct": str(gvf[0]), "minimum_vegetation_pct": str(gvf[1]),
        "snow_depth": "0.0", "snow_water_equivalent": "0.0", "canopy_water": "0.0",
        "have_relative_humidity": ".false.", "temperature_offset": "0.0", "temperature_scale": "1.0",
        "pressure_scale": "1.0", "precipitation_scale": "1.0",
    }
    for i, line in enumerate(header):
        key = line.split("=", 1)[0].strip()
        if key in updates:
            header[i] = f" {key} = {updates[key]}"
    rows = []
    # Row i forces the step ending at START + (i + 1) dt; the first stamp is the start.
    for i in range(-1, N_STEPS):
        j = max(i, 0)
        when = START + timedelta(seconds=(i + 1) * STEP_S)
        rows.append(f"{when:%Y %m %d %H %M} {f['wind'][j]:.12g} {f['tk'][j]:.12g} {f['q'][j]:.12g} "
                    f"{f['p'][j]:.12g} {f['sw'][j]:.12g} {f['lw'][j]:.12g} {f['rain'][j]:.12g}")
    (forcing_dir / "bondville.dat").write_text("\n".join(header + rows) + "\n")
    subprocess.run([str(SHIPPED / "create_point_data.exe")], cwd=forcing_dir, check=True, capture_output=True)
    for path in forcing_dir.glob("????????????.LDASIN_DOMAIN1"):
        path.rename(path.with_name(path.name[:10] + path.name[12:]))
    run = work / "hrldas"
    (run / "output").mkdir(parents=True)
    shutil.copy(SHIPPED / "NoahmpTable.TBL", run / "NoahmpTable.TBL")
    (run / "namelist.hrldas").write_text(namelist_set((SHIPPED / "namelist.hrldas_example").read_text(), {
        "HRLDAS_SETUP_FILE": f'"{forcing_dir / "hrldas_setup_single_point.nc"}"',
        "INDIR": f'"{forcing_dir}"', "OUTDIR": f'"{run / "output"}"',
        "START_YEAR": str(START.year), "START_MONTH": str(START.month), "START_DAY": str(START.day),
        "START_HOUR": str(START.hour), "START_MIN": str(START.minute), "KDAY": str(N_STEPS // 24),
        "FORCING_TIMESTEP": str(STEP_S), "NOAH_TIMESTEP": str(STEP_S), "SOIL_TIMESTEP": "0.0",
        "OUTPUT_TIMESTEP": str(STEP_S), "SPLIT_OUTPUT_COUNT": "0", "RESTART_FREQUENCY_HOURS": "0",
        "NOAHMP_OUTPUT": "1",
    }))
    subprocess.run([str(SHIPPED / "hrldas.exe")], cwd=run, check=True, capture_output=True)
    return netCDF4.Dataset(next((run / "output").glob("*LDASOUT_DOMAIN1")))


def run_driver(work: Path, f: dict[str, np.ndarray], cosz: np.ndarray, veg: int, gvf: tuple[float, float]) -> np.ndarray:
    run = work / "driver"
    run.mkdir(parents=True)
    shutil.copy(SHIPPED / "NoahmpTable.TBL", run / "NoahmpTable.TBL")
    (run / "namelist.hrldas").write_text(namelist_set((SHIPPED / "namelist.hrldas_example").read_text(), {
        "HRLDAS_SETUP_FILE": '"none"', "INDIR": '"none"', "OUTDIR": '"none"', "KDAY": "1",
        "FORCING_TIMESTEP": str(STEP_S), "NOAH_TIMESTEP": str(STEP_S), "SOIL_TIMESTEP": "0.0",
        "OUTPUT_TIMESTEP": str(STEP_S), "RESTART_FREQUENCY_HOURS": "0",
    }))
    (run / "control.nml").write_text(
        "&control\n"
        f" nrows={N_STEPS}, vegtype={veg}, soiltype=8, latitude=40.01, longitude=-88.37, tmn=285.0, skin={SKIN!r},\n"
        f" gvfmax={gvf[0]!r}, gvfmin={gvf[1]!r}, lai_initial=2.0,\n"
        f" soil_t={','.join(repr(v) for v in SOIL_T)}, soil_m={','.join(repr(v) for v in SOIL_M)}\n"
        "/\n")
    ends = [START + timedelta(seconds=(i + 1) * STEP_S) for i in range(N_STEPS)]
    julian = [(e.date() - datetime(e.year, 1, 1).date()).days + e.hour / 24.0 for e in ends]
    table = np.column_stack([f["sw"], f["lw"], f["tk"], f["wind"], f["q"], f["p"], f["rain"], cosz, julian,
                             [e.year for e in ends]])
    table.astype(np.float64).tofile(run / "forcing.bin")
    subprocess.run([str(DRIVER)], cwd=run, check=True, capture_output=True)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from ht_adapter import OUTPUT_NAMES
    raw = np.fromfile(run / "output.bin", dtype=np.float64).reshape(N_STEPS, len(OUTPUT_NAMES))
    return dict(zip(OUTPUT_NAMES, raw.T))


def compare(name: str, work: Path, f: dict[str, np.ndarray], veg: int, gvf: tuple[float, float]) -> dict:
    with run_hrldas(work, f, veg, gvf) as nc:
        nc.set_auto_mask(False)
        def series(var: str) -> np.ndarray:
            data = np.asarray(nc.variables[var][1:], dtype=np.float64)
            return data.reshape(data.shape[0], -1)
        cosz = series("COSZ")[:, 0]
        hrldas = {
            "FSA": series("FSA")[:, 0], "FIRA": series("FIRA")[:, 0], "HFX": series("HFX")[:, 0],
            "LH": series("LH")[:, 0], "GRDFLX": series("GRDFLX")[:, 0], "TRAD": series("TRAD")[:, 0],
            "SNEQV": series("SNEQV")[:, 0], "SNOWH": series("SNOWH")[:, 0], "CANLIQ": series("CANLIQ")[:, 0],
            "CANICE": series("CANICE")[:, 0], "ECAN": series("ECAN")[:, 0], "ETRAN": series("ETRAN")[:, 0],
            "EDIR": series("EDIR")[:, 0], "SOIL_T1": series("SOIL_T")[:, 0],
            "SOILWATER": (series("SOIL_M") * np.array(LAYERS) * 1000.0).sum(axis=1),
            "RUNOFF": np.diff(np.concatenate(([0.0], series("SFCRNOFF")[:, 0] + series("UGDRNOFF")[:, 0]))),
        }
    d = run_driver(work, f, cosz, veg, gvf)
    driver = {
        "FSA": d["RadSwAbsSfc"], "FIRA": d["RadLwNetSfc"], "HFX": d["HeatSensibleSfc"],
        "LH": d["HeatLatentCanopy"] + d["HeatLatentGrd"] + d["HeatLatentTransp"], "GRDFLX": d["HeatGroundTot"],
        "TRAD": d["TemperatureRadSfc"], "SNEQV": d["SnowWaterEquiv"], "SNOWH": d["SnowDepth"],
        "CANLIQ": d["CanopyLiqWater"], "CANICE": d["CanopyIce"], "ECAN": d["EvapCanopyNet"],
        "ETRAN": d["Transpiration"], "EDIR": d["EvapGroundNet"], "SOIL_T1": d["TemperatureSoil1"],
        "SOILWATER": d["SoilWater"],
        "RUNOFF": d["RunoffSurface"] + d["RunoffSubsurface"],
    }
    report = {"case": name, "steps": N_STEPS, "peak_swe_mm": float(hrldas["SNEQV"].max()),
              "total_runoff_mm": float(hrldas["RUNOFF"].sum()), "series": {}}
    for key in hrldas:
        delta = np.abs(hrldas[key] - driver[key])
        report["series"][key] = {"max_abs_difference": float(delta.max()), "max_abs_value": float(np.abs(hrldas[key]).max()),
                                 "identical_steps": int((delta == 0.0).sum())}
    return report


def main() -> int:
    root = Path(sys.argv[1]).resolve()
    f = forcing()
    reports = []
    for name, veg, gvf in (("croplands, class 12, as upstream's point example", 12, (96.0, 1.0)),
                           ("barren, class 16", 16, (0.0, 0.0))):
        work = root / f"veg{veg}"
        if work.exists():
            shutil.rmtree(work)
        reports.append(compare(name, work, f, veg, gvf))
    print(json.dumps(reports, indent=1))
    worst = max(s["max_abs_difference"] / max(s["max_abs_value"], 1.0) for r in reports for s in r["series"].values())
    print(f"worst difference relative to the series' largest value: {worst:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

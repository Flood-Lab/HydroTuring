#!/usr/bin/env python3
"""HydroTuring adapter for SHUD (SHUD-System/SHUD, commit 851fb8c).

The physics is upstream's: the `shud` executable compiled unmodified, reading the
text inputs it always reads and writing the binary outputs it always writes. This
file writes those inputs and turns those outputs into the contract's rows.
README.md states every choice and how much it moves results.

Mesh. Four triangles and one reach: a representative stretch of river with a
hillslope on each bank, sized and sloped from upstream's example catchment, North
Fork Cache Creek. Forcing is uniform, and the stated area enters only `dis`.

Forcing. SHUD computes its own potential evaporation by Penman-Monteith from
temperature, humidity, wind and net radiation; it has no PET input. Those three
are mocked row by row so that its PET comes out at the case's `pet` (see
`mock_meteorology`). Labelled in run.json.

Reading. Upstream's writer averages what it prints, states too, over the output
interval, so a store's end-of-row value is only available by printing at the
solver's own step. That is 144 records per day per variable, which would not fit a
run's tmpfs, so the output files are named pipes and each is read and reduced to
one value per row as it is written.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
import threading
from datetime import datetime
from pathlib import Path

import numpy as np

COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "mrso", "snw", "canopy", "gw", "channel", "stage"]
MODEL = {"name": "shud", "version": "2.2-851fb8c.1", "revision": "851fb8c7e146a533a9545f786e8f68fd8a8d817f"}

MINUTES = {"PT1D": 1440, "PT1H": 60}
UPSTREAM = Path(os.environ.get("HT_SHUD_UPSTREAM", "/model/upstream/ccw"))
PROJECT = "ht"
DIAG = bool(os.environ.get("HT_SHUD_DIAG"))
# For the sensitivity runs README.md reports. A case never sets these.
OPT_VEGFRAC = os.environ.get("HT_SHUD_VEGFRAC", "")
OPT_AQUIFER = os.environ.get("HT_SHUD_AQUIFER", "")

# Upstream's solver step and tolerances in ccw.cfg.para.
SOLVER_STEP_MIN = int(os.environ.get("HT_SHUD_SOLVER_STEP", 10))
ABSTOL = float(os.environ.get("HT_SHUD_TOL", 1e-4))
RELTOL = float(os.environ.get("HT_SHUD_TOL", 1e-4))

# The mesh: statistics of ccw.sp.*, see README.md.
REACH_LENGTH_M = 1167.0     # median length of a type 1 reach
HILLSLOPE_M = 776.0         # catchment area / (2 * total river length)
HILL_SLOPE = float(os.environ.get("HT_SHUD_HILL_SLOPE", 0.25))    # median element slope
RIVER_SLOPE = 0.076         # median reach slope
Z_OUTLET_M = 450.0          # the example's outlet elevation
LAND_COVER_ROW = 11         # the class of 1020 of its 1147 elements

# ccw.sp.att and ccw.cfg.calib give the river these. Their width term belongs to a
# network draining 196 km2; here the reach has its own type's 2 m.
RIVER_WIDTH_ADD_M = float(os.environ.get("HT_SHUD_WIDTH_ADD", 0.0))
# The river type's depth relative to the aquifer's, as in the example: 5.5 m in 30 m.
RIVER_DEPTH_OVER_AQUIFER = float(os.environ.get("HT_SHUD_DEPTH_RATIO", 5.5 / 30.0))
# For a case that states no capacity: the depth a 320 mm capacity, the suite's usual, gives.
DEFAULT_AQUIFER_M = 0.8
# Upstream refuses an element whose roots reach below its aquifer (RootReachLevel <= 0),
# which a small stated capacity would ask for. Roots keep the example's depth until they
# no longer fit, then this share of the aquifer.
MAX_ROOT_SHARE_OF_AQUIFER = 0.9

# Mock meteorology. Wind is constant at the mean of the example's forcing; humidity
# and radiation come from each row's own `pet` and `tas`, so no row knows the calendar.
WIND_MS = 1.82
WIND_MEASURE_M = 10.0       # HeightWindMeasure in Macros.hpp


def read_table(path: Path) -> np.ndarray:
    """One of upstream's tables: a size line, a header line, then rows."""
    return np.loadtxt(path, skiprows=2, ndmin=2)


def representative_parameters() -> dict:
    """Upstream's example parameters, the median over its elements for each column."""
    att = read_table(UPSTREAM / "ccw.sp.att")
    soil = read_table(UPSTREAM / "ccw.para.soil")[att[:, 1].astype(int) - 1]
    geol = read_table(UPSTREAM / "ccw.para.geol")[att[:, 2].astype(int) - 1]
    land = read_table(UPSTREAM / "ccw.para.lc")[LAND_COVER_ROW - 1]
    if OPT_VEGFRAC:
        land[2] = float(OPT_VEGFRAC)
    lai = read_table(UPSTREAM / "ccw.tsd.lai")[:, LAND_COVER_ROW]
    melt = read_table(UPSTREAM / "ccw.tsd.mf")[:, 1]
    return {
        "soil": np.median(soil, axis=0)[1:],
        "geol": np.median(geol, axis=0)[1:],
        "land": land[1:],
        "lai": float(lai.mean()),
        "mf": float(melt.mean()),
    }


def calibration() -> dict:
    values = {}
    for line in (UPSTREAM / "ccw.cfg.calib").read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            values[parts[0]] = float(parts[1])
    values["RIV_WDTH+"] = RIVER_WIDTH_ADD_M
    return values


def mock_meteorology(pet_mm_day: np.ndarray, tas: np.ndarray, z_m: float, land: np.ndarray, lai: float) -> tuple[np.ndarray, np.ndarray]:
    """Humidity and net radiation under which SHUD's two demands both equal `pet`.

    tReadForcing computes an open-water demand and a Penman-Monteith canopy demand,
    each linear in available energy and in the vapour pressure deficit, so a row's
    `pet` and `tas` determine both unknowns. The same expressions as MD_ET.cpp and
    is_sm_et.hpp, in mm/day and MJ/m2/day, solved. README.md shows what comes out and
    why two equal demands: one saturated, radiation-only sky is what was tried first,
    and its canopy demand was a third of its open-water demand.
    """
    albedo = land[0]
    lam = 2.501 - 0.002361 * tas
    pressure = 101.325 * ((293.0 - 0.0065 * z_m) / 293.0) ** 5.26
    gamma = 0.0016286 * pressure / lam
    es = 0.6108 * np.exp(17.27 * tas / (tas + 237.3))
    delta = 4098.0 * es / (tas + 237.3) ** 2
    rho = 3.486 * pressure / (275.0 + tas)
    uz = abs(WIND_MS) + 0.001
    u2 = uz * math.log(2.0 / 0.00137) / math.log(WIND_MEASURE_M / 0.00137)
    hc = lai * 0.5
    zm = hc * 1.3333
    d = 0.67 * hc
    ra = math.log(abs(zm - d) / (0.123 * hc)) * math.log(abs(zm - d) / (0.0123 * hc)) / (0.4 * 0.4 * uz)
    rsa = (200.0 / lai) / ra
    # pet * lam * (delta + gamma)          = delta * RG + gamma * 6.43 * (1 + 0.536 U2) * D
    # pet * lam * (delta + gamma (1 + r))  = delta * RG + rho Cp 86400 / ra * D
    a1 = gamma * 6.43 * (1.0 + 0.536 * u2)
    a2 = rho * 1.013e-3 * 86400.0 / ra
    b1 = pet_mm_day * lam * (delta + gamma)
    b2 = pet_mm_day * lam * (delta + gamma * (1.0 + rsa))
    deficit = (b2 - b1) / (a2 - a1)
    rh = np.clip(1.0 - deficit / es, 0.01, 1.0)
    deficit = es * (1.0 - rh)
    # What is left of the open-water demand is radiation's, W/m2 at the top of the
    # albedo and ground heat terms tReadForcing applies.
    available = np.maximum(b1 - a1 * deficit, 0.0) / delta
    ground = 0.4 * math.exp(-0.5 * lai)
    rn = available / (1.0 - ground) / (1.0 - albedo) / 0.0864
    return rh, rn


def write_inputs(root: Path, forcing: list[dict], step_min: int, aquifer_m: float, params: dict) -> dict:
    folder = root / "input" / PROJECT
    folder.mkdir(parents=True)
    (root / "output").mkdir()
    length, hill = REACH_LENGTH_M, HILLSLOPE_M

    # Nodes: the river line, then each bank's top. Elevation falls towards the river
    # and towards the outlet at x = length.
    def z(x: float, y: float) -> float:
        return Z_OUTLET_M + RIVER_SLOPE * (length - x) + HILL_SLOPE * abs(y)
    nodes = [(0.0, 0.0), (length, 0.0), (length, hill), (0.0, hill), (length, -hill), (0.0, -hill)]
    # Counter-clockwise triangles. A neighbour in column j lies across the edge
    # opposite node j. Elements 1 and 4 hold the river on one edge, 2 and 3 drain to
    # them. As in upstream's meshes the river lies along an edge between neighbours.
    triangles = [(1, 2, 3), (1, 3, 4), (1, 6, 5), (1, 5, 2)]
    neighbours = []
    for i, tri in enumerate(triangles):
        row = []
        for j in range(3):
            edge = {tri[(j + 1) % 3], tri[(j + 2) % 3]}
            other = [k + 1 for k, t in enumerate(triangles) if k != i and edge <= set(t)]
            row.append(other[0] if other else 0)
        neighbours.append(row)
    area = 0.5 * length * hill
    mesh_area = 4.0 * area

    lines = [f"{len(triangles)}\t8", "ID\tNode1\tNode2\tNode3\tNabr1\tNabr2\tNabr3\tZmax"]
    zc = []
    for i, tri in enumerate(triangles):
        zs = sum(z(*nodes[k - 1]) for k in tri) / 3.0
        zc.append(zs)
        lines.append(f"{i+1}\t{tri[0]}\t{tri[1]}\t{tri[2]}\t{neighbours[i][0]}\t{neighbours[i][1]}\t{neighbours[i][2]}\t{zs:.4f}")
    lines += [f"{len(nodes)}\t5", "ID\tX\tY\tAqDepth\tElevation"]
    for k, (x, y) in enumerate(nodes):
        lines.append(f"{k+1}\t{x:.3f}\t{y:.3f}\t{aquifer_m:.6f}\t{z(x, y):.4f}")
    (folder / f"{PROJECT}.sp.mesh").write_text("\n".join(lines) + "\n")

    lines = ["4\t9", "INDEX\tSOIL\tGEOL\tLC\tFORC\tMF\tBC\tSS\tilake"]
    lines += [f"{i+1}\t1\t1\t1\t1\t1\t0\t0\t0" for i in range(4)]
    (folder / f"{PROJECT}.sp.att").write_text("\n".join(lines) + "\n")

    # Upstream's river type 1 with its depth carried in proportion to the aquifer.
    depth = RIVER_DEPTH_OVER_AQUIFER * aquifer_m
    lines = ["1\t6", "Index\tDown\tType\tSlope\tLength\tBC", f"1\t-3\t1\t{RIVER_SLOPE}\t{length}\t0"]
    lines += ["1\t9", "Index\tDepth\tBankSlope\tWidth\tSinuosity\tManning\tCwr\tKsatH\tBedThick",
              f"1\t{depth:.6f}\t0\t2\t1\t0.04\t0.6\t0.1\t0.1"]
    lines += ["1\t6", "From.x\tFrom.y\tFrom.z\tTo.x\tTo.y\tTo.z", f"0\t0\t{z(0, 0):.3f}\t{length}\t0\t{z(length, 0):.3f}"]
    (folder / f"{PROJECT}.sp.riv").write_text("\n".join(lines) + "\n")
    lines = ["2\t4", "Index\tiRiv\tiEle\tLength", f"1\t1\t1\t{length}", f"2\t1\t4\t{length}"]
    (folder / f"{PROJECT}.sp.rivseg").write_text("\n".join(lines) + "\n")

    def table(name: str, header: str, values) -> None:
        body = "\t".join(f"{v:.8g}" for v in values)
        (folder / f"{PROJECT}.{name}").write_text(f"1\t{len(values) + 1}\n{header}\n1\t{body}\n")
    table("para.soil", "INDEX\tKsatV\tThetaS\tThetaR\tInfD\tAlpha\tBeta\thAreaF\tmacKsatV", params["soil"])
    table("para.geol", "INDEX\tKsatH\tKsatV\tThetaS\tThetaR\tvAreaF\tmacKsatH\tDmac", params["geol"])
    table("para.lc", "INDEX\tALBEDO\tVEGFRAC\tROUGH\tRZD\tSOILDGRD\tIMPAF", params["land"])

    n = len(forcing)
    n_days = n * step_min / 1440.0
    # Flat over the case: a table can only hold a calendar if it is given one.
    for name, label, value in (("tsd.lai", "X1", params["lai"]), ("tsd.mf", "MF", params["mf"]), ("tsd.rl", "X1", 0.1)):
        (folder / f"{PROJECT}.{name}").write_text(
            f"2\t2\t20000101\nTime_Day\t{label}\n0\t{value}\n{n_days + 2}\t{value}\n"
        )

    calib = calibration()
    (folder / f"{PROJECT}.cfg.calib").write_text("\n".join(f"{k}\t{v}" for k, v in calib.items()) + "\n")
    outputs = {k: 0 for k in (
        "DT_YE_SNOW", "DT_YE_IC", "DT_YE_SURF", "DT_YE_UNSAT", "DT_YE_GW", "DT_QE_ET", "DT_QE_PRCP",
        "DT_QE_INFIL", "DT_QE_RECH", "DT_QE_SURF", "DT_QE_SUB", "DT_QE_RSURF", "DT_QE_RSUB", "DT_YR_STAGE",
        "DT_QR_DOWN", "DT_QR_SURF", "DT_QR_SUB", "DT_QR_UP")}
    for key in ("DT_YE_SNOW", "DT_YE_IC", "DT_YE_SURF", "DT_YE_UNSAT", "DT_YE_GW", "DT_QE_ET", "DT_QE_PRCP", "DT_YR_STAGE", "DT_QR_DOWN"):
        outputs[key] = SOLVER_STEP_MIN
    if DIAG:
        for key in ("DT_QE_INFIL", "DT_QE_RECH", "DT_QE_SURF", "DT_QE_SUB", "DT_QR_SURF", "DT_QR_SUB", "DT_QR_UP"):
            outputs[key] = SOLVER_STEP_MIN
    para = {
        "VERBOSE": 0, "INIT_MODE": 2, "ASCII_OUTPUT": 0, "BINARY_OUTPUT": 1, "SPINUPDAY": 0, "NUM_OPENMP": 1,
        "SCR_INTV": 10 ** 9, "ABSTOL": ABSTOL, "RELTOL": RELTOL, "INIT_SOLVER_STEP": 1,
        "MAX_SOLVER_STEP": SOLVER_STEP_MIN, "LSM_STEP": 60, "START": 0, "END": repr(days_not_before(n * step_min)), **outputs,
    }
    (folder / f"{PROJECT}.cfg.para").write_text("\n".join(f"{k} {v}" for k, v in para.items()) + "\n")

    (folder / f"{PROJECT}.tsd.forc").write_text(
        f"1 20000101\n./input/{PROJECT}/\nID\tLon\tLat\tX\tY\tZ\tFilename\n1\t0\t0\t0\t0\t-9999\tforcing.csv\n"
    )
    tas = np.array([float(row["tas"]) for row in forcing])
    pet = np.array([float(row["pet"]) for row in forcing])
    rh, rn = mock_meteorology(pet, tas, float(np.mean(zc)), params["land"], params["lai"])
    with open(folder / "forcing.csv", "w") as fh:
        fh.write(f"{n}\t6\t20000101\t20000102\nTime_Day\tAPCP\tTMP\tRH\tWIND\tRN\n")
        for i, row in enumerate(forcing):
            fh.write(f"{days_not_after(i * step_min)!r}\t{float(row['pr']):.8f}\t{tas[i]:.6f}\t{rh[i]:.6f}\t{WIND_MS}\t{rn[i]:.6f}\n")
    return {"mesh_area_m2": mesh_area, "river_depth_m": depth, "mesh_bed_m": (zc[0] + zc[3]) / 2.0 - depth, "mock_rn_max_w_m2": float(rn.max()), "mock_rn_mean_w_m2": float(rn.mean()), "mock_rh_mean": float(rh.mean()), "mock_rh_min": float(rh.min()), "element_z": zc}


def days_not_after(minutes: float) -> float:
    """Days as a double whose minutes do not pass `minutes`, so a row's weather starts on time.

    SHUD compares its clock with the file's time in minutes, exactly, so a boundary written
    even a little late holds the last row's weather for another solver step.
    """
    days = minutes / 1440.0
    return float(np.nextafter(days, 0.0)) if days * 1440.0 > minutes else days


def days_not_before(minutes: float) -> float:
    """Days as a double whose minutes reach `minutes`, as SHUD truncates END to whole minutes."""
    days = minutes / 1440.0
    return float(np.nextafter(days, np.inf)) if days * 1440.0 < minutes else days


class PipeReader(threading.Thread):
    """Reduces one of upstream's binary outputs to one array per row as it is written."""

    def __init__(self, path: Path, samples_per_row: int, n_rows: int, reduce: str, keep_mean: bool = False) -> None:
        super().__init__(daemon=True)
        os.mkfifo(path)
        self.path, self.samples, self.n_rows, self.reduce = path, samples_per_row, n_rows, reduce
        self.keep_mean = keep_mean
        self.result: np.ndarray | None = None
        self.means: np.ndarray | None = None
        self.error: BaseException | None = None

    def run(self) -> None:
        try:
            with open(self.path, "rb", buffering=1 << 20) as fh:
                fh.read(1024 + 8)
                n_var = int(np.frombuffer(fh.read(8), dtype=np.float64)[0])
                fh.read(8 * n_var)
                width = n_var + 1
                result = np.full((self.n_rows, n_var), np.nan)
                # `last` keeps the store at the row's end; a stage also wants the row's mean.
                means = np.full((self.n_rows, n_var), np.nan) if self.keep_mean else None
                for row in range(self.n_rows):
                    raw = fh.read(8 * width * self.samples)
                    block = np.frombuffer(raw, dtype=np.float64)
                    if block.size != width * self.samples:
                        break
                    block = block.reshape(self.samples, width)[:, 1:]
                    result[row] = block.mean(axis=0) if self.reduce == "mean" else block[-1]
                    if means is not None:
                        means[row] = block.mean(axis=0)
                while fh.read(1 << 20):
                    pass
            self.result, self.means = result, means
        except BaseException as exc:  # reported by the caller once the model has exited
            self.error = exc


def parse_time(text: str) -> datetime:
    return datetime.fromisoformat(text.strip())


def simulate(forcing: list[dict], static: dict, timestep: str) -> tuple[list[dict], dict]:
    step_min = MINUTES[timestep]
    days = step_min / 1440.0
    n = len(forcing)
    params = representative_parameters()
    calib = calibration()
    sy = calib["GEOL_THETAS"] * params["geol"][2] - calib["GEOL_THETAR"] * params["geol"][3]

    capacity = static.get("soil_capacity_mm")
    aquifer_m = float(capacity) / 1000.0 / sy if capacity is not None else DEFAULT_AQUIFER_M
    if OPT_AQUIFER:
        aquifer_m = float(OPT_AQUIFER)
    params["land"][3] = min(params["land"][3], MAX_ROOT_SHARE_OF_AQUIFER * aquifer_m)
    samples = step_min // SOLVER_STEP_MIN

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        info = write_inputs(root, forcing, step_min, aquifer_m, params)
        out = root / "output" / f"{PROJECT}.out"
        out.mkdir()
        readers = {
            name: PipeReader(out / f"{PROJECT}.{name}.dat", samples, n, reduce, keep_mean=name == "rivystage")
            for name, reduce in (
                ("eleysnow", "last"), ("eleyic", "last"), ("eleysurf", "last"), ("eleyunsat", "last"), ("eleygw", "last"),
                ("rivystage", "last"), ("eleveta", "mean"), ("elevetp", "mean"), ("elevetic", "mean"), ("elevettr", "mean"),
                ("elevetev", "mean"), ("elevprcp", "mean"), ("elevnetprcp", "mean"), ("rivqdown", "mean"),
            )
        }
        if DIAG:
            for name in ("elevinfil", "elevexfil", "elevrech", "eleqsurf", "eleqsub", "rivqsurf", "rivqsub", "rivqup"):
                readers[name] = PipeReader(out / f"{PROJECT}.{name}.dat", samples, n, "mean")
        for reader in readers.values():
            reader.start()
        with open(root / "shud.log", "w") as log:
            done = subprocess.run(["/model/shud", PROJECT], cwd=root, stdout=log, stderr=subprocess.STDOUT)
        tail = "\n".join((root / "shud.log").read_text(errors="replace").splitlines()[-25:])
        if os.environ.get("HT_SHUD_KEEP"):
            import shutil
            for sub in ("input", "shud.log"):
                src = root / sub
                dst = Path(os.environ["HT_SHUD_KEEP"]) / sub
                shutil.copytree(src, dst, dirs_exist_ok=True) if src.is_dir() else shutil.copy(src, dst)
        if done.returncode != 0:
            raise SystemExit(f"shud failed with {done.returncode}:\n{tail}")
        for reader in readers.values():
            reader.join(timeout=120)
        data = {}
        for name, reader in readers.items():
            if reader.is_alive() or reader.error is not None or reader.result is None:
                raise SystemExit(f"could not read {name}: {reader.error}\n{tail}")
            if np.isnan(reader.result).any():
                raise SystemExit(f"{name} ended early\n{tail}")
            data[name] = reader.result
            if reader.means is not None:
                stage_mean = reader.means

    if DIAG:
        print(json.dumps({k: [round(float(x), 2) for x in v.sum(axis=0) * days * (1000.0 if not k.startswith(("riv", "eleq")) else 1000.0 / info["mesh_area_m2"])]
                          for k, v in data.items()}), file=sys.stderr)
    # Depths over the mesh, in mm; every element has the same area.
    mesh_area = info["mesh_area_m2"]
    mean = lambda name: data[name].mean(axis=1)
    soil = mean("eleyunsat") * sy * 1000.0
    ground = mean("eleygw") * sy * 1000.0
    snow = mean("eleysnow") * 1000.0
    canopy = mean("eleyic") * 1000.0
    pond = mean("eleysurf") * 1000.0
    width = 2.0 + calib["RIV_WDTH+"]
    river_mm = data["rivystage"][:, 0] * width * REACH_LENGTH_M / mesh_area * 1000.0
    # Upstream prints fluxes in m/day and river flow in m3/day.
    et = mean("eleveta") * 1000.0
    out_mm_day = data["rivqdown"][:, 0] / mesh_area * 1000.0

    # The datum a case declares for the bed, or the mesh's own: the river bed lies the
    # river's depth below the banks' mean surface.
    bed = float(static["bed_elevation_m"]) if "bed_elevation_m" in static else info["mesh_bed_m"]
    area_km2 = float(static["area_km2"])
    rows = [
        {
            "time": forcing[i]["time"],
            "pr": forcing[i]["pr"],
            "evspsbl": et[i],
            "mrro": out_mm_day[i],
            "dis": out_mm_day[i] * area_km2 / 86.4,
            "mrso": soil[i],
            "snw": snow[i],
            "canopy": canopy[i],
            "gw": ground[i],
            "channel": pond[i] + river_mm[i],
            "stage": bed + stage_mean[i, 0],
        }
        for i in range(n)
    ]

    pr = np.array([float(row["pr"]) for row in forcing])
    storage = soil + ground + snow + canopy + pond + river_mm
    residual = (pr[1:] - et[1:] - out_mm_day[1:]) * days - np.diff(storage)
    checks = {
        "aquifer_depth_m": aquifer_m,
        "specific_yield": sy,
        "root_depth_m": float(params["land"][3]),
        **{k: v for k, v in info.items() if k != "element_z"},
        "budget_residual_share_of_rain": float(residual.sum() / max((pr[1:] * days).sum(), 1e-9)),
        "max_row_budget_residual_mm": float(np.abs(residual).max()),
        "pet_printed_over_forcing": float(mean("elevetp").sum() * 1000.0 * days / max(sum(float(r["pet"]) for r in forcing) * days, 1e-9)),
        "prcp_printed_over_forcing": float(mean("elevprcp").sum() * 1000.0 / max(pr.sum(), 1e-9)),
    }
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
    if timestep not in MINUTES:
        raise SystemExit(f"unsupported timestep {timestep!r}")
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
        "mocked_forcing": {
            "WIND": "1.82 m/s every row, the example's mean",
            "RH, RN": "from that row's pet and tas, solving SHUD's own two demand equations so both equal pet",
        },
        "not_reported": {
            "snm": "qEleNetPrep is net of canopy capture, and the snow module's own output is not printed",
            "hfls, hfss, hfg and other energy terms": "SHUD computes no surface energy budget",
        },
        "notes": {
            "mesh": "four triangles and one reach, sized from upstream's example; the stated area enters only dis",
            "mrso": "unsaturated zone water, Yus * Sy, over the mesh",
            "gw": "saturated zone water, Ygw * Sy, over the mesh",
            "channel": "ponded surface water and river channel water, over the mesh",
            "mrro": "the reach's downstream flow, as a depth over the mesh",
            "stage": "bed_elevation_m, or the mesh's own bed, plus the reach's depth; the reach is the mesh's and dis is scaled to the stated area",
            "seed": "unused; the model is deterministic",
        },
        "checks": checks,
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

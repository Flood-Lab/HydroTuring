#!/usr/bin/env python3
"""HydroTuring adapter for HYPE 5.36.0, standard library only.

The physics is upstream's hype executable from the release archive, unmodified. The
configuration is upstream's demo2 from the same release: its par.txt, GeoClass.txt and
the GeoData.txt row of subbasin 3672, the demo's largest subbasin with neither a lake
nor a wetland. Those files are baked into the image unchanged, and every difference
for a case is named below.

HYPE's own basin output is what gets reported, one file for the one subbasin, at the
step of the case. Nothing is derived from a residual.

Step: HYPE's soil, snow and routing rates are per timestep, so a daily parameter file
is a statement about a day. The evaluation's native step is PT1D, where every value is
upstream's. At any other step (the resolution probe asks for the next finer one
whatever the manifest says) the rates are rescaled to the step here, because HYPE does
not do it itself:

    cmlt, cevp, mperc1, mperc2, mactrinf, rrcs3    amounts per step    x dt_days
    rrcs1, rrcs2, trrcs, srrcs                     fractions per step  1 - (1 - k)^dt_days

Forcing goes in as depth per row, as HYPE wants it, and fluxes come back as rates.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

MODEL = {"name": "hype", "version": "5.36.0.1"}
HYPE = os.environ.get("HYPE_EXE", "/model/hype")
TEMPLATE = Path(os.environ.get("HYPE_TEMPLATE", "/model/template"))
SUBID = "3672"

COLUMNS = ["time", "pr", "snm", "evspsbl", "mrro", "dis", "mrso", "snw", "canopy", "channel"]
SECONDS_PER_DAY = 86400.0
STEP_SECONDS = {"PT1D": 86400, "PT1H": 3600, "PT15M": 900, "PT5M": 300, "PT1M": 60}
STEP_CODE = {"PT1D": "1d", "PT1H": "1h", "PT15M": "15min", "PT5M": "5min", "PT1M": "1min"}

# HYPE's basin output codes this adapter reads, in the order it asks for them.
OUTPUTS = ["prec", "cprc", "cprf", "melt", "evap", "epot", "icpe", "crun", "cout", "soim", "snow", "clrv", "cmrv"]

AMOUNT_PER_STEP = {"cmlt", "cevp", "mperc1", "mperc2", "mactrinf", "rrcs3"}
FRACTION_PER_STEP = {"rrcs1", "rrcs2", "trrcs", "srrcs"}


def is_comment(line: str) -> bool:
    return line.startswith("!")


def par_txt(dt_days: float) -> tuple[str, list]:
    """Demo2's par.txt, with the per-step rates rescaled when the step is not a day."""
    out, changed = [], []
    for line in (TEMPLATE / "par.txt").read_text(encoding="latin-1").splitlines():
        parts = line.split()
        name = parts[0].lower() if parts and not is_comment(line) else ""
        if dt_days != 1.0 and name in AMOUNT_PER_STEP | FRACTION_PER_STEP:
            values = [float(v) for v in parts[1:]]
            if name in AMOUNT_PER_STEP:
                new = [v * dt_days for v in values]
            else:
                new = [1.0 - (1.0 - v) ** dt_days for v in values]
            changed.append({"parameter": name, "upstream": values, "used": new})
            line = "\t".join([parts[0]] + [f"{v:.9g}" for v in new])
        out.append(line)
    return "\n".join(out) + "\n", changed


def soil_capacity(rows: list[list[str]], fractions: dict[int, float], par: str) -> float:
    """What the class mixture can hold, mm: porosity share x soil depth, area weighted."""
    wc = {}
    for line in par.splitlines():
        parts = line.split()
        if parts and parts[0].lower() in ("wcwp", "wcfc", "wcep"):
            wc[parts[0].lower()] = [float(v) for v in parts[1:]]
    total = 0.0
    for row in rows:
        slc, soil, layers = int(row[0]), int(row[2]), int(row[10])
        if fractions.get(slc, 0.0) <= 0.0:
            continue
        depth = float(row[10 + layers])
        pore = sum(wc[key][soil - 1] for key in wc)
        total += fractions[slc] * pore * depth * 1000.0
    return total


def geoclass_txt(scale: float) -> str:
    """Demo2's GeoClass.txt, every soil depth multiplied by one factor."""
    out = []
    for line in (TEMPLATE / "GeoClass.txt").read_text(encoding="latin-1").splitlines():
        if not line.strip() or is_comment(line):
            out.append(line)
            continue
        parts = line.split()
        for k in (8, 9, 11, 12, 13):
            if k < len(parts):
                parts[k] = f"{float(parts[k]) * scale:.9g}"
        out.append("\t".join(parts))
    return "\n".join(out) + "\n"


def geodata_row() -> tuple[list[str], list[str]]:
    lines = (TEMPLATE / "GeoData.txt").read_text(encoding="latin-1").splitlines()
    header = lines[0].split("\t")
    for line in lines[1:]:
        row = line.split("\t")
        if row[header.index("SUBID")] == SUBID:
            return header, row
    raise SystemExit(f"subbasin {SUBID} is not in demo2's GeoData.txt")


def write_series(path: Path, stamps: list[str], values: list[float]) -> None:
    lines = [f"DATE\t{SUBID}"]
    lines += [f"{s}\t{v:.9g}" for s, v in zip(stamps, values)]
    path.write_text("\n".join(lines) + "\n")


def configure(tmp: Path, forcing: list[dict], static: dict, timestep: str, dt_days: float) -> dict:
    times = [datetime.fromisoformat(r["time"]) for r in forcing]
    subdaily = dt_days != 1.0
    fmt = "%Y-%m-%d %H:%M" if subdaily else "%Y-%m-%d"
    stamps = [t.strftime(fmt) for t in times]

    par, par_changes = par_txt(dt_days)
    (tmp / "par.txt").write_text(par)

    header, row = geodata_row()
    values = dict(zip(header, row))
    fractions = {int(k[4:]): float(v) for k, v in values.items() if k.startswith("SLC_")}
    area_km2 = float(static["area_km2"])
    values["MAINDOWN"] = "0"
    values["AREA"] = f"{area_km2 * 1.0e6:.6f}"
    # The case's channel lengths when it has them. Otherwise demo2's main river, and
    # HYPE's own default for the local stream, the square root of the land area.
    if static.get("main_channel_length_km") is not None:
        values["RIVLEN"] = f"{float(static['main_channel_length_km']) * 1000.0:.6f}"
    if static.get("centroid_channel_length_km") is not None:
        header = header + ["LOC_RIVLEN"]
        values["LOC_RIVLEN"] = f"{float(static['centroid_channel_length_km']) * 1000.0:.6f}"
    (tmp / "GeoData.txt").write_text("\t".join(header) + "\n" + "\t".join(values[h] for h in header) + "\n")

    classes = [
        line.split()
        for line in (TEMPLATE / "GeoClass.txt").read_text(encoding="latin-1").splitlines()
        if line.strip() and not is_comment(line)
    ]
    native = soil_capacity(classes, fractions, par)
    wanted = static.get("soil_capacity_mm")
    scale = 1.0 if wanted is None else float(wanted) / native
    (tmp / "GeoClass.txt").write_text(geoclass_txt(scale))

    write_series(tmp / "Pobs.txt", stamps, [float(r["pr"]) * dt_days for r in forcing])
    write_series(tmp / "Tobs.txt", stamps, [float(r["tas"]) for r in forcing])
    # petmodel 0's documented replacement: an Xobs series of repo takes the place of
    # HYPE's temperature formula for potential evaporation, in mm per timestep.
    xobs = ["!! Potential evaporation prescribed by the case", "name\trepo", f"0\t{SUBID}"]
    xobs += [f"{s}\t{float(r['pet']) * dt_days:.9g}" for s, r in zip(stamps, forcing)]
    (tmp / "Xobs.txt").write_text("\n".join(xobs) + "\n")

    info = [
        f"bdate\t{stamps[0]}",
        f"cdate\t{stamps[0]}",
        f"edate\t{stamps[-1]}",
        f"steplength\t{STEP_CODE[timestep]}",
        "readdaily\tyes",
        "readobsid\tno",
        f"basinoutput variable\t{' '.join(OUTPUTS)}",
        # 0 is every timestep. Demo2 asks for 1, a daily mean, which at a daily step is
        # the same row.
        "basinoutput meanperiod\t0",
        "basinoutput decimals\t6",
        f"basinoutput subbasin\t{SUBID}",
    ]
    (tmp / "info.txt").write_text("\n".join(info) + "\n")
    return {
        "subbasin": SUBID,
        "river_length_m": values["RIVLEN"],
        "soil_depth_scale": scale,
        "soil_capacity_demo2_mm": native,
        "soil_capacity_used_mm": native * scale,
        "par_changes": par_changes,
    }


def run_hype(tmp: Path) -> str:
    done = subprocess.run(
        [HYPE, f"{tmp}/"], cwd=tmp,
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, errors="replace",
    )
    logs = sorted(tmp.glob("hyss*.log"))
    log = logs[-1].read_text(errors="replace") if logs else ""
    # HYPE reports a bad setup in its log, and a run that never wrote its output has failed.
    if done.returncode != 0 or not (tmp / f"{int(SUBID):07d}.txt").exists():
        raise SystemExit(f"hype failed ({done.returncode}):\n{done.stdout[-1500:]}\n{log[-3000:]}")
    return log


def read_output(tmp: Path, n: int) -> dict[str, list[float]]:
    """HYPE's basin file, by column name. A sub-daily date is two words, the rest are one."""
    with open(tmp / f"{int(SUBID):07d}.txt", newline="") as fh:
        rows = [line.split("\t") for line in fh.read().splitlines() if line.strip()]
    header = [h.strip().lower() for h in rows[0]]
    body = rows[2:]
    if len(body) != n:
        raise SystemExit(f"hype wrote {len(body)} rows for {n} forcing rows")
    missing = [name for name in OUTPUTS if name not in header]
    if missing:
        raise SystemExit(f"hype did not write {missing}: {header}")
    return {name: [float(r[header.index(name)]) for r in body] for name in OUTPUTS}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    io_dir = request_path.parent
    request = json.loads(request_path.read_text())
    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        forcing = list(csv.DictReader(fh))
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    timestep = request.get("timestep", "PT1D")
    if timestep not in STEP_SECONDS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    step_seconds = STEP_SECONDS[timestep]
    dt_days = step_seconds / SECONDS_PER_DAY
    area_km2 = float(static["area_km2"])
    n = len(forcing)

    with tempfile.TemporaryDirectory() as name:
        tmp = Path(name)
        config = configure(tmp, forcing, static, timestep, dt_days)
        log = run_hype(tmp)
        out = read_output(tmp, n)

    to_mm = 1000.0 / (area_km2 * 1.0e6)
    per_day = 1.0 / dt_days
    table = io_dir / request["output"]["table"]
    table.parent.mkdir(parents=True, exist_ok=True)
    with open(table, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(COLUMNS)
        for k, row in enumerate(forcing):
            q = out["cout"][k]
            writer.writerow([
                row["time"],
                row["pr"],
                f"{(out['cprf'][k] + out['melt'][k]) * per_day:.8g}",
                f"{out['evap'][k] * per_day:.8g}",
                f"{q * SECONDS_PER_DAY * to_mm:.8g}",
                f"{q:.8g}",
                f"{out['soim'][k]:.8g}",
                f"{out['snow'][k]:.8g}",
                0.0,
                f"{(out['clrv'][k] + out['cmrv'][k]) * to_mm:.8g}",
            ])

    rain_mm = sum(float(r["pr"]) * dt_days for r in forcing)
    et_mm = sum(out["evap"])
    q_mm = sum(out["cout"]) * step_seconds * to_mm
    held = [out["soim"][k] + out["snow"][k] + (out["clrv"][k] + out["cmrv"][k]) * to_mm for k in range(n)]
    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "n_steps": n,
        "step_seconds": step_seconds,
        "template": "demo2 of release 5.36.0: par.txt, GeoClass.txt, GeoData.txt subbasin 3672",
        "config": config,
        "zero_stores": {
            "canopy": "HYPE has no interception store. icpe is a flux set by a precipitation correction "
                      "(pcluse), which demo2 leaves at zero, so nothing is intercepted at all",
        },
        "not_reported": {
            "gw": "no store below the soil column: demo2 uses no aquifer, and the three soil layers are all under mrso",
        },
        "mapping": {
            "snm": "cpRF + melt: rain goes straight to the soil and meltwater leaves the pack; HYPE holds no liquid in snow",
            "evspsbl": "evap, the subbasin evaporation, mm per step over the step length",
            "mrro": "cout over area_km2",
            "dis": "cout",
            "mrso": "soim, every soil layer including water below wilting point",
            "snw": "snow",
            "channel": "clrv + cmrv over area_km2",
            "pet": "Xobs repo, petmodel 0",
        },
        "water_budget": {
            "precipitation_mm": rain_mm,
            "hype_prec_mm": sum(out["prec"]),
            "hype_corrected_prec_mm": sum(out["cprc"]),
            "hype_epot_mm": sum(out["epot"]),
            "hype_icpe_mm": sum(out["icpe"]),
            "hype_crun_mm": sum(out["crun"]),
            "evspsbl_mm": et_mm,
            "mrro_mm": q_mm,
            "store_change_first_to_last_row_mm": held[-1] - held[0],
            "residual_after_first_row_mm": (
                sum(float(r["pr"]) * dt_days for r in forcing[1:])
                - sum(out["evap"][1:])
                - sum(out["cout"][1:]) * step_seconds * to_mm
                - (held[-1] - held[0])
            ),
        },
        "hype_warnings": sorted({line.strip()[:200] for line in log.splitlines() if "WARNING" in line})[:20],
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

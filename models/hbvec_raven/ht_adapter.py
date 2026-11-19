#!/usr/bin/env python3
"""HydroTuring adapter for HBV-EC as emulated by Raven, standard library only.

The physics is upstream's Raven executable at the pinned commit. The configuration
is Raven's own HBV-EC tutorial model, Alouette2 from RavenTutorialFiles.zip: its
Alouette2.rvi, .rvp, .rvh and .rvc are baked into the image unchanged, and every
difference from them for a case is one entry in the edit lists below, each asserted
to match exactly once. The gauge file is written here, from the case's forcing.

What Raven writes, not what this file derives, is what gets reported:

    WatershedStorage.nc    every store, end of step, in mm, and Raven's own MB error
    Hydrographs.nc         the outflow, averaged over the step, in m3/s
    AET custom output      the water that left to the atmosphere in the step, in mm

Raven's CSV outputs print six significant figures, so a 300 mm store is known to
0.001 mm: enough to fail a closure on a 0.02 mm drizzle, and enough for a change of
area to flip a last digit in a probe that asks for one part in a million. Its NetCDF
outputs are doubles, and :WriteNetCDFFormat is the only output option changed for that.
ncdump, which comes with the library, reads them back.

Elevation: the tutorial's 21 non-glacier HRUs span seven elevation bands, 161 m to
1500 m, and OROCORR_HBV corrects each HRU's temperature, precipitation and PET by its
height above the gauge. All three corrections are linear in that height, so with the
gauge at the HRUs' area-weighted mean elevation the catchment receives exactly the
case's precipitation and PET, spread over the bands as the tutorial spreads them. The
precipitation correction steepens above HBVEC_LAPSE_ELEV, which the tutorial leaves at
Raven's 5000 m default; rvh_hrus asserts both that and that every HRU is below it.

Step: Raven's rates are all in mm/d and its solver takes :TimeStep directly, so no
parameter is rescaled here. Forcing is written at the case's own step, as rates.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

MODEL = {"name": "hbvec_raven", "version": "4.12-81011f0-alouette2.1"}
RAVEN = os.environ.get("RAVEN_EXE", "/model/Raven")
TEMPLATE = Path(os.environ.get("RAVEN_TEMPLATE", "/model/template"))
RUN = "Alouette2"
NCDUMP = os.environ.get("NCDUMP", "ncdump")

COLUMNS = ["time", "pr", "snm", "evspsbl", "mrro", "dis", "mrso", "snw", "canopy", "gw", "channel"]
SECONDS_PER_DAY = 86400.0
STEP_SECONDS = {"PT1D": 86400, "PT1H": 3600, "PT15M": 900, "PT5M": 300, "PT1M": 60}

# Alouette values the template carries and a case replaces only if it says so.
TEMPLATE_THICKNESS_M = 1.0
TEMPLATE_POROSITY = 0.2
TEMPLATE_LATITUDE = 49.3
# The tutorial's HRUs, and the two of them that are glacier.
TEMPLATE_HRUS = 23
GLACIER_HRUS = {"34", "35"}
# Raven's default HBVEC_LAPSE_ELEV (GlobalParams.cpp:324), above which OROCORR_HBV's
# precipitation correction is no longer linear in height.
HBVEC_LAPSE_ELEV_M = 5000.0

# WatershedStorage.nc variable -> reported store. Anything else it holds with water in
# it stops the run, so a store can never be dropped silently.
STORES = {
    "Soil Water[0]": "mrso",
    "Soil Water[1]": "gw",
    "Soil Water[2]": "gw",
    "Snow": "snw",
    "Snow Melt (Liquid)": "snw",
    "Canopy": "canopy",
    "Canopy Snow": "canopy",
    "channel_storage": "channel",
    "reservoir_storage": "channel",
    "rivulet_storage": "channel",
    "Surface Water": "channel",
    "Ponded Water": "channel",
}
# Compartments the template declares that no remaining HRU fills once the glacier
# HRUs are gone. Checked to stay zero rather than assumed to.
EMPTY = {"Glacier Liquid Storage", "Glacier Ice"}
# Not stores: cumulative counters and bookkeeping.
BOOKKEEPING = {"time", "rainfall", "snowfall", "total", "cum_input", "cum_outflow", "MB_error", "Cum. Losses to Atmosphere"}


def replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise SystemExit(f"template edit does not match exactly once: {old!r}")
    return text.replace(old, new)


def read_template(ext: str) -> str:
    # The tutorial files have Windows line endings; the edits are written without them.
    return (TEMPLATE / f"{RUN}.{ext}").read_text().replace("\r\n", "\n")


def rvi_edits(start: str, n_steps: int, step_days: float) -> list[tuple[str, str, str]]:
    return [
        (":StartDate          1991-10-01 00:00:00", f":StartDate          {start}", "the case's first row"),
        (":Duration           1096 00:00:00", f":Duration           {n_steps * step_days:.10f}", "the case's length"),
        (":TimeStep           24:00:00", f":TimeStep           {step_days:.12f}", "the case's step"),
        (":Evaporation         PET_FROMMONTHLY", ":Evaporation         PET_DATA",
         "the case prescribes PET; the template derives it from Alouette monthly normals"),
        (":OW_Evaporation      PET_FROMMONTHLY", ":OW_Evaporation      PET_DATA",
         "the same series; :LakeEvaporation on SLOW_RESERVOIR draws on it (README.md)"),
        (":DefineHRUGroups Open_N, Open_S, Forest, Glacier", ":DefineHRUGroups Open_N, Open_S, Forest",
         "the Glacier reporting group loses its HRUs with the glacier (rvh)"),
        (":EvaluationMetrics NASH_SUTCLIFFE", ":WriteNetCDFFormat\n:CustomOutput CONTINUOUS AVERAGE AET ENTIRE_WATERSHED\n"
         ":CustomOutput CONTINUOUS AVERAGE Between:ATMOS_PRECIP.And.PONDED_WATER ENTIRE_WATERSHED\n"
         ":CustomOutput CONTINUOUS AVERAGE Between:SNOW_LIQ.And.PONDED_WATER ENTIRE_WATERSHED",
         "output only: no observed flow to score; doubles rather than six figures, and the per-step "
         "fluxes this adapter reports"),
    ]


def rvp_edits(thickness_m: float, canopy_mm: float | None) -> list[tuple[str, str, str]]:
    edits = [(
        "DEFAULT_P, 3, TOPSOIL, 1.0, FAST_RES,", f"DEFAULT_P, 3, TOPSOIL, {thickness_m:.10f}, FAST_RES,",
        "soil_capacity_mm from the case, as thickness x 1000 x POROSITY",
    )]
    if canopy_mm is not None:
        edits.append((
            "FOREST,          0.0,        10000,            10000,", f"FOREST,          0.0,        {canopy_mm:.6f},            {canopy_mm:.6f},",
            "canopy_capacity_mm from the case; the template's 10000 mm means no limit",
        ))
    return edits


RVH_EDITS = [
    ("  31, 32, 33, 34, 35\n", "  31, 32, 33\n", "Band7 loses the glacier HRUs"),
    (":HRUGroup Glacier\n  34, 35\n:EndHRUGroup\n", "", "the Glacier reporting group, now empty"),
]

RVC_EDITS = [
    ("   1, 1\n", "   1, 0\n",
     "the template's 1 m3/s is Alouette's; Raven would release it in the first step without "
     "taking it from any store"),
    (", 6.50657, 0,0\n", ", 6.50657\n", "the two glacier HRUs' entries, gone with them"),
]


def rvh_hrus(text: str, area_km2: float, latitude: float) -> tuple[str, float, dict]:
    """The HRU table with the glacier HRUs removed and the rest scaled to the case's area.

    Every HRU keeps its elevation, land use, aspect and share of the area. Returns the
    edited text and the remaining HRUs' area-weighted mean elevation.
    """
    head, sep, rest = text.partition(":HRUs\n")
    body, end, tail = rest.partition(":EndHRUs")
    if not sep or not end:
        raise SystemExit("template rvh has no :HRUs table")
    lines = body.split("\n")
    rows = [(i, [f.strip() for f in line.split(",")]) for i, line in enumerate(lines)
            if line.split(",")[0].strip().isdigit()]
    if len(rows) != TEMPLATE_HRUS:
        raise SystemExit(f"template rvh has {len(rows)} HRUs, expected {TEMPLATE_HRUS}")
    glacier = [f[0] for _, f in rows if f[6] == "GLACIER"]
    if set(glacier) != GLACIER_HRUS:
        raise SystemExit(f"template rvh glacier HRUs are {glacier}, expected {sorted(GLACIER_HRUS)}")
    kept = [(i, f) for i, f in rows if f[0] not in GLACIER_HRUS]
    if max(float(f[2]) for _, f in kept) >= HBVEC_LAPSE_ELEV_M:
        raise SystemExit("an HRU reaches HBVEC_LAPSE_ELEV; the gauge elevation would not keep the precipitation total")
    template_area = sum(float(f[1]) for _, f in kept)
    scale = area_km2 / template_area
    elevation = sum(float(f[1]) * float(f[2]) for _, f in kept) / template_area
    drop = {i for i, f in rows if f[0] in GLACIER_HRUS}
    for i, f in kept:
        f = list(f)
        f[1] = repr(float(f[1]) * scale)
        f[3] = f"{latitude:.6f}"
        lines[i] = "   \t " + ",\t ".join(f)
    lines = [line for i, line in enumerate(lines) if i not in drop]
    record = {
        "template": f"{TEMPLATE_HRUS} HRUs, {template_area:.4f} km2 without the glacier, latitude {TEMPLATE_LATITUDE}",
        "used": f"{len(kept)} HRUs, each area x {scale:.10g}, latitude {latitude:.6f}",
        "why": "area and latitude from the case; the glacier HRUs 34 and 35 are removed because the "
               "case has no glacier. Elevations, land use, aspects and area shares are the template's",
    }
    return head + sep + "\n".join(lines) + end + tail, elevation, record


def write_rvt(path: Path, forcing: list[dict], start: str, step_days: float, latitude: float, elevation: float) -> None:
    """The template's gauge, with the case's series in place of Alouette's.

    The template's monthly normals fed PET_FROMMONTHLY and are dropped with it. Its
    :RainCorrection 1.25885 and :SnowCorrection 1.04283 are gauge-undercatch factors
    for the Alouette gauge. The case's rain is exact, so both are left at Raven's
    default of 1. The gauge stands at the HRUs' mean elevation (module docstring).
    """
    lines = [
        ":Gauge ALUcomp",
        f"  :Latitude               {latitude:.6f}",
        "  :Longitude             -100",
        f"  :Elevation              {elevation:.10f}",
        "  :RainCorrection         1.0",
        "  :SnowCorrection         1.0",
        ":MultiData",
        f" {start} {step_days:.12f} {len(forcing)}",
        ":Parameters TEMP_AVE PRECIP PET OW_PET",
        ":Units C mm/d mm/d mm/d",
    ]
    lines += [f"{r['tas']} {r['pr']} {r['pet']} {r['pet']}" for r in forcing]
    lines += [":EndMultiData", ":EndGauge", ""]
    path.write_text("\n".join(lines))


def configure(tmp: Path, forcing: list[dict], static: dict, step_days: float) -> dict:
    if "HBVEC_LAPSE" in read_template("rvp"):
        raise SystemExit("the template sets an HBVEC_LAPSE parameter; the gauge elevation assumes Raven's defaults")
    start = datetime.fromisoformat(forcing[0]["time"]).strftime("%Y-%m-%d %H:%M:%S")
    capacity = static.get("soil_capacity_mm")
    thickness = TEMPLATE_THICKNESS_M if capacity is None else float(capacity) / (1000.0 * TEMPLATE_POROSITY)
    canopy = static.get("canopy_capacity_mm")
    latitude = float(static.get("latitude_deg", TEMPLATE_LATITUDE))
    plan = {
        "rvi": rvi_edits(start, len(forcing), step_days),
        "rvp": rvp_edits(thickness, None if canopy is None else float(canopy)),
        "rvh": RVH_EDITS,
        "rvc": RVC_EDITS,
    }
    record = {}
    texts = {}
    for ext, edits in plan.items():
        text = read_template(ext)
        for old, new, _ in edits:
            text = replace_once(text, old, new)
        texts[ext] = text
        record[ext] = [{"template": o.strip(), "used": n.strip(), "why": w} for o, n, w in edits]
    texts["rvh"], elevation, hrus = rvh_hrus(texts["rvh"], float(static["area_km2"]), latitude)
    record["rvh"].insert(0, hrus)
    for ext, text in texts.items():
        (tmp / f"{RUN}.{ext}").write_text(text)
    write_rvt(tmp / f"{RUN}.rvt", forcing, start, step_days, latitude, elevation)
    record["rvt"] = [
        {"template": "gauge ALUcomp at 582.358 m, RainCorrection 1.25885, SnowCorrection 1.04283, "
                     "monthly normals, Alouette series in ALUGauge.rvt",
         "used": f"gauge ALUcomp at {elevation:.4f} m, RainCorrection 1, SnowCorrection 1, no normals, "
                 "the case's TEMP_AVE PRECIP PET",
         "why": "see write_rvt"},
    ]
    return record


def run_raven(tmp: Path) -> None:
    out = tmp / "out"
    out.mkdir()
    done = subprocess.run(
        [RAVEN, str(tmp / RUN), "-o", f"{out}/"], cwd=tmp,
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    errors = out / "Raven_errors.txt"
    log = errors.read_text() if errors.exists() else ""
    # Raven exits 0 after an error and says so only in this file.
    if done.returncode != 0 or "ERROR" in log:
        raise SystemExit(f"Raven failed ({done.returncode}):\n{log[-3000:]}\n{done.stdout[-1000:]}")


def ncdump(path: Path) -> dict[str, list[float]]:
    """The numeric variables of a NetCDF file, every value at full precision."""
    text = subprocess.run(
        [NCDUMP, "-p", "17,17", str(path)], check=True, stdout=subprocess.PIPE, text=True,
    ).stdout
    variables: dict[str, list[float]] = {}
    for chunk in text.split("data:", 1)[1].split(";"):
        if "=" not in chunk:
            continue
        name, _, body = chunk.partition("=")
        values: list[float] = []
        for token in body.split(","):
            token = token.strip()
            try:
                values.append(float("nan") if token == "_" else float(token))
            except ValueError:
                values = []
                break
        if values:
            variables[name.strip().replace("\\", "")] = values
    return variables


def custom_output(out: Path, pattern: str, n: int) -> list[float]:
    """One value a step from a :CustomOutput file, whose only series besides time."""
    matches = glob.glob(str(out / pattern))
    if len(matches) != 1:
        raise SystemExit(f"expected one custom output for {pattern}, found {sorted(os.listdir(out))}")
    series = {k: v for k, v in ncdump(Path(matches[0])).items() if k != "time"}
    if len(series) != 1:
        raise SystemExit(f"{matches[0]} holds {sorted(series)}")
    values = next(iter(series.values()))
    if len(values) != n:
        raise SystemExit(f"{matches[0]} has {len(values)} rows for {n} forcing rows")
    return values


def steps(cumulative: list[float]) -> list[float]:
    """Mm in each step, from a :CustomOutput series Raven keeps cumulatively.

    Raven's BETWEEN fluxes are labelled mm/d and are in fact the running sum in mm
    (CHydroUnit::GetCumulFluxBet).
    """
    return [b - a for a, b in zip([0.0] + cumulative[:-1], cumulative)]


def stores(out: Path, n: int) -> tuple[dict, dict]:
    """Every store at the end of every step, from Raven's own table, the initial state first."""
    table = ncdump(out / "WatershedStorage.nc")
    totals = {key: [0.0] * (n + 1) for key in ("mrso", "snw", "canopy", "gw", "channel")}
    for name, values in table.items():
        if name in BOOKKEEPING:
            continue
        if len(values) != n + 1:
            raise SystemExit(f"WatershedStorage {name} has {len(values)} rows for {n} forcing rows")
        if name in STORES:
            column = totals[STORES[name]]
            for k, v in enumerate(values):
                column[k] += v
        elif name in EMPTY:
            if any(abs(v) > 1.0e-9 for v in values):
                raise SystemExit(f"{name} holds water ({max(values)}); it has no reported home")
        else:
            raise SystemExit(f"WatershedStorage variable {name!r} is not mapped; found {sorted(table)}")
    return totals, {
        "precip_mm_per_day": [a + b for a, b in zip(table["rainfall"][1:], table["snowfall"][1:])],
        "raven_mb_error_mm": table["MB_error"][-1],
        "cum_inputs_mm": table["cum_input"][-1],
        "cum_outflow_mm": table["cum_outflow"][-1],
    }


def discharge(out: Path, n: int) -> list[float]:
    q = ncdump(out / "Hydrographs.nc")["q_sim"]
    if len(q) == n + 1:
        return q[1:]
    if len(q) != n:
        raise SystemExit(f"Hydrographs has {len(q)} rows for {n} forcing rows")
    return q


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
    step_days = step_seconds / SECONDS_PER_DAY
    area_km2 = float(static["area_km2"])
    n = len(forcing)

    with tempfile.TemporaryDirectory() as name:
        tmp = Path(name)
        edits = configure(tmp, forcing, static, step_days)
        run_raven(tmp)
        out = tmp / "out"
        state, check = stores(out, n)
        q = discharge(out, n)
        aet_mm = custom_output(out, "*AET*Continuous*", n)
        rain_to_ground = steps(custom_output(out, "*BETWEEN_ATMOS_PRECIP_AND_PONDED*Continuous*", n))
        melt_to_ground = steps(custom_output(out, "*BETWEEN_SNOW_LIQ_AND_PONDED*Continuous*", n))

    to_mm_per_day = SECONDS_PER_DAY * 1000.0 / (area_km2 * 1.0e6)
    table = io_dir / request["output"]["table"]
    table.parent.mkdir(parents=True, exist_ok=True)
    with open(table, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(COLUMNS)
        for k, row in enumerate(forcing):
            writer.writerow([
                row["time"],
                row["pr"],
                f"{(rain_to_ground[k] + melt_to_ground[k]) / step_days:.15g}",
                f"{aet_mm[k] / step_days:.15g}",
                f"{q[k] * to_mm_per_day:.15g}",
                f"{q[k]:.15g}",
                f"{state['mrso'][k + 1]:.15g}",
                f"{state['snw'][k + 1]:.15g}",
                f"{state['canopy'][k + 1]:.15g}",
                f"{state['gw'][k + 1]:.15g}",
                f"{state['channel'][k + 1]:.15g}",
            ])

    rain_mm = sum(float(r["pr"]) * step_days for r in forcing)
    et_mm = sum(aet_mm)
    q_mm = sum(q) * to_mm_per_day * step_days
    held = [sum(state[key][k] for key in state) for k in range(n + 1)]
    d_store = held[-1] - held[0]
    raven_rain = sum(check["precip_mm_per_day"]) * step_days
    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "n_steps": n,
        "step_seconds": step_seconds,
        "template": "RavenTutorialFiles.zip Alouette2/Alouette2.rv{i,p,h,c}, the HBV-EC tutorial model",
        "edits_from_template": edits,
        "zero_stores": {},
        "not_reported": {},
        "mapping": {
            "evspsbl": "AET custom output, mm per step over the step length",
            "mrro": "Hydrographs.nc q_sim over area_km2",
            "dis": "Hydrographs.nc q_sim",
            "snm": "BETWEEN ATMOS_PRECIP and PONDED_WATER + BETWEEN SNOW_LIQ and PONDED_WATER",
            **{f"{v}<-{k}": k for k, v in STORES.items()},
        },
        "water_budget": {
            "precipitation_mm": rain_mm,
            "raven_precipitation_mm": raven_rain,
            "evspsbl_mm": et_mm,
            "mrro_mm": q_mm,
            "store_change_mm": d_store,
            "residual_mm": rain_mm - et_mm - q_mm - d_store,
            "raven_mb_error_mm": check["raven_mb_error_mm"],
            "raven_cum_inputs_mm": check["cum_inputs_mm"],
            "raven_cum_outflow_mm": check["cum_outflow_mm"],
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

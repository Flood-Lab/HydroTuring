#!/usr/bin/env python3
"""HydroTuring adapter for NOAA-OWP CFE, standard library only.

The physics is upstream's libcfebmi at the pinned commit, run by ht_cfe_driver.
This file is the units layer: it writes CFE's config for the case, hands the
forcing over in CFE's units for the row, and reports every store CFE has.

Configuration: upstream's configs/cfe_config_cat_87_pass.txt, Schaake partition
and Nash-cascade surface runoff, as in issue #127. Soil capacity from the case.

Step: CFE is hourly. Its constants are per hour or per step, and its surface
cascade hard-codes a one-hour step. At the hourly step every value below is
upstream's, unscaled. At any other step (the resolution probe asks for the next
finer one whatever the manifest says) those constants are rescaled to the step
here, because the model has no way to do that itself:

    Cgw, K_lf                 amounts per hour             x dt_hours
    K_nash_subsurface         fraction per hour            1 - (1 - k)^dt_hours
    K_nash_surface, Kinf      per hour, against a dt_h=1   x dt_hours
    nsubsteps_nash_surface    keeps K/n below one          max(10, ceil(10 dt_hours))

satdk, the Schaake term and PET take the step from time_step_size themselves.
"""

from __future__ import annotations

import argparse
import array
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

MODEL = {"name": "cfe", "version": "f9182df-ht.1"}
DRIVER = os.environ.get("CFE_DRIVER", "/model/ht_cfe_driver")
COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "mrso", "gw", "channel"]
SECONDS_PER_DAY = 86400.0
STEP_SECONDS = {"PT1D": 86400, "PT1H": 3600, "PT15M": 900, "PT5M": 300, "PT1M": 60}

# The order ht_cfe_driver.c writes them in.
OUTPUTS = [
    "q_out", "aet", "direct", "nash_lateral", "baseflow", "excess",
    "soil", "gw", "nash_surface", "nash_subsurface",
    "et_rain", "et_retention", "et_soil", "infiltration", "runon",
    "percolation", "lateral", "volin", "volout",
]

# Upstream's configs/cfe_config_cat_87_pass.txt at the pinned commit, in its order,
# which matters: gw_storage is read as a fraction of the max_gw_storage already seen.
UPSTREAM = {
    "soil_params.depth": 2.0,
    "soil_params.b": 4.05,
    "soil_params.satdk": 3.38e-06,
    "soil_params.satpsi": 0.355,
    "soil_params.slop": 0.01,
    "soil_params.smcmax": 0.439,
    "soil_params.wltsmc": 0.066,
    "soil_params.expon": 1.0,
    "soil_params.expon_secondary": 1.0,
    "max_gw_storage": 0.25,
    "Cgw": 1.8e-05,
    "expon": 6.0,
    "gw_storage": 0.125,
    "alpha_fc": 0.33,
    "soil_storage": 0.585626,
    "K_nash_subsurface": 0.03,
    "K_lf": 0.01,
    "nash_storage_subsurface": "0.0,0.0",
    "surface_runoff_scheme": "NASH_CASCADE",
    "N_nash_surface": 2,
    "K_nash_surface": 0.83089,
    "nsubsteps_nash_surface": 10,
    # Not in the shipped config: the defaults its parser falls back on.
    "Kinf_nash_surface": 0.001,
    "retention_depth_nash_surface": 0.001,
    "nash_storage_surface": "0.0,0.0",
    "surface_water_partitioning_scheme": "Schaake",
}


def configure(static: dict, step_seconds: int) -> tuple[dict, dict]:
    """CFE's config for this case, and a record of what differs from upstream."""
    config = dict(UPSTREAM)
    changed: dict[str, dict] = {}

    def change(key, value, why):
        changed[key] = {"upstream": config[key], "used": value, "why": why}
        config[key] = value

    capacity_mm = static.get("soil_capacity_mm")
    if capacity_mm is not None:
        # CFE's capacity is porosity x depth; the depth is the free one.
        change(
            "soil_params.depth",
            float(capacity_mm) / 1000.0 / float(config["soil_params.smcmax"]),
            "soil_capacity_mm from static.json as smcmax x depth",
        )

    dt_hours = step_seconds / 3600.0
    if dt_hours != 1.0:
        why = f"rescaled from one hour to the case's step of {dt_hours:g} h"
        change("Cgw", config["Cgw"] * dt_hours, why)
        change("K_lf", config["K_lf"] * dt_hours, why)
        change("K_nash_subsurface", 1.0 - (1.0 - config["K_nash_subsurface"]) ** dt_hours, why)
        change("K_nash_surface", config["K_nash_surface"] * dt_hours, why)
        change("Kinf_nash_surface", config["Kinf_nash_surface"] * dt_hours, why)
        change(
            "nsubsteps_nash_surface",
            max(config["nsubsteps_nash_surface"], math.ceil(10 * dt_hours)),
            why,
        )
    return config, changed


def write_config(path: Path, config: dict) -> None:
    lines = ["forcing_file=BMI"]
    for key, value in config.items():
        lines.append(f"{key}={value!r}" if isinstance(value, float) else f"{key}={value}")
    path.write_text("\n".join(lines) + "\n")


def run_cfe(config: dict, step_seconds: int, forcing: list[dict]) -> tuple[list[list[float]], str]:
    """Rows of OUTPUTS in metres per row, the initial state first, and CFE's stdout."""
    dt_days = step_seconds / SECONDS_PER_DAY
    feed = array.array("d")
    for row in forcing:
        # Depth for the row in mm (Update divides by 1000), and m/s.
        feed.append(float(row["pr"]) * dt_days)
        feed.append(float(row["pet"]) / 1000.0 / SECONDS_PER_DAY)
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        write_config(tmp_path / "cfe.cfg", config)
        with open(tmp_path / "forcing.bin", "wb") as fh:
            feed.tofile(fh)
        done = subprocess.run(
            [DRIVER, str(tmp_path / "cfe.cfg"), str(step_seconds),
             str(tmp_path / "forcing.bin"), str(tmp_path / "output.bin")],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        if done.returncode != 0:
            raise SystemExit(f"ht_cfe_driver exited {done.returncode}: {done.stderr[-2000:]}")
        values = array.array("d")
        size = os.path.getsize(tmp_path / "output.bin")
        with open(tmp_path / "output.bin", "rb") as fh:
            values.fromfile(fh, size // values.itemsize)
    width = len(OUTPUTS)
    rows = [values[i:i + width].tolist() for i in range(0, len(values), width)]
    if len(rows) != len(forcing) + 1:
        raise SystemExit(f"ht_cfe_driver wrote {len(rows)} rows for {len(forcing)} forcing rows")
    return rows, done.stdout


def budget(rows: list[list[float]], forcing: list[dict], step_seconds: int, stdout: str) -> dict:
    """What CFE's reported columns do with the water, and what CFE says it did."""
    index = {name: i for i, name in enumerate(OUTPUTS)}
    dt_days = step_seconds / SECONDS_PER_DAY
    rain_mm = sum(float(r["pr"]) * dt_days for r in forcing)

    def stored(row):
        return sum(row[index[k]] for k in ("soil", "gw", "nash_surface", "nash_subsurface"))

    et_mm = sum(r[index["aet"]] for r in rows[1:]) * 1000.0
    q_mm = sum(r[index["q_out"]] for r in rows[1:]) * 1000.0
    d_store_mm = (stored(rows[-1]) - stored(rows[0])) * 1000.0
    worst = 0.0
    steps_off = 0
    for before, after, row in zip(rows[:-1], rows[1:], forcing):
        rain = float(row["pr"]) * dt_days / 1000.0
        off = rain - after[index["aet"]] - after[index["q_out"]] - (stored(after) - stored(before))
        worst = max(worst, abs(off))
        steps_off += abs(off) > 1.0e-9
    return {
        "reported_columns": {
            "precipitation_mm": rain_mm,
            "evspsbl_mm": et_mm,
            "mrro_mm": q_mm,
            "store_change_mm": d_store_mm,
            "residual_mm": rain_mm - et_mm - q_mm - d_store_mm,
            "worst_step_residual_mm": worst * 1000.0,
            "steps_off_by_more_than_1e-6_mm": steps_off,
        },
        "cfe_stdout_warnings": sorted({
            line.strip() for line in stdout.splitlines() if "WARNING" in line or "problem" in line
        }),
        "cfe_mass_balance_check": [
            line.strip() for line in stdout.splitlines() if "residual" in line.lower()
        ],
    }


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

    timestep = request.get("timestep", "PT1H")
    if timestep not in STEP_SECONDS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    step_seconds = STEP_SECONDS[timestep]
    per_day = SECONDS_PER_DAY / step_seconds
    area_km2 = float(static["area_km2"])

    config, changed = configure(static, step_seconds)
    rows, stdout = run_cfe(config, step_seconds, forcing)

    index = {name: i for i, name in enumerate(OUTPUTS)}
    table = io_dir / request["output"]["table"]
    table.parent.mkdir(parents=True, exist_ok=True)
    with open(table, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(COLUMNS)
        for source, row in zip(forcing, rows[1:]):
            mrro = row[index["q_out"]] * 1000.0 * per_day
            writer.writerow([
                source["time"],
                source["pr"],
                f"{row[index['aet']] * 1000.0 * per_day:.15g}",
                f"{mrro:.15g}",
                f"{mrro * 1.0e-3 * area_km2 * 1.0e6 / SECONDS_PER_DAY:.15g}",
                f"{row[index['soil']] * 1000.0:.15g}",
                f"{row[index['gw']] * 1000.0:.15g}",
                f"{(row[index['nash_surface']] + row[index['nash_subsurface']]) * 1000.0:.15g}",
            ])

    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "n_steps": len(forcing),
        "step_seconds": step_seconds,
        "config": config,
        "changed_from_upstream_config": changed,
        "not_reported": {
            "snw": "CFE has no snow module; all precipitation enters as rain",
            "canopy": "CFE has no interception store; et_from_rainfall takes PET from the row's rain",
        },
        "mapping": {
            "mrro": "Q_OUT",
            "dis": "Q_OUT over area_km2",
            "evspsbl": "ACTUAL_ET",
            "mrso": "SOIL_STORAGE",
            "gw": "GW_STORAGE",
            "channel": "nash_surface_params.nash_storage + nash_storage_subsurface",
        },
        "water_budget": budget(rows, forcing, step_seconds, stdout),
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

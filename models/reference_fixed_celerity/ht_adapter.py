#!/usr/bin/env python3
"""Fixed-celerity negative control for the wave-celerity probe."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path

import numpy as np

MODEL = {"name": "reference_fixed_celerity", "version": "1.0.0"}
SECONDS_PER_DAY = 86400.0
CELERITY_M_S = 1.0


def _dt_seconds(forcing: list[dict]) -> float:
    if len(forcing) < 2:
        raise ValueError("fixed-celerity reference needs at least two forcing rows")
    t0 = datetime.fromisoformat(forcing[0]["time"])
    t1 = datetime.fromisoformat(forcing[1]["time"])
    dt = (t1 - t0).total_seconds()
    if dt <= 0:
        raise ValueError("forcing timestamps must increase")
    return dt


def simulate(forcing: list[dict], static: dict) -> list[dict]:
    area_km2 = float(static["area_km2"])
    reach_length_m = float(static["reach_length_m"])
    # Read and validate the remaining declared hydraulic geometry. It is
    # intentionally not allowed to affect the broken model's routing speed.
    if str(static["cross_section_shape"]).strip().lower() != "rectangular":
        raise ValueError("reference_fixed_celerity requires rectangular geometry")
    for key in ("width_m", "slope", "manning_n"):
        if float(static[key]) <= 0:
            raise ValueError(f"{key} must be positive")
    float(static["bed_elevation_m"])

    dt_s = _dt_seconds(forcing)
    q = np.asarray([
        float(row["q_in"])
        if "q_in" in row and str(row["q_in"]).strip() != ""
        else (
            max(float(row["pr"]) - float(row["pet"]), 0.0)
            * 1.0e-3 * area_km2 * 1.0e6 / SECONDS_PER_DAY
        )
        for row in forcing
    ], dtype=float)
    delay_s = reach_length_m / CELERITY_M_S
    delay_steps = delay_s / dt_s
    x = np.arange(len(q), dtype=float)
    # Fractional-step causal delay. Before the record starts, extend the first
    # base discharge so spinup does not manufacture an artificial front.
    shifted = x - delay_steps
    routed = np.interp(shifted, x, q, left=float(q[0]), right=float(q[-1]))
    return [{"time": row["time"], "dis": float(value)}
            for row, value in zip(forcing, routed)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent
    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        forcing = list(csv.DictReader(fh))
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    rows = simulate(forcing, static)

    output = io_dir / request["output"]["table"]
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["time", "dis"])
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "routing": {"celerity_m_s": CELERITY_M_S, "state_dependent": False},
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

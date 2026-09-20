#!/usr/bin/env python3
"""Finite-volume one-dimensional Saint-Venant reference adapter.

The state is water depth ``h`` and discharge per unit width ``q``.  The
adapter advances the conservative shallow-water equations

    d_t h + d_x q = 0
    d_t q + d_x(q^2 / h + g h^2 / 2) = g h (S_0 - S_f)

on a rectangular prismatic reach.  Interface fluxes use the local
Lax-Friedrichs (Rusanov) approximate Riemann solver.  Bed slope is explicit
and Manning friction is relaxed semi-implicitly, so the long-time solution
is selected by the momentum equation rather than by a normal-depth lookup.

The upstream ghost cell prescribes only discharge; both depth and discharge
are extrapolated at the downstream boundary.  No function in this module
evaluates or inverts the Manning normal-depth formula.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

MODEL = {"name": "reference_saint_venant", "version": "1.0.0"}
SECONDS_PER_DAY = 86400.0
GRAVITY = 9.80665
N_CELLS = 64
CFL = 0.45
MIN_DEPTH_M = 1.0e-4
MAX_STEPS = 120_000
CHECK_EVERY = 100
RELATIVE_STEADY_TOLERANCE = 2.0e-9
REQUIRED_STEADY_CHECKS = 4


def _rusanov_flux(
    h_left: np.ndarray,
    q_left: np.ndarray,
    h_right: np.ndarray,
    q_right: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return mass and momentum flux at every cell interface."""
    u_left = q_left / h_left
    u_right = q_right / h_right
    c_left = np.sqrt(GRAVITY * h_left)
    c_right = np.sqrt(GRAVITY * h_right)
    speed = np.maximum(np.abs(u_left) + c_left, np.abs(u_right) + c_right)

    mass = 0.5 * (q_left + q_right) - 0.5 * speed * (h_right - h_left)
    momentum_left = q_left * u_left + 0.5 * GRAVITY * h_left * h_left
    momentum_right = q_right * u_right + 0.5 * GRAVITY * h_right * h_right
    momentum = (
        0.5 * (momentum_left + momentum_right)
        - 0.5 * speed * (q_right - q_left)
    )
    return mass, momentum


def _advance(
    depth: np.ndarray,
    unit_discharge: np.ndarray,
    prescribed_unit_discharge: float,
    width_m: float,
    slope: float,
    manning_n: float,
    dx_m: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Advance one CFL-limited finite-volume step."""
    speed = np.max(
        np.abs(unit_discharge / depth) + np.sqrt(GRAVITY * depth)
    )
    dt_s = CFL * dx_m / max(float(speed), 1.0e-6)

    h_ext = np.empty(len(depth) + 2, dtype=float)
    q_ext = np.empty(len(depth) + 2, dtype=float)
    h_ext[1:-1] = depth
    q_ext[1:-1] = unit_discharge

    # One condition enters a subcritical reach at each end.  Upstream
    # discharge is prescribed while depth is extrapolated.  The downstream
    # boundary is open (zero gradient), so no Manning-derived depth is
    # supplied anywhere in the solve.
    h_ext[0] = depth[0]
    q_ext[0] = prescribed_unit_discharge
    h_ext[-1] = depth[-1]
    q_ext[-1] = unit_discharge[-1]

    mass_flux, momentum_flux = _rusanov_flux(
        h_ext[:-1], q_ext[:-1], h_ext[1:], q_ext[1:]
    )
    ratio = dt_s / dx_m
    new_depth = depth - ratio * (mass_flux[1:] - mass_flux[:-1])
    new_depth = np.maximum(new_depth, MIN_DEPTH_M)

    provisional_q = (
        unit_discharge
        - ratio * (momentum_flux[1:] - momentum_flux[:-1])
        + dt_s * GRAVITY * new_depth * slope
    )

    area = width_m * new_depth
    hydraulic_radius = area / (width_m + 2.0 * new_depth)
    # Backward Euler on dq/dt = -K q |q| has a closed positive root.  Using
    # that root, instead of linearising the denominator with provisional_q,
    # preserves the uniform-flow equilibrium for any CFL-limited dt.
    friction_k = (
        GRAVITY
        * manning_n**2
        / (new_depth * hydraulic_radius ** (4.0 / 3.0))
    )
    magnitude = np.abs(provisional_q)
    new_magnitude = (
        2.0 * magnitude
        / (1.0 + np.sqrt(1.0 + 4.0 * dt_s * friction_k * magnitude))
    )
    new_q = np.sign(provisional_q) * new_magnitude
    return new_depth, new_q, dt_s


def solve_steady_reach(
    discharge_m3s: float,
    width_m: float,
    slope: float,
    manning_n: float,
    reach_length_m: float,
    initial_state: tuple[np.ndarray, np.ndarray] | None = None,
    n_cells: int = N_CELLS,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Integrate the full dynamic equations until the reach is steady."""
    if discharge_m3s <= 0.0:
        raise ValueError("the Saint-Venant reference requires positive discharge")
    if min(width_m, slope, manning_n, reach_length_m) <= 0.0:
        raise ValueError("reach width, slope, roughness and length must be positive")

    if n_cells < 8:
        raise ValueError("the Saint-Venant reference needs at least eight cells")
    dx_m = reach_length_m / n_cells
    prescribed_q = discharge_m3s / width_m
    if initial_state is None:
        # Deliberately not a Manning depth: every seed starts from the same
        # weakly nonuniform, motionless pool and has to find its own balance.
        x = (np.arange(n_cells, dtype=float) + 0.5) / n_cells
        depth = 0.35 * (1.0 + 0.08 * np.sin(2.0 * math.pi * x))
        unit_discharge = np.zeros(n_cells, dtype=float)
    else:
        depth = np.asarray(initial_state[0], dtype=float).copy()
        unit_discharge = np.asarray(initial_state[1], dtype=float).copy()

    stable_checks = 0
    elapsed_s = 0.0
    last_change = float("inf")
    for step in range(1, MAX_STEPS + 1):
        old_depth = depth
        old_q = unit_discharge
        depth, unit_discharge, dt_s = _advance(
            old_depth,
            old_q,
            prescribed_q,
            width_m,
            slope,
            manning_n,
            dx_m,
        )
        elapsed_s += dt_s

        if step % CHECK_EVERY != 0:
            continue
        depth_change = float(
            np.max(np.abs(depth - old_depth)) / max(float(np.mean(depth)), MIN_DEPTH_M)
        )
        q_scale = max(abs(prescribed_q), 1.0e-9)
        q_change = float(np.max(np.abs(unit_discharge - old_q)) / q_scale)
        discharge_mismatch = abs(float(unit_discharge[-1]) - prescribed_q) / q_scale
        last_change = max(depth_change, q_change, discharge_mismatch)
        if last_change <= RELATIVE_STEADY_TOLERANCE:
            stable_checks += 1
            if stable_checks >= REQUIRED_STEADY_CHECKS:
                break
        else:
            stable_checks = 0
    else:
        raise RuntimeError(
            "Saint-Venant solve did not converge after "
            f"{MAX_STEPS} steps (relative change {last_change:.3e})"
        )

    center = n_cells // 2
    diagnostics = {
        "steps": step,
        "pseudo_time_s": elapsed_s,
        "relative_change": last_change,
        "gauge_depth_m": float(depth[center]),
        "gauge_discharge_m3s": float(width_m * unit_discharge[center]),
        "outlet_discharge_m3s": float(width_m * unit_discharge[-1]),
        "cells": n_cells,
        "cfl": CFL,
    }
    return depth, unit_discharge, diagnostics


def simulate(forcing: list[dict], static: dict) -> tuple[list[dict], list[dict]]:
    area_km2 = float(static["area_km2"])
    width_m = float(static["width_m"])
    bed_m = float(static["bed_elevation_m"])
    slope = float(static["slope"])
    manning_n = float(static["manning_n"])
    reach_length_m = float(static["reach_length_m"])
    shape = str(static.get("cross_section_shape", "rectangular")).strip().lower()
    if shape != "rectangular":
        raise ValueError("reference_saint_venant requires a rectangular section")

    state = None
    last_inflow = None
    gauge_discharge = None
    gauge_stage = None
    solves: list[dict] = []
    rows: list[dict] = []
    for item in forcing:
        effective_mm_day = max(float(item["pr"]) - float(item["pet"]), 0.0)
        inflow = effective_mm_day * 1.0e-3 * area_km2 * 1.0e6 / SECONDS_PER_DAY
        if inflow <= 0.0:
            raise ValueError("the generated Saint-Venant reach must stay wet")

        if last_inflow is None or not math.isclose(inflow, last_inflow, rel_tol=1e-12):
            depth, unit_discharge, diagnostics = solve_steady_reach(
                inflow,
                width_m,
                slope,
                manning_n,
                reach_length_m,
                state,
            )
            state = (depth, unit_discharge)
            center = N_CELLS // 2
            gauge_discharge = float(width_m * unit_discharge[center])
            gauge_stage = float(bed_m + depth[center])
            diagnostics["prescribed_inflow_m3s"] = inflow
            solves.append(diagnostics)
            last_inflow = inflow

        rows.append({
            "time": item["time"],
            "dis": gauge_discharge,
            "stage": gauge_stage,
        })
    return rows, solves


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
    rows, solves = simulate(forcing, static)

    output = io_dir / request["output"]["table"]
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["time", "dis", "stage"])
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "n_steps": len(rows),
        "solver": {
            "equations": "one-dimensional Saint-Venant continuity and momentum",
            "flux": "Rusanov finite volume",
            "friction": "semi-implicit Manning source relaxation",
            "normal_depth_lookup": False,
            "solves": solves,
        },
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

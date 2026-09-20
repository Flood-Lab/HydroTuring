"""Verification of the finite-volume Saint-Venant must-pass baseline."""

from __future__ import annotations

import importlib.util
import math

import numpy as np
import pytest

from hydroturing import registry


@pytest.fixture(scope="module")
def solver():
    path = registry.find_model("reference_saint_venant").path / "ht_adapter.py"
    spec = importlib.util.spec_from_file_location("saint_venant_reference", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _friction_ratio(module, depth, unit_discharge, width=100.0, slope=0.0015, n=0.035):
    area = width * depth
    radius = area / (width + 2.0 * depth)
    discharge = width * unit_discharge
    friction = (n * discharge / (area * radius ** (2.0 / 3.0))) ** 2
    return friction / slope


@pytest.mark.parametrize("discharge", [3.0, 10.0, 25.0])
def test_full_equations_converge_to_uniform_friction_balance(solver, discharge):
    assert not hasattr(solver, "normal_depth")
    depth, unit_q, diagnostics = solver.solve_steady_reach(
        discharge, 100.0, 0.0015, 0.035, 6000.0
    )
    ratio = _friction_ratio(solver, depth, unit_q)
    center = len(depth) // 2
    assert ratio[center] == pytest.approx(1.0, rel=2e-6)
    assert diagnostics["relative_change"] <= solver.RELATIVE_STEADY_TOLERANCE
    assert diagnostics["gauge_discharge_m3s"] == pytest.approx(discharge, rel=2e-6)


def test_equilibrium_is_independent_of_initial_state(solver):
    base_depth, base_q, _ = solver.solve_steady_reach(
        10.0, 100.0, 0.0015, 0.035, 6000.0
    )
    x = np.linspace(0.0, 2.0 * math.pi, len(base_depth), endpoint=False)
    initial = (1.25 * base_depth * (1.0 + 0.05 * np.sin(x)), 0.7 * base_q)
    depth, unit_q, _ = solver.solve_steady_reach(
        10.0, 100.0, 0.0015, 0.035, 6000.0, initial_state=initial
    )
    center = len(depth) // 2
    assert depth[center] == pytest.approx(base_depth[center], rel=2e-5)
    assert unit_q[center] == pytest.approx(base_q[center], rel=2e-5)


def test_gauge_solution_is_grid_converged(solver):
    coarse_h, coarse_q, _ = solver.solve_steady_reach(
        10.0, 100.0, 0.0015, 0.035, 6000.0, n_cells=32
    )
    fine_h, fine_q, _ = solver.solve_steady_reach(
        10.0, 100.0, 0.0015, 0.035, 6000.0, n_cells=64
    )
    assert coarse_h[len(coarse_h) // 2] == pytest.approx(
        fine_h[len(fine_h) // 2], rel=2e-5
    )
    assert coarse_q[len(coarse_q) // 2] == pytest.approx(
        fine_q[len(fine_q) // 2], rel=2e-5
    )


def test_simulation_solves_each_distinct_plateau(solver):
    forcing = []
    for day, effective in enumerate((1.0, 1.0, 2.0, 2.0, 4.0, 4.0)):
        forcing.append({
            "time": f"2001-01-{day + 1:02d}",
            "pr": 2.0 + effective,
            "pet": 2.0,
            "tas": 15.0,
        })
    static = {
        "area_km2": 150.0,
        "width_m": 100.0,
        "cross_section_shape": "rectangular",
        "bed_elevation_m": 80.0,
        "slope": 0.0015,
        "manning_n": 0.035,
        "reach_length_m": 6000.0,
    }
    rows, solves = solver.simulate(forcing, static)
    assert len(rows) == len(forcing)
    assert len(solves) == 3
    assert len({round(row["dis"], 8) for row in rows}) == 3
    assert rows[0]["stage"] < rows[2]["stage"] < rows[4]["stage"]

"""Sensible heat stored in a fixed layer must match its boundary heat input."""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL, PASS, CriterionResult, criterion, make_window, segments,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@criterion("soil_heat_storage")
def soil_heat_storage(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Compare interval-mean boundary fluxes with endpoint layer temperatures.

    The fixed areal heat capacity belongs to the configured control volume,
    not to a fit of the reported fluxes. This first case excludes phase change,
    heat advection and internal sources. A passing budget does not establish
    the accuracy of the temperature response or the thermal conductivity.
    """
    top = str(params.get("top", "hfg"))
    bottom = str(params.get("bottom", "hfg_bottom"))
    temperature = str(params.get("temperature", "tsoil_layer"))
    capacity_key = str(params.get("capacity", "soil_heat_capacity_areal"))
    initial_key = str(params.get("initial_temperature", "soil_temperature_initial"))
    segment_column = str(params.get("segment_column", "_phase"))
    threshold = float(params.get("threshold", 0.05))
    floor = float(params.get("floor", 1.0))
    if not np.isfinite(threshold) or threshold < 0 or not np.isfinite(floor) or floor <= 0:
        raise ValueError("soil_heat_storage needs a non-negative threshold and positive floor")
    if capacity_key not in run.case.static:
        raise ValueError(f"soil_heat_storage needs static '{capacity_key}' in J m-2 K-1")
    capacity = float(run.case.static[capacity_key])
    if not np.isfinite(capacity) or capacity <= 0:
        raise ValueError("soil_heat_storage needs a finite positive areal heat capacity")

    w = make_window(run, probe)
    for var in (top, bottom, temperature):
        if var not in w.table.columns:
            raise ValueError(f"soil_heat_storage needs '{var}' in the model result")
    phases = segments(w, segment_column)
    if run.case.spinup_steps:
        prior = float(w.state0[temperature])
    else:
        # With no spinup the first row is already an interval END. Using it
        # as its own starting value would silently discard the first change.
        if initial_key not in run.case.static:
            raise ValueError(f"soil_heat_storage needs static '{initial_key}' without spinup")
        prior = float(run.case.static[initial_key])

    temperatures = w.table[temperature].to_numpy(dtype=float)
    previous = np.r_[prior, temperatures[:-1]]
    top_flux = w.table[top].to_numpy(dtype=float)
    bottom_flux = w.table[bottom].to_numpy(dtype=float)
    finite = (
        np.isfinite(temperatures) & np.isfinite(previous)
        & np.isfinite(top_flux) & np.isfinite(bottom_flux)
    )
    if not finite.all():
        n_bad = int((~finite).sum())
        return CriterionResult(
            name="soil_heat_storage", status=FAIL,
            message=f"non-finite soil heat values on {n_bad} scored intervals",
            diagnostics={"non_finite_steps": n_bad},
        )

    net_input = top_flux - bottom_flux
    storage_rate = capacity * (temperatures - previous) / (w.dt_days * 86400.0)
    residual = net_input - storage_rate
    blocks = []
    for label, start, stop in phases:
        mean_abs = float(np.abs(residual[start:stop]).mean())
        mean_net = float(np.abs(net_input[start:stop]).mean())
        allowance = max(threshold * mean_net, floor)
        blocks.append({
            "label": label, "start": start, "stop": stop,
            "duration_hours": (stop - start) * w.dt_days * 24.0,
            "mean_abs_w_m2": mean_abs,
            "mean_abs_net_input_w_m2": mean_net,
            "mean_abs_storage_rate_w_m2": float(np.abs(storage_rate[start:stop]).mean()),
            "allowance_w_m2": allowance,
            "floor_controls": floor >= threshold * mean_net,
            "slack": mean_abs / allowance,
            "passed": mean_abs <= allowance,
        })
    failed = sum(not block["passed"] for block in blocks)
    worst = max(blocks, key=lambda block: block["slack"])
    return CriterionResult(
        name="soil_heat_storage", status=PASS if failed == 0 else FAIL,
        value=worst["slack"], threshold=1.0,
        message=(
            f"{failed} of {len(blocks)} soil heat phases fail; worst {worst['label']} "
            f"[{worst['start']}:{worst['stop']}] has mean absolute residual "
            f"{worst['mean_abs_w_m2']:.3g} W m-2 against "
            f"{worst['allowance_w_m2']:.3g} W m-2 allowed"
        ),
        diagnostics={
            "blocks": blocks, "failed_blocks": failed, "worst_block": worst,
            "relative_threshold": threshold, "floor_w_m2": floor,
            "heat_capacity_j_m2_k": capacity,
            "max_abs_residual_w_m2": float(np.abs(residual).max()),
        },
    )

"""Inspect a warm single-point HRLDAS run; see docs/noahmp-soil-heat-validation.md.

This is a supplemental native heat-equation diagnostic. Variable heat capacity
does not satisfy the fixed-capacity probe's assumptions, even when a residual
calculated with an approximate constant happens to be small.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
import pandas as pd

from hydroturing.criteria import get
from hydroturing.protocol import Case, RunResult


def phase_summary(hourly: pd.DataFrame, residual: np.ndarray) -> list[dict]:
    blocks = []
    for label, start, stop in (("heating", 24, 36), ("recovery", 36, 96)):
        mean_abs = float(np.abs(residual[start:stop]).mean())
        mean_net = float(np.abs(hourly.net_input_w_m2.iloc[start:stop]).mean())
        allowance = max(0.05 * mean_net, 1.0)
        blocks.append({
            "phase": label, "hours": stop - start,
            "mean_abs_residual_w_m2": mean_abs,
            "mean_abs_net_input_w_m2": mean_net,
            "comparison_allowance_w_m2": allowance,
            "within_comparison_allowance": mean_abs <= allowance,
        })
    return blocks


def evaluate(workdir: Path) -> dict:
    thermal = pd.read_csv(workdir / "native_thermal_steps.csv")
    column = pd.read_csv(workdir / "native_column_steps.csv")
    # These checks catch missing steps or mismatched output files before a
    # temperature difference is incorrectly interpreted as a model residual.
    for name, data in (("thermal", thermal), ("column", column)):
        if not np.isfinite(data.to_numpy(dtype=float)).all():
            raise ValueError(f"{name}: non-finite diagnostic output")
        if not np.array_equal(data.step, np.arange(1, len(data) + 1)):
            raise ValueError(f"{name}: expected consecutive single-point steps")
        if len(data[["grid_i", "grid_j"]].drop_duplicates()) != 1:
            raise ValueError(f"{name}: this experiment requires a single point")
        if not np.all(data.dt_s == 300) or len(data) != 1152:
            raise ValueError(f"{name}: expected 96 hours at 300 seconds")
        if not np.array_equal(data.time_end_s, data.step * data.dt_s):
            raise ValueError(f"{name}: diagnostic timestamps do not match the steps")
    if not np.array_equal(thermal.step, column.step):
        raise ValueError("thermal and full-column records do not align")
    if not np.all(column.soil_step_executed == 1) or not np.all(thermal.opt_stc == 1):
        raise ValueError("expected one prescribed-flux thermal solve per full model step")

    capacity = thermal.heat_capacity_area_j_m2_k.to_numpy()
    if not np.all(capacity > 0) or not np.all(thermal.z1_m == thermal.z1_m.iloc[0]):
        raise ValueError("layer thickness must stay fixed and heat capacity positive")
    previous = np.r_[column.t1_start_k.iloc[0], column.t1_final_k.iloc[:-1]]
    dt = thermal.dt_s.to_numpy()
    delta_thermal = thermal.t1_thermal_end_k - thermal.t1_start_k
    delta_reported = column.t1_final_k.to_numpy() - previous
    native_storage = capacity * delta_thermal / dt
    reported_storage = capacity * delta_reported / dt
    # The source term is logged from the native equation, not inferred by
    # balancing the temperature change. It should be zero in this bare case.
    net_input = thermal.g_top_w_m2 + thermal.penetrated_sw_w_m2 - thermal.g_bottom_w_m2
    step_data = pd.DataFrame({
        "hour": (thermal.step - 1) // 12,
        "g_top_mean_w_m2": thermal.g_top_w_m2,
        "g_bottom_mean_w_m2": thermal.g_bottom_w_m2,
        "net_input_w_m2": net_input,
        "native_storage_rate_w_m2": native_storage,
        "reported_storage_rate_w_m2": reported_storage,
        "native_residual_w_m2": net_input - native_storage,
        "reported_residual_w_m2": net_input - reported_storage,
        "t_end_k": column.t1_final_k,
        "capacity_j_m2_k": capacity,
    })
    hourly = step_data.groupby("hour", sort=True).mean()
    hourly["t_end_k"] = step_data.groupby("hour").t_end_k.last()
    hourly.to_csv(workdir / "native_hourly_budget.csv", index=True)

    # The comparison uses the probe's existing 5% / 1 W m-2 allowance. It
    # checks the native time-discrete C_n * delta(T_n) equation, not the
    # enthalpy change of a soil-water mixture with changing water mass.
    residuals = {
        "native_thermal_equation": hourly.native_residual_w_m2.to_numpy(),
        "full_column_temperature": hourly.reported_residual_w_m2.to_numpy(),
        "frozen_reported_temperature": hourly.net_input_w_m2.to_numpy(),
        "half_reported_temperature_change": (
            hourly.net_input_w_m2 - 0.5 * hourly.reported_storage_rate_w_m2
        ).to_numpy(),
    }
    reasons = []
    if not np.all(capacity == capacity[0]):
        reasons.append("native heat capacity changes with soil water content")
    if np.any(thermal.snow_layers != 0) or np.any(column.snow_swe_mm != 0):
        reasons.append("snow is present")
    if np.any(thermal.soil_ice_max_fraction != 0) or np.any(column.tsoil_min_k <= 273.16):
        reasons.append("ice or freezing conditions are present")
    if np.any(thermal.penetrated_sw_w_m2 != 0):
        reasons.append("penetrating radiation supplies a layer source")
    if np.any(column.heat_precip_advected_w_m2 != 0):
        reasons.append("precipitation heat advection is present")
    if np.any(column.veg_fraction != 0):
        reasons.append("the column is not entirely bare ground")
    post_thermal = column.t1_final_k.to_numpy() - thermal.t1_thermal_end_k.to_numpy()
    if np.any(post_thermal != 0):
        reasons.append("a later process changes the thermally solved temperature")

    # Show the actual fixed-C implementation's arithmetic separately. This
    # replay is explicitly not a physical verdict when its assumptions fail.
    times = pd.date_range("2020-06-20T06:00:00", periods=96, freq="h")
    phases = ["spinup"] * 24 + ["heating"] * 12 + ["recovery"] * 60
    case = Case(
        probe_id="validation/hrldas-noahmp-full-column", seed=0,
        forcing=pd.DataFrame({"time": times, "_phase": phases}),
        static={"soil_heat_capacity_areal": float(capacity[0])},
        spinup_steps=24, timestep="PT1H",
    )
    table = pd.DataFrame({
        "time": times, "hfg": hourly.g_top_mean_w_m2.to_numpy(),
        "hfg_bottom": hourly.g_bottom_mean_w_m2.to_numpy(),
        "tsoil_layer": hourly.t_end_k.to_numpy(),
    })
    result = get("soil_heat_storage")(
        RunResult(case, table, {}, 0.0), None, {"threshold": 0.05, "floor": 1.0}
    )
    report = {
        "scope": "supplemental native heat-equation check during a full HRLDAS/Noah-MP run",
        "formal_probe_gate": False,
        "steps": len(thermal), "native_dt_seconds": 300, "hours": 96,
        "layer_depth_m": float(-thermal.z1_m.iloc[0]),
        "native_heat_capacity_range_j_m2_k": [float(capacity.min()), float(capacity.max())],
        "native_heat_capacity_relative_range": float(np.ptp(capacity) / capacity[0]),
        "soil_moisture_range": [float(column.soil_moisture_1_end.min()),
                                float(column.soil_moisture_1_end.max())],
        "soil_temperature_range_k": [float(column.t1_final_k.min()), float(column.t1_final_k.max())],
        "max_post_thermal_temperature_change_k": float(np.abs(post_thermal).max()),
        "max_native_step_residual_w_m2": float(np.abs(net_input - native_storage).max()),
        "max_hrldas_energy_balance_error_w_m2": float(np.abs(column.energy_balance_error_w_m2).max()),
        "max_hrldas_water_balance_error_mm": float(np.abs(column.water_balance_error_mm).max()),
        "warm_bare_case_checks": {
            "max_snow_swe_mm": float(column.snow_swe_mm.max()),
            "max_soil_ice_fraction": float(column.soil_ice_max_fraction.max()),
            "max_vegetation_fraction": float(column.veg_fraction.max()),
            "max_abs_penetrated_sw_w_m2": float(np.abs(thermal.penetrated_sw_w_m2).max()),
            "max_abs_precip_advected_heat_w_m2": float(np.abs(column.heat_precip_advected_w_m2).max()),
        },
        "native_equation_comparisons": {
            name: phase_summary(hourly, values) for name, values in residuals.items()
        },
        "fixed_capacity_criterion_replay": {
            "constant_capacity_used_j_m2_k": float(capacity[0]),
            "applicable_to_full_run": not reasons,
            "inapplicability_reasons": reasons,
            "interpretation": "arithmetic replay only when assumptions fail; not a model physical verdict",
            "result": asdict(result),
        },
    }
    (workdir / "full_model_validation.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.workdir)
    for name, blocks in report["native_equation_comparisons"].items():
        print(name + ": " + "; ".join(
            f"{b['phase']} {b['mean_abs_residual_w_m2']:.6g} W m-2 "
            f"(allowance {b['comparison_allowance_w_m2']:.6g})" for b in blocks
        ))
    fixed = report["fixed_capacity_criterion_replay"]
    print("Fixed-capacity probe applicable:", fixed["applicable_to_full_run"])
    for reason in fixed["inapplicability_reasons"]:
        print("  " + reason)
    print("Report:", args.workdir / "full_model_validation.json")


if __name__ == "__main__":
    main()

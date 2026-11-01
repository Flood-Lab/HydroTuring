"""Head-derived storage bounds for a passive, source-free confined aquifer.

Positive transmissivity and a connected river boundary keep head between its
uniform initial value and the prescribed river head. Uniform storativity
converts that maximum principle into a bound on aggregate storage change.
Neither a particular recession shape nor a finite response time is imposed.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL, PASS, CriterionIncompatibleError, CriterionResult, criterion, make_window,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

NAME = "aquifer_equilibrium_bounds"
VARIANTS = ("control", "raised", "lowered")
STATIC_KEYS = (
    "area_km2", "aquifer_storage_coefficient", "aquifer_specific_yield",
    "aquifer_initial_head_m", "aquifer_top_m", "aquifer_bottom_m",
    "river_conductance_m2_per_day", "river_bottom_offset_m",
)
OUTPUTS = ("gw", "gw_sw_exchange", "gw_to_sw", "sw_to_gw")


def _experiment(runs: dict[str, RunResult]) -> dict[str, float]:
    """Validate the judge-side case before judging a model's physics."""
    for variant in VARIANTS:
        if variant not in runs:
            raise ValueError(f"{NAME} needs variant {variant!r}")
    case0 = runs["control"].case
    static = {key: float(case0.static[key]) for key in STATIC_KEYS}
    if not np.isfinite(list(static.values())).all():
        raise ValueError("aquifer static inputs must be finite")
    positive = ("area_km2", "aquifer_storage_coefficient", "aquifer_specific_yield",
                "river_conductance_m2_per_day", "river_bottom_offset_m")
    if any(static[key] <= 0.0 for key in positive):
        raise ValueError("aquifer area, storage coefficients and conductance must be positive")
    if static["aquifer_bottom_m"] >= static["aquifer_top_m"]:
        raise ValueError("aquifer bottom must be below its top")
    if case0.spinup_steps < 1 or case0.spinup_steps >= len(case0.forcing):
        raise ValueError("the experiment needs common spinup and a scored interval")
    h0 = static["aquifer_initial_head_m"]
    stages = {}
    for variant in VARIANTS:
        run = runs[variant]
        case = run.case
        if case.static != case0.static:
            raise ValueError("paired aquifer static attributes must be identical")
        if (case.spinup_steps != case0.spinup_steps or case.timestep != case0.timestep
                or len(case.forcing) != len(case0.forcing)
                or len(run.table) != len(case.forcing)
                or not case.forcing["time"].equals(case0.forcing["time"])):
            raise ValueError("paired aquifer runs must have aligned times, lengths and spinup")
        recharge = case.forcing["gw_recharge"].to_numpy(dtype=float)
        stage = case.forcing["sw_stage_m"].to_numpy(dtype=float)
        if not (np.isfinite(recharge).all() and np.isfinite(stage).all()):
            raise ValueError("prescribed aquifer forcing must be finite")
        if np.any(recharge != 0.0):
            raise ValueError("the equilibrium experiment requires zero recharge")
        start = case.spinup_steps
        if not np.allclose(stage[:start], h0, rtol=0.0, atol=1.0e-9):
            raise ValueError("all variants must spin up at the initial head")
        if not np.all(stage[start:] == stage[start]):
            raise ValueError("the scored river head must be constant")
        stages[variant] = float(stage[start])
    if (abs(stages["control"] - h0) > 1.0e-9
            or stages["raised"] <= h0 or stages["lowered"] >= h0):
        raise ValueError("control, raised and lowered heads must bracket the initial head")
    minimum_head = min(h0, *stages.values())
    maximum_bed = max(h0, *stages.values()) - static["river_bottom_offset_m"]
    if minimum_head <= static["aquifer_top_m"]:
        raise ValueError("the experiment must stay fully confined")
    if minimum_head - maximum_bed < 1.0:
        raise ValueError("riverbed-to-head clearance must be at least 1 m")
    with np.errstate(over="ignore", invalid="ignore"):
        capacity = np.float64(1000.0) * static["aquifer_storage_coefficient"]
        bounds = {v: float(capacity * (np.float64(h) - h0)) for v, h in stages.items()}
    if not np.isfinite(list(bounds.values())).all():
        raise ValueError("prescribed hydraulic storage bounds overflow")
    return bounds


@criterion(NAME, paired=True)
def aquifer_equilibrium_bounds(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict,
) -> CriterionResult:
    """Score every interval against head-derived bounds, with neutral guards."""
    bounds = _experiment(runs)
    tolerances = {
        key: float(params.get(key, default)) for key, default in (
            ("rel_tol", 1.0e-4), ("abs_tol_mm", 1.0e-4),
            ("initial_abs_tol_mm", 1.0e-4), ("neutral_abs_tol_mm", 1.0e-4),
            ("boundary_abs_tol_mm", 1.0e-6),
        )
    }
    if (not np.isfinite(list(tolerances.values())).all()
            or tolerances["rel_tol"] < 0.0
            or any(value <= 0.0 for key, value in tolerances.items() if key != "rel_tol")):
        raise ValueError("equilibrium allowances must be finite and positive")

    def invalid(message: str) -> CriterionResult:
        return CriterionResult(NAME, FAIL, message, threshold=1.0)

    for variant in VARIANTS:
        run = runs[variant]
        for column in OUTPUTS:
            if column not in run.table:
                raise ValueError(f"{NAME}: {variant} needs result column {column!r}")
        boundary_declared = run.model is not None and "gw_boundary" in run.model.emitted
        if boundary_declared and "gw_boundary" not in run.table:
            raise ValueError(f"{variant}: declared gw_boundary is absent from the result")
        columns = (*OUTPUTS, "gw_boundary") if "gw_boundary" in run.table else OUTPUTS
        if not np.isfinite(run.table[list(columns)].to_numpy(dtype=float)).all():
            return invalid(f"{variant}: non-finite aquifer storage or exchange")
        if "gw_boundary" in run.table:
            with np.errstate(over="ignore", invalid="ignore"):
                gross = float(np.abs(run.table["gw_boundary"].to_numpy(dtype=float)).sum() * run.case.dt_days)
            if not np.isfinite(gross):
                return invalid(f"{variant}: gross aquifer boundary activity overflows")
            if gross > tolerances["boundary_abs_tol_mm"]:
                raise CriterionIncompatibleError(
                    f"{variant}: active additional aquifer boundary ({gross:g} mm gross); "
                    "the equilibrium envelope requires no other aquifer source or sink"
                )

    control = make_window(runs["control"], probe)
    start = runs["control"].case.spinup_steps
    full_control = runs["control"].table
    inventory0 = float(full_control["gw"].iloc[0])
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        drift = float(np.max(np.abs(full_control["gw"].to_numpy(dtype=float) - inventory0)))
        gross = float(
            np.abs(full_control[["gw_to_sw", "sw_to_gw"]].to_numpy(dtype=float)).sum()
            * control.dt_days
        )
        gross = max(gross, float(np.abs(full_control["gw_sw_exchange"].to_numpy(dtype=float)).sum() * control.dt_days))
        initial = max(float(np.max(np.abs(
            runs[v].table["gw"].iloc[:start].to_numpy(dtype=float)
            - full_control["gw"].iloc[:start].to_numpy(dtype=float)
        ))) for v in ("raised", "lowered"))
    if not np.isfinite((drift, gross, initial)).all():
        return invalid("non-finite neutral or initial storage difference (arithmetic overflow)")
    diagnostics = {
        "control": {"max_storage_drift_mm": drift, "gross_exchange_mm": gross},
        "max_initial_storage_difference_mm": initial, "variants": {},
    }
    ratios = {
        "initial storage": initial / tolerances["initial_abs_tol_mm"],
        "control storage": drift / tolerances["neutral_abs_tol_mm"],
        "control exchange": gross / tolerances["neutral_abs_tol_mm"],
    }
    responses = []
    for variant in ("raised", "lowered"):
        other = make_window(runs[variant], probe)
        lower, upper = min(0.0, bounds[variant]), max(0.0, bounds[variant])
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            allowed = max(tolerances["rel_tol"] * abs(bounds[variant]), tolerances["abs_tol_mm"])
            delta = other.table["gw"].to_numpy(dtype=float) - control.table["gw"].to_numpy(dtype=float)
            excess = np.maximum(np.maximum(lower - delta, delta - upper), 0.0)
            error = float(np.max(excess))
            ratio = error / allowed
        if not (np.isfinite(delta).all() and np.isfinite(excess).all()
                and np.isfinite((allowed, ratio)).all()):
            return invalid(f"{variant}: non-finite storage envelope arithmetic (overflow)")
        diagnostics["variants"][variant] = {
            "lower_bound_mm": lower, "upper_bound_mm": upper,
            "allowed_mm": allowed, "max_excess_mm": error,
            "min_delta_storage_mm": float(np.min(delta)),
            "max_delta_storage_mm": float(np.max(delta)),
            "worst_step": int(np.argmax(excess)),
        }
        ratios[variant] = ratio
        responses.append(float(np.max(np.abs(delta))))
    if not np.isfinite(list(ratios.values())).all():
        return invalid("non-finite normalized storage envelope arithmetic")
    worst = max(ratios, key=ratios.get)
    value = ratios[worst]
    unexercised = max(responses) <= tolerances["abs_tol_mm"]
    diagnostics["no_hydraulic_response_exercised"] = unexercised
    diagnostics["ratios"] = ratios
    passed = value <= 1.0
    message = (
        "confined storage stays inside both prescribed head envelopes"
        if passed else f"equilibrium bounds fail {worst}: residual / allowance = {value:.6g}"
    )
    if passed and unexercised:
        message += "; no hydraulic response exercised"
    return CriterionResult(NAME, PASS if passed else FAIL, message, value, 1.0, diagnostics)

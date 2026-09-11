"""Check instantaneous surface temperature against total upward longwave.

For a uniform, opaque, snow-free gray surface with a fixed emissivity, the
upward longwave radiation at an instant is fixed by the surface temperature
and the sky:

    rlus = eps * sigma * ts**4 + (1 - eps) * rlds

Emission and reflected sky must agree with the reported flux at every step;
interval means and cancellation between steps are not allowed. The default
tolerance is max(0.005 * abs(rlus), 0.5 W m-2). Emissivity comes only from
static.json. This checks output consistency, not temperature accuracy.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

# Stefan-Boltzmann constant, W m-2 K-4 (CODATA 2018).
STEFAN_BOLTZMANN = 5.670374419e-8


@criterion("radiative_identity")
def radiative_identity(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Upward longwave must equal what the reported surface temperature emits
    plus the reflected share of the downward longwave, at every step."""
    upward = str(params.get("upward", "rlus"))
    temperature = str(params.get("temperature", "ts"))
    downward = str(params.get("downward", "rlds"))
    emissivity_key = str(params.get("emissivity", "eps"))
    sigma = float(params.get("sigma", STEFAN_BOLTZMANN))
    rel_tol = float(params.get("rel_tol", 0.005))
    abs_floor = float(params.get("abs_floor", 0.5))
    # NaN comparisons and infinite allowances can silently accept bad output.
    if not (np.isfinite(sigma) and sigma > 0.0):
        raise ValueError(f"radiative_identity needs a finite sigma > 0, not {sigma!r}")
    if not (np.isfinite(rel_tol) and rel_tol >= 0.0):
        raise ValueError(f"radiative_identity needs a finite rel_tol >= 0, not {rel_tol!r}")
    if not (np.isfinite(abs_floor) and abs_floor > 0.0):
        raise ValueError(f"radiative_identity needs a finite abs_floor > 0, not {abs_floor!r}")

    w = make_window(run, probe)
    for var in (upward, temperature):
        if var not in w.table.columns:
            raise ValueError(f"radiative_identity needs '{var}' in the model result")
    if downward not in w.forcing.columns:
        raise ValueError(
            f"radiative_identity needs forcing column '{downward}'; this probe's "
            "generator does not produce it"
        )

    # Output-derived emissivity would let a model force its own consistency.
    static = run.case.static
    if emissivity_key not in static:
        raise ValueError(
            f"radiative_identity needs '{emissivity_key}' in the case's static "
            "attributes; this probe's generator does not supply it"
        )
    eps = float(static[emissivity_key])
    if not 0.0 < eps <= 1.0:
        raise ValueError(f"emissivity '{emissivity_key}' is {eps!r}, outside (0, 1]")

    rlds = w.forcing[downward].to_numpy(dtype=float)
    if not np.isfinite(rlds).all():
        raise ValueError(f"forcing column '{downward}' has non-finite values")

    rlus = w.table[upward].to_numpy(dtype=float)
    ts = w.table[temperature].to_numpy(dtype=float)
    finite = np.isfinite(rlus) & np.isfinite(ts)
    if not finite.all():
        n_bad = int((~finite).sum())
        return CriterionResult(
            name="radiative_identity",
            status=FAIL,
            message=f"non-finite surface temperature or upward longwave on {n_bad} scored steps",
            diagnostics={"non_finite_steps": n_bad},
        )
    # The fourth power hides the sign: -290 K would otherwise emit as +290 K.
    below_zero = ts <= 0.0
    if below_zero.any():
        n_bad = int(below_zero.sum())
        return CriterionResult(
            name="radiative_identity",
            status=FAIL,
            message=(
                f"'{temperature}' is at or below 0 K on {n_bad} scored steps; the surface "
                "temperature must be reported in kelvin"
            ),
            diagnostics={"non_positive_kelvin_steps": n_bad},
        )

    # Finite inputs can still overflow; never turn inf/inf into a passing NaN.
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        expected = eps * sigma * ts**4 + (1.0 - eps) * rlds
        residual = rlus - expected
        allowance = np.maximum(rel_tol * np.abs(rlus), abs_floor)
        slack = np.abs(residual) / allowance
    if not np.isfinite(allowance).all():
        raise ValueError("radiative_identity needs a finite allowance; rel_tol overflows")
    if not np.isfinite(slack).all():
        return CriterionResult(
            name="radiative_identity",
            status=FAIL,
            message="surface radiation calculation overflowed on scored steps",
            diagnostics={"non_finite_calculation_steps": int((~np.isfinite(slack)).sum())},
        )
    violating = slack > 1.0

    n_bad = int(violating.sum())
    worst = int(slack.argmax())
    times = w.forcing["time"].astype(str).to_numpy()

    detail = (
        f"worst step {slack[worst]:.2f} of tolerance at {times[worst]} "
        f"(residual {residual[worst]:+.3g} W m-2 against {allowance[worst]:.3g} allowed; "
        f"eps {eps:.4g})"
    )
    if n_bad == 0:
        message = (
            f"upward longwave agrees with the reported surface temperature on all "
            f"{len(rlus)} steps: {detail}"
        )
    else:
        first = int(violating.argmax())
        message = (
            f"upward longwave contradicts the reported surface temperature on {n_bad} of "
            f"{len(rlus)} steps, first at index {first} ({times[first]}): "
            f"{rlus[first]:.3f} W m-2 reported, {expected[first]:.3f} expected from "
            f"{temperature} {ts[first]:.2f} K and {downward} {rlds[first]:.1f} W m-2; {detail}"
        )

    return CriterionResult(
        name="radiative_identity",
        status=PASS if n_bad == 0 else FAIL,
        value=float(slack[worst]),
        threshold=1.0,
        message=message,
        diagnostics={
            "violating_steps": n_bad,
            "scored_steps": len(rlus),
            "worst_slack": float(slack[worst]),
            "worst_step": {
                "index": worst,
                "time": times[worst],
                "upward_w_m2": float(rlus[worst]),
                "expected_w_m2": float(expected[worst]),
                "residual_w_m2": float(residual[worst]),
                "allowance_w_m2": float(allowance[worst]),
                "temperature_k": float(ts[worst]),
                "downward_w_m2": float(rlds[worst]),
            },
            "max_abs_residual_w_m2": float(np.abs(residual).max()),
            "mean_abs_residual_w_m2": float(np.abs(residual).mean()),
            # Identify weak fluxes where the absolute floor sets the bound.
            "floor_steps": int((rel_tol * np.abs(rlus) < abs_floor).sum()),
            "emissivity": eps,
            "sigma_w_m2_k4": sigma,
            "rel_tol": rel_tol,
            "abs_floor_w_m2": abs_floor,
        },
    )

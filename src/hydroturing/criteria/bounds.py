"""Physical bounds on states and fluxes.

This is the criterion that catches closure by construction, which is the
attack randomised forcing cannot touch. A model that solves for storage as
whatever balances the budget will pass any closure check, forever, on every
seed. What it cannot do is keep that storage physical: the invented storage
drifts, and over a long random sequence it leaves the range a real soil
column or snowpack can occupy.

That is why the model contract requires absolute storage states rather than
tendencies. Reporting dS/dt would hide exactly the quantity this checks.
"""

from __future__ import annotations

import math

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

EPS = 1e-6  # mm, absorbs floating point noise without hiding a real excursion


def _resolve(bound, static: dict) -> float:
    """Bounds may be a number, null (unbounded), or a key into static.json."""
    if bound is None:
        return math.inf
    if isinstance(bound, (int, float)):
        return float(bound)
    if isinstance(bound, str):
        if bound not in static:
            raise ValueError(f"state bound '{bound}' is not in the case attributes")
        return float(static[bound])
    raise ValueError(f"cannot interpret state bound {bound!r}")


@criterion("state_bounds")
def state_bounds(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Every reported storage must stay inside its physical range."""
    w = make_window(run, probe)
    static = run.case.static
    violations = []
    worst_excess, worst_var = 0.0, None

    for var, spec in params.items():
        if var not in w.table.columns:
            continue
        lo = _resolve(spec[0], static) if isinstance(spec, (list, tuple)) else 0.0
        hi = _resolve(spec[1], static) if isinstance(spec, (list, tuple)) else math.inf
        values = np.asarray(w.table[var], dtype=float)

        below = float(np.maximum(lo - values.min(), 0.0))
        above = float(np.maximum(values.max() - hi, 0.0)) if math.isfinite(hi) else 0.0
        excess = max(below, above)

        if excess > EPS:
            n_bad = int(((values < lo - EPS) | (values > hi + EPS)).sum())
            violations.append(
                f"{var} leaves [{lo:g}, {hi:g}] on {n_bad} steps "
                f"(range {values.min():.3g} to {values.max():.3g})"
            )
            if excess > worst_excess:
                worst_excess, worst_var = excess, var

    ok = not violations
    return CriterionResult(
        name="state_bounds",
        status=PASS if ok else FAIL,
        value=worst_excess,
        threshold=EPS,
        message="all storages stay physical" if ok else "; ".join(violations),
        diagnostics={"worst_variable": worst_var, "n_violations": len(violations)},
    )


@criterion("et_plausible")
def et_plausible(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Evapotranspiration must be non-negative and bounded by potential ET.

    Cheap, and it catches the degenerate escape of setting ET equal to
    precipitation so that the water balance closes with no runoff at all.
    """
    max_ratio = float(params.get("max_ratio_to_pet", 1.0))
    pet_var = params.get("pet_variable", "pet")
    w = make_window(run, probe)

    if "evspsbl" not in w.table.columns:
        raise ValueError("et_plausible needs 'evspsbl' in the model result")
    et = np.asarray(w.table["evspsbl"], dtype=float)

    if et.min() < -EPS:
        return CriterionResult(
            name="et_plausible", status=FAIL, value=float(et.min()), threshold=0.0,
            message=f"evapotranspiration goes negative (min {et.min():.4g} mm/day)",
        )

    if pet_var not in w.forcing.columns:
        return CriterionResult(
            name="et_plausible", status=PASS,
            message="ET non-negative; no potential ET in the forcing to bound it against",
        )

    pet = np.asarray(w.forcing[pet_var], dtype=float)
    cum_ratio = float(et.sum() / max(pet.sum(), 1e-12))
    ok = cum_ratio <= max_ratio + 1e-9

    return CriterionResult(
        name="et_plausible",
        status=PASS if ok else FAIL,
        value=cum_ratio,
        threshold=max_ratio,
        message=(
            f"cumulative ET is {cum_ratio:.3f} of potential ET "
            f"(limit {max_ratio:g})"
        ),
        diagnostics={
            "n_steps_over_pet": int((et > pet + EPS).sum()),
            "cumulative_et_mm": float(et.sum()),
            "cumulative_pet_mm": float(pet.sum()),
        },
    )

"""Cross-budget coherence: two budgets that close, describing one evaporation.

A model with a water head and an energy head can close both budgets and still
report a latent heat flux implying a different evaporation than the one it
reported as water. No single-budget criterion can see that, because within
either budget the numbers are consistent. This module holds the criteria that
compare across the two.

The identity is exact, not a parameterisation:

    LE = lambda * E

with LE in W m-2, E in kg m-2 s-1, and lambda the latent heat of the phase
change the water actually undergoes. Two things make it discriminating rather
than arithmetic.

First, lambda depends on temperature. `lambda_v(T) = 2.501e6 - 2361 T` J kg-1
(Brutsaert 1982), so a model using a constant is wrong by at most 3.4 percent
over a -20 to +35 degC record. That is invisible to the suite's 5 percent rule,
which is why this criterion is scored per step against a much tighter bound.

Second, lambda jumps at the phase change. Water leaving a snowpack sublimates
and takes `lambda_s = lambda_v(0) + lambda_f`, 13.3 percent more energy per
kilogram. A model that converts every kilogram at the vaporisation rate is
short of that on exactly the steps where snow is disappearing.

The sublimated mass is never inferred. A model that reports `sbl`, the
sublimating share of its evaporation, is held to the equality at every step,
because it has said which kilograms left as ice. A model that does not report
it is held only to the interval the two latent heats span, wherever a pack is
present or could arrive during the step.

An earlier version inferred the split instead, reading any loss from the snow
store on a dry sub-freezing day as sublimation. That is wrong, and review
caught it: a pack also loses water at its base, which is what Snow-17's DAYGM
term does. An honest model running a constant ground melt of 0.3 mm/day failed
on about 400 steps of 3650 while its latent heat was exactly right, and at
0.1 mm/day it passed with only 20 percent of the tolerance to spare. Asking
the model rather than guessing costs one optional variable and removes the
whole class of error.

Related work. Enforcing conservation inside a neural emulator, and measuring
the "physical inconsistency" left over when it is not enforced, is established
for atmospheric convection (Beucler, Pritchard, Rasp, Ott, Baldi and Gentine,
2021, Phys. Rev. Lett. 126, 098302). What is different here is the form and
the domain: a binary gate rather than a penalty term, over the surface water
and energy budgets, evaluated as a residual the model cannot fit away.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

# Latent heat of vaporisation as A + B*T, T in degC, J kg-1.
LAMBDA_V = (2.501e6, -2361.0)
# Latent heat of fusion, J kg-1. lambda_s = A + LAMBDA_F.
LAMBDA_F = 3.337e5
SECONDS_PER_DAY = 86400.0


@criterion("flux_identity")
def flux_identity(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Latent heat must equal the reported evaporation times its latent heat.

    Scored at every step. The tolerance is relative to the reported flux with
    an absolute floor, because the flux passes through zero and a percentage
    of nothing is not a bound.
    """
    flux = str(params.get("flux", "hfls"))
    water = str(params.get("water", "evspsbl"))
    temperature = str(params.get("temperature", "tas"))
    rel_tol = float(params.get("rel_tol", 0.005))
    abs_floor = float(params.get("abs_floor", 0.5))
    lam_a, lam_b = (float(v) for v in params.get("lambda_vapour", LAMBDA_V))
    lam_f = float(params.get("lambda_fusion", LAMBDA_F))
    phase_state = params.get("phase_state", "snw")
    snowfall_var = str(params.get("snowfall_from", "pr"))
    snow_threshold = float(params.get("snow_threshold_degC", 0.0))
    sublimation_var = params.get("sublimation", "sbl")

    w = make_window(run, probe)
    for var in (flux, water):
        if var not in w.table.columns:
            raise ValueError(f"flux_identity needs '{var}' in the model result")
    if temperature not in w.forcing.columns:
        raise ValueError(
            f"flux_identity needs forcing column '{temperature}'; this probe's "
            "generator does not produce it"
        )

    dt = w.dt_days
    tas = w.forcing[temperature].to_numpy(dtype=float)
    le = w.table[flux].to_numpy(dtype=float)
    # mm day-1 is kg m-2 day-1, so a per-step depth is a per-step mass.
    et_mass = w.table[water].to_numpy(dtype=float) * dt

    lam_v = lam_a + lam_b * tas
    lam_s = lam_a + lam_f
    seconds = dt * SECONDS_PER_DAY

    # Two regimes, and which one a step is in depends on what the model chose
    # to report rather than on anything this criterion assumes about snow.
    #
    #   reports `sbl`      it has said which kilograms left as ice, so the
    #                      equality holds at every step, pack or no pack:
    #                          LE = lambda_v (E - sbl) + lambda_s sbl
    #
    #   reports no `sbl`   snow-free steps where none could fall are still an
    #                      equality at lambda_v; anywhere a pack is or could
    #                      be present, only the interval:
    #                          lambda_v(T) E <= LE <= lambda_s E
    #
    # The interval is the honest bound on a model that has not told us, and it
    # is weak on purpose. A pack loses water to more than sublimation --
    # Snow-17's DAYGM puts it into the soil -- so a criterion that read the
    # pack's own mass loss as sublimation would fail an honest model.
    subl = np.zeros_like(et_mass)
    reported_split = bool(sublimation_var) and sublimation_var in w.table.columns
    if reported_split:
        column = pd.to_numeric(w.table[sublimation_var], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(column).all():
            raise ValueError(f"column {sublimation_var!r} has non-finite values")
        subl = column * dt

    pack_possible = np.ones_like(et_mass, dtype=bool)
    if phase_state and phase_state in w.table.columns and snowfall_var in w.forcing.columns:
        snw = w.table[phase_state].to_numpy(dtype=float)
        prior = float(w.state0[phase_state]) if phase_state in w.state0.index else snw[0]
        previous = np.concatenate(([prior], snw[:-1]))
        # A step is snow-free only if none could arrive during it: a pack that
        # falls and disappears inside one step begins and ends at zero.
        snowfall = (w.forcing[snowfall_var].to_numpy(dtype=float) > 0.0) & (tas < snow_threshold)
        pack_possible = (previous > 0.0) | (snw > 0.0) | snowfall

    failures: list[str] = []
    if reported_split:
        # The split is a claim about the model's own evaporation, so it has to
        # be one: never negative, never more than what evaporated, and zero
        # where the model itself reports no ice to lose. Without the last of
        # these, a model could report a fictitious sublimating share to bend
        # its effective lambda upwards on a warm day.
        checks = (
            (subl < -1e-9, "negative"),
            (subl > np.maximum(et_mass, 0.0) + 1e-9,
             "larger than the evaporation it is a share of"),
            ((~pack_possible) & (subl > 1e-9),
             "non-zero where the model reports no snow and none could fall"),
        )
        for mask, what in checks:
            if mask.any():
                failures.append(f"{sublimation_var!r} is {what} on {int(mask.sum())} steps")
        equality = (lam_v * (et_mass - subl) + lam_s * subl) / seconds
        lower = upper = equality
    else:
        liquid_only = lam_v * et_mass / seconds
        spanned = np.stack([liquid_only, lam_s * et_mass / seconds])
        lower = np.where(pack_possible, spanned.min(axis=0), liquid_only)
        upper = np.where(pack_possible, spanned.max(axis=0), liquid_only)

    tolerance = np.maximum(rel_tol * np.abs(le), abs_floor)
    residual = np.where(le < lower, le - lower, np.where(le > upper, le - upper, 0.0))
    slack = np.abs(residual) / tolerance
    violating = slack > 1.0
    required = np.clip(le, lower, upper)

    # The implied bulk latent heat is the single most readable diagnostic: a
    # model using a constant lambda reports the constant back, and one that
    # ignores sublimation reports a number below lambda_v.
    total_mass = float(et_mass.sum())
    implied = float((le * dt * SECONDS_PER_DAY).sum() / total_mass) if total_mass > 0 else float("nan")

    bounded = pack_possible & (not reported_split)
    worst = float(slack.max()) if slack.size else 0.0
    n_bad = int(violating.sum())
    ok = n_bad == 0 and not failures

    def worst_in(mask):
        return float(slack[mask].max()) if mask.any() else 0.0

    split_steps = int((subl > 0).sum())
    detail = f"implied lambda {implied:.4g} J/kg; worst step {worst:.2f} of tolerance"
    if reported_split:
        detail += (
            f" ({split_steps} steps report sublimation, worst there "
            f"{worst_in(subl > 0):.2f})"
        )
    else:
        detail += (
            f" (equality on {int((~bounded).sum())} snow-free steps, interval on "
            f"{int(bounded.sum())} where a pack is or could be present)"
        )

    if failures:
        message = "; ".join(failures) + f" [{detail}]"
    elif ok:
        message = f"latent heat matches the reported evaporation at every step: {detail}"
    else:
        first = int(violating.argmax())
        message = (
            f"latent heat contradicts the reported evaporation on {n_bad} of "
            f"{len(le)} steps, first at index {first} "
            f"({le[first]:.3f} W m-2 reported, {required[first]:.3f} required); {detail}"
        )

    return CriterionResult(
        name="flux_identity",
        status=PASS if ok else FAIL,
        value=worst,
        threshold=1.0,
        message=message,
        diagnostics={
            "implied_lambda_j_per_kg": implied,
            "worst_slack": worst,
            "worst_slack_sublimating": worst_in(subl > 0.0),
            "worst_slack_bounded": worst_in(bounded),
            "bounded_steps": int(bounded.sum()),
            "reported_split": bool(reported_split),
            "violating_steps": n_bad,
            "sublimating_steps": split_steps,
            "rel_tol": rel_tol,
            "abs_floor_w_m2": abs_floor,
        },
    )


@criterion("energy_closure")
def energy_closure(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Net radiation must be accounted for by the fluxes leaving the surface.

    Separate from `closure` rather than a denominator of it, because `closure`
    differences the probe's water states and this budget has none: `hfg` is
    the storage term, and a model that reports it has already said where the
    energy went. Keeping them apart also lets one probe score both budgets
    without two criteria of the same name.

    Net radiation changes sign every night, so the relative test carries an
    absolute floor as the module docstring of `closure` prescribes.
    """
    driver = str(params.get("driver", "rn"))
    sinks = list(params.get("sinks", ["hfls", "hfss", "hfg"]))
    threshold = float(params.get("threshold", 0.05))
    floor = float(params.get("floor", 2.0))

    w = make_window(run, probe)
    if driver not in w.forcing.columns:
        raise ValueError(
            f"energy_closure needs forcing column '{driver}'; this probe's "
            "generator does not produce it"
        )

    # The driver is taken from the forcing, never from what the model echoed
    # back, for the same reason `closure` does it: a model must not be able to
    # move its own denominator.
    drive = w.volume(w.forcing[driver].to_numpy(dtype=float))
    outflow = np.zeros(len(w.table))
    for var in sinks:
        if var not in w.table.columns:
            raise ValueError(f"energy_closure needs '{var}' in the model result")
        outflow += w.volume(w.table[var].to_numpy(dtype=float))

    step_residual = drive - outflow
    cumulative = float(step_residual.sum())
    total = float(np.abs(drive).sum())
    if total <= 0:
        return CriterionResult(
            name="energy_closure",
            status=FAIL,
            message=f"accumulated |{driver}| is zero; the case is degenerate",
        )

    relative = abs(cumulative) / total
    mean_abs = float(np.abs(step_residual).mean() / w.dt_days)
    ok = relative <= threshold or mean_abs <= floor

    return CriterionResult(
        name="energy_closure",
        status=PASS if ok else FAIL,
        value=relative,
        threshold=threshold,
        message=(
            f"cumulative residual {relative:.4%} of accumulated |{driver}| "
            f"(limit {threshold:.1%}, floor {floor:g} W m-2; "
            f"mean step residual {mean_abs:.3g} W m-2)"
        ),
        diagnostics={
            "cumulative_residual": cumulative,
            "denominator_total": total,
            "mean_step_residual_w_m2": mean_abs,
            "max_step_residual_w_m2": float(np.abs(step_residual).max() / w.dt_days),
            "floor_w_m2": floor,
        },
    )

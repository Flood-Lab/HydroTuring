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

The sublimated mass is not assumed. On a step with no precipitation and air
below freezing, no snow can fall and none can melt, so any decrease in the
snow water equivalent the model itself reported is sublimation. Where the
split cannot be inferred that way, because a pack is present and melting is
possible, the criterion asserts only the interval that the two latent heats
span rather than an equality. A criterion that failed an honest model whose
snow physics differs from the reference's would be worse than no criterion,
and a model that declines to report a snow state simply has the split
disabled rather than being judged against an assumption about it.

Related work. Enforcing conservation inside a neural emulator, and measuring
the "physical inconsistency" left over when it is not enforced, is established
for atmospheric convection (Beucler, Pritchard, Rasp, Ott, Baldi and Gentine,
2021, Phys. Rev. Lett. 126, 098302). What is different here is the form and
the domain: a binary gate rather than a penalty term, over the surface water
and energy budgets, evaluated as a residual the model cannot fit away.
"""

from __future__ import annotations

import numpy as np

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

    # Every step falls into one of three cases, and only two of them are an
    # equality. Which case a step is in is decided from the forcing and the
    # model's own snow state, never from an assumption about how the model
    # partitions evaporation.
    #
    #   no pack           all of it is liquid          LE = lambda_v(T) E
    #   pack, dry, frozen no snowfall and no melt is possible, so the pack's
    #                     own mass loss is the sublimated mass
    #                                                  LE = lambda_v L + lambda_s S
    #   pack, otherwise   melt and sublimation are not separable from what the
    #                     model reports, so the criterion asserts only the
    #                     interval the two latent heats span
    #                                     lambda_v(T) E <= LE <= lambda_s E
    #
    # The interval is weaker on purpose. Asserting a split the model never
    # reported would fail an honest model whose snow physics differs from the
    # reference's, and a criterion that does that is worse than no criterion.
    subl = np.zeros_like(et_mass)
    inferable = np.ones_like(et_mass, dtype=bool)
    split_steps = 0
    if phase_state and phase_state in w.table.columns and snowfall_var in w.forcing.columns:
        snw = w.table[phase_state].to_numpy(dtype=float)
        prior = float(w.state0[phase_state]) if phase_state in w.state0.index else snw[0]
        previous = np.concatenate(([prior], snw[:-1]))
        frozen_and_dry = (w.forcing[snowfall_var].to_numpy(dtype=float) <= 0.0) & (
            tas < snow_threshold
        )
        # A step is snow-free only if none could have arrived during it. A
        # pack that falls and sublimates away inside one step begins and ends
        # at zero, and asserting the liquid equality there would fail a model
        # that handled it correctly.
        snowfall = (w.forcing[snowfall_var].to_numpy(dtype=float) > 0.0) & (
            tas < snow_threshold
        )
        no_pack = (previous <= 0.0) & (snw <= 0.0) & ~snowfall
        pack_loss = np.clip(previous - snw, 0.0, None)
        subl = np.where(frozen_and_dry, np.minimum(pack_loss, np.clip(et_mass, 0.0, None)), 0.0)
        inferable = no_pack | frozen_and_dry
        split_steps = int((subl > 0).sum())

    equality = (lam_v * (et_mass - subl) + lam_s * subl) / seconds
    spanned = np.stack([lam_v * et_mass / seconds, lam_s * et_mass / seconds])
    lower = np.where(inferable, equality, spanned.min(axis=0))
    upper = np.where(inferable, equality, spanned.max(axis=0))

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

    liquid = inferable & (subl <= 0.0)
    worst = float(slack.max()) if slack.size else 0.0
    n_bad = int(violating.sum())
    ok = n_bad == 0

    def worst_in(mask):
        return float(slack[mask].max()) if mask.any() else 0.0

    detail = (
        f"implied lambda {implied:.4g} J/kg; worst step {worst:.2f} of tolerance "
        f"(liquid {worst_in(liquid):.2f}, sublimating {worst_in(subl > 0.0):.2f} "
        f"over {split_steps} steps, bounded {worst_in(~inferable):.2f} over "
        f"{int((~inferable).sum())})"
    )
    if ok:
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
            "worst_slack_liquid": worst_in(liquid),
            "worst_slack_sublimating": worst_in(subl > 0.0),
            "worst_slack_bounded": worst_in(~inferable),
            "bounded_steps": int((~inferable).sum()),
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

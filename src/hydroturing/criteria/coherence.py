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
and energy budgets, and in `partition_shift` a counterfactual rather than a
same-instant residual.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from hydroturing.criteria.base import (
    FAIL, PASS, CriterionResult, criterion, make_window, segments,
)
from hydroturing.criteria.response import pick
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
    the downward flux at the actual soil surface, already corrected by the
    adapter for storage above a deeper flux boundary if needed. Subsurface
    storage is not subtracted again. Keeping the budgets apart also lets one
    probe score both without two criteria of the same name.

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


@criterion("energy_closure_by_phase")
def energy_closure_by_phase(
    run: RunResult, probe: ProbeSpec, params: dict
) -> CriterionResult:
    """The mean absolute skin-budget residual must be small in every block.

    The generator labels contiguous day/night blocks independently of the
    sign of net radiation. Taking the absolute residual before integration
    prevents opposite errors from cancelling both within and across blocks.
    H and LE may have either sign; G is at the actual soil surface, so soil
    heat storage is not subtracted again from this zero-capacity skin budget.
    """
    driver = str(params.get("driver", "rn"))
    sinks = list(params.get("sinks", ["hfls", "hfss", "hfg"]))
    threshold = float(params.get("threshold", 0.05))
    floor = float(params.get("floor", 2.0))
    segment_column = str(params.get("segment_column", "_regime"))

    w = make_window(run, probe)
    if driver not in w.forcing.columns:
        raise ValueError(
            f"energy_closure_by_phase needs forcing column '{driver}'; "
            "this probe's generator does not produce it"
        )
    for var in sinks:
        if var not in w.table.columns:
            raise ValueError(f"energy_closure_by_phase needs '{var}' in the model result")
    phases = segments(w, segment_column)
    drive = w.forcing[driver].to_numpy(dtype=float)
    fluxes = w.table[sinks].to_numpy(dtype=float)
    finite = np.isfinite(drive) & np.isfinite(fluxes).all(axis=1)
    if not finite.all():
        n_bad = int((~finite).sum())
        return CriterionResult(
            name="energy_closure_by_phase",
            status=FAIL,
            message=f"non-finite surface energy values on {n_bad} scored steps",
            diagnostics={"non_finite_steps": n_bad},
        )

    # A Case has one fixed dt, so it cancels from sum(abs(r) * dt) / sum(dt).
    # Taking the mean directly also avoids rounding an exact tolerance-boundary
    # value upward through an unnecessary multiplication and division by dt.
    residual = drive - fluxes.sum(axis=1)
    blocks = []
    for label, start, stop in phases:
        duration_days = (stop - start) * w.dt_days
        mean_abs = float(np.abs(residual[start:stop]).mean())
        mean_driver = float(np.abs(drive[start:stop]).mean())
        allowance = max(threshold * mean_driver, floor)
        slack = (
            mean_abs / allowance if allowance > 0
            else (0.0 if mean_abs == 0 else float("inf"))
        )
        blocks.append({
            "label": label,
            "start": start,
            "stop": stop,
            "duration_hours": duration_days * 24.0,
            "mean_abs_w_m2": mean_abs,
            "mean_abs_driver_w_m2": mean_driver,
            "allowance_w_m2": allowance,
            "slack": slack,
            "passed": mean_abs <= allowance,
        })

    failed = sum(not block["passed"] for block in blocks)
    worst = max(blocks, key=lambda block: block["slack"])
    return CriterionResult(
        name="energy_closure_by_phase",
        status=PASS if failed == 0 else FAIL,
        value=worst["slack"],
        threshold=1.0,
        message=(
            f"{failed} of {len(blocks)} phase blocks fail; worst {worst['label']} "
            f"[{worst['start']}:{worst['stop']}] has mean absolute residual "
            f"{worst['mean_abs_w_m2']:.3g} W m-2 against "
            f"{worst['allowance_w_m2']:.3g} W m-2 allowed; "
            "assumes ground heat is mapped to the actual soil surface "
            "(an uncorrected deeper-boundary flux can also cause a residual)"
        ),
        diagnostics={
            "blocks": blocks,
            "failed_blocks": failed,
            "worst_block": worst,
            "relative_threshold": threshold,
            "floor_w_m2": floor,
        },
    )


@criterion("partition_shift", paired=True)
def partition_shift(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    """Energy a drying surface stops evaporating must warm the air.

    The two variants carry byte-identical net radiation and differ only in
    precipitation, so over the perturbed window

        sum (dLE + dH + dG) = sum dRn = 0

    is an identity rather than an expectation. There is no parameterisation in
    it and no tolerance worth arguing about: the latent heat a surface gives
    up when its soil dries has nowhere to go but the sensible and ground
    fluxes.

    Three things are asserted, and the second and third are what make it a
    test rather than an arithmetic check on a budget already scored.

    Responded    the latent flux must actually fall. A model whose energy side
                 never reads its water side does not notice the drought at
                 all, and fails here before the sum is ever examined.

    Conserved    the three changes must sum to zero.

    In software-testing terms this is a metamorphic relation rather than a
    budget check: it asserts how the output must change when the input is
    changed in a known way, which is what makes it reach a property that no
    amount of accuracy on either run can establish. Enforcing conservation
    inside a network is well studied (Beucler et al. 2021, Phys. Rev. Lett.
    126, 098302, which quantifies "physical inconsistency" across mass,
    momentum, radiation and energy budgets for neural emulators). Asserting
    the cross-budget response to a counterfactual is a different question, and
    the one a model cannot fit its way past.

    Went to air  most of it must appear in the sensible flux. A model can
                 satisfy the sum by letting the ground flux absorb whatever
                 the latent flux gives up, which over a season would require a
                 soil heat store that does not exist. Bounding the ground
                 share is what closes that dodge, and it is the half that
                 corresponds to the land-surface amplification of heat
                 extremes: a model that fails it cannot represent a hot
                 drought however well it fits a hydrograph.
    """
    latent = str(params.get("latent", "hfls"))
    sensible = str(params.get("sensible", "hfss"))
    ground = str(params.get("ground", "hfg"))
    driver = str(params.get("driver", "rn"))
    label = str(params.get("window_label", "_perturbed"))
    report_label = params.get("intensification_label", "_intensification")
    min_shift = float(params.get("min_shift", 0.05))
    ground_share = float(params.get("max_ground_share", 0.20))
    tolerance = float(params.get("tolerance", 0.05))

    control = make_window(pick(runs, params, "control", "control"), probe)
    perturbed = make_window(pick(runs, params, "perturbed", "perturbed"), probe)

    # Everything below indexes both variants with one mask, so the two have to
    # be the same shape and marked the same way. A generator that breaks either
    # invariant should be told so here rather than crashing inside numpy or,
    # worse, silently comparing two different stretches of record.
    if len(control.table) != len(perturbed.table):
        raise ValueError(
            f"the variants have different scored lengths ({len(control.table)} and "
            f"{len(perturbed.table)}); partition_shift compares them step for step"
        )
    for name, w in (("control", control), ("perturbed", perturbed)):
        if label not in w.forcing.columns:
            raise ValueError(
                f"partition_shift needs the generator to mark the perturbed window "
                f"with a '{label}' forcing column; the {name} variant has none"
            )
    mask = control.forcing[label].to_numpy(dtype=float) > 0
    if not np.array_equal(mask, perturbed.forcing[label].to_numpy(dtype=float) > 0):
        raise ValueError(
            f"the variants mark different steps as '{label}'; the perturbed window "
            "has to be the same stretch of record in both"
        )
    if not mask.any():
        raise ValueError(f"the '{label}' column marks no steps; the case is degenerate")

    # The identity only holds if the generator really did hold radiation
    # fixed. Checking it here means a generator that drifts is diagnosed as a
    # generator fault rather than reported as a model failure.
    rn_c = control.volume(control.forcing[driver].to_numpy(dtype=float))[mask]
    rn_p = perturbed.volume(perturbed.forcing[driver].to_numpy(dtype=float))[mask]
    available = float(rn_c.sum())
    drift = float(np.abs(rn_p - rn_c).sum())
    if available <= 0:
        raise ValueError(f"accumulated {driver} over the window is not positive")
    if drift > 1e-9 * max(available, 1.0):
        raise ValueError(
            f"the variants differ in {driver} by {drift:.4g} over the window; "
            "this criterion is an identity only while the driver is held fixed"
        )

    def total(w, var: str, over) -> float:
        if var not in w.table.columns:
            raise ValueError(f"partition_shift needs '{var}' in the model result")
        return float(w.volume(w.table[var].to_numpy(dtype=float))[over].sum())

    def change(var: str, over=None) -> float:
        over = mask if over is None else over
        return total(perturbed, var, over) - total(control, var, over)

    def evaporative_fraction(w, over) -> float:
        """LE / (LE + H) over a window, the form this literature reads in."""
        le, h = total(w, latent, over), total(w, sensible, over)
        turbulent = le + h
        return le / turbulent if abs(turbulent) > 1e-12 else float("nan")

    d_le, d_h, d_g = change(latent), change(sensible), change(ground)
    imbalance = d_le + d_h + d_g
    d_ef = evaporative_fraction(perturbed, mask) - evaporative_fraction(control, mask)

    # The stretch where the drying actually bites, reported and never gated.
    # The scoring window has to be long enough for the shift to be
    # unambiguous, which makes it longer than the timescale on which this
    # repartitioning is usually studied. Reporting the intensification stretch
    # separately turns a four-month total into a rate, without hanging a
    # threshold on a stretch whose position comes from the reference model.
    onset: dict[str, float] = {}
    if report_label and report_label in control.forcing.columns:
        onset_mask = control.forcing[report_label].to_numpy(dtype=float) > 0
        if onset_mask.any():
            onset = {
                "steps": int(onset_mask.sum()),
                "d_latent": change(latent, onset_mask),
                "d_sensible": change(sensible, onset_mask),
                "d_ground": change(ground, onset_mask),
                "d_evaporative_fraction": (
                    evaporative_fraction(perturbed, onset_mask)
                    - evaporative_fraction(control, onset_mask)
                ),
            }

    failures = []
    if d_le > -min_shift * available:
        failures.append(
            f"{latent} barely moved: {d_le:+.4g} against a required fall of "
            f"{min_shift * available:.4g} (W m-2 day, {min_shift:.0%} of the "
            f"{available:.4g} the window received)"
        )
    scale = max(abs(d_le), min_shift * available)
    if abs(imbalance) > tolerance * scale:
        failures.append(
            f"the three fluxes changed by {imbalance:+.4g} in total, not zero, "
            f"under net radiation that did not change (limit {tolerance:.0%} of "
            f"{scale:.4g})"
        )
    if abs(d_g) > ground_share * scale:
        failures.append(
            f"{ground} absorbed {abs(d_g) / scale:.2f} of the shift (limit "
            f"{ground_share:.2f}); the energy went into the soil instead of the air"
        )
    elif d_h <= 0:
        failures.append(f"{sensible} did not rise ({d_h:+.4g})")

    ok = not failures
    detail = (
        f"d{latent} {d_le:+.4g}, d{sensible} {d_h:+.4g}, d{ground} {d_g:+.4g} "
        f"(W m-2 day), sum {imbalance:+.3g}, dEF {d_ef:+.3f}"
    )
    if onset:
        detail += (
            f"; intensification stretch ({onset['steps']} steps) d{latent} "
            f"{onset['d_latent']:+.4g}, "
            f"d{sensible} {onset['d_sensible']:+.4g}, "
            f"dEF {onset['d_evaporative_fraction']:+.3f}"
        )
    return CriterionResult(
        name="partition_shift",
        status=PASS if ok else FAIL,
        value=abs(imbalance) / scale,
        threshold=tolerance,
        message=(
            f"the latent heat the drying surface gave up warmed the air: {detail}"
            if ok
            else "; ".join(failures) + f" [{detail}]"
        ),
        diagnostics={
            "d_latent": d_le,
            "d_sensible": d_h,
            "d_ground": d_g,
            "imbalance": imbalance,
            "available": available,
            "ground_fraction": abs(d_g) / scale,
            "d_evaporative_fraction": d_ef,
            "window_steps": int(mask.sum()),
            "onset": onset,
        },
    )


def _as_bool(value, default=True):
    """`bool("false")` is True, and an empty YAML value is None, not False."""
    if value is None:
        return default
    if isinstance(value, str):
        text = value.strip().lower()
        if text == "":
            return default
        return text not in ("0", "false", "no", "off")
    return bool(value)


@criterion("melt_energy")
def melt_energy(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Melt reported as water must equal the melt the energy budget paid for.

    The third place the two ledgers meet, after vaporisation in
    `flux_identity` and the drydown partition in `partition_shift`. A model
    can melt a degree-day depth and close a surface energy budget that never
    mentions fusion; both budgets balance and no single-budget criterion sees
    it. Here the surface residual is not required to vanish. It is required to
    be the melt energy:

        mean(rn - hfls - hfss - hfg)  ==  [lambda_f * M + dU] / (N * dt * 86400)

    with `M` the ice the model says it lost. Two properties of the case make
    `M` observable through a contract that carries neither a melt flux nor a
    snowfall flux:

    Neither term is inferred from `snw` alone. `snw` is the pack's total water
    in this suite's own adapters -- Snow-17 reports WE + LIQW + lagged excess
    and SUMMA reports scalarSWE, ice plus liquid -- so a fall in `snw` is net
    water leaving the pack and not evidence of a phase change. Melt retained
    as liquid moves no `snw` at all, and drainage of water that melted days
    earlier moves it without any fusion happening now. The ice is therefore
    taken as `snw - lwsnl`, and only its change is charged at the latent heat
    of fusion:

        M = -d(snw - lwsnl) - sum(sbl dt)

    This criterion is deliberately narrow: it is written for a block that opens
    cold and closes warm. Refreezing is not scored and no credit for it exists:
    a block with net ice gain fails `M > 0`, and refreezing inside a block that
    then re-melts leaves both endpoints unchanged and is invisible to an
    integrated balance. Two of the
    checks below -- no liquid on the opening row, and specific cold content not
    rising across the block -- are consequences of that same precondition, and
    each is a parameter a probe with a different block can relax.

    Cold content is a term, not an assumption. A pack below freezing spends
    energy warming towards zero that does no melting, and

        dU = -( csnow_end - csnow_start )

    is what that costs, with `csnow` the pack's cold content: the energy still
    needed to bring its ice to 0 C. Requiring it is what makes warming and
    fusion distinguishable; without it the two are the same number to any
    criterion, whatever the case is engineered to do. It is asked for as the
    energy rather than as a pack temperature because that is the quantity a
    budget spends -- Snow-17 already carries it as `NEGHS` -- and because
    converting a temperature back through an assumed heat capacity is wrong for
    any model whose capacity differs from the assumed one.

    One property of the case is still load-bearing: the scored block carries
    no precipitation, so no snowfall adds ice the contract cannot see. The
    criterion refuses a block that is not dry rather than measuring it.

    Scored over each block the generator labels, so a stretch where the pack
    is building cannot dilute the stretch where it is going.
    """
    driver = str(params.get("driver", "rn"))
    sinks = list(params.get("sinks", ["hfls", "hfss", "hfg"]))
    pack = str(params.get("pack", "snw"))
    liquid_var = str(params.get("liquid", "lwsnl"))
    cold_content = str(params.get("cold_content", "csnow"))
    # The cap below is the cold content the reported ice would hold at
    # `cap_margin_k` beneath the coldest air in the record, applied only to the
    # rows the identity reads. It bounds how much of its own fusion a model can
    # cancel by declaring the pack cold on a boundary row; it is not a claim
    # that a pack can never be colder than the air.
    c_ice = float(params.get("specific_heat_ice", 2100.0))  # J kg-1 K-1
    cap_margin_k = float(params.get("cold_content_cap_margin_k", 10.0))
    air = str(params.get("air_temperature", "tas"))
    # Absolute, and set where round-off cannot reach them. `bounds.py` uses
    # 1e-6 mm for masses; 1 J m-2 of cold content moves the block balance by
    # under 1e-6 W m-2, and the largest excess an honest reference shows on the
    # specific-cold-content check below is 2.2e-8 J m-2 over a 20-seed sweep.
    mass_tol = float(params.get("mass_tolerance_mm", 1.0e-6))
    energy_tol = float(params.get("cold_content_tolerance_j", 1.0))
    # This case opens its block frozen, so an honest model carries no liquid
    # there. That is a property of this generated case and not of snowpacks in
    # general, which is why it is a parameter a different probe can relax.
    # 1 mm rather than a hair: a model with a smooth freezing curve holds trace
    # liquid in a cold pack, and an attacker buying room here gains at most
    # lambda_f * 1 mm / (60 * 86400) = 0.064 W m-2 against a 2 W m-2 floor.
    max_opening_liquid = float(params.get("max_opening_liquid_mm", 1.0))
    check_specific_cold = _as_bool(params.get("specific_cold_content_cannot_rise"), True)
    # A dry block's pack cannot gain water. This bounds gains in the pack alone
    # and is not a per-step form of `closure`, which is cumulative and spans
    # every store: a one-step excursion in one store that returns nets to
    # nothing there, which is why `closure` cannot see any of this.
    #
    # The tolerance is coupled to `floor`, and to the step. A model that lowers
    # the opening pack and closes the offset gradually at `g` mm per step
    # understates the demand by
    #
    #     lambda_f * g / (dt * 86400)   W m-2
    #
    # independent of block length, because the total gain divides by the same
    # number of steps. Keeping that inside the floor needs
    # g <= floor * dt * 86400 / lambda_f, which is 0.518 mm at a daily step and
    # the 2 W m-2 floor -- but only 0.0216 mm at PT1H, where this same 0.5 mm
    # would let about 46 W m-2 through. A sub-daily probe reusing this
    # criterion has to set `max_pack_gain_mm_per_step` from its own step.
    #
    # That figure is also the *least* the tolerance lets through, not a ceiling:
    # each step's gain is net of the drainage the model reports, so a pack
    # reported not to drain while its runoff carries the water away gets more
    # room. That is a multi-row fiction, named in the probe's Scope rather than
    # guarded; closing it needs a pack outflow flux.
    #
    # Deposition a model reports in `sbl` is netted out below and needs no
    # allowance; what is left is unreported gain. Change `floor`, the step or
    # this tolerance and check the other two.
    max_pack_gain = float(params.get("max_pack_gain_mm_per_step", 0.5))
    sublimation_var = str(params.get("sublimation", "sbl"))
    threshold = float(params.get("threshold", 0.05))
    floor = float(params.get("floor", 2.0))
    lam_f = float(params.get("lambda_fusion", LAMBDA_F))
    segment_column = str(params.get("segment_column", "_regime"))
    scored_label = str(params.get("scored_label", "melt"))
    # The pack the case builds has to actually appear, or a model that reports
    # no snow at all would satisfy the identity with two zeroes. A share, not a
    # depth, so it holds for any seed.
    #
    # How much of that pack a model then melts is deliberately not scored. A
    # model that melts part of it and pays for exactly that part is coupled
    # correctly, which is the only thing this criterion is entitled to judge;
    # failing it for melting slowly would score a calibration choice as a
    # conservation violation. See the probe README, "Scope".
    min_peak_share = float(params.get("min_peak_share_of_snowfall", 0.5))
    snowfall_label = str(params.get("snowfall_label", "accumulation"))
    precipitation = str(params.get("precipitation", "pr"))
    dry_tolerance_mm = float(params.get("dry_tolerance_mm", 1e-9))

    w = make_window(run, probe)
    if driver not in w.forcing.columns:
        raise ValueError(
            f"melt_energy needs forcing column '{driver}'; this probe's "
            "generator does not produce it"
        )
    for var in (*sinks, pack, liquid_var, cold_content):
        if var not in w.table.columns:
            raise ValueError(f"melt_energy needs '{var}' in the model result")

    # Beside the other column checks, not further down: a table that trips a
    # contract check while `tas` is missing should report the missing column
    # rather than whichever check it happened to reach first.
    if air not in w.forcing.columns:
        raise ValueError(
            f"melt_energy needs forcing column '{air}', in degrees Celsius, to "
            f"bound '{cold_content}' on the rows the identity reads. The cold "
            "content is a model's own statement about its own pack; without it "
            "the bound silently disappears and a model can cancel its own fusion "
            "by declaring the pack cold on one row"
        )

    coldest = float(w.forcing[air].to_numpy(dtype=float).min())

    blocks_all = segments(w, segment_column)
    scored = [b for b in blocks_all if b[0] == scored_label]
    if not scored:
        raise ValueError(
            f"melt_energy found no '{scored_label}' block in column "
            f"'{segment_column}'; the generator has to label the stretch the "
            "pack is melting over"
        )

    # The melt inference `M = -d(snw - lwsnl) - sum(sbl dt)` is only melt where nothing
    # is being added to the pack. This generator builds a dry block for exactly
    # that reason, but a probe reusing the criterion might not, and snowfall
    # inside the block nets out of the pack change and quietly shrinks the melt.
    # Refuse rather than measure the wrong thing.
    if precipitation not in w.forcing.columns:
        raise ValueError(
            f"melt_energy needs forcing column '{precipitation}' to check that "
            "the scored block is dry. The melt inference rests on nothing being "
            "added to the pack, and a criterion that cannot check its own "
            "precondition must not measure past it"
        )
    else:
        rain = w.forcing[precipitation].to_numpy(dtype=float)
        for label, start, stop in scored:
            fell = float(w.volume(rain[start:stop]).sum())
            if fell > dry_tolerance_mm:
                raise ValueError(
                    f"melt_energy scores '{label}' as a melt block, but "
                    f"{fell:.3f} mm of '{precipitation}' falls within it. Melt "
                    "is inferred from the pack, which holds only where nothing "
                    "is added to it: the scored block has to be dry"
                )

    drive = w.forcing[driver].to_numpy(dtype=float)
    fluxes = w.table[sinks].to_numpy(dtype=float)
    snw = w.table[pack].to_numpy(dtype=float)
    lwsnl = w.table[liquid_var].to_numpy(dtype=float)
    # The pack's cold content as the model reports it, J m-2 and positive. No
    # round trip through a temperature and an assumed heat capacity: this is
    # the quantity the budget spends, and the one Snow-17 already carries.
    csnow = w.table[cold_content].to_numpy(dtype=float)
    # The ice, which is the only thing a phase change moves.
    ice = snw - lwsnl
    subl = (
        w.table[sublimation_var].to_numpy(dtype=float)
        if sublimation_var in w.table.columns
        else np.zeros(len(w.table))
    )

    # `subl` belongs here because it feeds `melted` below. Without it a NaN in
    # one sbl row is not reported as the contract violation it is: it
    # propagates into the melt and the demand, and the criterion fails with a
    # message reading "melt of nan mm demands nan W/m2".
    finite = (
        np.isfinite(drive)
        & np.isfinite(fluxes).all(axis=1)
        & np.isfinite(snw)
        & np.isfinite(lwsnl)
        & np.isfinite(csnow)
        & np.isfinite(subl)
    )
    if not finite.all():
        n_bad = int((~finite).sum())
        return CriterionResult(
            name="melt_energy",
            status=FAIL,
            message=(
                f"non-finite snow, sublimation or energy values on {n_bad} "
                "scored steps"
            ),
            diagnostics={"non_finite_steps": n_bad},
        )

    # --- what the model reported has to be a pack at all -------------------
    #
    # Each of these is a contract violation rather than physics, and says so,
    # because each is a route by which a self-reported diagnostic on a single
    # boundary row could otherwise decide a verdict.
    def contract(message, value, **diag):
        return CriterionResult(
            name="melt_energy", status=FAIL, value=value,
            message="contract: " + message, diagnostics=diag,
        )

    for column, values, tol in (
        (pack, snw, mass_tol), (liquid_var, lwsnl, mass_tol),
        (cold_content, csnow, energy_tol),
    ):
        if np.any(values < -tol):
            worst = int(np.argmin(values))
            return contract(
                f"'{column}' is {float(values[worst]):.3g} on {int((values < -tol).sum())} "
                "scored steps. A stored mass and a cold content are both "
                "non-negative; a negative one on the step before the block "
                "would move the measured melt by its whole value",
                float(values[worst]), column=column, min_value=float(values[worst]),
            )
    if np.any(ice < -mass_tol):
        worst = int(np.argmin(ice))
        return contract(
            f"'{liquid_var}' exceeds '{pack}' on {int((ice < -mass_tol).sum())} scored "
            f"steps, worst {float(lwsnl[worst]):.1f} mm of liquid in a "
            f"{float(snw[worst]):.1f} mm pack. The liquid is held inside the pack "
            "and cannot be more than all of it",
            float(ice[worst]), min_ice_mm=float(ice[worst]),
        )
    empty_but_cold = (ice <= mass_tol) & (csnow > energy_tol)
    if empty_but_cold.any():
        worst = int(np.argmax(empty_but_cold))
        return contract(
            f"'{cold_content}' is {float(csnow[worst]):.3g} J/m2 on "
            f"{int(empty_but_cold.sum())} steps carrying no ice. Cold content is "
            "the energy needed to bring ice to 0 C, and there is none to bring",
            float(csnow[worst]), steps=int(empty_but_cold.sum()),
        )

    # The rows whose self-reports can actually move the verdict: the step before
    # each block, where `ice_before` and `csnow_before` are read, and its last
    # step, where `ice_after` and `csnow_after` are. Everything between them is
    # never read by the identity, and asserting a physical bound on states the
    # criterion does not evaluate would be a stronger claim than this probe is
    # entitled to make.
    def boundary(start, stop):
        if start <= 0:
            before = (float(w.state0[pack]) - float(w.state0[liquid_var]),
                      float(w.state0[liquid_var]), float(w.state0[cold_content]))
        else:
            before = (float(ice[start - 1]), float(lwsnl[start - 1]),
                      float(csnow[start - 1]))
        after = (float(ice[stop - 1]), float(lwsnl[stop - 1]), float(csnow[stop - 1]))
        return before, after

    # This case opens its block frozen -- the air is at or below -12 C through
    # accumulation and at or below -14 C for the block's first twenty days -- so
    # an honest model carries no liquid there. Without this, `ice_before` is a
    # self-report on a single row that no other term constrains, and moving it
    # moves the melt attributed to the model. The share guard above cannot
    # oppose that, because it is evaluated on the same row and moves with it.
    #
    # This is a property of the generated case, not of snowpacks: a pack may
    # legitimately hold liquid at the start of a melt block elsewhere, which is
    # why the tolerance is a parameter.
    for label, start, stop in scored:
        (_, liquid_before, _), _ = boundary(start, stop)
        if liquid_before > max_opening_liquid:
            return contract(
                f"'{liquid_var}' is {liquid_before:.3g} mm on the step before "
                f"the '{label}' block, over the {max_opening_liquid:g} mm "
                f"allowed. This probe opens its block on a frozen pack -- the "
                f"coldest air in the record is {coldest:.1f} C -- so an honest "
                "model carries no liquid there; liquid declared on that row "
                "moves the ice the melt is measured from, and no other term "
                "constrains it",
                liquid_before, block=label, opening_liquid_mm=liquid_before,
            )

    # A dry block's pack cannot gain water. It bounds gains in `snw` alone and
    # is not a per-step form of `closure`, which is cumulative and spans every
    # store: `ice_before` is read only to be differenced, so a dip on that row
    # moves the melt attributed to the model and nothing else notices.
    #
    # Deposition a model reports in `sbl` is netted out here, so an honest pack
    # that gains water by deposition passes whatever the tolerance.
    # `snw` alone, which is the store the identity differences. Summing the
    # canopy with it would leave the same dip reachable by moving the water
    # into `canopy` on the opening row: the sum does not move, and nothing
    # else looks at either store.
    for label, start, stop in scored:
        prior = float(w.state0[pack]) if start <= 0 else float(snw[start - 1])
        series = np.concatenate(([prior], snw[start:stop]))
        gain = np.diff(series) + w.volume(subl[start:stop])
        if gain.size and float(gain.max()) > max_pack_gain:
            worst = int(np.argmax(gain))
            when = w.forcing["time"].iloc[start + worst] if "time" in w.forcing else (
                start + worst
            )
            return contract(
                f"'{pack}' gains {float(gain[worst]):.3g} mm of water in one step "
                f"of the '{label}' block, at {when}, over the {max_pack_gain:g} mm "
                "allowed. The block carries no precipitation, so a pack that "
                "grows was not reported consistently, and a one-step excursion "
                "that returns nets to nothing in any cumulative check",
                float(gain[worst]), block=label, step=str(when),
                max_gain_mm=float(gain.max()),
            )

    # Cold content per unit ice cannot rise over a block that opens cold and
    # ends warm. Inflating `csnow_before` to buy room raises the demand by at
    # least as much. This closes one *construction* of the opening-row flip --
    # the one that copies the last row's cold content -- and not the route:
    # scaled to the declared ice instead, the excess is round-off and this
    # check is silent, leaving only the opening-liquid check. It does close the
    # variant that sets the last row's cold content at the cap.
    #
    # Tied to the same case precondition: specific cold content may legitimately
    # rise over a block that ends colder than it began, through refreezing or
    # renewed cooling, so a probe whose block does not end warm turns this off.
    if check_specific_cold:
        for label, start, stop in scored:
            (ice_b, _, csnow_b), (ice_a, _, csnow_a) = boundary(start, stop)
            if ice_b > mass_tol and ice_a > mass_tol:
                allowed = ice_a * csnow_b / ice_b
                if csnow_a > allowed + energy_tol:
                    return contract(
                        f"'{cold_content}' per unit ice rises over the '{label}' "
                        f"block: {csnow_a:.4g} J/m2 on {ice_a:.1f} mm of ice "
                        f"against the {allowed:.4g} J/m2 that the opening state "
                        "allows. This check assumes a block that opens cold and "
                        "ends warm, where the pack left cannot be colder per "
                        "kilogram than the pack it started from; a probe whose "
                        "block ends colder sets "
                        "specific_cold_content_cannot_rise false",
                        csnow_a, block=label, allowed_j=allowed,
                    )

    # The reported pack state has to stay physically plausible on the rows this
    # identity reads, which is narrower than saying a pack can never be colder
    # than the air. Snow-17 bounds its own deficit as a mass fraction
    # (NEGHS <= 0.33 * WE, worth 52 K) rather than as a temperature, and a thin
    # pack early in accumulation is legitimately colder than a temperature
    # bound allows. Those rows are never read by the identity, so they are not
    # bounded here.
    for label, start, stop in scored:
        (ice_b, _, csnow_b), (ice_a, _, csnow_a) = boundary(start, stop)
        for where, ice_row, csnow_row in (
            ("the step before", ice_b, csnow_b), ("the last step of", ice_a, csnow_a),
        ):
            with np.errstate(over="ignore", invalid="ignore"):
                cap = c_ice * ice_row * max(0.0, -(coldest - cap_margin_k))
            if np.isfinite(cap) and csnow_row > cap + energy_tol:
                return contract(
                    f"'{cold_content}' is {csnow_row:.3g} J/m2 on {where} the "
                    f"'{label}' block, carrying {ice_row:.1f} mm of ice, above "
                    f"the {cap:.3g} J/m2 that ice would hold at "
                    f"{coldest - cap_margin_k:.1f} C, {cap_margin_k:.0f} K below "
                    "the coldest air in the record. The reported pack state has "
                    "to stay physically plausible on the rows this identity "
                    "reads",
                    csnow_row, block=label, row=where, cap_j=float(cap),
                )

    # Measured on the ice standing at the step before the block opens, not on
    # the peak anywhere in the window. A peak lets a model that melts and drains
    # its pack before the block still satisfy the guard, and lets one boundary
    # row decide it. The winter here stays at or below -12 C, so an honest model
    # still holds all of its snowfall when the block opens.
    opening = [
        (float(w.state0[pack]) - float(w.state0[liquid_var])) if start <= 0
        else float(ice[start - 1])
        for _, start, _ in scored
    ]
    peak = float(max(opening)) if opening else 0.0
    snowfall = sum(
        float(w.volume(w.forcing[precipitation].to_numpy(dtype=float)[start:stop]).sum())
        for label, start, stop in blocks_all
        if label == snowfall_label
    )
    if snowfall > 0.0 and peak < min_peak_share * snowfall:
        return CriterionResult(
            name="melt_energy",
            status=FAIL,
            value=peak,
            threshold=min_peak_share * snowfall,
            message=(
                f"no pack to melt: {peak:.1f} mm of ice where the block opens, "
                f"against {snowfall:.1f} mm of snowfall, under the "
                f"{min_peak_share:.0%} the case requires"
            ),
            diagnostics={"peak_swe_mm": peak, "snowfall_mm": snowfall},
        )

    # Finite inputs can still overflow: a melt large enough sends `demanded` to
    # infinity, and with it the allowance, so `gap <= allowance` becomes
    # `inf <= inf` and a meaningless run reports PASS. Same defence the
    # radiation criterion carries, and for the same reason.
    with np.errstate(over="ignore", invalid="ignore"):
        residual = drive - fluxes.sum(axis=1)
    seconds = w.dt_days * SECONDS_PER_DAY
    blocks = []
    for label, start, stop in scored:
        # Ice at the step before the block, through the shared helper that owns
        # the off-by-one, and at its last step.
        if start <= 0:
            ice_before = float(w.state0[pack]) - float(w.state0[liquid_var])
            csnow_before = float(w.state0[cold_content])
        else:
            ice_before = float(ice[start - 1])
            csnow_before = float(csnow[start - 1])
        ice_after, csnow_after = float(ice[stop - 1]), float(csnow[stop - 1])

        with np.errstate(over="ignore", invalid="ignore"):
            # Fusion: only the ice that changed phase, sublimation removed.
            melted = -(ice_after - ice_before) - float(w.volume(subl[start:stop]).sum())
            # Cold content: what warming a sub-freezing pack cost, which does
            # no melting and which no case can be engineered to rule out.
            # Spending it reduces the deficit, so the energy is minus the change.
            d_internal = -(csnow_after - csnow_before)
            n = stop - start
            demanded = (lam_f * melted + d_internal) / (n * seconds)
            available = float(residual[start:stop].mean())
            gap = abs(available - demanded)
            allowance = max(threshold * abs(demanded), floor)
        if not all(np.isfinite(v) for v in
                   (melted, d_internal, demanded, available, gap, allowance)):
            return CriterionResult(
                name="melt_energy",
                status=FAIL,
                message=(
                    f"the melt-energy calculation overflowed on the "
                    f"'{label}' block; the reported pack or surface fluxes are "
                    "outside the range the identity can be evaluated in"
                ),
                diagnostics={
                    "block": label,
                    "melt_mm": melted,
                    "demanded_w_m2": demanded,
                    "available_w_m2": available,
                },
            )
        left = float(snw[stop - 1])
        # Reported so a reviewer can see how far the melt got, never scored.
        blocks.append({
            "label": label,
            "start": start,
            "stop": stop,
            "melt_mm": melted,
            "fusion_w_m2": lam_f * melted / (n * seconds),
            "cold_content_w_m2": d_internal / (n * seconds),
            "demanded_w_m2": demanded,
            "available_w_m2": available,
            "gap_w_m2": gap,
            "allowance_w_m2": allowance,
            "slack": gap / allowance if allowance > 0 else float("inf"),
            "snw_end_mm": left,
            "share_of_peak_left": left / peak if peak > 0 else 0.0,
            "passed": gap <= allowance and melted > 0.0,
        })

    stalled = [b for b in blocks if b["melt_mm"] <= 0.0]
    if stalled:
        return CriterionResult(
            name="melt_energy",
            status=FAIL,
            value=stalled[0]["melt_mm"],
            threshold=0.0,
            message=(
                f"the pack did not melt over the {scored_label} block: "
                f"{stalled[0]['melt_mm']:.2f} mm lost from a {peak:.1f} mm pack "
                "under forcing well above freezing"
            ),
            diagnostics={"blocks": blocks, "peak_swe_mm": peak},
        )

    failed = sum(not b["passed"] for b in blocks)
    worst = max(blocks, key=lambda b: b["slack"])
    return CriterionResult(
        name="melt_energy",
        status=PASS if failed == 0 else FAIL,
        value=worst["slack"],
        threshold=1.0,
        message=(
            f"melt of {worst['melt_mm']:.1f} mm demands "
            f"{worst['demanded_w_m2']:.2f} W/m2 and the surface budget has "
            f"{worst['available_w_m2']:.2f} W/m2"
            + (
                f"; agree within {worst['allowance_w_m2']:.2f} W/m2"
                if failed == 0
                # Signed, because the failing branch fires in both directions
                # and "paid for melt that never appeared as water" is the
                # opposite defect from "melted ice it never paid for".
                else (
                    f"; {'short by' if worst['available_w_m2'] < worst['demanded_w_m2'] else 'over by'}"
                    f" {worst['gap_w_m2']:.2f} W/m2, past the "
                    f"{worst['allowance_w_m2']:.2f} W/m2 allowed"
                )
            )
        ),
        diagnostics={"blocks": blocks, "peak_swe_mm": peak, "n_blocks": len(blocks)},
    )

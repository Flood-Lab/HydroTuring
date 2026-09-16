"""A declared head-driven exchange must answer the head it was given.

`gwex` is the contract's confession channel: a model whose catchment gains or
loses water across its boundary — regional groundwater, an inter-basin
transfer, a prescribed withdrawal — declares that flux, and its budget closes
over a term it admits to rather than hiding it in the residual. `AGENTS.md`
puts it plainly: declared, it is a source in the budget and the budget can
close; hidden, it is the residual. What the contract does not yet say is what a
declared exchange may look like, and that is the escape this criterion closes.

Write the budget without the exchange,

    R0 = pr - evspsbl - mrro - d(stores)/dt          [mm per step]

so that the budget the suite checks is R = R0 + gwex. A model that closes has
R ~ 0 and therefore gwex ~ -R0 identically: on a CREST implementation over ten
years, max |R| per step is 2e-13 mm/day. So comparing the declared flux with
the residual it closes is circular — it is the same number twice — and nothing
about the *magnitude* of gwex, or its timing on its own, can separate an honest
exchange from an invented one. A discriminator has to use information that is
not algebraically fixed by closure. This criterion uses a counterfactual: the
same weather, the same model, the same seed, and one change to the one driver
the case prescribes.

The case supplies an external hydraulic head, `gwh`, as a visible forcing
column, and runs the record three times: with the head as given, raised by a
constant, and lowered by the same constant. The shift begins with the first
scored step; the spinup is byte-identical across the runs. That timing is not
incidental. A responsive aquifer answers a shift by admitting the water that
raises its own head to meet the new external one, and then by nothing: an
early draft applied the shift from the start of spinup, and a native MODFLOW
aquifer with a five-day time constant had finished answering before scoring
began, so its paired difference read as zero and it failed for responding
correctly. With the shift at the scored start, the answer falls inside the
window that is scored.

`gwh` is the prescribed external hydraulic head associated with the model's
declared `gwex`, and `gwex` is positive into the catchment. The probe requires
the column, so a model that does not list it in `needs_forcing` or
`uses_forcing` is INCOMPATIBLE and never reaches this criterion: declaring it
is a semantic opt-in — the model asserts that its external exchange responds
monotonically to this potential, as a general-head boundary does — and it is
what the criterion holds the model to. With

    G       = sum_t gwex_t * dt          [mm, over the scored record]
    G_gross = sum_t |gwex_t| * dt        [mm, control run]

the assertion is

    G(raised)  - G(control)  >=  +max(s * G_gross, eps)
    G(lowered) - G(control)  <=  -max(s * G_gross, eps)

with s a small share and eps a floating-point tolerance. Raising the external
head must bring more water in, and lowering it less, by a share of the
exchange the model itself declared.

Why a share and not an absolute minimum. An absolute floor would fail a
genuine boundary of small conductance, whose response is small; that was the
argument an earlier draft made for having no floor at all. It does not carry
over to a share. For a general-head boundary Q = C (H - h) the response to a
shift and the gross exchange both scale with C, so their ratio does not:
against a constant internal head it is dh * T / sum|H - h| dt, about 2 for this
head series; against an internal head that moves it is smaller, bounded below
by about dh / TV(H), the shift over the head's total variation — about 2e-3
here — because in the fast limit both the response (S * dh) and the gross
(S * TV(H)) scale with storativity and neither with conductance. The floor is
set at a share below that physical minimum, which no honest boundary of any
conductance can fall under. What falls under it is the model the criterion
exists to catch: one whose declared exchange is the day's accounting error
with a token head-proportional term added. Its response scales with the token
conductance; its gross does not, because the gross is the error. The token
control in this repository sits six orders of magnitude below the floor.

What the share does not close. A model whose declared `gwex` mixes a genuine
head-driven part with a large unrelated one — a deep loss, a withdrawal — has a
gross that the unrelated part inflates and a share that is honest but small.
Such a model is indistinguishable from the token cheat from outside, and the
floor is placed low, at a thousandth, so as to tolerate a mixture of that order
rather than to catch every conceivable cheat. The residual escape is a model
that reads the head, keeps the sink, and sizes its token term to a thousandth
of the sink: disclosed, and the price of not asking the contract to split the
head-driven part of `gwex` into its own variable.

Two earlier drafts gated on how often the exchange reversed direction, on its
own and then relative to the head's turning points. Both were wrong — the
second by a false theorem, since between two turning points of the external
head the flux crosses zero as often as the internal head moves — and every
reversal statistic is now reported and none of it is gated.

Two things are decided regardless. A declared loss cannot exceed the water
that was there to lose, in every variant. And an exchange below rounding error
is not an exchange: under a fixed share of the record's precipitation the
criterion passes with the reason named, never an unscored state.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (FAIL, PASS, CriterionResult, Window,
                                       criterion, make_window, segments)
from hydroturing.criteria.response import pick
from hydroturing.protocol import RunResult
from hydroturing.spec import STATE_VARS, ProbeSpec


def reversal_fraction(g: np.ndarray) -> tuple[float, int]:
    """Magnitude-weighted share of a series that is spent reversing direction.

    Each sign change is charged the smaller of the two magnitudes that
    straddle it, so passing slowly through zero is almost free and flipping at
    full amplitude is not. Returns the share and the raw number of sign changes.
    """
    g = np.asarray(g, dtype=float)
    gross = float(np.abs(g).sum())
    if len(g) < 2 or gross <= 0.0:
        return 0.0, 0
    flip = np.sign(g[1:]) * np.sign(g[:-1]) < 0
    charged = np.minimum(np.abs(g[1:]), np.abs(g[:-1])) * flip
    return float(charged.sum() / gross), int(flip.sum())


def significant_reversals(g: np.ndarray, deadband: float) -> int:
    """Sign changes whose straddling magnitudes both clear a share of the mean."""
    g = np.asarray(g, dtype=float)
    if len(g) < 2:
        return 0
    scale = float(np.abs(g).mean())
    if scale <= 0.0:
        return 0
    flip = np.sign(g[1:]) * np.sign(g[:-1]) < 0
    big = np.minimum(np.abs(g[1:]), np.abs(g[:-1])) >= deadband * scale
    return int((flip & big).sum())


def turning_points(h: np.ndarray, deadband: float) -> int:
    """Local extrema of a series: sign changes of its step-to-step change."""
    return significant_reversals(np.diff(np.asarray(h, dtype=float)), deadband)


def _reversal_diagnostics(w: Window, var: str, driver: str, label: str,
                          column: str, deadband: float) -> dict:
    """Everything the earlier drafts gated on, reported and not judged."""
    g = w.volume(w.table[var])
    h = np.asarray(w.forcing[driver], dtype=float)
    frac, flips = reversal_fraction(g)
    out = {"record_reversal_fraction": frac, "record_sign_changes": flips}
    blocks = [(a, b) for lab, a, b in segments(w, column) if lab == label]
    if blocks:
        rev_g = rev_h = 0
        num = den = 0.0
        for a, b in blocks:
            seg = g[a:b]
            rev_g += significant_reversals(seg, deadband)
            rev_h += turning_points(h[a:b], deadband)
            f, _ = reversal_fraction(seg)
            d = float(np.abs(seg).sum())
            num += f * d
            den += d
        out.update({
            "quiescent_exchange_reversals": rev_g,
            "quiescent_driver_turning_points": rev_h,
            "quiescent_reversal_fraction": (num / den) if den > 0 else 0.0,
            "quiescent_steps": int(sum(b - a for a, b in blocks)),
        })
    return out


def _availability_shortfall(w: Window, var: str) -> tuple[int, float]:
    """Steps on which a declared loss exceeds the water there was to lose."""
    g = w.volume(w.table[var])
    present = [v for v in STATE_VARS if v in w.table.columns]
    stored = (w.table[present].to_numpy(dtype=float).sum(axis=1) if present
              else np.zeros(len(w.table)))
    opening = np.concatenate([[w.storage_initial(tuple(present))], stored[:-1]])
    supply = opening + (w.volume(w.forcing["pr"]) if "pr" in w.forcing.columns
                        else np.zeros(len(g)))
    short = np.maximum(-g - supply, 0.0)
    return int((short > 1e-6).sum()), (float(short.max()) if len(short) else 0.0)


@criterion("exchange_response", paired=True)
def exchange_response(runs: dict[str, RunResult], probe: ProbeSpec, params: dict) -> CriterionResult:
    """Raising the prescribed head must bring more water in; lowering it, less."""
    var = str(params.get("variable", "gwex"))
    driver = str(params.get("driver", "gwh"))
    share = float(params.get("min_response_share", 1.0e-3))
    epsilon = float(params.get("epsilon", 1.0e-8))
    label = str(params.get("quiescent_label", "dry"))
    column = str(params.get("segment_column", "_regime"))
    deadband = float(params.get("deadband", 0.05))
    negligible = float(params.get("negligible_share_of_pr", 1.0e-4))
    check_avail = bool(params.get("check_availability", True))
    weather = list(params.get("weather", ["pr", "tas", "pet"]))

    control_run = pick(runs, params, "control", "control")
    raised_run = pick(runs, params, "raised", "raised")
    lowered_run = pick(runs, params, "lowered", "lowered")
    control = make_window(control_run, probe)
    raised = make_window(raised_run, probe)
    lowered = make_window(lowered_run, probe)

    for w in (control, raised, lowered):
        if var not in w.table.columns:
            raise ValueError(
                f"exchange_response needs '{var}' in the model result; the probe "
                "has to require it so a model without one is N/A rather than judged"
            )
        if driver not in w.forcing.columns:
            raise ValueError(
                f"exchange_response needs the prescribed driver '{driver}' in the "
                "forcing; the generator has to supply the head the exchange answers"
            )
    if not (len(control.table) == len(raised.table) == len(lowered.table)):
        raise ValueError("the variants produced windows of different lengths")

    # The probe requires the driver, so the harness only sends models that
    # declared it. A manifest that reaches this point without the driver is a
    # contract violation and is said so, not silently excused; a RunResult with
    # no manifest at all — built by hand, as the tests do — is judged, because
    # the default of this criterion is to judge.
    model = control_run.model
    if model is not None:
        declared = tuple(model.needs_forcing) + tuple(model.uses_forcing)
        if driver not in declared:
            raise ValueError(
                f"exchange_response reached a model that does not declare '{driver}' "
                "in needs_forcing or uses_forcing; the probe's requires.forcing should "
                "have made it INCOMPATIBLE before any run"
            )

    # Only the head may differ between the variants, and only from the first
    # scored step: the spinup has to be identical so that every variant enters
    # the scored record in the same state and the response falls inside it.
    full = {"control": control_run.case.forcing, "raised": raised_run.case.forcing,
            "lowered": lowered_run.case.forcing}
    spin = control_run.case.spinup_steps
    for col in weather:
        if col not in full["control"].columns:
            continue
        a = np.asarray(full["control"][col], dtype=float)
        for name in ("raised", "lowered"):
            b = np.asarray(full[name][col], dtype=float)
            if len(a) != len(b) or not np.allclose(a, b, rtol=0, atol=1e-9):
                raise ValueError(
                    f"'{col}' differs between control and {name}; only '{driver}' may "
                    "change, or the response cannot be attributed to it"
                )
    h_ctl = np.asarray(full["control"][driver], dtype=float)
    for name in ("raised", "lowered"):
        h = np.asarray(full[name][driver], dtype=float)
        if not np.allclose(h[:spin], h_ctl[:spin], rtol=0, atol=1e-9):
            raise ValueError(
                f"'{driver}' differs between control and {name} during spinup; the "
                "shift has to begin with the scored record, or a fast aquifer answers "
                "it before scoring begins and the response is never seen"
            )
    shift_up = float(np.mean(np.asarray(raised.forcing[driver], dtype=float)
                             - np.asarray(control.forcing[driver], dtype=float)))
    shift_down = float(np.mean(np.asarray(lowered.forcing[driver], dtype=float)
                               - np.asarray(control.forcing[driver], dtype=float)))
    if not (shift_up > 0 and shift_down < 0):
        raise ValueError(
            f"the raised variant must lift '{driver}' and the lowered one drop it "
            f"over the scored record (got {shift_up:+.3g} and {shift_down:+.3g})"
        )

    g0 = control.volume(control.table[var])
    gross = float(np.abs(g0).sum())
    net = float(g0.sum())
    rain = float(control.volume(control.forcing["pr"]).sum()) if "pr" in control.forcing.columns else 0.0

    floor = negligible * rain if rain > 0 else 0.0
    if gross <= floor or gross <= 0.0:
        return CriterionResult(
            name="exchange_response", status=PASS, value=None,
            message=(f"negligible exchange: {gross:.3g} mm gross over the control "
                     f"record, under {negligible:g} of its {rain:.0f} mm of rain"
                     if gross > 0 else "no exchange declared"),
            diagnostics={"gross_mm": gross, "net_mm": net, "outcome": "negligible",
                         "gross_share_of_pr": (gross / rain) if rain else None},
        )

    up = float(raised.volume(raised.table[var]).sum() - net)
    down = float(lowered.volume(lowered.table[var]).sum() - net)
    # A share of the model's own gross exchange, floored by a floating-point
    # tolerance. The share is conductance-invariant for a head-driven boundary
    # (see the module docstring), so a weak boundary is not penalised for
    # being weak; a token head term on top of an accounting sink is.
    required = max(share * gross, epsilon * max(gross, 1.0))
    response_share = min(up, -down) / gross

    diagnostics = {
        "driver": driver,
        "head_shift_m": {"raised": shift_up, "lowered": shift_down},
        "response_mm": {"raised": up, "lowered": down},
        "response_share_of_gross": response_share,
        "required_share": share,
        "required_mm": required,
        "gross_mm": gross,
        "net_mm": net,
        "directionality": abs(net) / gross,
        "gross_share_of_pr": (gross / rain) if rain else None,
        "net_share_of_pr": (abs(net) / rain) if rain else None,
        "net_direction": "into the catchment" if net > 0 else "out of the catchment",
    }
    diagnostics.update(_reversal_diagnostics(control, var, driver, label, column, deadband))

    failures = []
    if abs(up) < required and abs(down) < required:
        failures.append(
            f"the model declares it consumes '{driver}' but its exchange is the same "
            f"with the head raised by {shift_up:+.2f} m and lowered by {shift_down:+.2f} m "
            f"({up:+.3g} mm and {down:+.3g} mm against a required {required:.3g} mm, "
            f"{share:g} of its {gross:.0f} mm gross exchange); it does not answer the "
            f"driver it claims to follow"
        )
    else:
        if up < required:
            failures.append(
                f"raising the prescribed head by {shift_up:+.2f} m changed the "
                f"integrated exchange by {up:+.3g} mm where at least {required:.3g} mm "
                f"more inflow was required ({share:g} of the {gross:.0f} mm gross exchange)"
            )
        if down > -required:
            failures.append(
                f"lowering it by {shift_down:+.2f} m changed the integrated exchange "
                f"by {down:+.3g} mm where at least {required:.3g} mm less inflow was "
                f"required ({share:g} of the {gross:.0f} mm gross exchange)"
            )
    diagnostics["outcome"] = "failed" if failures else "judged"

    n_short = None
    if check_avail:
        n_short, worst = 0, 0.0
        for w in (control, raised, lowered):
            n, m = _availability_shortfall(w, var)
            n_short += n
            worst = max(worst, m)
        if n_short:
            failures.append(
                f"a declared loss exceeds the water available to lose on "
                f"{n_short} steps across the variants (worst by {worst:.3g} mm)"
            )
    diagnostics["steps_short_of_availability"] = n_short

    ok = not failures
    return CriterionResult(
        name="exchange_response",
        status=PASS if ok else FAIL,
        value=response_share,
        threshold=share,
        message=(
            f"the declared exchange answers the prescribed head: {up:+.4g} mm over "
            f"the record with it raised {shift_up:+.1f} m, {down:+.4g} mm with it "
            f"lowered, {response_share:.3g} of its {gross:.0f} mm gross exchange "
            f"(required {share:g}); stays available"
            if ok else "; ".join(failures)
        ),
        diagnostics=diagnostics,
    )


__all__ = ["exchange_response", "reversal_fraction", "significant_reversals", "turning_points"]

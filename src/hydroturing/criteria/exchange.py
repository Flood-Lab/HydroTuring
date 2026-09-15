"""A declared exchange must answer the driver it was given.

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

The case supplies an external head, `gwh`, as a visible forcing column, and runs
the record three times: with the head as given, raised by a constant, and
lowered by the same constant. `gwh` is the prescribed external hydraulic head
associated with the model's declared `gwex`, and `gwex` is positive into the
catchment. A model that lists `gwh` in `needs_forcing` or `uses_forcing` is
thereby asserting that its external exchange responds monotonically to this
potential — the semantic of a general-head boundary, `Q = C (H_ext - h)`. That
is an opt-in, not a consequence of reading the column, and it is what the
criterion holds the model to. With

    G = sum_t gwex_t * dt          [mm over the scored record]

the assertion is

    G(raised) > G(control) > G(lowered)

up to a numerical tolerance of 1e-8 * max(G_gross, 1 mm), where G_gross is
sum_t |gwex_t| dt over the control. The tolerance is for floating point and
nothing else. No physical minimum is asked for: a boundary with a very small
conductance responds by a very small amount and is entirely physical, so the
only thing the gate excludes is a response that is absent or wrongly signed. A
model that declares the head and never reads it produces the same exchange in
all three runs — the harness gives every variant of a seed the same model seed
— and answers with exactly zero. That is the cheat: the declared-exchange
channel used as a sink for the day's accounting error, with a claim of physical
provenance it cannot back.

The paired design isolates the intervention. `gwh` is the only quantity varied
between the runs; the model's internal states and heads of course respond and
differ, and that is the response being measured.

Two earlier drafts of this criterion gated on how often the exchange reversed
direction, first on its own inside rainless windows, then relative to the
prescribed head's turning points. Both were wrong, and the second was wrong by
a false theorem: between two turning points of the external head the flux
`h_ext - h_catchment` can cross zero any number of times if the catchment's own
head moves, and in any real aquifer it does — a native MODFLOW run of a confined
aquifer behind an oscillating general-head boundary reverses its boundary flux
with no rain, no ET and no pumping. Reversal frequency is not a law of
head-driven exchange unless the internal head is constrained too, and the case
does not prescribe the internal head. So every reversal statistic this
criterion computes — inside the rainless windows, over the record, and against
the head's turning points — is reported and none of it is gated.

Two things are decided for every model regardless of what it declares. A
declared loss cannot exceed the water that was there to lose. And an exchange
below rounding error is not an exchange: under a fixed share of the record's
precipitation the criterion passes with the reason named, never an unscored
state. The magnitude of the exchange is reported and never gated.
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


@criterion("exchange_response", paired=True)
def exchange_response(runs: dict[str, RunResult], probe: ProbeSpec, params: dict) -> CriterionResult:
    """Raising the prescribed head must bring more water in; lowering it, less."""
    var = str(params.get("variable", "gwex"))
    driver = str(params.get("driver", "gwh"))
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

    # Only the head may differ between the variants, or the response cannot be
    # attributed to it.
    for col in weather:
        if col not in control.forcing.columns:
            continue
        a = np.asarray(control.forcing[col], dtype=float)
        for other in (raised, lowered):
            b = np.asarray(other.forcing[col], dtype=float)
            if not np.allclose(a, b, rtol=0, atol=1e-9):
                raise ValueError(
                    f"'{col}' differs between the variants; only '{driver}' may "
                    "change, or the response cannot be attributed to it"
                )
    h0 = np.asarray(control.forcing[driver], dtype=float)
    shift_up = float(np.mean(np.asarray(raised.forcing[driver], dtype=float) - h0))
    shift_down = float(np.mean(np.asarray(lowered.forcing[driver], dtype=float) - h0))
    if not (shift_up > 0 and shift_down < 0):
        raise ValueError(
            f"the raised variant must lift '{driver}' and the lowered one drop it "
            f"(got {shift_up:+.3g} and {shift_down:+.3g})"
        )

    g0 = control.volume(control.table[var])
    gross = float(np.abs(g0).sum())
    net = float(g0.sum())
    rain = float(control.volume(control.forcing["pr"]).sum()) if "pr" in control.forcing.columns else 0.0

    declared = ()
    if control_run.model is not None:
        declared = tuple(control_run.model.needs_forcing) + tuple(control_run.model.uses_forcing)
    consumes = driver in declared

    floor = negligible * rain if rain > 0 else 0.0
    if gross <= floor or gross <= 0.0:
        return CriterionResult(
            name="exchange_response", status=PASS, value=None,
            message=(f"negligible exchange: {gross:.3g} mm gross over the control "
                     f"record, under {negligible:g} of its {rain:.0f} mm of rain"
                     if gross > 0 else "no exchange declared"),
            diagnostics={"gross_mm": gross, "net_mm": net, "outcome": "negligible",
                         "driver_consumed": consumes,
                         "gross_share_of_pr": (gross / rain) if rain else None},
        )

    up = float(raised.volume(raised.table[var]).sum() - net)
    down = float(lowered.volume(lowered.table[var]).sum() - net)
    # Floating-point tolerance only. A tiny conductance gives a tiny response
    # and is physical; what is excluded is no response, or the wrong sign.
    required = epsilon * max(gross, 1.0)

    diagnostics = {
        "driver": driver,
        "driver_consumed": consumes,
        "head_shift_m": {"raised": shift_up, "lowered": shift_down},
        "response_mm": {"raised": up, "lowered": down},
        "tolerance_mm": required,
        "gross_mm": gross,
        "net_mm": net,
        "directionality": abs(net) / gross,
        "gross_share_of_pr": (gross / rain) if rain else None,
        "net_share_of_pr": (abs(net) / rain) if rain else None,
        "net_direction": "into the catchment" if net > 0 else "out of the catchment",
    }
    diagnostics.update(_reversal_diagnostics(control, var, driver, label, column, deadband))

    failures = []
    if consumes:
        if up < required:
            failures.append(
                f"raising the prescribed head by {shift_up:+.2f} m changed the "
                f"integrated exchange by {up:+.3g} mm; more inflow was required "
                f"(tolerance {required:.2g} mm, floating point only)"
            )
        if down > -required:
            failures.append(
                f"lowering it by {shift_down:+.2f} m changed the integrated exchange "
                f"by {down:+.3g} mm; less inflow was required "
                f"(tolerance {required:.2g} mm, floating point only)"
            )
        if failures and abs(up) <= required and abs(down) <= required:
            failures = [
                f"the model declares it consumes '{driver}' but its exchange is "
                f"identical with the head raised by {shift_up:+.2f} m and lowered by "
                f"{shift_down:+.2f} m ({up:+.3g} mm and {down:+.3g} mm); it does not "
                f"read the driver it claims to follow"
            ]
        diagnostics["outcome"] = "judged" if not failures else "failed"
    else:
        diagnostics["outcome"] = "driver not consumed; response reported, not judged"

    n_short = None
    if check_avail:
        worst = 0.0
        n_short = 0
        for w in (control, raised, lowered):
            g = w.volume(w.table[var])
            present = [v for v in STATE_VARS if v in w.table.columns]
            stored = (w.table[present].to_numpy(dtype=float).sum(axis=1) if present
                      else np.zeros(len(w.table)))
            opening = np.concatenate([[w.storage_initial(tuple(present))], stored[:-1]])
            supply = opening + (w.volume(w.forcing["pr"]) if "pr" in w.forcing.columns
                                else np.zeros(len(g)))
            short = np.maximum(-g - supply, 0.0)
            n_short += int((short > 1e-6).sum())
            worst = max(worst, float(short.max()) if len(short) else 0.0)
        if n_short:
            failures.append(
                f"a declared loss exceeds the water available to lose on "
                f"{n_short} steps across the variants (worst by {worst:.3g} mm)"
            )
    diagnostics["steps_short_of_availability"] = n_short

    ok = not failures
    if ok:
        if consumes:
            message = (f"the declared exchange answers the prescribed head: "
                       f"{up:+.4g} mm over the record with it raised {shift_up:+.1f} m, "
                       f"{down:+.4g} mm with it lowered; stays available")
        else:
            message = (f"the model does not declare it consumes '{driver}', so its "
                       f"response ({up:+.1f} / {down:+.1f} mm) is reported and not judged; "
                       f"stays available")
    return CriterionResult(
        name="exchange_response",
        status=PASS if ok else FAIL,
        value=min(up, -down) if consumes else None,
        threshold=required if consumes else None,
        message=message if ok else "; ".join(failures),
        diagnostics=diagnostics,
    )


__all__ = ["exchange_response", "reversal_fraction", "significant_reversals", "turning_points"]

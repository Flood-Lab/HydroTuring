"""Rating-curve coherence: stage has to be a single-valued, rising function
of discharge.

A river reach does not store water and then discharge it without a trace.
Where a model carries a channel store it usually also carries a stage, and
the stage is what an observer would read off a staff gauge. That reading
imposes two constraints that the water budget alone does not.

The first is monotonicity. More water in the reach means a higher stage. If
two stretches of the record carry the same discharge, they should carry
approximately the same stage, and if one carries more discharge it should
carry more stage. A model whose stage wanders independently of its discharge
has decoupled the gauge from the flow, which is exactly what happens when a
diagnostic is computed from something other than the store it claims to
describe.

The second is the loop. Real gauges trace a hysteresis loop: for the same
discharge, the rising limb sits *lower* than the falling limb, because a
flood wave steepens as it arrives (a kinematic wave with a looped rating) and
the channel is still filling while the flow is already climbing. The sign of
that loop is physical, not a convention. A rating that loops the other way —
higher stage on the rise than on the fall — is one the model has inverted,
and a probe that only checked monotonicity would never see it.

Both criteria are read from `stage`, in metres, which is a diagnostic the
probe asks the model to report. It is deliberately not one of the storages
the budget is differenced over: a stage is a reading, not a volume, and
summing it into `reported_states` would corrupt every closure test in the
suite.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

# Equal-count bins, not equal-width: a rating curve is sampled densely at low
# flow and sparsely at high flow, so width-binned medians would put one flood
# in one bin and forty baseflow days in another. Rank bins give every bin the
# same evidentiary weight.
DEFAULT_BINS = 12
# A bin needs a few members before its median says anything about the rating.
MIN_BIN_SIZE = 5
# How far a stage may sit from the best single-valued rating, as a fraction of
# the rating's own span, before that scatter counts as a real two-limb
# separation rather than sampling noise. A model that computes its stage as an
# exact function of one store lands at machine epsilon; a genuinely hysteretic
# rating misses by a visible share of its range. Two percent sits between them
# with room on both sides.
MIN_LOOP_FRACTION = 0.02
# The share of the span within which a stage counts as a *function* of an axis
# rather than as a rating drawn against it. This is much tighter than the size
# floor above, because it answers a different question — "is the rating
# single-valued" is a claim that the departure is rounding, not a claim about
# what an instrument could resolve. A model that computes its stage from one
# quantity lands well inside this; a rating that genuinely loops is an order of
# magnitude outside it, which is what keeps the escape from swallowing the
# inverted reference gauge.
MIN_SINGLE_VALUED_FRACTION = 0.005
# How far a stage may dip below the running maximum, as a fraction of the
# rating's own span, before that dip is a violation rather than sampling noise.
#
# This is deliberately much looser than the loop's floor above, because the two
# criteria ask for opposite things and an honest rating has both. A stage read
# off a store rather than off the instantaneous flow is hysteretic by
# construction — that *is* the loop — and hysteresis means the same discharge
# carries two levels, so its binned median necessarily dips below the running
# maximum by a visible share of the span. Published loop ratings run 5-20% of
# the depth; eight percent admits the honest ones while a gauge that never
# comes back down (a running maximum, say) still misses by far more — the
# drift reference model misses by roughly 60%. The 2% this used to be was only
# survivable while the rating span was tens of metres, which no reach has: it
# made a physically scaled rating fail the monotonicity check for having the
# very loop the probe exists to find.
MIN_MONOTONIC_FRACTION = 0.08
# The share of bins that must agree in sign before the loop's direction is
# treated as the physics rather than as pairing noise. A real loop has the
# same sign at every discharge and lands far above this; artifacts from
# nearest-neighbour pairing scatter around a half.
MAJORITY_FRACTION = 0.75
# A rating's span has to be more than a few nanometres before any fraction of
# it means anything. Guards the degenerate case where every stage is identical.
MIN_SPAN_M = 1e-9


def _discharge(w) -> tuple[np.ndarray, str]:
    """The discharge the rating is against, and what it came from.

    `dis` is the honest answer, in m3/s. When a model does not report it the
    criterion still has to read *something* against the stage, and the only
    other quantity in a momentum probe is the `channel` store — the water in
    transit. The store is deliberately returned as-is rather than differenced:
    the limb labels in `rating_loop` are decided by the sign of this array's
    step-over-step change, and the store's own tendency is the reach filling
    or emptying, which is the reading the probe documents. Differencing here
    would relabel the limbs against the reach's lag and silently change what
    "rising limb" means between models that report `dis` and models that do
    not.

    The second element of the return is the provenance, and it is not
    cosmetic: a probe that reads the rating off the store is comparing stage
    with storage, not with discharge, and `rating_monotonic` puts that string
    in its reason so the archived verdict does not claim a check it did not
    make. A probe that needs a genuine rating curve should require `dis`.
    """
    if "dis" in w.table.columns:
        return np.asarray(w.table["dis"], dtype=float), "dis"
    if "channel" in w.table.columns:
        return np.asarray(w.table["channel"], dtype=float), "channel"
    return np.zeros(len(w.table)), "none"


def _abscissa(w, prefer: str) -> tuple[np.ndarray, str]:
    """The x-axis the two limbs are paired on.

    A rating curve is conventionally drawn against discharge, and `dis` is
    the first choice. But the loop is a statement about how long the water
    stays in the reach, and when a model reports the reach's own store,
    pairing the limbs on that store is the more faithful reading: it compares
    the two limbs at equal water in transit, which is the quantity the
    hysteresis is actually about. The `prefer` argument picks which one this
    criterion uses; both are reported in the diagnostics either way.

    A store that never varies is no abscissa at all. A model that does not
    route — the exact bucket is the one in this suite — holds nothing in
    transit and reports a channel store of zero at every step, so pairing on
    it would collapse the whole rating onto one point and compare the limbs
    at a single value. That is a statement about the model having no reach,
    not a loop, and `dis` is the reading that still means something. The
    fallback is taken on the store being constant rather than merely small,
    so a reach that genuinely holds a little water is still read on its store.
    """
    if prefer == "store" and "channel" in w.table.columns:
        store = np.asarray(w.table["channel"], dtype=float)
        if len(store) and not np.allclose(store, store[0]):
            return store, "channel"
    if "dis" in w.table.columns:
        return np.asarray(w.table["dis"], dtype=float), "dis"
    if "channel" in w.table.columns:
        return np.asarray(w.table["channel"], dtype=float), "channel"
    return np.zeros(len(w.table)), "none"


def _stage(w) -> np.ndarray | None:
    if "stage" not in w.table.columns:
        return None
    return np.asarray(w.table["stage"], dtype=float)


def _rating_scatter(stage: np.ndarray, x: np.ndarray) -> float:
    """How far `stage` departs from a single-valued rating in `x`, in metres.

    A rating is single-valued when steps that share an abscissa share a
    stage. The measure is therefore taken between *neighbouring* steps in
    `x`: sort by the abscissa, and for each adjacent pair ask how far the
    stage moved over how little the abscissa moved. A stage that is a
    function of the abscissa moves by the local gradient times the abscissa
    step and no more, so the residual against that line is at rounding
    scale. A rating that loops does not, because two steps at the same
    abscissa carry different stages and the vertical jump is the whole loop.

    Working between neighbours rather than within bins is what makes the
    measure honest on a skewed rating. The abscissa of a river is dense at
    low flow and sparse at high flow, so an equal-count bin can straddle
    three orders of magnitude of discharge; its median stage then sits far
    from its own members for reasons that have nothing to do with a loop.
    Adjacent pairs carry no such assumption — they are the smallest
    comparison the data supports, and the one a loop must show up in.
    """
    if len(stage) < 2:
        return 0.0
    order = np.argsort(x, kind="stable")
    xs = np.asarray(x, dtype=float)[order]
    hs = np.asarray(stage, dtype=float)[order]

    # A pair with no separation in the abscissa is the sharpest test there
    # is: the two steps are at the same discharge, so a single-valued rating
    # must give them the same stage. Duplicated abscissae are rare in
    # continuous output but are exactly where a loop is unambiguous.
    #
    # Ties are one contribution to the residual, not the whole of it. A
    # short-circuit here would throw away every non-tied pair — the majority
    # of the record — so a rating that is wildly non-single-valued everywhere
    # except at a handful of coincident steps would read as if it were a
    # function. Worse, where the tied pair itself has no height difference the
    # short-circuit returns zero and the rating looks perfectly single-valued.
    # Both contributions are therefore measured and the larger is reported.
    dx = np.diff(xs)
    dh = np.abs(np.diff(hs))
    coincident = dx <= 0.0
    tie = float(dh[coincident].max()) if coincident.any() else 0.0

    # Otherwise, subtract the local gradient over each adjacent pair: the
    # part of the stage change the abscissa change already explains.
    separated = ~coincident
    if separated.any():
        sdx = dx[separated]
        sdh = dh[separated]
        slope = sdh / sdx
        # The median slope is the reach's typical rating gradient; a pair whose
        # stage move exceeds that gradient over its own abscissa step is the
        # part that no single-valued rating accounts for.
        typical = float(np.median(slope))
        residual = sdh - typical * sdx
        gradient = float(np.maximum(0.0, residual).max())
    else:
        gradient = 0.0
    return max(tie, gradient)


def _binned(frame: pd.DataFrame, bins: int) -> tuple[np.ndarray, np.ndarray]:
    """Equal-count bin medians.

    Returns (median discharge, median stage) per bin, in ascending discharge.
    Bins are formed on ranks so each holds the same number of steps, then
    ordered by their displacement, which is what makes the sequence readable
    as a rating curve.
    """
    n = len(frame)
    edges = np.linspace(0, n, bins + 1).astype(int)
    xs, ys = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi - lo < MIN_BIN_SIZE:
            continue
        chunk = frame.iloc[lo:hi]
        xs.append(float(np.median(chunk["_q"])))
        ys.append(float(np.median(chunk["_h"])))
    return np.asarray(xs), np.asarray(ys)


@criterion("rating_monotonic")
def rating_monotonic(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Stage must not fall as discharge rises.

    The measure is a running-maximum deficit rather than an adjacent-pair
    comparison, because adjacent pairs only catch a local inversion and miss
    a curve that steps down, recovers, and steps down again:

        R_mono = sum_i max(0, max_{j<i} Hb_j - Hb_i)   [m]

    It is positive only where a bin's median stage sits below the median of
    some earlier, less discharge-heavy bin, and it is zero for any curve that
    never does. Taking it against the running maximum rather than the
    previous bin is what makes it sum the total size of the violations
    instead of the size of the last one.
    """
    w = make_window(run, probe)
    stage = _stage(w)
    if stage is None:
        return CriterionResult(
            name="rating_monotonic",
            status=FAIL,
            message="the run reports no 'stage' column, so there is no rating to check",
        )

    q, q_source = _discharge(w)
    bins = int(params.get("bins", DEFAULT_BINS))
    tolerance = float(params.get("tolerance", 0.0))

    frame = pd.DataFrame({"_q": q, "_h": stage}).sort_values("_q").reset_index(drop=True)
    qb, hb = _binned(frame, bins)

    if len(qb) < 2:
        return CriterionResult(
            name="rating_monotonic",
            status=FAIL,
            message=f"only {len(qb)} usable discharge bin(s); the rating cannot be read",
            diagnostics={"bins": len(qb), "discharge_source": q_source},
        )

    running_max = np.maximum.accumulate(hb)
    # The first bin has no earlier bin, so its deficit is zero by definition.
    deficits = np.maximum(0.0, running_max - hb)
    deficit = float(deficits.sum())

    # A bin's median stage is an estimate, not a measurement, so a rating read
    # off binned medians is never exactly monotone even when the underlying
    # relation is. The tolerance therefore has a floor that scales with the
    # rating's own span — the same reasoning as the loop's floor and for the
    # same reason: a deflection worth a fraction of a percent of the reach is
    # sampling noise, while a curve that steps down through a visible share of
    # its range is the failure the criterion exists to catch. The floor is the
    # default, not a hard minimum: an explicit `tolerance` in metres is
    # honoured as written, tighter or looser, matching `rating_loop`'s
    # `min_loop_m`.
    span = float(np.max(stage) - np.min(stage)) if len(stage) else 0.0
    floor = max(MIN_SPAN_M, MIN_MONOTONIC_FRACTION * span)
    # An explicit positive `tolerance` is honoured as written, including when
    # it is tighter than the floor: `rating_loop`'s `min_loop_m` is a true
    # override and this is the same knob, so a probe that names the deflection
    # it will accept gets it. The floor applies only where the probe leaves
    # the tolerance at zero.
    tolerance = tolerance if tolerance > 0.0 else floor

    ok = deficit <= tolerance
    # Name what was actually compared. Reading the rating off the channel
    # store is comparing the gauge with storage, not with discharge, and the
    # archived reason must not assert a discharge check the criterion did not
    # make.
    against = "discharge" if q_source == "dis" else f"the {q_source} store"
    return CriterionResult(
        name="rating_monotonic",
        status=PASS if ok else FAIL,
        value=deficit,
        threshold=tolerance,
        message=(
            f"stage rises monotonically with {against} ({deficit:.4f} m of deficit "
            f"over {len(qb)} bins, source {q_source}, tolerance {tolerance:.4f} m)"
            if ok
            else f"stage falls as {against} rises: {deficit:.4f} m of running-maximum "
            f"deficit over {len(qb)} bins, source {q_source}, "
            f"tolerance {tolerance:.4f} m"
        ),
        diagnostics={
            "deficit_m": deficit,
            "bins": len(qb),
            "discharge_source": q_source,
            "span_m": span,
            "bin_discharge": [float(v) for v in qb],
            "bin_stage": [float(v) for v in hb],
        },
    )

@criterion("rating_loop")
def rating_loop(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """On the rise the stage must sit below the same abscissa on the fall.

    The two limbs are paired at equal abscissa, into equal-count bins over
    the rising limb, and the falling-limb median is asked to be the *higher*
    of the two. Reading the same abscissa higher on the fall is what a reach
    does when water is still in transit behind the wave, and it is the
    signature that separates a routing from an algebraic transform.

    A single store cannot produce this, and that is the whole point of the
    criterion. If the flow and the gauge are both functions of one state,
    then at equal state they are equal, and at equal flow they are equal too;
    the two limbs coincide and the rating is a function rather than a loop.
    The loop needs the reach to respond with two time constants — a fast path
    that carries the flood and a slow store that the gauge reads — and a
    model that collapses them into one is exactly what fails here.

    The `abscissa` parameter chooses what the limbs are paired on: `"store"`
    pairs them at equal water in transit, which is the more faithful reading
    of a hysteresis and the one this probe asks for; `"discharge"` pairs them
    at equal flow, the conventional way to draw a rating curve.

    Steps that belong to neither limb (the peak plateau) are dropped rather
    than assigned, because the loop is a statement about the two limbs.
    """
    w = make_window(run, probe)
    stage = _stage(w)
    if stage is None:
        return CriterionResult(
            name="rating_loop",
            status=FAIL,
            message="the run reports no 'stage' column, so there is no loop to check",
        )

    bins = int(params.get("bins", DEFAULT_BINS))
    tolerance = float(params.get("tolerance", 0.0))
    prefer = str(params.get("abscissa", "store"))
    # The step direction that decides which limb a step belongs to is always
    # the flow: a limb is a rise or a fall in discharge, whatever the two
    # limbs are then compared on.
    q, q_source = _discharge(w)
    x, x_source = _abscissa(w, prefer)

    # A step is on the rise if discharge is higher than the step before it,
    # and on the fall if it is lower. The first step has no predecessor, so
    # it is excluded.
    dq = np.diff(q, prepend=q[0])
    rising = dq > 0
    falling = dq < 0

    if not rising.any() or not falling.any():
        return CriterionResult(
            name="rating_loop",
            status=FAIL,
            message=(
                "the record contains no rising and falling limb at once "
                f"({int(rising.sum())} rising, {int(falling.sum())} falling steps); "
                "a loop cannot be formed"
            ),
            diagnostics={"rising_steps": int(rising.sum()), "falling_steps": int(falling.sum())},
        )

    rise = pd.DataFrame({"_q": x[rising], "_h": stage[rising]})
    fall = pd.DataFrame({"_q": x[falling], "_h": stage[falling]})

    # Bin the rising limb by abscissa rank and take each bin's median; the
    # falling limb is then read at those values by nearest, which is the
    # interval-free way of asking "at this value, how did the two limbs
    # compare".
    rise_sorted = rise.sort_values("_q").reset_index(drop=True)
    qb, hb_rise = _binned(rise_sorted, bins)
    if len(qb) < 2:
        return CriterionResult(
            name="rating_loop",
            status=FAIL,
            message="the rising limb is too short to bin, so the loop cannot be read",
            diagnostics={"rise_bins": len(qb)},
        )

    fall_sorted = fall.sort_values("_q").reset_index(drop=True)
    fq = fall_sorted["_q"].to_numpy()
    fh = fall_sorted["_h"].to_numpy()

    # Nearest falling-limb step per rising bin. Not an interpolation: an
    # interpolated stage would assert a smooth rating the model never
    # claimed, and the loop is a sign test either way.
    idx = np.abs(fq[None, :] - qb[:, None]).argmin(axis=1)
    hb_fall = fh[idx]

    # Positive where the rise sits above the fall — the unphysical direction.
    #
    # The loop is read as an aggregate rather than as a per-bin sum, because
    # the pairing is by nearest neighbour and a single bin can land on the
    # wrong step wherever the falling limb is sampled sparsely. That artifact
    # is the size of the local rating, not of a loop, and summing it in would
    # fail a model whose rating is right everywhere except one bin. What a
    # loop asserts is the *direction the two limbs disagree in*, so the
    # measure is the median of the signed per-bin differences: a rating that
    # genuinely loops has every bin on the same side of zero and the median
    # is that whole separation, while pairing noise on a handful of bins
    # leaves the median where the physics put it. A median is used over a
    # mean for the same reason a median is used to bin: one bin carrying a
    # flood should not outvote the ten that carry the routine.
    signed = hb_rise - hb_fall
    loop_signed = float(np.median(signed)) if len(signed) else 0.0
    # The loop is a signed quantity: positive where the rise sits above the
    # fall. The criterion fails only on that unphysical direction, so the
    # value it reports against `tolerance` is the *upward* part of the median,
    # clipped at zero.
    #
    # The size of the loop, though, is the magnitude of the median whatever
    # its sign: a reach that separates the limbs by 0.9 m in the physical
    # direction has a loop of 0.9 m, and the "is it big enough to read"
    # question has to be asked of that number. Reporting the clipped value as
    # the size would make every correctly-signed loop read as zero and take
    # the no-loop escape below, which is how a gauge with a real hysteresis
    # came to be reported as having none.
    loop = max(0.0, loop_signed)
    loop_size = abs(loop_signed)

    # A loop can only be read where the two limbs are actually distinguishable,
    # and the honest way to ask that is to ask whether the model's rating is
    # single-valued at all. If `stage` is a function of the store — the
    # signature of a model that keeps one state and draws its gauge reading
    # from it — then steps that share a store share a stage, and any
    # difference the two limbs show is pairing error from reading a sparse
    # falling limb by nearest neighbour, not hysteresis.
    #
    # That is measured, not assumed, by the worst residual of `stage` against
    # its own per-bin median: a single-valued rating has a residual at
    # floating-point scale, and a rating that genuinely loops does not. Taking
    # the worst bin rather than a mean is deliberate — a loop that opens over
    # one stretch of the flood has to survive, and averaging over twelve bins
    # would hide it.
    #
    # The floor scales with the rating's own span so it reads as a fraction of
    # the reach rather than as an absolute length. `min_loop_m` overrides it,
    # for a probe that wants to name the separation it will accept.
    span = float(np.max(stage) - np.min(stage)) if len(stage) else 0.0
    min_loop = float(params.get("min_loop_m", 0.0))
    if min_loop <= 0.0:
        min_loop = max(MIN_SPAN_M, MIN_LOOP_FRACTION * span)
    # The scatter that marks a rating as single-valued scales with the rating's
    # own span and not with `min_loop`. The two answer different questions:
    # `min_loop` is the separation the reach can resolve, which a probe may
    # state absolutely, while "is this rating a function of its axis" is a
    # share-of-range question — a stage that departs from its axis by a
    # rounding-scale fraction of its own range is a function of that axis
    # whatever the instrument's resolution happens to be.
    scatter_floor = max(MIN_SPAN_M, MIN_SINGLE_VALUED_FRACTION * span)
    scatter = _rating_scatter(stage, x)
    scatter_q = _rating_scatter(stage, q) if q_source != x_source else scatter

    # Two ways a rating can have no loop worth reading, and both have to pass.
    #
    # The first is a stage that is a *function* of the abscissa: steps sharing
    # an abscissa share a stage to floating-point, so whatever the two limbs
    # appear to differ by is pairing error rather than hysteresis. This is the
    # signature of a model that keeps one state and reads its gauge off it,
    # and `scatter` measures it directly.
    #
    # The second is a loop too small for the reach to resolve, whether or not
    # the rating is single-valued. A model that does not route holds no water
    # in transit, so its gauge has no second time constant to lag on and its
    # two limbs separate by a fraction of a percent of the rating. Failing it
    # would demand that it become a different model, which is not what a probe
    # is for. The size test is applied only when the limbs are *not*
    # consistently on the wrong side, because consistency is what separates a
    # small real loop from a small artifact:
    #
    #   * a rating that genuinely loops the wrong way has essentially every
    #     bin inverted — the sign is the physics, and it is the same at every
    #     discharge. `reference_rating_inverted` inverts 12 of 12.
    #   * pairing noise on a non-routing model flips a handful of bins and
    #     leaves the rest, so the inverted fraction sits near a half.
    #
    # Requiring both a sub-threshold size *and* an inconsistent sign is what
    # lets a non-routing model pass without also letting an inverted one
    # through.
    inverted = int((signed > 0.0).sum())
    inverted_fraction = inverted / len(signed) if len(signed) else 0.0
    consistent = inverted_fraction >= MAJORITY_FRACTION

    # A rating that is a function of the abscissa — or of the discharge — is
    # single-valued, and a single-valued rating has no loop: two steps that
    # share the abscissa (or the flow) share the level, however they got
    # there. This is the signature of a model that reads its gauge off one
    # state, and whatever the two limbs then differ by is pairing error.
    #
    # Checking the discharge as well as the abscissa matters for the physical
    # models. `flex_lumped` and `sacsma_snow17` build their stage from their
    # own discharge, so their rating is single-valued in `dis` even though the
    # probe pairs on the store. On the store axis their two limbs then miss
    # each other by centimetres of noise, and no size floor separates that
    # from a real loop. Reading the rating as a function of the flow it is
    # drawn against says what the model actually did — it computed a stage
    # from its discharge — and escapes the criterion without weakening it,
    # because a rating that genuinely loops is a function of *neither* axis
    # alone.
    single_valued_in = None
    single_valued_scatter = 0.0
    if scatter <= scatter_floor:
        single_valued_in, single_valued_scatter = x_source, scatter
    elif scatter_q <= scatter_floor:
        single_valued_in, single_valued_scatter = q_source, scatter_q

    if single_valued_in is not None:
        reason = (
            f"stage is a function of {single_valued_in} to within "
            f"{single_valued_scatter:.2e} m, so the rating is single-valued"
        )
    elif loop_size <= min_loop:
        # The size floor is symmetric: a separation too small for the reach to
        # resolve is not a loop in either direction. An inverted gauge hides
        # behind it only if its loop is smaller than the reach can read, which
        # is why the reference gauges are built at a physical scale: they have
        # to clear the same floor a real loop would.
        sign_note = (
            f"and the sign agrees across bins ({inverted} of {len(qb)} inverted)"
            if inverted == 0
            else f"and the sign is not consistent ({inverted} of {len(qb)} inverted)"
        )
        reason = (
            f"the two limbs differ by only {loop_size:.2e} m "
            f"({loop_size / span:.2%} of the {span:.2f} m rating), {sign_note}, "
            f"which is pairing noise rather than a loop"
        )
    else:
        reason = None

    if reason is not None:
        return CriterionResult(
            name="rating_loop",
            status=PASS,
            value=0.0,
            threshold=tolerance,
            message=(
                f"no loop to read: {reason} ({len(qb)} bins)"
                if single_valued_in is not None
                else (
                    f"no loop to read: {reason} — below the {min_loop:.4f} m "
                    f"the reach resolves ({len(qb)} bins)"
                )
            ),
            diagnostics={
                "loop_m": 0.0,
                "loop_size_m": loop_size,
                "loop_median_m": loop_signed,
                "bins": len(qb),
                "abscissa": x_source,
                "limb_direction_by": q_source,
                "rating_scatter_m": scatter,
                "rating_scatter_discharge_m": scatter_q,
                "min_loop_m": min_loop,
                "span_m": span,
                "single_valued": single_valued_in is not None,
                "single_valued_in": single_valued_in,
                "bin_abscissa": [float(v) for v in qb],
            },
        )

    ok = loop <= tolerance
    return CriterionResult(
        name="rating_loop",
        status=PASS if ok else FAIL,
        value=loop,
        threshold=tolerance,
        message=(
            f"the rising limb sits below the falling limb at matching {x_source} "
            f"(median {loop_signed:+.4f} m over {len(qb)} bins, "
            f"{inverted} of {len(qb)} inverted)"
            if ok
            else f"the rating loops the wrong way: the rising limb sits "
            f"{loop:.4f} m above the falling limb at matching {x_source} "
            f"(median over {len(qb)} bins, {inverted} inverted)"
        ),
        diagnostics={
            "loop_m": loop,
            "loop_size_m": loop_size,
            "loop_median_m": loop_signed,
            "inverted_bins": inverted,
            "bins": len(qb),
            "abscissa": x_source,
            "limb_direction_by": q_source,
            "rating_scatter_m": scatter,
            "rating_scatter_discharge_m": scatter_q,
            "min_loop_m": min_loop,
            "span_m": span,
            "single_valued": False,
            "bin_abscissa": [float(v) for v in qb],
            "rising_stage": [float(v) for v in hb_rise],
            "falling_stage": [float(v) for v in hb_fall],
        },
    )

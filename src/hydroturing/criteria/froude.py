"""Froude regime: a mild-sloped reach cannot carry supercritical flow.

Momentum is the one conservation law the existing suite reads only through
its shadow. `momentum/routing-conservation` bounds the water a reach holds
and `momentum/stage-discharge-monotonic` asks whether the gauge moves with
the flow, but neither asks whether the *pair* the model reports is one the
reach could physically carry. A stage and a discharge are two readings of
the same cross-section, and open-channel flow puts a hard constraint on how
fast that cross-section can deliver water:

    Fr = v / sqrt(g * D) <= 1

with `v` the section-averaged velocity, `D` the hydraulic depth and `g`
gravity. Fr above one is not a small error in a budget — it is a different
regime. Supercritical flow cannot be sustained on a mild slope: the reach
would have to be steep enough, or short enough, that the wave outruns the
water it is made of. A model that reports it is not solving anything
resembling momentum, however well its mass closes.

For a wide rectangular reach of width `w` carrying `Q` at depth `d` the
velocity is `Q / (w * d)` and the hydraulic depth is `d`, so

    Fr = Q / (w * d**1.5 * sqrt(g))

which is the form scored here. Every quantity on the right is something the
model has already declared: `Q` from `dis`, or from `mrro` over the
catchment area when the model reports only a depth rate; `d` from the `stage`
diagnostic; `w` from the case's static attributes, so all models are judged
in the same channel rather than in one of their own choosing.

`stage` is a water-surface elevation on the case's fixed vertical datum, and
the flow depth is `stage - bed_elevation_m`: the case declares the datum, the
probe requires it, and `depth_series` does the subtraction. A fixed datum
cancels from a difference but not from a level, so a criterion that forms a
depth — or any ratio or power of the level, which `d**(3/2)` is — may not infer
one. That is why the bed is in `requires.static` rather than assumed: a model
that never read the datum is judged against a number it never saw, and the
suite's rule for that is N/A (INCOMPATIBLE), not a conservation violation. The
helper raises `CriterionIncompatibleError` when the reported column is a depth
where the contract asks for a level, and this criterion lets that exception
reach the harness rather than turning it into a local failure: the
classification is the harness's to make, and a seed it records as incompatible
must not arrive in the archive as a violation. A level reported against some
other zero is not merely a deep reach — an offset only ever adds to the depth,
and the depth enters at the three-halves power, so an offset is the cheapest
possible way to look subcritical.

Froude is a hard constraint, so the tolerance is absolute and small: the
5% margin is for the rectangular-cross-section approximation and for the
rounding in a stage that an adapter derived from a normal-depth relation,
not for a model being a little wrong about its own hydraulics.

The measure is nevertheless a **frequency**, not a worst case. Froude is
computed per step and every step is a separate claim about the reach, so
what is scored is the share of scored steps on which the reach went
supercritical, with a default limit of zero. Reporting the frequency keeps
the failure legible — a model that spikes once at a numerical edge looks
different from one that is supercritical through every flood — and leaves
the limit a probe parameter for a future case where a small share is the
honest reading.

A step is left out when nothing is moving through it, and for no other reason —
and "nothing" is read from every flow the model reports, not from the one this
criterion prefers. `dis` is preferred, so a model reporting both columns could
otherwise drop a step from the scored set with a zero in that column while its
`mrro` still carried the water, which under the contract is the same outflow
(`mrro` is the total runoff and `channel` holds what the routing has not
released yet). Where the preferred column is silent and the other is not, the
step is scored on the other. Only a step both readings call still is excused,
and `min_scored_fraction` bounds how much of a record that can be.
The depth carries a floor of one centimetre, but the floor is a floor on the
*division* rather than a dryness test: the depth is clamped to it, so a receding
flow still has a number to divide by, and a step with water in it is scored at
the clamped depth however shallow the model says it is. Deciding dryness on the
depth itself would be backwards here. At a given discharge the shallowest step is
the *fastest* one in the record — `Fr = Q / (w sqrt(g) d**1.5)` grows without
bound as `d` falls — so a depth test drops exactly the steps this criterion exists
to catch, and the direction of the error is the wrong way round: a model carrying
its flood-peak discharge at a millimetre of depth would have those peaks excused
while every step it reported honestly was still scored, so making its worst steps
*more* wrong would move it from a failure to a pass.

Leaving a step out is for steps the *reach* makes unscorable, not for values the
criterion cannot read. A non-finite stage or discharge is refused rather than
dropped: `NaN` compares false against every bound, so the mask would quietly
score one step fewer and say nothing about it — an invalid input leaving the
assessment by the same door an honest dry step uses.

The depth has a ceiling, and a probe gives it a value: `max_depth_m`.
It is a section sanity bound and not the defence against a datum mismatch — that
job belongs to the declared bed and the subtraction — but a depth deeper than
the section could hold is not a depth in it, and a model that ignored the
declared datum and reported a level is exactly what it catches once the offset
is large. Unlike the floor this one *refuses* the step rather than clamping it:
clamping puts the reading back inside the range the bound exists to exclude, and
an offset only ever adds to the depth. A probe that knows its reach sets the
ceiling; without one the criterion takes the reading as it stands, because how
deep is too deep for a section is a statement about the case rather than about
`stage`.

What the criterion cannot see is stated in the probe's README: it reads the
pair, not either reading against the truth. A model that derives its stage
from its own discharge through *any* normal-depth relation is consistent by
construction — `Fr` collapses to `sqrt(S)/n * d**(1/6) / sqrt(g)`, which
depends on the depth and not on the flow — so a model whose discharge is
wrong but internally consistent stays subcritical whatever the error. The
absolute readings belong to the mass probes; this one asks only whether the
two numbers the model reports could both be true of the declared section.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL,
    PASS,
    CriterionIncompatibleError,
    CriterionResult,
    criterion,
    depth_series,
    make_window,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

# Gravity, m/s2.
GRAVITY = 9.81
SECONDS_PER_DAY = 86400.0

# The diagnostic read as the water surface, and the flux read as the flow.
DEFAULT_STAGE = "stage"
DEFAULT_DISCHARGE = "dis"
# The runoff flux, in mm/day, used when the model reports no `dis`.
DEFAULT_RUNOFF = "mrro"
# Static keys for the reach geometry. All of them are supplied by the case, so
# every model is judged in the same channel.
DEFAULT_WIDTH_KEY = "width_m"
DEFAULT_AREA_KEY = "area_km2"
# A depth floor of one centimetre, applied to the division and not to the mask:
# enough that a receding flow does not divide by zero, small enough that it never
# rescues a real cross-section. A step is excused for carrying no flow, never for
# being shallow.
DEFAULT_MIN_DEPTH_M = 0.01
# No ceiling unless a probe sets one. A depth deeper than the declared section
# could hold is not a depth in it, but how deep is too deep depends on the reach
# the probe declared — which is why it is the probe's to give rather than a
# constant here. It is a section sanity bound: the datum itself is fixed by the
# contract and subtracted, not guarded by a value test.
DEFAULT_MAX_DEPTH_M = None
# Fr may exceed one by this much. It covers the rectangular approximation and
# the rounding in a normal-depth stage, not a wrong velocity.
DEFAULT_TOLERANCE = 0.05
# The share of scored steps that may be supercritical. Zero by default: the
# constraint is per step, and a single supercritical step on a mild reach is
# already the finding.
DEFAULT_MAX_EXCEED_FRACTION = 0.0
# Below this share of the record being scored, the case is degenerate for this
# criterion: a reach that is dry or flat throughout has no regime to read.
DEFAULT_MIN_SCORED_FRACTION = 0.05


def _depth(depth: np.ndarray, min_depth_m: float,
           max_depth_m: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """The depth the division uses, and the mask of depths this section holds.

    The depth arrives already formed: `depth_series` subtracts the case's
    declared `bed_elevation_m` from the reported `stage`, so nothing here has to
    know about the datum and nothing here may assume one.

    Two bounds doing two different jobs. `min_depth_m` is a floor on the
    division and not a dryness test: the returned depth is clamped to it so a
    receding flow has a number to divide by, and the caller decides dryness from
    the flow instead. `max_depth_m` is the mask, and it refuses rather than
    clamps: a step whose depth is deeper than the declared section could hold is
    not a reading of that section at all, and clamping it back into range would
    re-admit the reading the bound exists to exclude, because an offset only ever
    adds to the depth.
    """
    within = np.ones(depth.shape, dtype=bool)
    if max_depth_m is not None:
        within = depth <= max_depth_m
    return np.maximum(depth, min_depth_m), within


def _discharge(run: RunResult, w, params: dict) -> tuple[np.ndarray | None, str | None, str | None]:
    """The flow through the reach, in m3/s, where it came from, and why not.

    `dis` is the reach's own discharge and is preferred whenever the model
    reports it. A model that reports only `mrro` is still scorable: a depth
    rate over the catchment area is a volume rate, and the conversion is the
    same one the suite's other flux criteria use, so the two sources are
    interchangeable rather than merely comparable.

    Returns `(q, source, problem)`. `problem` is `None` when the flow could be
    read, and a sentence naming what is missing when it could not — the caller
    turns that into a failure rather than a crash, because a criterion that
    cannot make its statement must not be reported as satisfied.
    """
    table = w.table
    if DEFAULT_DISCHARGE in table.columns:
        return np.asarray(table[DEFAULT_DISCHARGE], dtype=float), "dis", None

    runoff_var = str(params.get("runoff", DEFAULT_RUNOFF))
    if runoff_var not in table.columns:
        return None, None, (
            f"the run reports neither '{DEFAULT_DISCHARGE}' nor '{runoff_var}', "
            "so there is no flow to read the water surface against"
        )
    static = run.case.static
    area_key = str(params.get("area_key", DEFAULT_AREA_KEY))
    if area_key not in static:
        return None, None, (
            f"the run reports '{runoff_var}' but the case declares no "
            f"'{area_key}', so a depth rate cannot be turned into a discharge"
        )
    area_km2 = float(static[area_key])
    q = _from_runoff(np.asarray(table[runoff_var], dtype=float), area_km2)
    return q, runoff_var, None


def _from_runoff(rate: np.ndarray, area_km2: float) -> np.ndarray:
    """A depth rate in mm/day over the catchment, as a volume rate in m3/s."""
    return np.maximum(rate, 0.0) * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY


def _second_reading(run: RunResult, w, params: dict,
                    source: str | None) -> tuple[np.ndarray | None, str | None]:
    """The other reading of the same outflow, where the model reports both.

    A step is excused when nothing is moving through the reach, and the model is
    the one reporting whether anything is. `dis` is preferred, so a model that
    reports both can take a step out of the scored set by writing a zero in that
    one column while its `mrro` still carries the water — and the steps worth
    zeroing are exactly the ones a shallow gauge fails on. The two columns are
    two readings of the same outflow under the contract (`AGENTS.md`: `mrro` is
    total runoff and `channel` holds what the routing has not released yet), so
    where `dis` is zero and the runoff is not, the runoff is what the model says
    moved and the step is scored on it.

    Returns `(q, source)` for the reading not already in use, or `(None, None)`
    when the model reports only one of them, or the case declares no area to
    convert a depth rate with.
    """
    if source != "dis":
        return None, None
    runoff_var = str(params.get("runoff", DEFAULT_RUNOFF))
    if runoff_var not in w.table.columns:
        return None, None
    area_key = str(params.get("area_key", DEFAULT_AREA_KEY))
    static = run.case.static
    if area_key not in static:
        return None, None
    rate = np.asarray(w.table[runoff_var], dtype=float)
    return _from_runoff(rate, float(static[area_key])), runoff_var


@criterion("froude_subcritical")
def froude_subcritical(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Flow in a mild-sloped reach must stay subcritical.

    Reports the share of scored steps on which Fr exceeded `1 + tolerance`,
    and fails when that share passes `max_exceed_fraction`. A step is scored
    when water is moving through it — a non-zero reading in either flow column
    the model reports — and its depth is one the declared section could hold, so
    a step is excused for being dry or for being deeper than the reach and never
    for being shallow. A record with too few scored steps is degenerate rather
    than passing.
    """
    stage_var = str(params.get("stage", DEFAULT_STAGE))
    width_key = str(params.get("width_key", DEFAULT_WIDTH_KEY))
    min_depth_m = float(params.get("min_depth_m", DEFAULT_MIN_DEPTH_M))
    max_depth_m = params.get("max_depth_m", DEFAULT_MAX_DEPTH_M)
    max_depth_m = None if max_depth_m is None else float(max_depth_m)
    tolerance = float(params.get("tolerance", DEFAULT_TOLERANCE))
    max_fraction = float(params.get("max_exceed_fraction", DEFAULT_MAX_EXCEED_FRACTION))
    min_scored_fraction = float(params.get("min_scored_fraction", DEFAULT_MIN_SCORED_FRACTION))

    w = make_window(run, probe)
    # A missing column is a failure, not a crash and not a pass: the criterion
    # cannot make its statement, and one that cannot be scored must not be
    # reported as satisfied. The harness normally short-circuits a model that
    # does not report the diagnostic before any criterion runs, so this is a
    # second line of defence rather than the first.
    if stage_var not in w.table.columns:
        return CriterionResult(
            name="froude_subcritical",
            status=FAIL,
            message=(
                f"the run reports no '{stage_var}' column, so there is no water "
                "surface to read against the flow"
            ),
        )

    static = run.case.static
    # The case declares the section and the verdict is a statement about a pair
    # of readings in it, so a case that declares no width is a statement about
    # the case rather than about the model: refuse, as the `area_key` branch
    # below does. Falling back to a width would score every model in a channel
    # no case ever declared, and the message would name a number that exists
    # nowhere in the run.
    if "width_m" not in params and width_key not in static:
        return CriterionResult(
            name="froude_subcritical",
            status=FAIL,
            message=(
                f"the case declares no '{width_key}', so there is no section to "
                "read a stage against; the criterion cannot make its statement"
            ),
            diagnostics={"width_key": width_key},
        )
    width_m = float(params.get("width_m", static.get(width_key)))
    if width_m <= 0.0:
        return CriterionResult(
            name="froude_subcritical",
            status=FAIL,
            message=(
                f"the reach width is {width_m!r} m, so the section the case "
                "declares has no width to carry a flow through"
            ),
        )

    # The depth is `stage - bed_elevation_m`, formed by the shared helper: the
    # datum is the case's to declare and the probe's to require, and this
    # criterion may not infer one. `depth_series` raises rather than defaulting
    # the bed to zero — a default would make a model that never read the datum
    # score exactly as one that did — so a case that declares none is refused
    # here, the way a case with no width is refused above.
    try:
        depth_all = depth_series(run)
    except CriterionIncompatibleError:
        # A convention mismatch, not a violation. `depth_series` raises this when
        # the reported column is a depth where the contract asks for a level on
        # the case's datum, and the classification belongs to the harness: it
        # turns the exception into N/A (INCOMPATIBLE) for that seed, which is
        # what the suite says about a model judged against a number it never
        # read. Returning a failed criterion here instead would put a
        # conservation violation in the archive for a model that reported a
        # different quantity, and would take the seed out of the harness's own
        # accounting of how many seeds were scored.
        raise
    except ValueError as exc:
        return CriterionResult(
            name="froude_subcritical",
            status=FAIL,
            message=(
                "the case declares no usable bed elevation, so the flow depth "
                f"cannot be formed from the reported stage ({exc})"
            ),
        )

    q, q_source, problem = _discharge(run, w, params)
    if problem is not None:
        return CriterionResult(name="froude_subcritical", status=FAIL, message=problem)
    q_second, second_source = _second_reading(run, w, params, q_source)
    n = len(w.table)

    # `depth_series` returns the whole record; the criterion scores the window
    # `make_window` took, which is the same table with the spinup dropped from
    # the front, so the two line up on the same slice.
    depth_window = depth_all[run.case.spinup_steps:]

    # The flow this criterion divides by, before anything is masked on it. The
    # magnitude rather than the signed value, because a discharge is a magnitude
    # here: a step that reports a negative flow is water moving, and reading the
    # sign as "nothing to score" would let a model take its worst steps out of
    # the sample by flipping them — the same move as making them shallow.
    #
    # Where the preferred column reports exactly no flow and the model's other
    # reading of the same outflow does, the step is scored on that reading. Only
    # a step both readings report as still is excused. Without this a model that
    # reports both columns picks which steps are scored: a zero in `dis` costs it
    # nothing that this criterion can see, and the steps worth zeroing are the
    # ones a shallow gauge fails on. Measured on `reference_shallow_rating`, the
    # probe's own must-fail: zeroing `dis` on its supercritical steps left 144 of
    # 1460 steps scored on the first gate seed and turned the verdict into a
    # pass, with the runoff column still carrying every drop.
    flow = np.abs(q)
    n_from_second = 0
    if q_second is not None:
        silent = flow == 0.0
        flow = np.where(silent, np.abs(q_second), flow)
        n_from_second = int((silent & (flow > 0.0)).sum())

    # A non-finite value is not a step to skip. The mask below keeps only what
    # it can compare — `NaN` is neither above zero nor below the ceiling, so both
    # halves of it are false — and a single bad value in the peak of a flood would
    # leave the scored set quietly: one
    # step fewer, a slightly smaller share, and a pass. That is an invalid input
    # leaving the assessment by the door an honest dry step uses, and the two
    # series being scored are the model's own readings, so the pair is refused
    # rather than thinned — the same reading `radiative_identity` takes of a
    # non-finite surface temperature. Read on the window the criterion scores,
    # so a spinup value it never looks at cannot fail it, and on the flow as it
    # will be divided by, so a bad value in the reading a silent step falls back
    # to is refused rather than read as stillness.
    non_finite = ~np.isfinite(depth_window) | ~np.isfinite(flow)
    if non_finite.any():
        n_bad = int(non_finite.sum())
        return CriterionResult(
            name="froude_subcritical",
            status=FAIL,
            message=(
                f"non-finite stage or flow on {n_bad} of {n} scored steps, so "
                "the pair this criterion reads is not defined on the steps it "
                "scores and no share of them can be reported"
            ),
            diagnostics={"non_finite_steps": n_bad, "discharge_source": q_source},
        )

    depth, within_section = _depth(depth_window, min_depth_m, max_depth_m)
    # Dryness is read from the flow and not from the depth. A step is excused when
    # nothing is moving through it; a *shallow* depth is not dryness but the
    # fastest flow the record can hold, and the depth has already been clamped to
    # the floor above so that every flowing step has a number to divide by.
    scored = within_section & (flow > 0.0)
    n_scored = int(scored.sum())
    # Steps the ceiling refused: a depth the section could not hold. Counted
    # apart from the dry ones so that the message below can say which of the two
    # bounds emptied the record.
    n_too_deep = int((depth > max_depth_m).sum()) if max_depth_m is not None else 0

    # `n_scored == 0` is checked rather than left to the fraction: a probe that
    # sets `min_scored_fraction` to zero to turn the floor off makes `0/x < 0`
    # false, and the shares below would then divide by zero and raise instead of
    # returning a result — an ERROR row for what is really a dry record.
    if n == 0 or n_scored == 0 or n_scored / n < min_scored_fraction:
        refused = (
            f"; {n_too_deep} steps report a depth above the {max_depth_m:g} m "
            "ceiling, which is not a depth in the section the case declares"
            if n_too_deep
            else ""
        )
        return CriterionResult(
            name="froude_subcritical",
            status=FAIL,
            message=(
                "the reach is dry or flat through the record "
                f"({n_scored} of {n} steps carry flow within the section, "
                f"needs {min_scored_fraction:.0%}){refused}; there is no flow "
                "regime to read and the case is degenerate for this criterion"
            ),
            diagnostics={
                "scored_steps": n_scored,
                "total_steps": n,
                "too_deep_steps": n_too_deep,
                "max_depth_m": max_depth_m,
            },
        )

    # Fr = |Q| / (w * d^1.5 * sqrt(g)). Every factor is something the model
    # declared, so a model cannot escape by reporting a stage and a discharge
    # that disagree: the disagreement is exactly what is being scored. The
    # magnitude of the flow, because the number is a ratio of speeds and a
    # signed numerator would come out negative, which is below the limit at
    # every step — a reading that cannot fail rather than one that passes.
    fr = np.zeros(n)
    fr[scored] = flow[scored] / (width_m * depth[scored] ** 1.5 * GRAVITY**0.5)

    limit = 1.0 + tolerance
    exceed = np.zeros(n, dtype=bool)
    exceed[scored] = fr[scored] > limit
    n_exceed = int(exceed.sum())
    fraction = float(n_exceed / n_scored)
    worst = float(fr[scored].max())
    median_fr = float(np.median(fr[scored]))

    ok = fraction <= max_fraction
    # The share is reported in both branches, and the "every scored step" claim
    # is only made when it is true: `max_exceed_fraction` is a knob a future case
    # can raise, and a message archived in the result column must not state the
    # opposite of what was measured.
    if n_exceed == 0:
        detail = (
            f"the reach stays subcritical on every scored step "
            f"({n_scored} steps, worst Fr {worst:.3f}, median {median_fr:.3f}, "
            f"limit {limit:.2f})"
        )
    elif ok:
        detail = (
            f"the reach goes supercritical on {fraction:.1%} of scored steps "
            f"({n_exceed} of {n_scored}, worst Fr {worst:.3f}, median "
            f"{median_fr:.3f}, limit {limit:.2f}), within the "
            f"{max_fraction:.1%} this probe allows"
        )
    else:
        detail = (
            f"the reach goes supercritical on {fraction:.1%} of scored steps "
            f"({n_exceed} of {n_scored}, worst Fr {worst:.3f}, median "
            f"{median_fr:.3f}, limit {limit:.2f}, allowed {max_fraction:.1%}); a "
            "mild-sloped reach cannot carry Fr > 1, so the stage and the "
            "discharge it reports are not two readings of the same cross-section"
        )
    # Steps scored on the model's other reading of the same outflow are named in
    # the message as well as counted, because the archive keeps the sentence: a
    # verdict that rests on the runoff column where the discharge column read
    # zero should say so rather than leave it to a diagnostic.
    if n_from_second:
        detail += (
            f"; {n_from_second} of the scored steps report no '{DEFAULT_DISCHARGE}' "
            f"and are read from '{second_source}', which is the same outflow"
        )
    return CriterionResult(
        name="froude_subcritical",
        status=PASS if ok else FAIL,
        value=fraction,
        threshold=max_fraction,
        message=detail,
        diagnostics={
            "discharge_source": q_source,
            "second_reading": second_source,
            "steps_read_from_second_reading": n_from_second,
            "width_m": width_m,
            "scored_steps": n_scored,
            "exceeding_steps": n_exceed,
            "exceed_fraction": fraction,
            "max_froude": worst,
            "median_froude": median_fr,
            "froude_limit": limit,
            "min_depth_m": min_depth_m,
            "max_depth_m": max_depth_m,
            "too_deep_steps": n_too_deep,
        },
    )

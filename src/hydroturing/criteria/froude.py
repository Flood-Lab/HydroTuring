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

`stage` is a depth above the reach bed by contract, which is what makes `d`
readable off it with no datum to subtract. A model that reports an absolute
level, or a level above any other zero, has reported a different quantity, and
the criterion has no way to convert it: the one thing worse than refusing such
a reading is scoring it, because an offset only ever adds to the depth and a
depth enters at the three-halves power, so an offset is the cheapest possible
way to look subcritical.

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

Steps where the reach is effectively dry are not scored. The depth carries a
floor of one centimetre so that a receding flow does not divide by zero, but
a step that only clears the floor because of the floor says nothing about
velocity, and scoring it would let a dry reach fail on arithmetic.

The floor has a mirror, and a probe can give it a value: `max_depth_m`.
A depth deeper than the section could be is not a depth in it, and scoring one
is how a datum offset passes — a stage reported from some other zero is not a
depth at all, and since only the depth is read, the offset makes every step
look slower than it is. A probe that knows its reach sets the ceiling; without
one the criterion takes the reading as it stands, because how deep is too deep
for a section is a statement about the case rather than about `stage`.

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

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
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
# A depth floor of one centimetre: enough that a receding flow does not divide
# by zero, small enough that it never rescues a real cross-section. Steps that
# only clear it because of it are not scored.
DEFAULT_MIN_DEPTH_M = 0.01
# No ceiling unless a probe sets one. A reading deep enough to be an elevation
# rather than a depth is one the criterion cannot use, but how deep is too deep
# depends on the reach the probe declared — which is why it is the probe's to
# give rather than a constant here.
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


def _depth(stage: np.ndarray, min_depth_m: float,
           max_depth_m: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """The reported depth, and the mask of steps that are depths in this section.

    `stage` is a depth above the reach bed by contract, so it is used as it
    stands: there is no datum to subtract, and a reading taken from another
    zero is not convertible into one.

    The mask is the point of both bounds. A step whose depth only clears
    `min_depth_m` because the floor lifted it is a dry step wearing a number,
    and the velocity computed on it is arithmetic, not hydraulics. A step whose
    depth passes `max_depth_m` is not a depth in the declared section at all —
    it is what a stage reported from another zero looks like once it is read as
    one — and scoring it would let exactly that pass, comfortably, because an
    offset only ever adds to the depth.
    """
    deep_enough = stage > min_depth_m
    if max_depth_m is not None:
        deep_enough &= stage <= max_depth_m
    return np.maximum(stage, min_depth_m), deep_enough


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
    rate = np.asarray(table[runoff_var], dtype=float)
    q = np.maximum(rate, 0.0) * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY
    return q, runoff_var, None


@criterion("froude_subcritical")
def froude_subcritical(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Flow in a mild-sloped reach must stay subcritical.

    Reports the share of scored steps on which Fr exceeded `1 + tolerance`,
    and fails when that share passes `max_exceed_fraction`. A step is scored
    when it reports a depth the declared section could hold — above the floor
    and, when the probe sets `max_depth_m`, below it — and a record with too
    few such steps is degenerate rather than passing.
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

    stage = np.asarray(w.table[stage_var], dtype=float)
    q, q_source, problem = _discharge(run, w, params)
    if problem is not None:
        return CriterionResult(name="froude_subcritical", status=FAIL, message=problem)
    n = len(stage)

    depth, deep_enough = _depth(stage, min_depth_m, max_depth_m)
    scored = deep_enough & (q > 0.0)
    n_scored = int(scored.sum())
    # Steps the ceiling refused: a depth the section could not hold, which is
    # what a stage reported from another zero looks like once it is read as a
    # depth. Counted apart from the dry ones so that the message below can say
    # which of the two bounds emptied the record.
    n_too_deep = int((stage > max_depth_m).sum()) if max_depth_m is not None else 0

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
                f"({n_scored} of {n} steps deep enough to score, "
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

    # Fr = Q / (w * d^1.5 * sqrt(g)). Every factor is something the model
    # declared, so a model cannot escape by reporting a stage and a discharge
    # that disagree: the disagreement is exactly what is being scored.
    fr = np.zeros(n)
    fr[scored] = q[scored] / (width_m * depth[scored] ** 1.5 * GRAVITY**0.5)

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
    return CriterionResult(
        name="froude_subcritical",
        status=PASS if ok else FAIL,
        value=fraction,
        threshold=max_fraction,
        message=detail,
        diagnostics={
            "discharge_source": q_source,
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

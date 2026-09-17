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
DEFAULT_BED_KEY = "bed_elevation_m"
DEFAULT_AREA_KEY = "area_km2"
# Fallbacks, used only when the case hands over no geometry.
FALLBACK_WIDTH_M = 18.0
FALLBACK_BED_M = 0.0
# A depth floor of one centimetre: enough that a receding flow does not divide
# by zero, small enough that it never rescues a real cross-section. Steps that
# only clear it because of it are not scored.
DEFAULT_MIN_DEPTH_M = 0.01
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


def _depth(stage: np.ndarray, bed_m: float, min_depth_m: float) -> tuple[np.ndarray, np.ndarray]:
    """Water depth above the bed, and the mask of steps deep enough to score.

    `stage` is whatever the model declared as its water surface. When a model
    reports an absolute level rather than a depth, the case's bed elevation is
    subtracted; the default bed is zero, which is the honest reading for a
    stage derived from a normal-depth relation, since that relation already
    returns a depth above the bed.

    The second mask is the point of the floor: a step whose depth only clears
    `min_depth_m` because the floor lifted it is a dry step wearing a number,
    and the velocity computed on it is arithmetic, not hydraulics.
    """
    above = stage - bed_m
    return np.maximum(above, min_depth_m), above > min_depth_m


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
    and fails when that share passes `max_exceed_fraction`.
    """
    stage_var = str(params.get("stage", DEFAULT_STAGE))
    width_key = str(params.get("width_key", DEFAULT_WIDTH_KEY))
    bed_key = str(params.get("bed_key", DEFAULT_BED_KEY))
    min_depth_m = float(params.get("min_depth_m", DEFAULT_MIN_DEPTH_M))
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
    width_m = float(params.get("width_m", static.get(width_key, FALLBACK_WIDTH_M)))
    bed_m = float(params.get("bed_elevation_m", static.get(bed_key, FALLBACK_BED_M)))
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

    depth, deep_enough = _depth(stage, bed_m, min_depth_m)
    scored = deep_enough & (q > 0.0)
    n_scored = int(scored.sum())

    if n == 0 or n_scored / n < min_scored_fraction:
        return CriterionResult(
            name="froude_subcritical",
            status=FAIL,
            message=(
                "the reach is dry or flat through the record "
                f"({n_scored} of {n} steps deep enough to score, "
                f"needs {min_scored_fraction:.0%}); there is no flow regime "
                "to read and the case is degenerate for this criterion"
            ),
            diagnostics={"scored_steps": n_scored, "total_steps": n},
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
    return CriterionResult(
        name="froude_subcritical",
        status=PASS if ok else FAIL,
        value=fraction,
        threshold=max_fraction,
        message=(
            f"the reach stays subcritical on every scored step "
            f"({n_scored} steps, worst Fr {worst:.3f}, median {median_fr:.3f}, "
            f"limit {limit:.2f})"
            if ok
            else (
                f"the reach goes supercritical on {fraction:.1%} of scored steps "
                f"({n_exceed} of {n_scored}, worst Fr {worst:.3f}, limit "
                f"{limit:.2f}); a mild-sloped reach cannot carry Fr > 1, so the "
                f"stage and the discharge it reports are not two readings of "
                f"the same cross-section"
            )
        ),
        diagnostics={
            "discharge_source": q_source,
            "width_m": width_m,
            "bed_elevation_m": bed_m,
            "scored_steps": n_scored,
            "exceeding_steps": n_exceed,
            "exceed_fraction": fraction,
            "max_froude": worst,
            "median_froude": median_fr,
            "froude_limit": limit,
            "min_depth_m": min_depth_m,
        },
    )

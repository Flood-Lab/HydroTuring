# momentum/froude-regime

A reach is a cross-section, and a stage and a discharge are two readings of
it. This probe asks whether the pair a model reports is one that section
could carry: on a mild slope, **flow must stay subcritical**.

    Fr = v / sqrt(g * D) = Q / (w * d**1.5 * sqrt(g)) <= 1

Every quantity on the right is something the model has already declared.
`Q` comes from `dis`, or from `mrro` over the catchment area when the model
reports only a depth rate; `d` comes from the `stage` diagnostic; `w` is the
width of the reach the case declares, so all models are judged in the same
channel rather than in one of their own choosing. Nothing is fitted and
nothing is calibrated — the section is given, and the question is whether
the two numbers agree about how deep the flow in it is.

Fr above one is not a small error in a budget. It is a different regime:
the water outruns the wave that carries it. Supercritical flow cannot be
sustained on a mild slope — the reach would have to be steep enough, or
short enough, that gravity stops setting the speed — and the case declares
the slope mild exactly so that seeing it is a violation rather than the
expected behaviour.

## What no budget can see

The mass closes perfectly when the gauge lies. A stage is not a store:
nothing is differenced over it, and there is no second reading to close it
against, so a model can conserve every millimetre of water and still report
a depth that has nothing to do with the flow it is carrying.
`momentum/routing-conservation` bounds the water the reach holds and
`momentum/stage-discharge-monotonic` asks whether the gauge *moves with* the
discharge, but neither asks whether the depth reported is one the section
could deliver that flow through. That is the gap this fills.

The nearest neighbour is `momentum/uniform-flow-friction-consistency` (#114),
and it asks a different question: whether the depth, the slope and the roughness
a model reports are mutually consistent with each other at *steady* flow,
through a friction balance. A model can satisfy that balance at every plateau
and still be scored here, because a pair can be internally consistent and be a
pair from the wrong regime — the friction balance says what depth this flow
needs, and this criterion says which regime that depth puts the flow in. Neither
subsumes the other, and neither is a duplicate of the other's mechanism.

This is not an assumption about what the other criteria do; it is what the
must-fail model is built to demonstrate. Its gauge applies the right law —
Manning's — to the wrong geometry: the depth a section five times wider would
need, which is what a rating curve reused from another reach looks like. Its
water is the exact bucket's, its depth varies (CV(depth) **0.56 to 0.94**
across the gate seeds read over the whole record, **identical** to the honest
gauge's — a constant factor cancels from a coefficient of variation, so no
variability threshold can separate them) and rises monotonically with the flow,
and its rating is single-valued. Driven through every probe in the suite that can ask
it anything, it **passes all of them** — including all three criteria of
`momentum/stage-discharge-monotonic` — and fails only `froude_subcritical`.

## Why a frequency, and why the limit is zero

Froude is a hard constraint, so the tolerance is absolute and small: the 5%
margin covers the rectangular-cross-section approximation and the rounding
in a stage an adapter derived from a normal-depth relation, not a model
being a little wrong about its own hydraulics.

The measure is nevertheless the **share of scored steps** that went
supercritical rather than the worst step, with a default limit of zero.
Every step is a separate claim about the reach, so the share is the honest
reading, and it keeps two very different failures apart: a model that spikes
once at the edge of the record reports a small fraction, while a reach that
is supercritical through every flood reports most of the record. The limit
is a parameter, so a future case where a small share is the honest reading
can say so without touching the criterion.

A step is **not scored** only when nothing is moving through it. The depth
carries a one-centimetre floor so a receding flow does not divide by zero, but
the floor is a floor on the division and not a dryness test: a step carrying
water is scored at the clamped depth however shallow the model says it is.
Deciding dryness on the depth would be backwards here — at a given discharge
the shallowest step is the *fastest* in the record, `Fr = Q / (w sqrt(g) d**1.5)`
growing without bound as `d` falls — so it would excuse exactly the steps the
probe exists to catch, and in the wrong direction: a model carrying its flood
peaks at a millimetre of depth would be failed only while it was *nearly*
honest, and would pass once it was wrong enough for the peaks to leave the
sample. Dryness is read from the flow, and the flow is read as a magnitude, so a
step reported with a negative discharge is water moving rather than an excuse.

## Where the separation comes from

A gauge that reads the flow it is carrying — Manning normal depth in the
declared section — is subcritical by construction: substituting the
normal-depth relation into Fr leaves `Fr = sqrt(S) / n * d**(1/6) / sqrt(g)`,
which on a mild slope with a realistic roughness sits near 0.4 at every flow
the case produces. Across the four must-pass baselines the worst Fr runs between
**0.376 and 0.453** over the gate's three seeds — and between **0.369 and
0.466** over seeds 0–99 — and no scored step exceeds one on any seed.

The must-fail puts the same flow through a section five times wider, so its
gauge reports a depth `5**0.6` (~2.5) times shallower and the Froude number
the pair implies is `5**0.9` (~4.3) times the honest one. On the gate seeds it is
supercritical on **83.3–90.1%** of scored steps, with a worst Fr of 1.77 to 1.93
and a median above one; over seeds 0–99 the same figures are **82.2–99.3%**,
1.75 to 1.98, and a median of 1.17 to 1.24. The populations do not touch: the
honest gauges' worst step (0.453 on the gate seeds) is 3.9 times below the
failing model's worst, and the failing model sits above the limit on most of the
record rather than at one marginal step.

| Criterion | Asserts |
| --- | --- |
| `froude_subcritical` | the reach stays subcritical on every scored step |
| `non_degenerate` | the runoff varies and responds to the weather |

Must-fail: `reference_shallow_rating`, the exact bucket with its gauge drawn
for a section five times wider than the one the case declares. Its water is
conserved to the floating point and its runoff is the bucket's, so the defect
is confined to the gauge — and the gauge is not obviously broken. It rises
with the flow, monotonically, and the whole rest of the suite accepts it. The
single change is the width the rating is evaluated at, which is what isolates
the one thing this probe asserts.

## Limitations

Five things this probe does not see — or refuses to score — so that a green
result is not read as more than it is.

**It reads the pair, not either reading against the truth.** A model that
derives its stage from its own discharge through a normal-depth relation is
consistent by construction, whatever that discharge is: substituting the
relation into `Fr` leaves `sqrt(S) / n * d**(1/6) / sqrt(g)`, which depends on
the depth and not on the flow, so a discharge that is wrong by a constant
factor still yields a subcritical pair. The absolute readings are the mass
probes' business; this one asks only whether the two numbers a model reports
could both be true of the declared section.

**The cross-section is rectangular, and it is the case's.** The relation
scored is the wide-rectangular one; on a compound or vegetated section the
depth a given discharge needs differs, and the 5% margin covers that
approximation rather than a wrong velocity. The width is the case's and not
the model's, so every model is judged in one reach — which is the intent, but
it also means the probe never checks that a model chose a sensible geometry
of its own.

Because that section decides the verdict, the probe declares
`requires.static: [width_m, slope, manning_n, bed_elevation_m]`. A model that
never read them — one carrying its own river width, or its own roughness — is
recorded **N/A (INCOMPATIBLE)** rather than failed, because the suite's rule is
that a verdict may not rest on an input the model never saw. The roughness is in
that list for the same reason the width is and not as an afterthought: a reach
is mild only relative to its roughness, so a gauge drawn with `n = 0.012` inside
the declared width is consistent with itself and still implies a velocity this
reach could not deliver.

**The datum is the contract's, and the criterion subtracts it.** `stage` is a
water-surface elevation on the fixed vertical datum the case declares
(`AGENTS.md`), and the flow depth is `stage - bed_elevation_m`. The probe
therefore requires `bed_elevation_m` as well as the section: a model that never
read the datum is judged against a number it never saw and is **N/A
(INCOMPATIBLE)**, not failed. That is how a gauge reporting the depth itself —
the convention this probe used to carry — is recorded: `depth_series` recognises
the mismatch and raises, the criterion lets the exception through, and the
harness turns it into N/A for that seed rather than into a conservation
violation.

What the subtraction does **not** catch is an offset added to a correctly formed
level, and that window is measured rather than hypothetical. Adding it to the
elevation the must-fail model reports: **+0.5 m** turns one of the three gate
seeds into a pass, **+1 m** turns all three, and **+5 m** leaves all three
comfortable. Only at **+10 m** does the ceiling below start refusing every step.
A metre-scale datum error is therefore indistinguishable from deep water, which
is why the datum is fixed by the contract and subtracted, not inferred from the
numbers a model reports.

`max_depth_m` is a **section sanity bound**, not that refusal: it is set here to
10 m, four times the case's declared bankfull depth and 2.3 times the deepest
reading a must-pass baseline gives over the gate seeds — **4.44 m**, from
`reference_bucket` on the second gate seed, whose worst Fr there is 0.453 — so it
never comes near honest flow. Over seeds 0–99 the deepest honest reading is
**5.27 m**, on seed 17 at a worst Fr of 0.466, which still leaves the ceiling
1.9x of room.
What it catches is the other end: a model that
ignored the declared datum and reported a level offset far enough to put its
depth past what the section can hold, which is refused rather than scored. The
ceiling refuses and does not clamp, deliberately: an offset only ever adds to the
depth, so clamping the reading back into range would re-admit exactly what the
bound exists to exclude.

**A non-finite reading is refused rather than skipped.** The mask compares
against bounds, and `NaN` compares false against every one of them, so a step
carrying a bad value would leave the scored set silently — one step fewer, and a
share measured on what is left. The criterion refuses the record instead and
reports how many steps were unreadable. Leaving a step out is for steps the
*reach* makes unscorable, not for values the arithmetic cannot carry.

**A shallow depth is not an excuse.** A step with no water in it is skipped; a
step with a millimetre of water in it is the *fastest* step in the record and is
scored at the one-centimetre floor. An earlier version of this criterion decided
dryness from the depth instead, and a model could buy a pass with it: reporting
its flood peaks at a millimetre left the 438 steps it would have failed on out of
the sample, while the remaining 70% of the record was still scored, comfortably
above the 5% `min_scored_fraction`, so nothing reported the third of the flood
that went unexamined. Reading dryness from the flow closes that door.

What is still open, and is the one number worth revisiting: a step is excused
when the model reports *exactly* no flow through it, and the model is the one
reporting. A model that zeroed its worst steps would still be bounded only by
`min_scored_fraction` — at 5%, a fifth of the way down. After the change the
honest baselines score 95.9% of the record or more, so the floor has room to be
raised towards them; that is a decision for the case rather than for the
criterion, and it is stated here rather than left to be discovered.

## What would make it sharper

`Fr = v / sqrt(g * D)` is formed here from `Q` and `d` through the
wide-rectangular relation. Once a model can report the section-averaged
velocity — `vel`, proposed in #95 and not in `main` yet — the numerator can be
the model's own velocity instead, which removes the rectangular assumption from
one of the two terms and lets the probe score a section it did not have to
assume.

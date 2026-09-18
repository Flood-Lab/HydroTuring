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

This is not an assumption about what the other criteria do; it is what the
must-fail model is built to demonstrate. Its gauge applies the right law —
Manning's — to the wrong geometry: the depth a section five times wider would
need, which is what a rating curve reused from another reach looks like. Its
water is the exact bucket's, its stage varies (cv 0.56 to 0.94 across the
gate seeds, the same variability the honest gauge has, so no variability
threshold can separate them) and rises monotonically with the flow, and its
rating is single-valued. Driven through every probe in the suite that can ask
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

Steps where the reach is effectively dry are **not scored**. The depth
carries a one-centimetre floor so a receding flow does not divide by zero,
but a step that only clears the floor because of the floor says nothing
about velocity — scoring it would let a dry reach fail on arithmetic.

## Where the separation comes from

A gauge that reads the flow it is carrying — Manning normal depth in the
declared section — is subcritical by construction: substituting the
normal-depth relation into Fr leaves `Fr = sqrt(S) / n * d**(1/6) / sqrt(g)`,
which on a mild slope with a realistic roughness sits near 0.4 at every flow
the case produces. Across the four must-pass baselines the worst Fr over the
gate's three seeds runs between **0.376 and 0.453**, and no scored step
exceeds one on any seed.

The must-fail puts the same flow through a section five times wider, so its
gauge reports a depth `5**0.6` (~2.5) times shallower and the Froude number
the pair implies is `5**0.9` (~4.3) times the honest one. It is supercritical
on **83–90%** of scored steps, with a worst Fr of 1.77 to 1.93 and a median
above one. The populations do not touch: the honest gauges' worst step (0.453)
is 3.9 times below the failing model's worst, and the failing model sits above
the limit on most of the record rather than at one marginal step.

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

Three things this probe does not see, so that a green result is not read as
more than it is.

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

Because that width is the denominator, the verdict rests on the model having
drawn its stage for the declared section, and the probe declares
`requires.static: [width_m]` accordingly. A model that never read `width_m` —
one carrying its own river width, or reporting a water level in its own
datum — is recorded **N/A (INCOMPATIBLE)** rather than failed: the suite's
rule for a verdict that would rest on an input the model never saw. Failing
it instead would put a conservation violation in the archive for a model
judged against a number it never read.

**Steps at the depth floor are skipped.** A reach that is dry or nearly dry
is not scored, so a model that only misbehaves at baseflow — reporting a
velocity it could not have when there is almost no water — passes. The
opposite case is covered: a gauge drawn for a wider section is too shallow at
every flow, and is caught on most of the record rather than only in floods.

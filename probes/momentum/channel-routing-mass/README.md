# momentum/channel-routing-mass

A reach holds the water it has generated as runoff and not yet released.
Call that store `channel`. This probe asks a model that reports it to obey
the one thing a store with nothing entering it can do: **drain**.

The physical statement is deliberately weak and one-sided. On a step where
no new water enters the reach, the store may fall — indeed it must, while it
releases what it holds — but it may not rise. Water that fell earlier is
still draining through, and that is exactly what the criterion allows. What
it forbids is a store that fills from nothing.

That is the gap `momentum/routing-conservation` leaves open. That probe bounds
the store: never negative, never more than a fifteen-day hydrograph of the
recent peak could have put there. This one asks the complementary question —
not *how much* is in the store, but *which way it moved* — and the answer does
not depend on the size of the error, only on whether it is applied where
nothing enters to mask it.

The two are not two readings of one bound, and the difference is measurable
rather than rhetorical. A leak that accumulates is caught by both: the store
walks away from its own hydrograph. An inflow that does **not** accumulate is
caught only here — `reference_unreported_inflow` routes exactly, closes its
budget to the floating point, and keeps the store inside the bound with an
order of magnitude to spare, yet rises on half of its recession steps.

No whole-catchment budget can see this. The error lives entirely inside the
channel: the water is subtracted from the reach, not from the soil or the
snowpack, so `closure` closes and `runoff_bounds` holds. The contradiction
is between the store and the weather, and it takes a recession to expose it.

## Why a frequency, not a worst case

An honest model can show an isolated small rise. A discrete unit
hydrograph renormalised over a partial history moves the store by about
`1e-5` of its own scale; the first steps of a record, where the hydrograph
is still filling from a zero state, move it by more. Scoring the worst step
would fail such a model on one startup sample, which is arithmetic and not
physics.

So the measure is the **share of recession steps on which the store rose**
above a relative floor. A reach that receives water it never generated rises on
nearly every recession step where nothing else moves the store:
`reference_unreported_inflow` sits at 0.48 to 0.52, and
`reference_stuck_router`, whose leak accumulates, at 0.57 to 0.67 over the gate
seeds. The honest side has to be read off the submitted models rather than the
reference ones — the four reference baselines sit at 0.00, but
`sacsma_snow17` rises on 8.7% of its recession steps and `summa` on 1.5%, and
the section below says what actually moves the store there. Against 8.7% the
threshold of 0.25 is about three times the honest share and about a third of
the failing one.

The floor a rise must clear scales with the store's own maximum, so it reads
as a fraction of the reach rather than as an absolute depth, and the first
steps of the scored window are excluded so that a hydrograph filling from
zero is not scored as a reach filling from nothing.

## What "no water enters" means

It is read off the **forcing**, not off the model: a step is in recession
when it and the three days before it carried no rain. Three days lets the
fastest surface response of an event clear before the step is counted, so a
small rise from a storm that has not quite finished draining is not scored.
A model cannot declare itself in recession to dodge the test, because the
label comes from the weather it was given.

The label is about **rain**, and a rainless step is not necessarily a step on
which nothing moves the store. On SAC-SMA the rises that survive the threshold
are not a slow hillslope still delivering water: the reviewer's analysis of this
archive puts 93% to 98% of them on a day the evaporative demand drops, against
half of all recession steps, which points to riparian evaporation drawn from
channel inflow — water leaving the reach when the atmosphere asks for it and
returning when it stops. `settle_days` clears the fastest surface response and
`max_rising_fraction` absorbs what is left, which is why the threshold is a
share rather than zero.

| Criterion | Asserts |
| --- | --- |
| `recession_drainage` | on rainless stretches the channel store does not rise |
| `non_degenerate` | the runoff varies and responds to the weather |

Must-fail, two of them, because they fail for different reasons and the second
exists to show the first's reason is not the only one:

- `reference_stuck_router`, a three-day triangular kernel summing to 0.9, so a
  tenth of every day's runoff never leaves and the shortfall accumulates across
  every recession step. `momentum/routing-conservation` catches it too, by the
  size of the error.
- `reference_unreported_inflow`, exact routing and a bounded inflow the model
  never reports — 0.05 mm over the record, applied on rainless steps. That is
  far inside the store bound: over 23 seeds, `routing-conservation` passes it
  with an order of magnitude to spare on this probe's weather and on its own.
  It fails here because the store rises on half of its recession steps. This is
  the baseline that shows the direction test is its own question.

## Where this departs from the accepted proposal

#71 asked for a residual: runoff leaving the reach over a recession against the
fall in `channel`, within 5%. The contract cannot form it. There is no variable
for water arriving at the reach, and `dis` is `mrro` times area, so the paired
quantities a model reports are its own runoff and its own store — such a
residual would compare a model against itself and pass by construction.

What ships instead is the one-sided version of the same physics: on steps where
no rain fell, the store may not rise. It is weaker as a statement and stronger
as a test, because it needs no inflow term, and it is a different question from
the store bound rather than a second way to read it — the two must-fail
baselines above are the demonstration. The ROADMAP entry keeps its original line
for the same reason: what this probe claims is what its criteria measure, and
the departure is recorded here rather than by editing the claim.

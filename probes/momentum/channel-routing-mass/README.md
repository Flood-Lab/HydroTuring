# momentum/channel-routing-mass

A reach holds the water it has generated as runoff and not yet released.
Call that store `channel`. This probe asks a model that reports it to obey
the one thing a store with nothing entering it can do: **drain**.

The physical statement is deliberately weak and one-sided. On a step where
no new water enters the reach, the store may fall — indeed it must, while it
releases what it holds — but it may not rise. Water that fell earlier is
still draining through, and that is exactly what the criterion allows. What
it forbids is a store that fills from nothing.

That is the gap `momentum/routing-conservation` leaves open. That probe
bounds the store: never negative, never more than a fifteen-day hydrograph
of the recent peak could have put there. A router that loses a fixed share
of every step stays inside those bounds for a long time, because the lost
water accumulates slowly and the bound scales with the flow. This probe asks
the complementary question — not *how much* is in the store, but *whether it
moved for a reason* — and a leak applied on every step is a store that
rises on every step where nothing is coming in to mask it.

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
above a relative floor. An honest router sits near zero — under 0.04 across
the reference models — while a router that loses a fixed share of every step
rises on nearly every recession step, because the error is applied every
step and nothing enters to hide it. Measured across seeds, the honest
routers stay below 0.04 and `reference_leaky_router` sits above 0.85; the
threshold of 0.25 has an order of magnitude of room on each side.

The floor a rise must clear scales with the store's own maximum, so it reads
as a fraction of the reach rather than as an absolute depth, and the first
steps of the scored window are excluded so that a hydrograph filling from
zero is not scored as a reach filling from nothing.

## What "no water enters" means

It is read off the **forcing**, not off the model: a step is in recession
when it and the three steps before it carried no rain. Three days lets the
fastest surface response of an event clear before the step is counted, so a
small rise from a storm that has not quite finished draining is not scored.
A model cannot declare itself in recession to dodge the test, because the
label comes from the weather it was given.

| Criterion | Asserts |
| --- | --- |
| `recession_drainage` | on rainless stretches the channel store does not rise |
| `non_degenerate` | the runoff varies and responds to the weather |

Must-fail: `reference_leaky_router`, a triangular kernel summing to 0.9, so
a tenth of every day's runoff never leaves and the shortfall accumulates in
the channel across every recession step.

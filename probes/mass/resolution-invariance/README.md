# mass/resolution-invariance

**Integrated runoff must not depend on the step the weather was given at.**

## The physics

Water is conserved under temporal aggregation. Take a month of weather
recorded minute by minute and average it over each hour: every hour carries
exactly the water its sixty minutes did. A model that has learned the
physics returns the same volumes from either record, because the same water
went in and the same catchment dealt with it. Runoff volume, evaporation
volume and the storage left at the end of the month must agree.

What is *not* invariant matters as much, and is the reason this probe sits
under mass rather than energy or momentum. A peak discharge, a stage, a
velocity head v²/2g or a momentum flux ρQv is nonlinear in the flow, and the
mean of a square is not the square of a mean. By Jensen's inequality,
aggregation can only attenuate such quantities, never sharpen them: the
hourly-mean hydrograph's peak is at most the minute hydrograph's peak, and
so on down the chain. That ordering is a rule for the energy and momentum
tracks, testable once a probe routes flow through a channel: **aggregation
cannot increase a peak, an energy head, or a momentum flux.** This probe
judges the volumes only.

## What it catches

Models are built at a step. Train an LSTM on hourly forcing and it learns
that one row is one hour; the rate of recession, the size of a wet-day
depth, the number of rows a storm occupies are all baked in. Run it at the
minute step and each row is still an hour to it, so it drains its stores and
consumes its evaporative demand sixty times too fast. Nothing in a
single-resolution evaluation can see this. Running the same weather at two
steps does, and the volumes disagree by tens of percent of the rainfall.

The deliberately broken baseline, `reference_fixed_step`, is the reference
bucket with its step hard-wired to a day. On the three gate seeds its runoff
and evaporation volumes differ between the minute and hourly runs by 20 to
87 percent of the month's precipitation; the reference bucket, which turns
rates into depths with the step it is given, agrees to within 0.1 percent.

A model that only runs at one step is not run at all. It is scored
INCOMPATIBLE, with the steps it lacks named, because the honest statement
about such a model is that it cannot take the transform, not that it failed
it.

## How the case is generated

`generate.py` draws one month of minute-by-minute weather per seed, with ten
days of spinup in front: Poisson storm arrivals, lognormal durations, gamma
depths spread over five-minute bursts with lognormal weights so the
intensity varies inside the hour, a diurnal temperature cycle on a wandering
daily anomaly, and potential evaporation confined to daylight. The `hourly`
variant is the same draw averaged over each hour. The weather is drawn once
and then aggregated, never drawn twice: a second draw at the coarse step
would be a different month, and the comparison would no longer isolate the
step.

Forcing stays in the contract's units at every step. Precipitation and
potential evaporation are rates in mm per day, so a one-millimetre burst in
one minute is a rate of 1440 mm/day, and the per-step depth is always the
rate times the step. A model reads the step from `timestep` in the request.

## Criteria

| Criterion | What it asserts |
| --- | --- |
| `closure` | the budget closes on the minute record (invariance alone is satisfiable by reporting zeros twice) |
| `state_bounds` | storages stay physical |
| `forcing_fidelity` | the model echoes the precipitation it was given, at the minute step |
| `resolution_invariance` | runoff and evaporation volumes, and the storage at the end of the month, agree between the minute and hourly runs to within 5 percent of the precipitation that fell |

## Baselines

- `must_pass: reference_bucket`, step-aware, identical at a daily step to
  the model every mass probe already passes.
- `must_fail: reference_fixed_step` on `resolution_invariance`.

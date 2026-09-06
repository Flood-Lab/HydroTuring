# mass/resolution-invariance

**Integrated runoff must not depend on the step the weather was given at.**

## The physics

Water is conserved under temporal aggregation. Take a month of weather
recorded minute by minute and average it over each hour, or each day: every
hour and every day carries exactly the water its minutes did. A model that
has learned the physics returns the same volumes from any of the records,
because the same water went in and the same catchment dealt with it. Runoff
volume, evaporation volume and the storage left at the end of the month
must agree.

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

## What it measures

Models are built at a step. Train an LSTM on hourly forcing and it learns
that one row is one hour; the rate of recession, the size of a wet-day
depth, the number of rows a storm occupies are all baked in. Run it at the
minute step and each row is still an hour to it, so it drains its stores and
consumes its evaporative demand sixty times too fast. Nothing in a
single-resolution evaluation can see this. Running the same weather at two
steps does, and the size of the effect is a number: how far the month's
runoff moved because the step moved, in percent of the rain that fell.

Every model is measured, whatever its manifest says about the steps it
supports. The harness serves the weather at the minute, the hour and the
day, runs a model at its own step and at the next finer one (the finest
pairs with the next coarser), and integrates both. A daily model is therefore
compared daily against hourly; an hourly model hourly against minutes. The
declaration in the manifest is where the model lives, not a way out, and an
adapter must feed the rows at the step it is handed rather than resample
them, or the probe cannot see the model at all.

The deliberately broken baseline, `reference_fixed_step`, is the reference
bucket with its step hard-wired to a day. Run on hourly rows it treats each
hour as a day, and on the three gate seeds its runoff volume differs from
its own daily run by 33 to 49 percent of the month's precipitation. The
reference bucket, which turns rates into depths with the step it is given,
agrees to about 1 percent between daily and hourly, which is the
discretisation of its own equations at a daily step and the floor the 5
percent rule sits on, and to 0.1 percent between hourly and minute. The
five daily-only reference models that predate this probe are measured too,
and fail the same way, because they are the old bucket. The catchment has a
shallow soil (120 mm) so that storms actually produce saturation-excess
floods within the month; on the closure probe's deeper soil this weather
never leaves baseflow.

Only runoff is required, so a model that reports streamflow and nothing else
is measured here rather than stopped at INCOMPLETE. Closure is the closure
probe's business; this one asks a question every runoff model can answer.

## How the case is generated

`generate.py` draws one month of minute-by-minute weather per seed, with ten
days of spinup in front: Poisson storm arrivals, lognormal durations, gamma
depths spread over five-minute bursts with lognormal weights so the
intensity varies inside the hour, a diurnal temperature cycle on a wandering
daily anomaly, and potential evaporation confined to daylight. The `hourly`
and `daily` variants are the same draw averaged over each hour and each
day. The weather is drawn once and then aggregated, never drawn twice: a
second draw at a coarser step would be a different month, and the
comparison would no longer isolate the step.

Forcing stays in the contract's units at every step. Precipitation and
potential evaporation are rates in mm per day, so a one-millimetre burst in
one minute is a rate of 1440 mm/day, and the per-step depth is always the
rate times the step. A model reads the step from `timestep` in the request.

## Criteria

| Criterion | What it asserts |
| --- | --- |
| `resolution_invariance` | runoff volume, and evaporation volume and end-of-month storage where the model reports them, agree between the model's step and the next finer one to within 10 percent of the precipitation that fell; the worst disagreement is the reported value. The limit is calibrated on physical models: an exact bucket moves 1 percent, the two FLEX models 2.5 and 5.1, because a partition that is nonlinear in intensity legitimately answers hourly and daily rain differently; a fixed-step model moves 49 |
| `non_degenerate` | runoff actually varies, so the invariance cannot be satisfied by reporting nothing twice; the runoff ratio is reported rather than judged on a month and the rainfall-response check is off |

## Baselines

- `must_pass: reference_bucket`, step-aware, identical at a daily step to
  the model every mass probe already passes.
- `must_fail: reference_fixed_step` on `resolution_invariance`, and
  `reference_degenerate` on `non_degenerate`.

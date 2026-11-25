# mass/resolution-invariance

**Integrated runoff must not depend on the step the weather was given at.**

## The physics

Hand a model the same weather at two steps, a day as one daily row or as
twenty-four hourly rows each carrying the day's rates, and the same water at
the same intensity goes in. A model that integrates its equations with the
step it is given returns the same runoff volume, the same evaporation volume
and the same storage at the end of the month, up to the discretisation of
those equations at the coarser step.

That holds only when the intensity is held. Conservation fixes the total at
every step, not how a model splits it, and averaging rain over a longer row
lowers its intensity. A model whose split depends on intensity then answers
the averaged rain differently, and is right to: infiltration excess, rain
arriving faster than the soil takes it, is the plain case. A bucket with a
2 mm/h infiltration capacity, given storms of five-minute bursts averaged to
hours and to days, moves its month's runoff between the hourly and the daily
run by more than 10 percent of the rain on 87 of 200 seeds, and by up to 59
percent, with no error anywhere in it. Version 2 of this probe served the weather that way and failed such
models for their physics. This version holds the intensity, so what is left
to move the volumes is the step.

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
its own daily run by 18 to 49 percent of the month's precipitation. The
reference bucket, which turns rates into depths with the step it is given,
agrees to about 1 percent between daily and hourly, which is the
discretisation of its own equations at a daily step and the floor the 10
percent limit sits on. The five daily-only reference models that predate
this probe are measured too, and fail the same way, because they are the
old bucket. The catchment has a shallow soil (120 mm) so that storms
actually produce saturation-excess floods within the month; on the closure
probe's deeper soil this weather never leaves baseflow.

Holding the intensity costs some power, and over 100 seeds it can be
counted. The fixed-step bucket moves more than 10 percent on 80 of them,
against 96 under version 2, whose bursts it also answered wrongly; on a dry
month a bucket that treats an hour as a day can still land close to its
daily self. Across the three seeds a model is run on it is still caught
with near certainty. The physical models go over the limit on 2 of the 100
seeds, both FLEX-Topo at 10.0 and 16.3 percent, against 4 under version 2
(FLEX-Topo and SAC-SMA twice each).

Only runoff is required, so a model that reports streamflow and nothing else
is measured here rather than stopped at INCOMPLETE. Closure is the closure
probe's business; this one asks a question every runoff model can answer.

## How the case is generated

`generate.py` draws one month of weather per seed, with ten days of spinup
in front, minute by minute: Poisson storm arrivals, lognormal durations,
gamma depths, a diurnal temperature cycle on a wandering daily anomaly, and
potential evaporation confined to daylight. It averages that draw over each
day, which is the daily record version 2 served, and holds each day flat:
the `hourly` variant is every day repeated over its 24 hours and the
`minute` variant over its 1440 minutes. Every step carries the same water at
the same rates, and only the length of the row changes. The weather is drawn
once, never redrawn per step: a second draw would be a different month, and
the comparison would no longer isolate the step.

Holding the day flat also removes the diurnal cycle from the finer steps.
That is deliberate. A diurnal cycle is weather inside the day that a daily
row cannot carry, and serving it at the hour alone would hand the finer run
information the coarser one never had.

Forcing stays in the contract's units at every step. Precipitation and
potential evaporation are rates in mm per day, and the per-step depth is
always the rate times the step. A model reads the step from `timestep` in
the request.

## Criteria

| Criterion | What it asserts |
| --- | --- |
| `resolution_invariance` | runoff volume, and evaporation volume and end-of-month storage where the model reports them, agree between the model's step and the next finer one to within 10 percent of the precipitation that fell; the worst disagreement is the reported value. The limit is calibrated on physical models: on the gate seeds the exact bucket moves 1.1 percent, FLEX-Lumped 1.7, FLEX-Topo 4.0 and SAC-SMA with Snow-17 0.7, the discretisation of their own equations at a daily step; a fixed-step model moves 18 to 49 |
| `non_degenerate` | runoff actually varies, so the invariance cannot be satisfied by reporting nothing twice; the runoff ratio is reported rather than judged on a month and the rainfall-response check is off |

## Baselines

- `must_pass: reference_bucket`, step-aware, identical at a daily step to
  the model every mass probe already passes, and the three other physical
  models.
- `must_fail: reference_fixed_step` on `resolution_invariance`, and
  `reference_degenerate` on `non_degenerate`.
- `tests/test_resolution.py` adds a Horton bucket, step-aware with a 2 mm/h
  infiltration capacity. It passes this version to about 1 percent and fails
  version 2's weather by 21 to 35 percent on the seeds it is checked on.

## Versions

- **3**: one daily record held flat over the hours and the minutes, so the
  intensity is the same at every step. Infiltration excess no longer moves
  the volumes.
- **2**: minute weather with storms in five-minute bursts, averaged to hours
  and days. Failed correct infiltration-excess models for their physics.

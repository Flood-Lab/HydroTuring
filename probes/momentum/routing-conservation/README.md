# momentum/routing-conservation

`channel` was added to the contract so that a model with a unit hydrograph
or a channel store could report the water it holds in transit and close
its budget. This probe asks whether that store behaves like one. Two
assertions: it is never negative, and at any step it holds at most
`max_lag_days` (15) times the largest runoff rate of the preceding month,
which is the most a hydrograph of that length can have accumulated from
the flow feeding it.

It sits under the momentum law because the channel is where momentum is
conserved, and conservation of the water in it is the shadow that every
routing scheme casts on the mass budget. A model that does not route reports
its channel store as identically zero, as the exact bucket does, and passes
trivially; the two FLEX models, which route through a triangular lag, are
the real test.

| Criterion | Asserts |
| --- | --- |
| `routing_conservation` | the channel store is non-negative and never exceeds 15 days of the largest recent runoff, plus `min_allowance_mm` |
| `non_degenerate` | runoff varies with the weather |

Must-fail: `reference_stuck_router`, the bucket with a three-day kernel that
sums to 0.9.

## The dry spells are what make a small leak visible

The allowance is a multiple of the largest runoff of the preceding month, so
the weather decides how strict it is — and a record that never goes dry keeps
the ceiling at a storm-time level all year. Over a rainless spell the runoff
falls to the baseflow and then to nothing, and thirty days later the peak has
left the window: the ceiling comes down with it.

That is the difference the four storms and the spells that follow them make. A
kernel summing to 0.999 releases all but a tenth of a percent of each day's
runoff, so it leaves a residue of `0.001 x the runoff so far`. On the gate seeds
that residue reaches **0.80 mm**, against a storm-time ceiling of **433 to 635
mm** — five hundred times smaller, and inside the bound everywhere the peak is
read from a storm. On this record the ceiling has decayed below it during the
drains, and the kernel fails on every gate seed.
`tests/test_routing_conservation_generator.py` pins both the weather's
properties and that separation.

## The absolute allowance

`min_allowance_mm` is the one part of the allowance that is not a multiple of
the recent peak, and it exists for the reach that has stopped being fed and
still holds the dead storage its own hydraulics keep. With the proportional term
decayed to nothing, a residue of a few thousandths of a millimetre would
otherwise be scored as a leak.

`wflow_sbm`'s river is the case it exists for, and its need is the largest
measured on this case: **0.0026 mm** through the droughts, against which 0.05 mm
leaves it about nineteen times. Over the three gate seeds the models that route
sit at **0.01 to 0.58 of the whole allowance**:

| Model | Worst ratio to the allowance |
| --- | --- |
| `summa` | 0.01 |
| `dhbv2` | 0.04 |
| `flex_topo` | 0.11 |
| `lisflood` | 0.12 |
| `wflow_sbm`, `sacsma_snow17` | 0.16 |
| `cwatm` | 0.47 |
| `flex_lumped` | 0.58 |

The tightest of those leaves about 1.7x of room, and the smallest fault this
probe exists to catch — a kernel keeping a tenth of a percent of each day's
runoff — is caught at **4.4x**. The value sits between them by measurement, not
by taste.

## Limitations

- **It needs the peak to decay.** The allowance is a multiple of the recent
  peak runoff, so a record whose runoff never falls — a glacier, a large
  regulated river, a case with perennial baseflow — keeps the ceiling high and
  hides the same leak. The generator's dry spells are what close that gap here;
  the criterion cannot close it on its own.
- **Timing is not read.** A kernel that sums to one but puts the water in the
  wrong day passes. This bounds what the store holds, not when it releases it;
  `momentum/routing-lag-consistency` is the probe that asks the timing question.
- **A small loss is not read directly.** What is compared is the store's level
  against its bound, so a kernel summing to 0.999 is caught only because the
  ceiling shrinks under it. Nothing here reads the fraction of the flow that was
  retained, and on a record with no dry spell that kernel passes.
- **The absolute allowance forgives a fixed depth.** Any model whose dead
  storage in the reach stays under `min_allowance_mm` is not distinguished from
  one that releases everything; the term trades that sensitivity for not
  failing models on channel geometry, and its size is set by the measurements
  above.

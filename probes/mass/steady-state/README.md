# mass/steady-state

**Under constant weather the catchment settles and the budget balances.**

## The physics

Hold rain, temperature and demand constant long enough and every store
finds the level at which what comes in equals what goes out. From then on
nothing changes: runoff and every storage hold still, runoff cannot exceed
the rain, and where evaporation is reported, runoff plus evaporation equal
the rain exactly, because storage is no longer changing. This is the one
place a budget can be checked without knowing the storage at all, which is
why the residual test here reaches models that report no storage.

## What it measures

A year of ordinary weather, then three years in which every day brings
2.5 mm of rain, 12 degC and 2.0 mm of potential evaporation. The three
seeds differ only in the spinup year, so the catchment enters the constant
stretch from three different states and must settle to the same place.
The last year is judged: every reported variable must vary by less than 1
percent of its level, runoff must not exceed the rain, and where
evaporation is reported the budget must balance to 1 percent of the rain.

What fails it: `reference_restless`, the reference bucket with a recession
that speeds up and slows down on its own thirty-day clock, whose runoff
varies by 60 percent in the last year of constant weather. Two of the older
baselines fail as well, for the right reasons: the closure cheater, whose
storage is the residual of the budget and drifts without end under
constant forcing, and the leaky bucket, whose 15 percent sink leaves the
budget 15 percent short at steady state.

## Results on the gate seeds

The exact bucket settles to 0.500 mm/day of runoff and 2.000 mm/day of
evaporation under 2.5 mm/day of rain, with zero variation in the last year
and the budget balanced to floating point.

## Criteria

| Criterion | What it asserts |
| --- | --- |
| `steady_state` | in the last year nothing reported varies by more than 1 percent, runoff stays below the rain, and runoff plus evaporation equal the rain where evaporation is reported |

## Baselines

- `must_pass: reference_bucket`
- `must_fail: reference_restless` on `steady_state`

# mass/steady-state

**Under constant weather the catchment settles and the budget balances.**

## The physics

Hold rain, temperature and demand constant long enough and every store
finds the level at which what comes in equals what goes out. From then on
nothing changes: runoff and every storage hold still, runoff cannot exceed
what enters, which is the rain plus any inflow the model declares as
`gwex`, and where evaporation is reported, rain plus the declared exchange
equal runoff plus evaporation exactly, because storage is no longer
changing. This is the one place a budget can be checked without knowing the
storage at all, which is why the residual test here reaches models that
report no storage.

## What it measures

A year of ordinary weather, then three years in which every day brings
2.5 mm of rain, 12 degC and 2.0 mm of potential evaporation. The seeds
differ only in the spinup year, so the catchment enters the constant
stretch from different states and must settle to the same place. The last
year is judged: every reported variable, every store included, must vary
by less than 1 percent of its level, runoff must not exceed the rain, the
declared inflow and whatever the stores released over the year, and where evaporation is reported the budget must
balance to 1 percent of the rain, with the storage change over the year
counted in it.

## Settling is not drift

How long a model takes to settle is its own business, and two years is not
always enough. The net input under this weather is only 0.5 mm/day, so a
catchment that leaves the spinup year dry spends a year refilling its
tension water before its slow stores start to fill. Version 1 of this
probe required stillness alone, and `sacsma_snow17`, a physical model
every probe must pass, failed it on 35 of seeds 0 to 399. Its worst seed,
139, entered the constant stretch with 75 mm of soil water against a
settled 274 mm, and its runoff was still rising by 2.3 percent in the
judged year, straight towards the same 0.512 mm/day every other seed
settles at. Calibrated slower baseflow makes it worse: with the lower
zone's primary recession `LZPK` at 0.002/day, version 1 failed 58 of 60
seeds, and at 0.001/day all 60.

So a variable that has not stopped is excused while it is still
approaching a level. Over the last two years, counting every step up or
down, its last year must travel at most 70 percent as far as the year
before. A geometric approach with an e-folding time tau gives
exp(-365 / tau), which passes up to tau = 1023 days, and so does a sum of
stores settling at different rates that rises past its level once on the
way, as SAC-SMA's two lower-zone stores do. A constant drift, a recession
on its own clock and noise all travel as far each year as the year before.

| Model | Seeds | Verdict | Distance ratio, last year to the one before |
| --- | --- | --- | --- |
| `sacsma_snow17` | 0–399 | 400 pass | at most 0.03 where excused |
| `sacsma_snow17`, `LZPK` 0.002, `LZFPM` 250 | 0–59 | 60 pass | at most 0.31 |
| `sacsma_snow17`, `LZPK` 0.001, `LZFPM` 300 | 0–59 | 60 pass | at most 0.60 |
| `reference_bucket`, `flex_lumped`, `flex_topo` | 0–99 | 100 pass each | nothing to excuse |
| `reference_calendar` | 0–99 | 100 fail | at least 0.81 |
| `reference_noise_sink` | 0–99 | 100 fail | at least 0.88 |
| `reference_restless`, `reference_cheater`, `reference_climatology`, `reference_slow_drift`, `reference_gw_slow_drift`, `reference_leaky_router`, `reference_stuck_router`, `reference_recharge_exchange` | 0–99 | 100 fail each | 1.00 |

The slow parameter sets were run in a scratch copy and are not packaged.
The limit is a statement about one record: a clock whose period is longer
than the last two years can look like an approach from inside them, and a
store that fills at a constant rate for more than a year, because its
deficit is that large against the net input, is indistinguishable from a
drift until it stops.

## What fails it

`reference_restless`, the reference bucket with a recession that speeds up
and slows down on its own thirty-day clock, whose runoff varies by 60
percent in the last year of constant weather and travels exactly as far as
it did the year before. The older baselines fail for the right reasons as
well: the closure cheater, whose storage is the residual of the budget and
drifts without end under constant forcing, and the leaky bucket, whose 15
percent sink leaves the budget 15 percent short at steady state.

## Results on the gate seeds

The exact bucket settles to 0.500 mm/day of runoff and 2.000 mm/day of
evaporation under 2.5 mm/day of rain, with zero variation in the last year
and the budget balanced to floating point.

## Criteria

| Criterion | What it asserts |
| --- | --- |
| `steady_state` | in the last year nothing reported varies by more than 1 percent unless it is still approaching a level, runoff stays below the rain plus declared inflow plus what the stores released, and rain plus declared exchange equal runoff, evaporation and the storage change where evaporation is reported |

## Baselines

- `must_pass: reference_bucket, flex_lumped, flex_topo, sacsma_snow17`
- `must_fail: reference_restless` on `steady_state`

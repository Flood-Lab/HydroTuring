# mass/dry-down

**Without rain, runoff can only fall, and only stored water can drain.**

## The physics

A year of ordinary weather fills the catchment. Then no rain falls for two
years, while temperature and evaporative demand keep their seasons. Two
things are then exact. Nothing arrives to raise the flow, so runoff, and
every storage the model reports, can only decrease. And the total that
runs off in the dry years cannot exceed what the catchment was holding
when the rain stopped, which is bounded above by its storage capacities:
here 320 mm of soil and 2 mm of canopy.

## What it measures

The criterion sums the runoff over the 730 rainless days and compares it
with the capacity bound, and checks that weekly means of runoff and of
every reported storage never rise by more than 2 percent of their level
after the first week, which allows for routing delay and for a model's own
sampling noise but not for a season. The exact model never rises at all.

What fails it is a model whose runoff comes from something other than the
water it has: `reference_climatology` emits the long-term mean for the day
of the year whatever falls, and over two dry years runs off 613 mm from a
catchment that could hold 322, rising with each winter. An LSTM that has
learned a baseflow floor from its training basins does the same thing more
quietly.

## Results on the gate seeds

The exact bucket drains 55 mm in two dry years against the 322 mm bound
and never rises. The climatology model runs off 613 mm and rises by 3.5
percent a week each autumn. The streamflow-only reference drains 42 mm and
passes: a model needs no evaporation to be right about this.

## Criteria

| Criterion | What it asserts |
| --- | --- |
| `dry_down` | runoff over the dry stretch stays under the storage bound; runoff and storages never rise |
| `non_degenerate` | runoff is there and varies as it drains, so a constant zero cannot pass |

## Baselines

- `must_pass: reference_bucket`
- `must_fail: reference_climatology` on `dry_down`

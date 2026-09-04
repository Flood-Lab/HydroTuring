# mass/warming-response

**Warmer air over the same rain must mean less runoff, and cooler air more.**

## The physics

Hold the rain fixed and warm the air. The demand for water from the
atmosphere rises with temperature, so more of the same water evaporates and
less of it runs off. Cool the air and the reverse happens. Integrated over
a long enough stretch, where the change in storage is small next to the
fluxes, the extra evaporation and the lost runoff are the same water: per
unit of demand added or removed, the two responses sum to about one, with
the split set by how water-limited the catchment is.

The probe pushes the driver both ways and scores each direction against
the control on its own. That is not redundancy. A model can carry the
relationship as a one-way rule, learned from the hot dry years in its
training record and never exercised the other way; it can saturate, so
that more demand does nothing once the soil is dry while less demand still
does; or it can clip an output at zero. Any of these gets one direction
right and the other wrong, and only the pair shows it.

This is not a closure statement. A model can close its budget exactly on
both runs and still get the relationship backwards, or ignore the
temperature altogether and close just as well. What is being tested is the
sign, and the size, of one internal relationship the model must carry if it
is doing hydrology rather than fitting a hydrograph: the way streamflow
responds to temperature at fixed precipitation.

## The case

Ten years of daily weather for a temperate, rain-dominated catchment: mean
air temperature around 15 degC, four sub-zero days in a decade, about 880
mm of rain and 1000 mm of potential evaporation a year, so the catchment
sits on the water-limited side and evaporation is free to rise with demand.
Snow is avoided on purpose. In a snowy catchment warming also moves the
melt, and the sign of the runoff response inside any one season depends on
timing rather than on evaporation.

The `warmer` and `cooler` variants are the same draw with the air three
degrees warmer and three degrees cooler over the scored years, potential
evaporation recomputed from each, which adds and removes about 140 mm a
year of demand. The spinup is untouched so every run enters the window
from the same state, and the precipitation is byte-identical, which the
criterion checks.

## Criteria

| Criterion | What it asserts |
| --- | --- |
| `response_sign` | per unit of demand added or removed, integrated runoff moves the opposite way and evaporation, where reported, the same way, by at least a tenth and by no more than all of it, under warming and under cooling separately; the reported value is the weakest runoff response per unit of demand |
| `non_degenerate` | runoff varies with the weather, so a constant cannot pass |

Only runoff is required, so a model that reports streamflow and nothing else
is measured here. A submitted model's evaluation window is widened to a full
year (`min_window_days: 365`), because inside a month the runoff response is
dominated by storage and can point either way.

## Baselines

- `must_pass: reference_bucket`. Its evaporation is limited by soil moisture
  and by demand, so warming raises it and lowers runoff by the same water,
  and cooling does the reverse: on the gate seeds, 0.20 of the demand
  change under warming and 0.24 under cooling, the asymmetry being the soil
  running short of water to give when demand rises.
- `must_fail: reference_degenerate`, which evaporates all precipitation
  whatever the temperature and so responds by exactly nothing, and
  `reference_streamflow_only`, which runs off a fixed share of a store that
  only ever sees rain and never reads the temperature at all.

## Adapting a model to this probe

The probe changes `tas` and `pet` together, in both directions, and holds
`pr`. An adapter has to carry that into the model's own inputs consistently: every input derived
from temperature or demand should move, and nothing derived from
precipitation may. Say in the model's README which inputs the probe's
temperature reaches. For the Google Flood Hub model, for instance, it
reaches the two temperature products, the mocked net longwave radiation and
the temperature and demand attributes, while the four precipitation
products and the radiation tied to cloudiness stay fixed.

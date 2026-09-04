# mass/warming-response

**Warmer air over the same rain must mean more evaporation and less runoff.**

## The physics

Hold the rain fixed and warm the air. The demand for water from the
atmosphere rises with temperature, so more of the same water evaporates and
less of it runs off. Integrated over a long enough stretch, where the
change in storage is small next to the fluxes, the extra evaporation and
the lost runoff are the same water: their shares of the added demand sum to
about one, with the split set by how water-limited the catchment is.

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

The `warmer` variant is the same draw with the air three degrees warmer
over the scored years and potential evaporation recomputed from the warmer
air, which adds about 140 mm a year of demand. The spinup is untouched so
both runs enter the window from the same state, and the precipitation is
byte-identical, which the criterion checks.

## Criteria

| Criterion | What it asserts |
| --- | --- |
| `response_sign` | integrated runoff falls, and evaporation rises where the model reports it, by at least a tenth of the added demand and by no more than all of it; the reported value is the runoff change as a share of the added demand |
| `non_degenerate` | runoff varies with the weather, so a constant cannot pass |

Only runoff is required, so a model that reports streamflow and nothing else
is measured here. A submitted model's evaluation window is widened to a full
year (`min_window_days: 365`), because inside a month the runoff response is
dominated by storage and can point either way.

## Baselines

- `must_pass: reference_bucket`. Its evaporation is limited by soil moisture
  and by demand, so warming raises it and lowers runoff by the same water.
- `must_fail: reference_degenerate`, which evaporates all precipitation
  whatever the temperature and so responds by exactly nothing, and
  `reference_streamflow_only`, which runs off a fixed share of a store that
  only ever sees rain and never reads the temperature at all.

## Adapting a model to this probe

The probe changes `tas` and `pet` together and holds `pr`. An adapter has
to carry that into the model's own inputs consistently: every input derived
from temperature or demand should move, and nothing derived from
precipitation may. Say in the model's README which inputs the probe's
temperature reaches. For the Google Flood Hub model, for instance, it
reaches the two temperature products, the mocked net longwave radiation and
the temperature and demand attributes, while the four precipitation
products and the radiation tied to cloudiness stay fixed.

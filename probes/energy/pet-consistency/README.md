# energy/pet-consistency

The first probe under the energy law, built without an energy-reporting
model, because the surface energy balance has a water shadow every model
casts: potential evaporation is the energy available for latent heat, and
what a model actually evaporates says whether it used that energy when it
could and only when it could.

Four assertions, the regimes taken from the model's own reported soil
store so that every model has both: on the days the store is in the
wettest fifth of its range, cumulative evaporation is at least 0.7 of
potential (energy-limited) and no more than potential; on the days it is
in the driest fifth, the ratio to potential is lower than in the wettest
fifth by at least 0.05 (water-limited: evaporation efficiency does not
fall as the soil gets wetter); and over the record evaporation is at most
potential (the ceiling `et_plausible` also holds).

| Criterion | Asserts |
| --- | --- |
| `demand_consistency` | the three ratios above |
| `non_degenerate` | runoff and evaporation vary with the weather |

Must-fail: `reference_thirsty`, the bucket whose evaporation is a fixed two
percent of the soil store per day and never consults the demand. It
conserves water exactly and every mass probe passes it; this one does not.

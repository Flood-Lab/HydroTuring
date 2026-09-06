# mass/runoff-bounds

Over a long record the water budget of any catchment reads `P = ET + Q + ΔS`
with `0 ≤ ET ≤ PET` and `|ΔS|` bounded by what the catchment can hold. So
runoff alone satisfies

    P − PET − S_max ≤ Q ≤ P + S_max

and a model that reports runoff and nothing else can be asked whether its
runoff is possible. Below the lower bound it has evaporated more than the
atmosphere asked for, or lost water; above the upper it has made some.
`S_max` is the sum of the catchment's stated capacities plus whatever the
model reports at the start in stores that have no capacity to name (snow,
groundwater, channel). Both bounds carry a two percent slack on the rain.

This is the probe that makes streamflow-only models falsifiable on mass.
On the closure probe such a model is INCOMPLETE by construction; here it
has to answer.

| Criterion | Asserts |
| --- | --- |
| `runoff_bounds` | integrated runoff over ten years lies inside the bounds above |
| `non_degenerate` | runoff varies with the weather |

Physical models: the exact bucket runs off about 0.4 of the rain on this
catchment with demand at 0.85 of rain, well inside `[0, 1.04]`. The
degenerate model (Q = 0) fails the lower bound; a bucket reporting 1.3 ×
its runoff fails the upper. A model with a learned source term that
returns more than the rain fails the upper bound with nothing but its
hydrograph.

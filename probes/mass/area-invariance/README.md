# mass/area-invariance

A lumped catchment's budget is written per unit area. Tell the same model
the same weather on the same catchment at ten times the area, and every
depth in millimetres, runoff, evaporation, every store, must come back
identical; only discharge in m³/s scales. There is no tolerance to tune:
the two runs agree to floating point or they do not.

The probe exists because a model evaluated this week did not: δHBV 2.0's
regional-groundwater term reads the upstream area, and its runoff plus
evaporation went from 1.00 to 2.03 of the rain as the area went from 8 to
2500 km². A network can also learn area as a proxy for climate or
regime, and then a physically meaningless change moves its answer.

| Criterion | Asserts |
| --- | --- |
| `invariance` | every reported depth is unchanged between 250 and 2500 km², to 1e-6 relative |
| `non_degenerate` | runoff varies with the weather |

Must-fail: `reference_area_leak`, the bucket losing a share of runoff that
grows with the stated area.

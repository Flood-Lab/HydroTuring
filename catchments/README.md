# Catchments

Real gauged basins whose HydroATLAS attributes every probe hands to every
model. Only the weather is generated; see `docs/catchments.md` for the
proposal and `src/hydroturing/weather.py` for the generator that matches
each basin's monthly climatology.

| Id | Basin | Source |
| --- | --- | --- |
| `alpine-snow` | South Fork Payette River at Lowman, ID | `camels_13235000` |
| `maritime-rain` | Tilton River near Cinebar, WA | `camels_14236200` |
| `humid-karst` | North Fork River near Tecumseh, MO | `camels_07057500` |
| `lowland-agricultural` | Fish Creek near Crystal, MI | `camels_04115265` |

Attributes are copied verbatim from Caravan (Kratzert et al. 2023,
CC-BY-4.0), which aggregates HydroATLAS v1.0 over each basin. Add a basin
from a local copy of Caravan with `scripts/catchment_from_caravan.py`.

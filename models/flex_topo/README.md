# flex_topo

FLEX-Topo, the semi-distributed model of Zhi Li's
[HydrologicModels](https://github.com/chrimerss/HydrologicModels)
repository (`semi-distributed_model/`, commit `cc0aa6f`), in the pool as a
**physical reference**: every probe must pass it, and a probe that fails it
is examined before the model is.

Three landscape units run on the same weather, plateau, hillslope and
wetland, with fractions 0.471 / 0.468 / 0.061 computed from the
repository's Wark catchment grids with its own rule (hillslope: slope
> 11; plateau: HAND > 5 and slope < 11; wetland: 0 < HAND ≤ 5; the 5 %
of channel cells at HAND 0 fall in no class and the three are normalised).
Each unit has an interception store, an unsaturated store with a beta
partition and a fast reservoir. The plateau percolates `Pmax · Su/Sumax`
to a shared slow reservoir; the hillslope sends a share `D` of its fast
partition there instead; the wetland percolates nothing and draws
capillary rise `Cmax · (1 − Su/Sumax)` from it. Total runoff is the slow
outflow plus the area-weighted fast outflows, through a triangular lag.
Parameters are `B_run_model.py`'s Wark sets.

Four stated deviations from the repository code: transpiration takes
FLEX's `min(1, ·)` limit; the hillslope's preferential share `D` of the
fast partition is moved once, to the slow reservoir (the original also
subtracts it from the unsaturated store, losing it); the wetland's capillary rise is limited by and
removed from the slow reservoir with the *wetland* fraction (the original
passes the plateau fraction into the wetland routine, which creates or
destroys water whenever the two differ, and here they differ eightfold);
and the catchment's capacities come from `static.json`, the three `Sumax`
scaled to the catchment's soil capacity keeping their ratios and every
`Imax` set to the canopy capacity. The step is a parameter, rescaled as in
`flex_lumped`.

Reports area-weighted `canopy`, `mrso`, `gw` (the shared slow reservoir)
and `channel` (fast reservoirs plus water in the lag). No snow module:
`snw` is identically zero.

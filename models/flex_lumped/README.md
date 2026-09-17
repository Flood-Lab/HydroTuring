# flex_lumped

The lumped FLEX/HBV model of Zhi Li's
[HydrologicModels](https://github.com/chrimerss/HydrologicModels)
repository (`lumped_model/HBVMod.py`, commit `cc0aa6f`), in the pool as a
**physical reference**: every probe must pass it, and a probe that fails it
is examined before the model is. A probe pull request is checked against
this model and against `flex_topo` first.

Structure: an interception store (`Imax`), an unsaturated store (`Sumax`,
`beta`) whose beta partition sends a share of effective rain to a fast
linear reservoir (`Kf`), transpiration limited by soil moisture (`Ce`),
percolation `Pmax · Su/Sumax` to a slow linear reservoir (`Ks`), and a
triangular lag (`Tlag`). Parameters are the repository's Wark values;
`Sumax` and `Imax` are taken from the catchment's `static.json`, because a
physical model is told its catchment.

Relative to the repository code, the adapter changes transpiration and
timestep handling as follows. Transpiration is FLEX's published
`Ep · min(1, Su/(Ce·Sumax))` (the original has no `min`, and a wet soil with
`Ce < 1` transpires above demand; on the closure probe the difference is one
percent of ET). The step is a parameter, with per-day fractions rescaled as
`1 − (1 − k)^dt` and the lag in days. The adapter also optionally configures
the repository's native `Tlag` from supplied channel geometry, as described
below.

When `area_km2`, `main_channel_length_km` and
`centroid_channel_length_km` are supplied together, the adapter configures
the model's native triangular lag with a constant-celerity travel-time
approximation that is independent of the routing-lag probe's Snyder criterion:

```text
channel_travel_time = centroid_channel_length_km / (1.0 m s-1)
Tlag = 2 * (channel_travel_time + dt / 2)
```

The fixed `1.0 m s-1` flood-wave celerity is a synthetic prior, not a universal
river constant or a value fitted to the probe. The choice follows
[Beven (2020, Appendix equations A18--A20 and Figure A3)](https://doi.org/10.5194/hess-24-2655-2020),
which distinguishes mean water velocity from kinematic-wave celerity, derives
`c = dQ/dA`, and uses `1.0 m s-1` for an upland-channel example. The underlying
generalised kinematic routing method is described by
[Beven (1979)](https://doi.org/10.1029/WR015i005p01238). The value also lies
within the independently observed regional range of
approximately 0.8--1.6 m s-1 reported by
[Le Mesnil et al. (2021)](https://doi.org/10.5194/hess-25-1259-2021).
The adapter extrapolates this first-order celerity across the probe's synthetic
catchment-size ladder and uses `Lc / c` as the triangle's mode timescale. It
does not claim a site-calibrated hydraulic relation or an exact peak-time
prediction.
The original FLEX `Tlag` is the full base of a symmetric triangle whose bins
start at the beginning of a generated-runoff row. The half-step aligns the
target travel time, which starts at the generated-runoff interval centre, with
that row-based kernel; it is 12 hours at `PT1D` and scales with `dt` at finer
supported steps. The factor two converts the target triangle mode to its full
base. The kernel remains causal, non-negative and unit-sum. Without both
channel lengths the calibrated Wark value
`Tlag = 1.1 d` remains the documented fallback; an area supplied alone is
validated but does not define a travel path.

Reports `pr, evspsbl, mrro` and the stores `canopy` (Si), `mrso` (Su),
`gw` (Ss) and `channel` (the fast reservoir plus water inside the lag).
No snow module: `snw` is identically zero and precipitation below freezing
is treated as rain, which conserves water and is wrong about timing.

On the closure probe's first gate seed the budget closes to 4e-16 of the
precipitation; runoff ratio 0.40, ET 0.66 of potential; runoff volume at
the hourly step within 2.4 % of the daily one.

With a prescribed human withdrawal in the forcing (`abstr`), the adapter
takes it from the unsaturated store first and from the day's runoff second,
and declares what it removed as a negative `gwex`; absent the column it is
unchanged.

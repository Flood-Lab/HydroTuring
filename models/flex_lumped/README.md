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

Two stated deviations from the repository code: transpiration is FLEX's
published `Ep · min(1, Su/(Ce·Sumax))` (the original has no `min`, and a
wet soil with `Ce < 1` transpires above demand; on the closure probe the
difference is one percent of ET), and the step is a parameter, with
per-day fractions rescaled as `1 − (1 − k)^dt` and the lag in days.

Reports `pr, evspsbl, mrro` and the stores `canopy` (Si), `mrso` (Su),
`gw` (Ss) and `channel` (the fast reservoir plus water inside the lag).
No snow module: `snw` is identically zero and precipitation below freezing
is treated as rain, which conserves water and is wrong about timing.

On the closure probe's first gate seed the budget closes to 4e-16 of the
precipitation; runoff ratio 0.40, ET 0.66 of potential; runoff volume at
the hourly step within 2.4 % of the daily one.

# momentum/uniform-flow-friction-consistency

A steady, uniform reach has no local acceleration and no water-surface
gradient left to balance: the Saint-Venant momentum equation reduces to
equality of bed slope and friction slope. For the rectangular reach supplied
by this probe,

```text
d   = stage - bed_elevation_m
A   = width_m * d
R_h = A / (width_m + 2 d)
S_f = (manning_n * dis / (A R_h^(2/3)))^2.
```

The criterion evaluates the final 365 days and requires
`mean(abs(S_f - slope)) / slope <= 0.05`. It also reports the 95th percentile
and maximum residual. Absolute residuals are averaged, so a rating that is too
deep on one part of the block and too shallow on another cannot pass by
cancellation.

Before the momentum balance is read, the final block must actually be steady:
the coefficients of variation of discharge and water depth must each be at
most 1%. Non-finite values, non-positive discharge and stage at or below the
bed fail rather than being silently skipped. The generator keeps the entire
record under constant wet, temperate forcing and gives the models a one-year
spin-up before three scored years, so all five seeds have a wet steady block.

`area_km2` is present because the physical hydrologic baselines use it to turn
their runoff depth into the `dis` they report. It is not in `requires.static`:
the criterion consumes `dis` directly and never converts `mrro`, so area is
not part of the asserted identity. The required static inputs are exactly the
four quantities in the residual: `width_m`, `bed_elevation_m`, `slope` and
`manning_n`.

The independent must-pass gate contains the bucket, both FLEX models and
SAC-SMA/Snow-17, as well as the exact `reference_uniform_flow`. The physical
models already compute their gauge from their own discharge and the supplied
Manning geometry. Their rating uses the wide-channel approximation; the
generated reach is wide enough that evaluating it with the exact rectangular
hydraulic radius remains inside the 5% allowance.

Two deliberately wrong references pin discrimination:

- `reference_wrong_roughness` reports a smooth, steady and monotone rating but
  builds it with `n = 0.08` instead of the case roughness;
- `reference_wrong_slope` uses `S_0 = 0.01` instead of the declared mild
  slope.

Both can conserve water, remain subcritical and look like ordinary ratings.
They fail because their discharge, stage and declared geometry cannot belong
to one steady momentum balance. Passing this probe makes only that necessary
steady-state claim; transient timing and hysteresis remain the work of
`routing-lag-consistency` and `stage-discharge-monotonic`.

## Source

The residual follows the Manning relation and rectangular-section definitions
in USGS Water-Resources Investigations Report 83-4251:
<https://doi.org/10.3133/wri834251>.

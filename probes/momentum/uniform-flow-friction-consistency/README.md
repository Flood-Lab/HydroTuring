# momentum/uniform-flow-friction-consistency

For a steady, uniform rectangular reach, the local Saint-Venant momentum
equation reduces to equality of bed slope and friction slope:

```text
d   = stage - bed_elevation_m
A   = width_m * d
R_h = A / (width_m + 2 d)
S_f = (manning_n * dis / (A R_h^(2/3)))^2.
```

`stage` is water-surface elevation measured from the same vertical datum as
`bed_elevation_m`; it is not water depth. The generator explicitly declares
`cross_section_shape: rectangular`, so the wetted perimeter and hydraulic
radius used above are unambiguous.

## Case and verdict

After a one-year low-flow spin-up, the forcing visits low, medium and high
constant-flow plateaus for one year each. The criterion scores the final 90
days of every plateau. Each scored block must satisfy all of the following:

- discharge and depth are finite and positive;
- the coefficients of variation of discharge and depth are at most 1%;
- the first-to-last-quarter change in either quantity is at most 1%;
- `mean(abs(S_f / S_0 - 1)) <= 0.05`.

Absolute residuals are averaged, so opposite-signed errors cannot cancel.
The report also records the 95th-percentile and maximum residual and the mean
and maximum Froude number at each plateau. The Froude number is diagnostic,
not an additional verdict in this probe.

The 5% allowance covers numerical convergence and the wide-channel stage
diagnostics of the existing physical baselines. Across the generated wide
rectangular sections, replacing the exact hydraulic radius with depth gives a
much smaller discrepancy than 5%. The two negative controls sit deliberately
just outside the boundary: a 3% high Manning coefficient produces
`1 - 1/1.03^2 = 5.74%` residual, while a 6% high slope produces 6% residual.
This pins discrimination near the stated tolerance instead of testing only
grossly broken ratings.

## Independent momentum baseline

`reference_saint_venant` is the must-pass model that tests the momentum
interpretation directly. It advances the conservative one-dimensional
shallow-water equations for depth and unit discharge with a Rusanov
finite-volume flux, explicit bed slope and an implicit Manning-friction
source. The upstream boundary prescribes discharge and the downstream
boundary is open. Every case starts from a weakly non-uniform motionless pool;
the solver neither evaluates nor inverts a normal-depth formula.

Tests require the numerical reach to:

- converge at low, medium and high discharge;
- reach the same equilibrium from different initial states;
- give the same interior gauge solution on 32- and 64-cell grids; and
- reproduce `S_f / S_0 = 1` independently after convergence.

The bucket, both FLEX models, SAC-SMA/Snow-17 and
`reference_uniform_flow` remain useful adapter/geometry checks, but the
Saint-Venant solver is what makes a passing gate depend on a model that
actually advances momentum.

The required `reach_length_m` sets the finite-volume cell length and therefore
the CFL time step and travel distance in that solver. `area_km2` is supplied
only so hydrologic baselines can convert runoff depth into their reported
discharge; the criterion consumes `dis` directly.

## Scope

Passing establishes one necessary steady-state momentum balance at three
flows. It does not establish transient wave timing or rating hysteresis; those
remain the jobs of `routing-lag-consistency` and
`stage-discharge-monotonic`.

## Sources

- MacDonald et al. (1997), *Analytic Benchmark Solutions for Open-Channel
  Flows*, <https://doi.org/10.1061/(ASCE)0733-9429(1997)123:11(1041)>.
- Delestre et al. (2013), *SWASHES: a compilation of shallow water analytic
  solutions for hydraulic and environmental studies*,
  <https://doi.org/10.1002/fld.1839>.
- Kurganov (2018), *Finite-volume schemes for shallow-water equations*,
  <https://doi.org/10.1017/S0962492918000028>.
- USGS Water-Resources Investigations Report 83-4251, rectangular-section
  Manning definitions, <https://doi.org/10.3133/wri834251>.

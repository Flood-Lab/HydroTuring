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
constant-flow plateaus for 924 days each. The slowest must-pass groundwater
store has a recession coefficient of 0.006 d^-1, so every plateau provides
`ceil(5 / 0.006) = 834` settling days before the final 90-day scored block.
Each scored block must satisfy all of the following:

- discharge and depth are finite and positive;
- the coefficients of variation of discharge and depth are at most 1.5%;
- the first-to-last-quarter change in either quantity is at most 1%;
- `mean(abs(S_f / S_0 - 1)) <= 0.05`.

The CV allowance includes a bounded daily numerical limit cycle without
mistaking it for a still-draining store. Across 500 generator seeds, the
largest such cycle in the packaged SAC-SMA/Snow-17 baseline was 1.47%, while
the separate quarter-shift gate still rejects a secular drift above 1%.

A plateau that fails the steadiness precondition is skipped rather than used
to judge momentum. Every steady plateau is still scored, so a non-steady
plateau cannot hide a friction violation measured on another plateau. The
probe is `N/A (INCOMPATIBLE)` only when none of the three plateaus is steady,
because only then has the experiment established no uniform-flow state in
which the friction balance can be evaluated. Reports retain the residual and
steadiness diagnostics for every skipped plateau.

All five requested seeds must produce a score (`min_scored_fraction: 1.0`).
Steadiness is measured on the model's own reported discharge and depth, so the
model decides which seeds it is scored on: any floor below 1 is a pass it can
buy by making the one seed it would lose on unscoreable. A ripple of 2.2% on a
single plateau is enough to do it. Over 260 generator seeds no packaged
baseline loses a plateau, so the floor costs them nothing.

The floor guards the pass and not the verdict. A friction violation measured on
a seed that was scored still decides the outcome however many other seeds went
`N/A`, because otherwise the same ripple would erase a violation instead of
buying a pass. Seed coverage is written into the criterion message and the
archive row either way.

Absolute residuals are averaged, so opposite-signed errors cannot cancel.
The report also records the 95th-percentile and maximum residual and the mean
and maximum Froude number at each plateau. The Froude number is diagnostic,
not an additional verdict in this probe.

The 5% allowance covers numerical convergence and the wide-channel stage
diagnostics of the existing physical baselines. Across 500,000 draws from the
generator, an honest wide-channel rating had a worst three-plateau residual
below 1.97%. The practical detection floor is therefore a Manning-roughness
error of about 2.5%: a 2% error can pass, whereas the exact-section negative
control at 3% high roughness has `1 - 1/1.03^2 = 5.74%` residual. The second
control uses a 6% high slope and produces 6% residual. These controls pin the
boundary without claiming sensitivity below what the hydraulic-radius
approximation permits.

Manning friction with the supplied `manning_n` is part of this probe's case
contract; the probe does not compare friction laws. A model using Chézy,
Darcy–Weisbach or a depth-dependent roughness should not declare that it
consumes this Manning parameter and is `N/A (INCOMPATIBLE)` here. A future
criterion using a model-reported section velocity could test a more
closure-independent momentum identity.

## Independent momentum baseline

`reference_saint_venant` supplies a numerically independent must-pass momentum
solution. It advances the conservative one-dimensional
shallow-water equations for depth and unit discharge with a Rusanov
finite-volume flux, explicit bed slope and an implicit Manning-friction
source. The upstream boundary prescribes discharge and the downstream
boundary is open. Every case starts from a weakly non-uniform motionless pool;
the solver neither evaluates nor inverts a normal-depth formula.

This is a pseudo-steady baseline, not a transient routing model: whenever the
inflow changes it integrates to equilibrium, memoizes that solution, and emits
the equilibrium gauge values over the plateau. Its Manning equilibrium follows
from the prescribed friction source. Its independence from the rating adapters
is in the finite-volume discretisation, boundary treatment and convergence from
a non-equilibrium state, not in using a different friction closure.

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

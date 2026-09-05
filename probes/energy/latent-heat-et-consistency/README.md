# energy/latent-heat-et-consistency

Ten years of generated daily weather and net radiation over a lumped
catchment. The question is not whether either budget closes. It is whether the
evaporation the model reports as water and the evaporation implied by the
latent heat flux it reports are the same evaporation.

```
LE = lambda(T, phase) * E
```

This is the first probe in the suite that no single-budget criterion can
substitute for. A model with a water head and an energy head can close the
water budget to floating point, close the energy budget to floating point, and
still be reporting two different evaporations, because the contradiction lives
strictly between the ledgers. `reference_two_head` is that model, and it is
what this probe exists to catch.

## Why the tolerance is not five percent

The suite's engineering rule is five percent of the driving flux. It is the
right rule for a budget and the wrong rule here, for a reason that is
arithmetic rather than a matter of taste.

`lambda_v(T) = 2.501e6 - 2361 T` J kg-1. Over a -20 to +35 degC record a model
that uses a constant instead is wrong by at most **3.4 percent**. Five percent
cannot see it. `flux_identity` is therefore scored at every step against
**0.5 percent** of the reported flux, with an absolute floor because the flux
passes through zero and a percentage of nothing is not a bound.

The floor is **0.5 W m-2**, not the 2 W m-2 that a first reading of the
nighttime problem suggests. At 2 W m-2 a constant-lambda model survives every
step whose evaporation is below roughly 3.5 mm/day, which on this record is
most of them: it trips 4 steps of 3650 rather than 33, and a different seed
could let it through entirely. That number came from running it.

Two remarks for anyone tempted to argue the tolerance is too tight. First, no
observation enters this probe: the reference budget is exact by construction
and the criterion measures a model against itself, so observational
uncertainty is not what is being bounded. Second, 0.5 percent is one to two
orders of magnitude below the spread between published ET products, which is
routinely 10 to 30 percent at annual scale. That gap is the point.

### The obvious objection, and why it argues the other way

*Observed surface energy budgets do not close. How can a benchmark demand
half a percent?*

They do not, and the size of the gap is well documented: energy balance
closure at eddy covariance sites varies systematically in space and time (Cui
& Chui 2019), and a recent correction across 172 flux towers moves the mean
imbalance from -14.99 to -0.65 W m-2 (Zhang et al. 2024) — that is, the raw
imbalance is of the order of 15 W m-2, roughly thirty times this probe's
floor.

That is precisely why the test belongs on the synthetic track and measures a
model against itself. If an observed budget were the reference, no tolerance
would be defensible, because the reference would be the thing that does not
close. ROADMAP already says this about the real-data track: "observed budgets
do not close, so observations can never be the reference". This probe is the
same argument applied to the energy budget, and it is the reason
`flux_identity` should never be ported to that track.

## The phase change, and what is not assumed

Water leaving a snowpack sublimates and takes `lambda_s = lambda_v(0) +
lambda_f`, **13.3 percent** more energy per kilogram. A model that converts
every kilogram at the vaporisation rate is short of that on exactly the steps
where a pack is losing mass to the atmosphere.

The sublimated mass is never assumed. Each step falls into one of three cases,
decided from the forcing and the model's own snow state:

| Case | What is asserted |
| --- | --- |
| No pack, and none can fall during the step | `LE = lambda_v(T) E` |
| Pack, no precipitation, air below freezing: nothing can fall and nothing can melt, so the pack's own reported mass loss is the sublimated mass | `LE = lambda_v L + lambda_s S` |
| Pack, anything else: melt and sublimation are not separable from what the model reports | `lambda_v(T) E <= LE <= lambda_s E` |

The third case is deliberately weak. Asserting a split the model never
reported would fail an honest model whose snow physics differs from the
reference's, and a criterion that does that is worse than no criterion. The
discrimination survives it: `reference_sublimation_blind` is caught on the
second case, which this record supplies about 400 steps of per seed.

A model that reports no snow state simply has the split disabled rather than
being judged against an assumption about it.

## The case

`generate.py`, from a recorded seed, ten years plus a year of spinup. The
water side is the record of `mass/catchment-closure`, unchanged, so a model
already adapted to that probe sees nothing unfamiliar. Two things are added.

Net radiation comes from a clear-sky cycle at 40 degN, attenuated by its own
AR(1) cloudiness process, with an albedo that switches below freezing and an
outgoing longwave term that follows air temperature. It crosses zero on winter
days, which is why the energy criteria carry a floor.

Potential evapotranspiration is then derived from that net radiation by
Priestley-Taylor with alpha = 1.26, rather than being an independent function
of temperature. A reviewer will ask why the demand and the available energy in
a coherence probe are consistent with one another, and this is the answer:
they are not two forcings, they are one. The consequence that matters is that
PET stays small but non-zero on cold clear days, without which no snowpack
would ever sublimate and half the criterion would have nothing to score.

Typical record: about 800 mm/yr of precipitation, 550 mm/yr of potential ET,
net radiation averaging 79 W m-2, and roughly 550 dry sub-zero days in 4015.

## Related work

Enforcing conservation inside a neural emulator, and quantifying the "physical
inconsistency" left when it is not enforced, is established for atmospheric
convection: Beucler et al. (2021) constrain neural networks emulating
convective processes either architecturally or through the loss, and measure
the violation across mass, momentum, radiation and energy budgets. That is the
nearest prior art to this probe and it should be read first.

Three things here are different. The verdict is binary rather than a penalty
term, because a benchmark has to be able to reject. The budgets are the
surface water and energy budgets of a catchment rather than a convective
column. And the sibling probe, `energy/evaporative-partition`, asserts a
cross-budget response to a counterfactual rather than a same-instant residual,
which is a property no amount of accuracy on either run establishes.

## What each reference model does here

| Model | Both budgets close? | Verdict |
| --- | --- | --- |
| `reference_coupled` | yes | PASS, every criterion |
| `reference_two_head` | yes, to 1e-15 | `flux_identity`, 3495 steps of 3650, implied lambda 3.13e6 J kg-1 |
| `reference_constant_lambda` | yes | `flux_identity` |
| `reference_sublimation_blind` | yes | `flux_identity`, sublimating steps only |
| `reference_energy_leak` | water only | `energy_closure` |
| `reference_bucket` | water only | INCOMPLETE: reports no energy fluxes |

That last row is the suite-level consequence of the first energy probe, and it
is the correct verdict rather than a regression. A water-only model has not
violated conservation of energy; it has declined to be falsifiable about it,
which is what INCOMPLETE has always meant here.

## References

- Beucler, T., Pritchard, M., Rasp, S., Ott, J., Baldi, P., & Gentine, P. (2021).
  Enforcing analytic constraints in neural networks emulating physical systems.
  *Physical Review Letters*, 126, 098302. https://doi.org/10.1103/PhysRevLett.126.098302
- Cui, W., & Chui, T. F. M. (2019). Temporal and spatial variations of energy balance
  closure across FLUXNET research sites. *Agricultural and Forest Meteorology*, 271, 12-21.
  https://doi.org/10.1016/j.agrformet.2019.02.026
- Zhang, W., Nelson, J. A., Miralles, D. G., Mauder, M., Migliavacca, M., Poyatos, R.,
  Reichstein, M., & Jung, M. (2024). A new post-hoc method to reduce the energy imbalance
  in eddy covariance measurements. *Geophysical Research Letters*, 51(2), e2023GL107084.
  https://doi.org/10.1029/2023GL107084
- Brutsaert, W. (1982). *Evaporation into the Atmosphere*. Reidel.

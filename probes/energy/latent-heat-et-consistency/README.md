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
most of them: on the first gate seed it trips 3 steps of 3650 rather than
993, and a different seed could let it through entirely. Those numbers came
from running it.

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

## The phase change, and why the criterion asks instead of guessing

Water leaving a snowpack sublimates and takes `lambda_s = lambda_v(0) +
lambda_f`, **13.3 percent** more energy per kilogram. A model that converts
every kilogram at the vaporisation rate is short of that on exactly the steps
where a pack is losing mass to the atmosphere.

Knowing which kilograms those were is the whole difficulty, and the first
version of this criterion got it wrong. It inferred the split: on a step with
no precipitation and air below freezing, nothing can fall and nothing can
melt, so any decrease in the reported snow store must be sublimation. That
reasoning has a hole. **A pack also loses water at its base.** Snow-17's
`DAYGM` term is exactly that, and this repository's own `sacsma_snow17` port
runs it at 0.1 mm/day.

Review reproduced the consequence. Taking `reference_coupled` unchanged except
for a constant ground melt into the soil, with the latent heat still exactly
right:

| Honest model | old criterion |
| --- | --- |
| no ground melt | pass, worst 0.00 of tolerance |
| ground melt 0.1 mm/day (this repository's Snow-17 setting) | pass, worst 0.77 |
| ground melt 0.3 mm/day | **fail on 375 to 429 steps of 3650**, worst 2.32 |

The probe's own stated principle is that a criterion which fails an honest
model whose snow physics differs from the reference's is worse than no
criterion. The inference did exactly that, with 23 percent of the tolerance to
spare at the repository's own default. So the criterion no longer infers.

**A model that reports `sbl`** — the sublimating share of its evaporation, a
component of `evspsbl` and never an addition to it — is held to the equality
at every step, pack or no pack, because it has said which kilograms left as
ice:

```
LE = lambda_v(T) * (E - sbl) + lambda_s * sbl
```

The claim has to be a real one: `sbl` may not be negative, may not exceed the
evaporation it is a share of, and must be zero where the model itself reports
no snow and none could fall. Without that last condition a warm-season model
could report a fictitious sublimating share to bend its effective lambda
upwards.

**A model that does not report `sbl`** is held to the equality at `lambda_v`
on steps where it reports no pack and none could arrive, and to the interval

```
lambda_v(T) * E  <=  LE  <=  lambda_s * E
```

wherever a pack is or could be present. Nothing is guessed.

The same ground-melt model under the criterion as it now stands, in both
postures:

| Ground melt | reports `sbl` | declines to report |
| --- | --- | --- |
| 0.0 mm/day | pass, worst 0.00 | pass |
| 0.1 mm/day | pass, worst 0.00 | pass |
| 0.3 mm/day | pass, worst 0.00 | pass |
| 1.0 mm/day | pass, worst 0.00 | pass |

**What this costs.** A model that declines to report `sbl` cannot be caught on
the phase change at all; the interval is wide enough to admit a
sublimation-blind model exactly on its lower bound. That is the honest price
of not guessing, and it is the right price: without the model saying which
kilograms were ice, blindness and a draining pack base are indistinguishable
from the outside. `reference_sublimation_blind` reports `sbl` and is caught on
every sublimating step whose flux is large enough for the missing 13.3 percent
to clear the floor, between 70 and 112 of roughly 400 sublimating steps on
each of ten seeds, rather than only the dry frozen ones, so the discrimination
is stronger where the information exists and absent where it does not.

## A choice the benchmark is making, stated plainly

`flux_identity` at 0.5 percent fails any model that uses a **constant** latent
heat of vaporisation, and `reference_constant_lambda` exists to make sure it
does. That is a position, not a conservation law, and it should be argued with
rather than discovered.

A constant lambda is a common convention in land-surface schemes, and a model
using one is internally consistent: its water and its energy describe the same
evaporation, at a lambda that is wrong by at most 3.4 percent over this
record. What the criterion asserts is a specific thermodynamic
parameterisation, `lambda_v(T) = 2.501e6 - 2361 T`, not merely that the two
budgets agree.

The benchmark chooses to call the constant a violation, for two reasons. The
identity is thermodynamics rather than a modelling choice, and 3.4 percent is
two orders of magnitude above what an exact model achieves here. And a
tolerance loose enough to admit a constant lambda is loose enough to admit a
model whose lambda is not a lambda at all, which is the failure the probe
exists for.

If the maintainers would rather this probe test coherence alone and leave the
temperature dependence to a separate criterion, the change is one parameter:
`lambda_vapour: [2.45e6, 0.0]` with `rel_tol` raised to the suite's 5 percent
turns it into that probe, and `reference_constant_lambda` would then belong to
the criterion that replaces it. This is worth settling before an LSM group
meets it in a report rather than after.

## The case

`generate.py`, from a recorded seed, ten years plus a year of spinup. The
water side is the record of `mass/catchment-closure`, unchanged, so a model
already adapted to that probe sees nothing unfamiliar. Two things are added.

Net radiation comes from a clear-sky cycle at 40 degN, attenuated by its own
AR(1) cloudiness process, with an albedo that switches below freezing and an
outgoing longwave term that follows air temperature. It averages under 20 W
m-2 in December and January and goes negative on a handful of overcast days,
which is why the energy criteria carry a floor.

Potential evapotranspiration is then derived from that net radiation by
Priestley-Taylor with alpha = 1.26, rather than being an independent function
of temperature. A reviewer will ask why the demand and the available energy in
a coherence probe are consistent with one another, and this is the answer:
they are not two forcings, they are one. The consequence that matters is that
PET stays small but non-zero on cold clear days, without which no snowpack
would ever sublimate and half the criterion would have nothing to score.

Typical record: about 850 mm/yr of precipitation, 775 mm/yr of potential ET,
net radiation averaging 90 W m-2, roughly 500 dry sub-zero days in 4015, and
about 400 sublimating steps in the 3650 that are scored.

### The seasonal cycle, corrected

The first release of the generator had the sign of the shortwave seasonal
term flipped: clear-sky radiation peaked at the winter solstice, at 378 W
m-2 in December against 112 in June, and Priestley-Taylor demand followed it,
270 mm over the five winter months against 154 over summer. Every annual
total in the paragraph above looked plausible, which is how it survived
review. The consequence was hydrological rather than energetic. The wet
season carried the year's strongest demand, the soil never reached capacity,
saturation excess all but vanished, and the exact model's runoff became 96
percent baseflow from a store that drains at 0.6 percent a day. Its
correlation with seven-day rainfall, which `non_degenerate` requires to be at
least 0.05, then depended on the draw: over 200 seeds it ranged from -0.017
upward and fell short on 13 percent of them, so a five-seed evaluation of the
exact model would have failed about half the time, and the gate was green only
because its five fixed seeds happened to pass.

The temptation was to switch the response test off for this probe, which the
criterion allows. That would have hidden a wrong case rather than fixed one.
With the sign corrected the water side is what it was designed to be, the
record of `mass/catchment-closure` under a demand of the same size and season,
and the exact model's rain response has the same distribution here that
`reference_bucket` has there: a minimum of 0.09 and a median of 0.20 over the
same 200 seeds, with no failures. Nothing about the tolerances above changed.
The step counts quoted for the reference models were re-measured on the
corrected case; the ground-melt table records the review of the earlier
criterion and stands as history.

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
| `reference_two_head` | yes, to 1e-15 | `flux_identity`, 3583 steps of 3650, implied lambda 3.31e6 J kg-1 |
| `reference_constant_lambda` | yes | `flux_identity` |
| `reference_sublimation_blind` | yes | `flux_identity`, sublimating steps only |
| `reference_energy_leak` | water only | `energy_closure` |
| `reference_bucket` | water only | INCOMPLETE: reports no energy fluxes |
| `flex_lumped`, `flex_topo`, `sacsma_snow17` | water only | INCOMPLETE, the same way |

Those last two rows are the suite-level consequence of the first probe that
needs energy fluxes, and they are the correct verdict rather than a
regression. A water-only model has not violated conservation of energy; it has
declined to be falsifiable about it, which is what INCOMPLETE has always meant
here. It is why `must_pass` names `reference_coupled` alone: every
hand-written physical model in the suite returns INCOMPLETE on this probe, so
none of them can guard its tolerance.

`min_window_days` is 3650 for the same family of reasons. A submitted model is
otherwise scored on a 30-day flood window, which on this record falls in the
spring thaw, where the phase change supplies between none and nine sublimating
steps instead of about four hundred, so the half of `flux_identity` that
exists to see it would have almost nothing to score. The expectations here
hold over the record, so the record is what a model is judged on.

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

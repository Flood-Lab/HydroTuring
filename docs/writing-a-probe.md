# Writing a probe

A probe is a conservation law made executable. It defines a generated case, a
set of binary criteria, and the reference models it must be able to separate.

## Propose it first

Open a [probe proposal](../../../issues/new?template=probe_proposal.yml)
before you build. The form asks how an unphysical model would pass your probe,
and it is much cheaper to discover that it would in a paragraph than in three
hundred lines. A maintainer labels the issue `accepted` and assigns it to you;
then fork the repository, work through the rest of this page, and open a pull
request from your fork when the gate is green.

One merged probe earns co-authorship on the benchmark paper, as do five
accepted model proposals. See [CONTRIBUTING.md](../CONTRIBUTING.md#credit).

## Start from a template

```bash
ht init-probe --list-templates           # what shapes are available
ht init-probe                            # writes probe-draft.yaml
#   ... fill in the fields ...
ht init-probe --from probe-draft.yaml    # creates probes/<law>/<slug>/
```

`--template <kind>` starts from a probe of that shape instead of the plain
one, and brings a matching generator skeleton with it.

| Template | The question it asks |
| --- | --- |
| `default` | does the budget close over one generated case |
| `extrapolation-space` | does it still close on catchments outside the range models are fitted to |
| `extrapolation-time` | does it still close under conditions outside anything earlier in the record |
| `counterfactual` | add water to the same case and ask where it went |
| `invariance` | change something the physics does not depend on and require nothing to move |

Each of these ships filled-in baselines and passes `ht gate` as scaffolded, so
you can run the gate before you have written a line of physics and see it
separate the reference models. Then replace the placeholder case with yours.

The draft is heavily commented. Lines beginning `#!` are guidance and are
stripped from the generated `probe.yaml`; your own `#` comments are kept.

What you get already validates and already honours the length contract, so
your first `ht gate` fails on your physics rather than on scaffolding.

```
probes/<law>/<slug>/
  probe.yaml     spec, schema-validated, carries your authorship
  generate.py    generate(seed) -> (DataFrame, dict)
  README.md      the physics in prose
```

`authors` in `probe.yaml` is the record the paper's author list is built
from, so fill it in as it should be printed: `name` and `affiliation`
(department and city, or "independent") for every author, and `orcid` where
you have one. The build checks that every merged probe's authors carry an
affiliation. See [CONTRIBUTING.md](../CONTRIBUTING.md#credit).

`probes/mass/catchment-closure/` is the reference implementation. Read it
before you start.

## The generator

```python
def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    ...
    return forcing, static
```

Requirements:

- deterministic given the seed, byte for byte
- a `time` column
- exactly the period plus the spinup, counted in rows at the probe's step:
  `period_years * 365 + spinup_days` rows for a daily probe, and for a
  sub-daily one `(period_days + spinup_days) / step`, where `period_days`
  may replace `period_years` in `case`
- fluxes as rates in mm per day whatever the step, so a minute of rain at
  one millimetre is a rate of 1440 mm/day
- no committed data files; CI rejects anything over 1 MB under `probes/`

Determinism is not a nicety. A benchmark whose failures cannot be reproduced
is unusable the first time a result is disputed.

## Criteria

Each entry in `criteria` names a registered criterion and its parameters.
Every one is binary.

| Criterion | Asserts | Scored over |
| --- | --- | --- |
| `closure` | the budget closes to within a share of the driving flux | one run |
| `event_water_closure` | every complete precipitation event satisfies `abs(R) <= max(threshold * P, absolute_tolerance_mm)`; defaults 0.05 and 0.001 mm; reports the worst residual / allowance against 1, with up to 20 failed events and summary percentiles | complete post-spinup wet events in one run, using supplied rain and all reported water stores |
| `state_bounds` | every reported storage stays physical | one run |
| `et_plausible` | ET is non-negative and bounded by potential ET | one run |
| `non_degenerate` | the partition and the response are non-trivial | one run |
| `forcing_fidelity` | the model reports back the forcing it was given | one run |
| `regime_transfer` | closure holds out of range as well as in range | labelled stretches |
| `counterfactual_response` | added or removed water is partitioned, not absorbed; `perturbed` may name one variant or a list, each scored against the control | paired runs |
| `invariance` | a transform the physics ignores changes nothing | paired runs |
| `resolution_invariance` | integrated volumes agree between the same weather at two steps | paired runs at different steps |
| `response_sign` | perturb one driver both ways, hold the rest: each response must point the way physics says, by a real share of the change in demand | paired runs |
| `causality` | nothing may change before an added storm, and runoff must answer it after | paired runs |
| `dry_down` | without rain, runoff and storages only fall, and no more drains than was held | one run, rainless record |
| `steady_state` | under constant weather everything settles, runoff stays below the rain, and the budget balances | one run, constant record |
| `monotone_response` | scaling a storm up a ladder cannot lower runoff, add more runoff than rain, or fail to run off most of an extreme | paired runs, a ladder |
| `runoff_bounds` | integrated runoff lies between rain minus demand minus storage and rain plus storage; needs runoff only | one run |
| `response_nonnegativity` | after an added storm the perturbed runoff is never below the control's, on any step | paired runs |
| `antecedent_monotonicity` | the same storm after more rain runs off more, and no more than the extra rain | paired runs |
| `phase_invariance` | the same water as rain instead of snow leaves the integrated volumes within a share of the rain | paired runs |
| `demand_consistency` | evaporation reaches demand when the model's own soil is wettest, never exceeds it, and falls when driest | one run |
| `energy_closure_by_phase` | the mean absolute surface-energy residual in each contiguous day or night stays within the larger of the relative radiation tolerance and the absolute flux floor | one run, contiguous day/night blocks |
| `radiative_identity` | upward longwave equals what the reported surface temperature emits plus the reflected downward longwave, at every step, within the larger of a relative tolerance and an absolute floor; emissivity comes from `static.json` | one run, instantaneous values |
| `soil_heat_storage` | interval boundary heat input agrees with fixed-layer temperature change and prescribed heat capacity | one run, separate heating/recovery phases |
| `melt_energy` | over a labelled melt block, the surface energy residual equals the fusion the reported ice change demanded plus what warming the pack cost; the ice is `snw - lwsnl` and the warming is `-d(csnow)`, both self-reported, so the two rows it differences carry seven contract checks, one of them read on every step | one run, labelled blocks |
| `routing_conservation` | the channel store is non-negative and never exceeds `max_lag_days` of the largest recent runoff | one run |
| `lag_time_bounds` | rainfall-to-runoff peak lag lies inside a broad duration-corrected Snyder envelope derived from public catchment geometry | paired runs, a geometry ladder |
| `scaling_monotonicity` | peak lag does not materially reverse and grows by a resolvable amount across catchment scales | paired runs, a geometry ladder |
| `rating_monotonic` | stage does not fall against its running maximum as the abscissa rises: equal-count bin medians are taken over the abscissa and the summed running-maximum deficit is compared with an explicit `tolerance` when one is given, and otherwise with 8% of the rating's span. The share is that large because a stage read off a store is hysteretic by construction, and its binned rating dips below its own running maximum by a visible fraction of the span for that reason alone | one run |
| `rating_loop` | where the gauge loops against the reach's store, the loop must be small enough to be noise or run the right way: at the same storage the rising limb sits lower than the falling one. A single-valued rating, or a loop below `min_loop_m` with an inconsistent sign across bins, is read as "no loop" and passes | one run |

Picking a denominator for `closure` and `regime_transfer`:

| Denominator | Use for | Floor |
| --- | --- | --- |
| `sum_pr` | water budgets | not needed, precipitation is non-negative |
| `sum_abs_rn` | energy budgets | required, net radiation crosses zero nightly |
| `sum_inflow` | routing | not needed |

### Labelled stretches

`regime_transfer` scores parts of one record separately and compares them,
which is how a failure confined to a rare stretch is kept from being diluted
by an ordinary one. The generator marks each step with a `_regime` column.

Any column whose name begins with `_` is an annotation: the criteria see it and
the harness strips it before staging, so the model never receives it. That is
not a detail. A column called `_regime` whose values are `ordinary` and
`anomaly` would otherwise tell a model exactly which part of the record it is
being judged on.

### Paired runs

`case.variants` names the cases the generator builds for one seed:

```yaml
case:
  generator: generate.py
  variants: [control, wetter]
```

The model is run once per variant and paired criteria are handed all the
results. The generator takes the variant as well as the seed:

```python
def generate(seed: int, variant: str = "control") -> tuple[pd.DataFrame, dict]:
```

A probe whose expectation only holds over a long enough stretch, such as the
sign of a response to warming, sets `case.min_window_days` and a submitted
model's evaluation window is widened to at least that.

`invariance` takes `unchanged` (must be reported and must not move),
`scaled` (must move by a factor) and `optional` (must not move *if the
model reports it*), so a store a model lacks is not an invariance failure.

The first variant is the control. Every non-paired criterion is scored against
it alone, so adding a variant to a probe never silently changes what its
existing criteria measure. Declaring variants without a paired criterion, or a
paired criterion without variants, is rejected at load time rather than
becoming a silent no-op.

Draw everything that comes from the seed before you branch on the variant. Two
variants that differ in the weather as well as in the perturbation cannot
isolate the perturbation, and the comparison means nothing.

### Variants at different steps

A resolution transform is the same weather at two or more steps. Declare
the step of any variant that does not run at `case.timestep`, and say how
the variants are selected for a model:

```yaml
case:
  timestep: PT1M
  period_days: 30
  spinup_days: 10
  variants: [minute, hourly, daily]
  timesteps:
    hourly: PT1H
    daily: PT1D
  variant_selection: native_and_finer
```

The control keeps `case.timestep`. Each variant's row count follows its own
step, the harness hands the model the step in `request.json`, and criteria
integrate each run with its own step. With `native_and_finer` a model is
run at the variant matching the step its manifest declares (or the nearest
coarser one) and at the next finer variant; the finest pairs with the next
coarser. The variant at the model's step is the control for single-run
criteria. Without it, every variant runs, control first.

Aggregate, do not redraw: a coarse variant must carry exactly the water of
the fine one, and `resolution_invariance` refuses a pair that does not. The
manifest's list of steps is not a gate on such a probe. How the model copes
with the other step is what the probe measures, and it reports the answer
as a share of the precipitation. See `probes/mass/resolution-invariance`.

## Baselines

```yaml
baselines:
  must_pass: [reference_bucket, flex_lumped, flex_topo, sacsma_snow17]
  must_fail:
    reference_leaky: closure
    reference_cheater: state_bounds
    reference_degenerate: non_degenerate
```

`must_pass` guards against tolerance drift: if a physical model ever fails
your probe, the probe is wrong until shown otherwise. Four are required:
the exact bucket, the two hand-written FLEX models from
chrimerss/HydrologicModels, which conserve water but partition it with the
nonlinearities a real conceptual model has, and the NWS's SAC-SMA with
Snow-17, the operational model, ported from its Fortran. They are what calibrates a
tolerance: the resolution probe's limit was moved from 5 to 10 percent when
FLEX-Topo, exactly conservative, moved 5.1 percent between an hourly and a
daily step because its partition answers intensity. A probe pull request is
run against these four before anything else. `must_fail` pins which criterion does
the catching, so a probe cannot appear to work while catching things for the
wrong reason.

A probe whose verdict rests on a forcing or static input that none of the four
consumes still needs a physical model in `must_pass`. An exact
domain-specific reference may sit beside it but cannot replace it: the
reference is written alongside the criterion and shares its assumptions, so a
gate that only the reference passes tests the criterion against a copy of
itself. Extend at least one
physical model's adapter to use the input, bump the model's version and
re-archive its rows, as `flex_lumped` maps channel geometry to its triangular
lag for `momentum/routing-lag-consistency`. Declare the input under `requires`;
the models that still cannot consume it are archived as N/A (INCOMPATIBLE). Do
not add an input declaration to a model whose adapter does not actually use it
merely to make the gate run.

On a probe with one case per seed, `must_fail` also decides what the report
says when a model passes: the `detail` column of `models/result.csv` names
each criterion it lists, with that criterion's own message, and no other. A
precondition no baseline is pinned to, such as `closure` declared ahead of
the identity it protects, is therefore never reported as the probe's result.
A probe with `variants` reports its paired criteria instead, since
everything beside them is a precondition checked on the control.
`tests/test_report_detail.py` pins what every merged probe reports.

The reference models available today:

| Model | What it does | Caught by |
| --- | --- | --- |
| `reference_bucket` | conserves water exactly by construction | nothing, it must pass |
| `flex_lumped` | lumped FLEX/HBV from chrimerss/HydrologicModels; conservative, nonlinear partition | nothing, it must pass |
| `flex_topo` | FLEX-Topo, three landscape units sharing a groundwater store | nothing, it must pass |
| `sacsma_snow17` | SAC-SMA + Snow-17 + gamma unit hydrograph, the NWS operational model | nothing, it must pass |
| `reference_leaky` | hides a silent 15% sink | `closure`, `counterfactual_response` |
| `reference_cheater` | solves for storage as whatever balances the budget; runoff is a fixed share of rain | `state_bounds`, `response_sign`, `counterfactual_response` |
| `reference_degenerate` | evaporates all precipitation, produces no runoff | `non_degenerate`, `counterfactual_response`, `response_sign` |
| `reference_in_sample` | exact in range, leaks outside it | `regime_transfer` |
| `reference_calendar` | recession drifts with the calendar year | `invariance` |
| `reference_fixed_step` | treats every row as a day whatever the step is | `resolution_invariance` |
| `reference_streamflow_only` | reports runoff only, from a store that never reads the temperature | `N/A (INCOMPLETE)` on budget probes; `response_sign` |
| `reference_anticipating` | reports runoff smoothed over a centred window, three days of the future in every value | `causality` |
| `reference_climatology` | the seasonal mean, whatever the weather; never reads the rain | `dry_down` |
| `reference_saturating` | daily runoff capped at 25 mm; flat beyond its training range | `monotone_response` |
| `reference_restless` | a recession with its own thirty-day clock; never settles | `steady_state` |
| `reference_overflowing` | reports its runoff plus 80% of the rain again | `runoff_bounds` |
| `reference_area_leak` | loses a share of runoff that grows with the stated area | `invariance` (area) |
| `reference_overshooting` | a derivative term sharpens its hydrograph | `response_nonnegativity` |
| `reference_sublimating` | loses 40% of every snowfall unreported | `phase_invariance` |
| `reference_thirsty` | evaporates a fixed share of its soil store whatever the demand | `demand_consistency` |
| `reference_stuck_router` | a routing kernel summing to 0.9 | `routing_conservation` |
| `reference_snyder_router` | consumes public catchment geometry and routes rain with a conservative triangular unit hydrograph whose peak follows the duration-corrected Snyder lag from the excess-rainfall centroid | nothing, it must pass the routing-lag probe |
| `reference_instant_router` | accepts the geometry but returns runoff in the rainfall row at every scale | `lag_time_bounds` |
| `reference_inverse_router` | uses individually plausible lags that reverse once as catchment scale grows | `scaling_monotonicity` |
| `reference_rating` | the bucket with a real rating curve: yield enters a shallow floodplain and a deep channel reservoir, and the stage is the depth the channel's volume makes in a fixed bed | must pass `momentum/stage-discharge-monotonic` |
| `reference_rating_drift` | derives its stage from a slowly decaying running maximum of discharge (`peak = max(q, 0.997 * peak)` per day), so the gauge ratchets up with each flood far faster than it relaxes | `rating_monotonic` |
| `reference_rating_inverted` | reads the loop backwards, high while the flood is arriving and low once it is leaving | `rating_loop` |
| `reference_flat_stage` | reports a constant stage, so there is no rating and no loop | `non_degenerate` |
| `reference_coupled` | the bucket with snow sublimation and a surface energy budget; every kilogram converted at the latent heat of the phase it actually underwent | nothing, it must pass the energy probes |
| `reference_soil_heat` | a synthetic fixed-layer fixture with conductive boundary fluxes and temperature integrated consistently | nothing, it must pass `soil_heat_storage` |
| `reference_frozen_soil` | keeps the conductive fluxes but reports a constant soil temperature | `soil_heat_storage` |
| `reference_half_soil` | keeps the conductive fluxes but halves the reported temperature change | `soil_heat_storage` |
| `reference_snow_energy` | an energy-balance snowpack: melt bought from its own surface budget, liquid held in the pore space and reported as `lwsnl`, cold content carried and reported as `csnow` | nothing, it must pass `energy/snowmelt-energy-water` |
| `reference_degree_day` | the same pack melted on air temperature, with an energy budget that closes around its evaporation alone | `melt_energy` |
| `reference_warming_free` | melts on the energy available and reports its cold content, but charges the budget for the fusion alone | `melt_energy` |
| `reference_two_head` | a water head and an energy head that never meet; both budgets close and the latent heat implies an evaporation it never reported | `flux_identity`, `partition_shift` |
| `reference_constant_lambda` | converts every kilogram at one latent heat of vaporisation | `flux_identity` |
| `reference_sublimation_blind` | converts snow sublimation at the latent heat of vaporisation instead of sublimation | `flux_identity` |
| `reference_energy_leak` | discards 15% of net radiation | `energy_closure` |
| `reference_ground_dodge` | coherent and closed, but its sensible flux never reads the soil, so the ground flux absorbs a drydown's whole shift | `partition_shift` |
| `reference_diurnal_bias` | shifts sensible heat to leave opposite day and night energy residuals that cancel over the full record | `energy_closure_by_phase` |
| `reference_radiative` | the coupled reference with a skin: temperature from the sensible flux through a fixed conductance, upward longwave from that skin at the given emissivity | nothing, it must pass the radiation probe |
| `reference_air_emitter` | reports the skin's temperature but emits at the air's, reflected sky unchanged | `radiative_identity` |
| `reference_no_reflection` | reports emission alone as the total upward longwave, the reflected sky left out | `radiative_identity` |

The three soil-heat references require incoming `rsds` and `rlds`, `tas`,
`pr`, and explicit layer depth, areal heat capacity and initial temperature.
They have no prescribed-`rn` input path or default thermal layer. Their
manifests therefore mark cases without these inputs `N/A (INCOMPATIBLE)`.
`reference_soil_heat` is a synthetic fixture; the optional
[native Noah-MP validation](noahmp-soil-heat-validation.md) supplies separate
physical-model evidence.

If your probe needs a broken model that does not exist yet, add it under
`models/` alongside the probe. A criterion with nothing that trips it is
untested.

## Before opening a PR

The pull request comes from your fork and closes the proposal issue. Title it
`[PROBE: <law>] <description>`. Review is two passes, one on the physics and
one on the implementation.

```bash
ht validate
ht gate --probe <law>/<slug>
pytest -q
```

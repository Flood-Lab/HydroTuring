<p align="center">
  <img src="res/HydroTuring.png" alt="HydroTuring Initiative" width="540">
</p>

<p align="center">
  <b><a href="https://flood-lab.github.io/HydroTuring/">flood-lab.github.io/HydroTuring</a></b>
  &middot; English, Español, 中文
</p>

<p align="center">
  <a href="#the-probes"><img alt="probes merged" src="https://img.shields.io/endpoint?url=https%3A%2F%2Fflood-lab.github.io%2FHydroTuring%2Fbadges%2Fprobes.json"></a>
  <a href="#models"><img alt="models evaluated" src="https://img.shields.io/endpoint?url=https%3A%2F%2Fflood-lab.github.io%2FHydroTuring%2Fbadges%2Fmodels.json"></a>
  <a href="CONTRIBUTORS.md"><img alt="contributors" src="https://img.shields.io/github/contributors/Flood-Lab/HydroTuring?color=1f6f8b&label=contributors"></a>
  <a href="LICENSE"><img alt="licence" src="https://img.shields.io/badge/licence-PolyForm%20Noncommercial%201.0.0-555"></a>
  <a href="https://discord.gg/7SQb6bUZD"><img alt="Discord" src="https://img.shields.io/badge/Discord-join%20the%20community-5865F2?logo=discord&logoColor=white"></a>
</p>

A benchmark that asks one question of any AI hydrologic model: **does it
conserve what physics says it must conserve?**

Not whether it fits a hydrograph. Whether its water budget closes, its energy
budget closes, and its routing conserves momentum. A model passes HydroTuring
only when every criterion of every probe passes.

> **Bring a probe or a model, join the paper.** The suite is only as good as
> the physics people bring to it, and only as interesting as what has been put
> through it. Anyone may propose either — you do not need to be invited,
> affiliated, or known to us.
>
> **One merged probe earns co-authorship on the HydroTuring paper. So do five
> accepted model proposals.** You do not have to be able to package a model to
> propose it; say so on the form and the work is assigned to someone who can,
> and it still counts as yours.
>
> Start with a [probe](ROADMAP.md#probes-we-want) or a
> [model](../../issues/new?template=model_submission.yml), or come talk it
> over on [Discord](https://discord.gg/7SQb6bUZD) first. Terms in
> [CONTRIBUTING.md](CONTRIBUTING.md#credit).

```
$ ht run --model reference_bucket
reference_bucket v1.0.0  ->  ✅ PASS (OK)  [1/1 probes passed]
  ✅ PASS  mass/catchment-closure
        ✅  closure            cumulative residual 0.0000% of sum_pr (limit 5.0%)
        ✅  state_bounds       all storages stay physical
        ✅  et_plausible       cumulative ET is 0.664 of potential ET (limit 1)
        ✅  non_degenerate     partition and variability are non-trivial
        ✅  forcing_fidelity   reported forcing matches the input
```

<p align="center">
  <img src="res/hydroturing-infra.svg" alt="Any AI hydrologic model enters through the /io contract, HydroTuring probes it against mass, energy and momentum, and a verdict with its reason comes out" width="800">
</p>

## The idea

A benchmark that only checks closure is trivially gamed. A set of deliberately
broken models lives in this repository to prove it, and no probe is merged
until it passes an exactly conservative model and catches the broken ones it
names.

`reference_cheater` is the one worth dwelling on. Its closure residual is
**exactly zero on every seed, forever**. Randomising the forcing cannot touch
it, because the cheat is in its internal wiring rather than its memory. It is
caught only because the model contract requires absolute storage states rather
than tendencies, so the storage it invents has to stay physical, and it does
not.

## The probes

Twenty-two: sixteen under mass, five under energy and one under momentum.
Each was merged only after the acceptance gate saw it pass its declared
exact reference and fail a purpose-built broken one on the named criterion.
Four physical models, a bucket that conserves water exactly, two
hand-written FLEX models and the NWS's SAC-SMA with Snow-17, must pass
every probe that can ask them anything; four of the five energy probes need
outputs they do not report and are not scored for them. A probe that fails a
physical model is examined before the model is; that is the first thing done
with any probe pull request. Eleven of the twenty-two can be scored on a model
that reports runoff and nothing else. `ht list` prints them;
[ROADMAP.md](ROADMAP.md#probes-we-want) has the ten more we want, all
unclaimed.

| Probe | Law | What it asks | The broken model it catches |
| --- | --- | --- | --- |
| [`mass/catchment-closure`](probes/mass/catchment-closure) | mass | Does the water budget close over ten generated years? | `reference_leaky`, `reference_cheater`, `reference_degenerate` |
| [`mass/resolution-invariance`](probes/mass/resolution-invariance) | mass | The same month at the minute, the hour and the day: do the integrated volumes agree? | `reference_fixed_step`, `reference_degenerate` |
| [`mass/warming-response`](probes/mass/warming-response) | mass | The same rain with the air 3 °C warmer and 3 °C cooler: does runoff move the way physics says, in both directions? | `reference_degenerate`, `reference_streamflow_only` |
| [`mass/causality`](probes/mass/causality) | mass | One storm added mid-record: nothing may change before it, and runoff must answer after it. | `reference_anticipating` |
| [`mass/dry-down`](probes/mass/dry-down) | mass | Two years without rain: runoff can only fall, and no more may drain than the catchment held. | `reference_climatology` |
| [`mass/steady-state`](probes/mass/steady-state) | mass | Three years of the same day: does everything settle, and does the budget balance once it has? | `reference_restless` |
| [`mass/extreme-rain`](probes/mass/extreme-rain) | mass | The largest storm scaled up to ten times: runoff may not fall, nor exceed the rain that was added. | `reference_saturating` |
| [`mass/extreme-event-closure`](probes/mass/extreme-event-closure) | mass | Overlap rainfall events toward fitted 100-year depths in one median-wet year of twenty; check each event's water budget so long-record averaging cannot hide a loss. | `reference_in_sample` |
| [`mass/runoff-bounds`](probes/mass/runoff-bounds) | mass | Over ten years, is the runoff possible at all: at least rain minus demand minus storage, at most rain plus storage? The mass question a runoff-only model has to answer. | `reference_degenerate`, `reference_overflowing` |
| [`mass/area-invariance`](probes/mass/area-invariance) | mass | The same weather on the same catchment told as ten times larger: every depth must be identical. | `reference_area_leak` |
| [`mass/response-nonnegativity`](probes/mass/response-nonnegativity) | mass | One 120 mm storm added: from that day on, runoff may never be lower than without it. | `reference_overshooting` |
| [`mass/antecedent-monotonicity`](probes/mass/antecedent-monotonicity) | mass | The same storm after a dry month and a wet one: the wetter catchment runs off more, and no more than the extra water. | `reference_cheater` |
| [`mass/phase-counterfactual`](probes/mass/phase-counterfactual) | mass | The same water falling as rain instead of snow: timing moves, the integrated volumes may not. | `reference_sublimating` |
| [`mass/time-origin-invariance`](probes/mass/time-origin-invariance) | mass | The same weather under a 28-year calendar shift that preserves seasons and leap days: evaporation, runoff and water stores must agree. | `reference_calendar`, `reference_degenerate`, `reference_leaky` |
| [`mass/precipitation-counterfactual`](probes/mass/precipitation-counterfactual) | mass | The same seed 20% wetter, 10% wetter and 20% drier: the water added or removed must be partitioned among evaporation, runoff and storage, and runoff must rise from drier to wetter. | `reference_cheater`, `reference_leaky`, `reference_degenerate` |
| [`mass/human-abstraction`](probes/mass/human-abstraction) | mass | A prescribed net irrigation withdrawal must leave the budget: the same weather run with and without it, and the difference must account for exactly the abstracted volume. | `reference_abstraction_blind`, `reference_leaky` |
| [`energy/pet-consistency`](probes/energy/pet-consistency) | energy | Evaporation reaches demand when the model's own soil is wettest, stays below it, and falls when the soil is driest. | `reference_thirsty` |
| [`energy/latent-heat-et-consistency`](probes/energy/latent-heat-et-consistency) | energy | The evaporation a model reports as water and the evaporation implied by the latent heat it reports: are they the same evaporation? | `reference_two_head`, `reference_constant_lambda`, `reference_sublimation_blind`, `reference_energy_leak` |
| [`energy/evaporative-partition`](probes/energy/evaporative-partition) | energy | One summer without rain under net radiation that did not change: the latent heat a drying surface gives up has to warm the air. | `reference_two_head`, `reference_ground_dodge` |
| [`energy/surface-energy-closure`](probes/energy/surface-energy-closure) | energy | Does each day and night close its hourly surface energy budget, without opposite errors cancelling? | `reference_diurnal_bias` |
| [`energy/radiation-consistency`](probes/energy/radiation-consistency) | energy | The surface temperature a model reports and the upward longwave it reports: do they describe one surface, hour by hour, at the emissivity it was given? | `reference_air_emitter`, `reference_no_reflection` |
| [`momentum/routing-conservation`](probes/momentum/routing-conservation) | momentum | The channel store is never negative and never holds more than its hydrograph can. | `reference_stuck_router` |

## Models

`models/` holds three kinds. A submitted model is there to be evaluated, and
every run of one is appended to [models/result.csv](models/result.csv); the
physical models' runs are archived there too. A
physical model is there to test the probes: the exact bucket, two
hand-written conceptual models from
[chrimerss/HydrologicModels](https://github.com/chrimerss/HydrologicModels)
and the NWS's SAC-SMA with Snow-17 must pass every probe that can ask them
anything, so a probe that fails one is wrong until shown otherwise. All four
report water and no energy, so all four are N/A, with reason INCOMPLETE, on
the four probes that need an energy output: the three that need the latent,
sensible and ground heat fluxes, `energy/latent-heat-et-consistency`,
`energy/evaporative-partition` and `energy/surface-energy-closure`, and
`energy/radiation-consistency`, which needs a surface temperature and its
upward longwave. Not scored rather than passing: they have not
violated conservation of energy, they have declined to be falsifiable about
it, exactly as `reference_streamflow_only` does on the budget probes. A broken model is broken in one specific way, so that no criterion
goes untested.

A standing counts the probes a model passed out of the probes that could score
it. A probe that needs a variable the model does not report, or that the model
cannot consume, is N/A for it and in neither number, so the totals differ
between models.

| Model | Kind | What it does | Standing |
| --- | --- | --- | --- |
| [`google_flood_forecast`](models/google_flood_forecast) | submitted | The mean-embedding forecast LSTM behind Google Flood Hub, at the published weights. Predicts discharge and nothing else. | **FAIL (VIOLATION)**, 5 of 11 probes passed. It passes the runoff bounds, memory, area, causality and steady state, and fails step, extreme rain, phase, warming, dry-down, and a 0.18 mm/day dip after an added storm |
| [`dhbv2`](models/dhbv2) | submitted | δHBV 2.0, the MHPI group's differentiable HBV: neural networks write the parameters of a bucket model that reports its stores and its evaporation. | **FAIL (VIOLATION)**, 11 of 18 probes passed. Its learned regional-groundwater term, declared as `gwex`, closes the budget to 1e-8; what remains is a learned field capacity twice the catchment's, a response to doubled rain above the rain added, a runoff depth that changes with the area it is told, a third more runoff when snow falls as rain, and on `mass/human-abstraction` it never reads the prescribed withdrawal, so it accounts for none of the 380 mm and trips that probe's `state_bounds` on the same learned field capacity |
| [`wflow_sbm`](models/wflow_sbm) | submitted | Deltares' Wflow.jl SBM: a soil column with unsaturated and saturated stores, interception, snow and kinematic-wave routing, adapted in Julia and run on one representative cell at the resolution of Wflow's Moselle model. | **FAIL (VIOLATION)**, 15 of 18 probes passed. Its water budget closes to 2e-4 of the rain and it passes the area, routing, step, calendar, causality, extreme-rain, phase, precipitation-counterfactual, warming and human-abstraction probes, the last by taking the prescribed withdrawal through Wflow's own water demand and allocation. That 2e-4 is water its river kinematic wave creates when open-water evaporation exceeds the rain on the river and the wave floors its discharge rather than drying the channel; on `mass/extreme-event-closure` it is 0.0016 mm on one 0.006 mm drizzle day on one seed, past that probe's 0.001 mm floor, the only one of the probe's events to fail. Its other two failures both come from its soil column: a storm that does not fill the column makes almost no runoff, and under the shipped mapping a wet month's extra water has already been evaporated down to the rooting depth when the storm arrives, so the wetter catchment runs off barely more than the dry one, about a thousandth of the storm where the probe asks for two hundredths; and once the root zone dries in a rainless spell, drainage from the unsaturated store raises a water table that sits below the roots, and lateral flow rises with no rain. Both move with how the stated soil capacity is mapped onto soil thickness and roots. |
| [`summa`](models/summa) | submitted | SUMMA 4.0.0, a land model that solves the water and energy balances of canopy, snow, soil and aquifer with one implicit solver, run as one lumped HRU from its shipped test-case setup, with the radiation and humidity it needs mocked from each forcing row. | **FAIL (VIOLATION)**, 8 of 21 probes passed. Its water budget closes to rounding and its surface energy budget to 0.2% of its own net radiation; what fails is a latent heat of vaporisation held at its 0 C value, a net radiation of its own that is not the probe's `rn` (every hourly phase, and a drought that warms the surface), a canopy that, as the shipped setup configures SUMMA, intercepts all rain and keeps what freezes near 0 C, holding up to 27 mm of ice against a 2 mm capacity on the probes that score it (also the only failure on the precipitation counterfactual), leaf area that keeps its seasons under constant weather, runoff that follows the melt rather than the rain on one seed, a 5.2% runoff response to turning snow into rain on another, Green-Ampt runoff that appears only at the hourly step, and no human water use, so the prescribed withdrawal on the human-abstraction probe is never taken |
| [`cwatm`](models/cwatm) | submitted | CWatM 1.11, IIASA's Community Water Model and an ISIMIP global hydrological model, run on one grid cell with every store it carries reported. | **FAIL (VIOLATION)**, 14 of 18 probes passed. Its budget closes to 0.03%, and its own water-demand module pumps a prescribed withdrawal out of groundwater to within 0.004%. That 0.03% is water its capillary rise creates, a few thousandths of a millimetre a day, and on `mass/extreme-event-closure` it fails 3 to 8 one-day drizzle events of under 0.07 mm per seed, where 0.001 to 0.005 mm more leaves or is stored than fell, against that probe's allowance of 5% of the rain or 0.001 mm, whichever is larger. It also fails on a groundwater reservoir with no dt that drains 24 times too fast at an hourly step (52% of the rain); on a 0.29 mm/day dip after an added storm, where preferential flow turns surface runoff into interflow that leaves through a slower runoff-concentration lag; and on evaporation on wet soil at 0.70 of demand, near the cap its crop coefficients set. The last two are packaging choices as much as results: this package follows the CWatM-Earth-30min template its parameters come from, and with `preferentialFlow = False`, the setting of the pinned model repository's own 30′ templates, both pass |
| [`lisflood`](models/lisflood) | submitted | LISFLOOD 5.0.0, the EC Joint Research Centre's distributed model behind EFAS and GloFAS, stepped through its own Python framework on one representative 5 km cell and reporting every store its own water balance module counts. | **FAIL (ERROR)**, 15 of 18 probes passed. Its budget closes to 1e-13 mm per step; it reports no heat fluxes and no surface temperature, so the four energy probes that need an energy output cannot ask it anything and are N/A. It fails the step probe because its potential infiltration is a pore-space storage multiplied by the step length, so rain falling within hours runs off at the hourly step (13.0% of the rain between PT1H and PT1D). The two ten-year probes take about 90 s per run under amd64 emulation against a 60 s budget, so their rows are ERROR from the host's speed, and those two ERROR rows alone make the verdict FAIL (ERROR). Run outside the limit, `mass/precipitation-counterfactual` passes every criterion. On `mass/human-abstraction`, LISFLOOD's own water-use rule takes the groundwater share in full and the rest only from channel water above an environmental-flow reserve, recording what the channel cannot give as shortage. On this one-cell water region it withdraws 17 to 33% of the prescription, from the sourced reserve to none, and leaves 67 to 83% where the probe allows 5% |
| `reference_bucket` | exact | conserves water exactly by construction | must pass every probe that can ask it anything; N/A on the four energy probes that need outputs it does not report |
| [`flex_lumped`](models/flex_lumped) | physical | lumped FLEX/HBV: interception, beta-partitioned unsaturated store, fast and slow reservoirs, triangular lag | must pass every probe that can ask it anything; **PASS**, 18 of 18 |
| [`flex_topo`](models/flex_topo) | physical | FLEX-Topo: plateau, hillslope and wetland units on real Wark fractions sharing one groundwater store | must pass every probe that can ask it anything; **PASS**, 18 of 18 |
| [`sacsma_snow17`](models/sacsma_snow17) | physical | the NWS's SAC-SMA with Snow-17 and a gamma unit hydrograph, ported from the legacy Fortran and checked against it | must pass every probe that can ask it anything; **PASS**, 18 of 18 |
| `reference_coupled` | exact | the bucket with snow sublimation and a surface energy budget: every kilogram converted at the latent heat of the phase it actually underwent | must pass every criterion of the three energy-flux probes; supports daily and hourly steps |
| `reference_abstraction_blind` | broken | the same bucket, blind to the prescribed withdrawal, so the two variants come out identical | caught by `human_abstraction` |
| `reference_two_head` | broken | a water head and an energy head that never meet: both budgets close to 1e-15 and the latent heat implies an evaporation it never reported | caught by `flux_identity` |
| `reference_constant_lambda` | broken | coherent, but converts every kilogram at the same latent heat of vaporisation | caught by `flux_identity` |
| `reference_sublimation_blind` | broken | coherent for liquid water, but converts snow sublimation at the latent heat of vaporisation instead of sublimation | caught by `flux_identity` |
| `reference_ground_dodge` | broken | coherent, both budgets close, and the sensible flux never reads the soil: when the soil dries the ground flux absorbs the whole shift | caught by `partition_shift` |
| `reference_diurnal_bias` | broken | shifts sensible heat so the energy residual is +20 W/m2 by day and -20 W/m2 by night; the full-record residual cancels | caught by `energy_closure_by_phase` |
| `reference_radiative` | exact | the coupled reference with a skin: temperature diagnosed from the sensible heat flux through a fixed conductance, upward longwave that skin's emission plus the reflected sky at the emissivity it was given | must pass every criterion of `energy/radiation-consistency` |
| `reference_air_emitter` | broken | reports the skin's temperature but emits at the air's, reflected sky unchanged | caught by `radiative_identity` |
| `reference_no_reflection` | broken | reports emission alone as the total upward longwave, the reflected sky left out: 3 to 23 W/m2 that a 5% tolerance would not reliably see | caught by `radiative_identity` |
| `reference_energy_leak` | broken | discards 15% of net radiation; the energy counterpart of `reference_leaky` | caught by `energy_closure` |
| `reference_leaky` | broken | hides a silent 15% sink | caught by `closure` |
| `reference_cheater` | broken | solves for storage as whatever balances the budget | caught by `state_bounds` |
| `reference_degenerate` | broken | evaporates all precipitation, produces no runoff | caught by `non_degenerate`, `response_sign` |
| `reference_fixed_step` | broken | treats every row as a day whatever the step is | caught by `resolution_invariance` |
| `reference_anticipating` | broken | smooths runoff over a centred window, so three days of the future are in every value | caught by `causality` |
| `reference_climatology` | broken | emits the seasonal mean whatever falls, and keeps flowing without rain | caught by `dry_down` |
| `reference_saturating` | broken | caps its daily runoff, so an extreme storm adds rain and no runoff | caught by `monotone_response` |
| `reference_restless` | broken | a recession on an internal thirty-day clock, so it never settles | caught by `steady_state` |
| `reference_overflowing` | broken | reports its runoff plus 80% of the rain again, from nowhere | caught by `runoff_bounds` |
| `reference_area_leak` | broken | loses a share of runoff that grows with the area it is told | caught by `invariance` (area) |
| `reference_overshooting` | broken | a derivative term sharpens its hydrograph, so an added storm lowers later flow | caught by `response_nonnegativity` |
| `reference_sublimating` | broken | loses 40% of every snowfall to an unreported sublimation | caught by `phase_invariance` |
| `reference_thirsty` | broken | evaporates a fixed share of its soil store, never reading demand; conserves water exactly | caught by `demand_consistency` |
| `reference_stuck_router` | broken | a routing kernel summing to 0.9, so a tenth of every day's runoff never leaves the channel | caught by `routing_conservation` |
| `reference_streamflow_only` | honest limit | reports discharge only, from a store that never reads the temperature | N/A (INCOMPLETE) on budget probes; caught by `response_sign` |
| `reference_in_sample` | broken | removes surface runoff above a fixed 55 mm daily precipitation cutoff | caught by `event_water_closure` |
| `reference_calendar` | broken | a recession that drifts with the calendar year | caught by `invariance` (time origin) |

## How a case is generated

Forcing is generated fresh at run time from a recorded seed. Nothing is
committed as data, so there is nothing to memorise, and a reviewer reviews a
short deterministic script instead of a binary blob. Each probe runs several
seeds and all of them must pass, so no model gets through on a lucky draw.
Any run reproduces exactly with `ht run --seed <n>`.

This is an **unseen-sample** guarantee, not a claim that the public generator's
distribution is secret. For evaluation against generators or data unavailable
during training, keep a second probe tree outside the repository and add it
with `--probe-root /secure/hidden-probes`. The model container receives only
opaque case metadata and the inputs needed for inference; the host retains the
probe identity, generator seed, annotations and scoring code. See
[docs/adapting-a-model.md](docs/adapting-a-model.md#private-evaluation-suites).

## Verdicts

A probe that can ask the model something passes or fails it, and the reason
is recorded separately from the verdict, because these mean different things:

- `VIOLATION` the model reported its budget and the budget did not close.
- `ERROR` the adapter or benchmark machinery failed. This is operational, not
  a scientific verdict, and is the one outcome that makes `ht run` exit 2.

A probe that cannot ask the model anything is `N/A`: not scored, and neither
a pass nor a fail. That happens two ways:

- `INCOMPLETE` the model never reported enough to be checked. Every
  streamflow-only model lands here on the budget probes. It has not violated
  conservation; it has declined to be falsifiable.
- `INCOMPATIBLE` the model and probe disagree on timestep, on a forcing or
  static input one needs and the other does not supply or declare, or on
  paired-perturbation support, so running them would not be meaningful.

A model passes when at least one probe could be put to it and every probe
that could be put to it passes. Otherwise it fails, with the worst reason
among the scored probes that did not pass; an unscored probe never supplies
that reason. A model that no probe could be put to at all is `N/A` too, and
`ht run` exits 1 for it as it does for a FAIL.

A scored verdict is one bit. Everything under it stays quantitative, so a
paper can show that one model leaks 6% and another 40% long before anyone
crosses the line.

## Quick start

```bash
pip install -e '.[dev]'

ht init-probe --list-templates            # probe shapes to start from
ht init-probe --template invariance      # start a new probe
ht list                                  # probes and models
ht validate                              # schema-check everything
ht gate                                  # the probe acceptance gate
ht run --model reference_bucket          # evaluate one model
ht run --model my-model --json out.json  # machine-readable report
```

Without installing, `./ht` (or `ht.cmd` in a Windows shell) runs the CLI
straight from `src/`, as does `python -m hydroturing` with `src` on
`PYTHONPATH`. Set `HT_ASCII=1`
for reports with the words and no marks, which is also what you get
automatically wherever the output stream cannot carry them.

## Submitting a model

Your model may be written in any language. It ships as a container plus a
thin adapter that reads `/io/request.json` and writes `/io/output/result.csv`.
The adapter is usually thirty lines. See [docs/adapting-a-model.md](docs/adapting-a-model.md),
and `AGENTS.md` if you are having a coding agent build the sandbox for you.

**We want more models, and we want the ones you think matter.** Every
published rainfall-runoff model, every LSTM, every foundation model with a
hydrologic claim is in scope.

Start with a [model proposal](../../issues/new?template=model_submission.yml).
The form accepts AI-based, AI+physics and physics models, with separate links
for the model code and optional pretrained weights, and asks whether there is
a time window you want the test to run over.

**You do not need to be able to package it yourself.** The *packaging status*
field decides who the issue is assigned to and nothing else: ready means it is
yours to finish, help needed means the maintainer takes it. Either way the
proposal counts as yours, and five accepted proposals earn co-authorship on
the benchmark paper the same as one merged probe. Knowing which models are
worth putting through the benchmark is a judgement about the field, and it is
not one the maintainer can make alone.

Once a proposal is accepted, fork and build; the pull request comes from the
fork and closes the issue.

A submitted model is scored on the largest flood event of the generated
record rather than on all ten years of it: by default a month for a daily
model and a week for an hourly one, with the full spinup in front, located
by the probe's own reference model. That is what keeps a model that takes
seconds per forecast inside the probe's time budget. `window_days` in
`model.yaml` changes it; `full` asks for the whole record.

The container runs with no network and never sees the probe code, so a model
cannot read the tolerance it is being judged against. Every evaluation is
appended to [models/result.csv](models/result.csv).

**A model that fails is worth proposing, and so is one most probes cannot
score.** `INCOMPLETE` is the current state of nearly every published
rainfall-runoff model on the budget probes — it has not violated
conservation, it has declined to be falsifiable — and recording that honestly
is a large part of what this is for.

## Contributing a probe

**Please submit one.** A benchmark with four probes tests four things; the
reason this repository is open is that the physics worth testing is wider
than any one group knows. If you have spent time with a conservation law that
AI models get wrong, that law is a probe, and we would rather have it from you
than approximate it ourselves.

**Contributors of merged probes are co-authors on the benchmark paper.** The
threshold is one probe, merged and passing the acceptance gate, or five
accepted model proposals. This is how model intercomparison projects have
always worked in this field: you contribute an experiment, you are an author
on the paper that reports it. [CONTRIBUTING.md](CONTRIBUTING.md#credit) has
the full terms — author order, the right to decline, and what happens before
anything is submitted.

**[ROADMAP.md](ROADMAP.md) lists the probes we want**, each with a difficulty
and a note on what it discriminates. Claiming one is easier than inventing
one, but inventing one is welcome too; propose it first so nobody builds it
twice. The two we most want are the cross-budget consistency probes: a model
can close its water budget and its energy budget while being incoherent
between them, and nothing in the suite currently notices.

Propose first, then fork, then open a pull request — the same three steps for
a probe or a model. The proposal issue is where we find out whether a probe
discriminates, which is cheaper to learn in a paragraph than in three hundred
lines. Once it is labelled `accepted` it is assigned to you and nobody else
will build it.

Then do not start from a blank page:

```bash
ht init-probe --list-templates           # the shapes available
ht init-probe --template <kind>          # writes probe-draft.yaml
#   ... fill in the fields ...
ht init-probe --from probe-draft.yaml    # creates probes/<law>/<slug>/
ht gate --probe <law>/<slug>             # prove it discriminates
```

The templates cover conservation over one case, extrapolation in space and in
time, counterfactual response, and invariance. Each scaffolds into a probe that
already passes the gate against a placeholder case, so you can watch it
separate the reference models before writing any physics, then replace the
case with yours. `docs/writing-a-probe.md` has the details.

The gate is the only bar that matters, and it is a technical one: your probe
must pass an exact physical model and catch the deliberately broken ones. See
[GOVERNANCE.md](GOVERNANCE.md) for how disagreements about tolerances get
settled.

## Community

Questions, probe ideas, a model you would like evaluated, a tolerance you
think is wrong: the [HydroTuring Discord](https://discord.gg/7SQb6bUZD) is
where that conversation happens, before and alongside the issues.

## Status

Suite `0.1.0`, pre-release. Twenty-two probes, sixteen mass, five energy, one momentum, synthetic track only. More
energy and momentum probes, and the real-data track, are next. The harness runs paired cases and
scores labelled regimes, so the generalisation probes on the roadmap —
extrapolation in space and time, counterfactual response, invariance — are
unblocked and unclaimed. Scores are only comparable within a suite version.

## Layout

```
src/hydroturing/   harness: protocol, runners, criteria, scoring
probes/<law>/<id>/ probe.yaml, generate.py
models/<name>/     model.yaml, Dockerfile, adapter
schemas/           JSON schemas for both spec files
```

# GR4J (airGR) under HydroTuring

GR4J is the four-parameter daily rainfall-runoff model of Perrin, Michel and Andréassian
(2003, *J. Hydrol.* 279, [10.1016/S0022-1694(03)00225-7](https://doi.org/10.1016/S0022-1694(03)00225-7)).
It is packaged here from INRAE's airGR 1.7.9, built from the CRAN mirror at
[`cran/airGR@4bb401c`](https://github.com/cran/airGR/tree/4bb401c6529e20116db357369c9f2bc34f959992).
The compiled kernel runs unchanged: the adapter calls `RunModel_GR4J` and reads what it returns.

Proposed as [#33](https://github.com/Flood-Lab/HydroTuring/issues/33) and again as
[#96](https://github.com/Flood-Lab/HydroTuring/issues/96). #96's proposer says a Dockerfile and
adapter exist, but no pull request was opened and none was available, so this package was
written independently. It takes #96's pin and its reading of the exchange as `gwex`.

Licence: GPL-2. The source the image was built from stays in the image, under `/opt/airGR`.

## What the adapter reports

| Column | What it is |
| --- | --- |
| `pr` | the forcing, echoed as text |
| `evspsbl` | `AE`, actual evapotranspiration |
| `mrro` | `Qsim`, total runoff |
| `dis` | `mrro × area_km2 / 86.4`, m3/s |
| `gwex` | `AExch`, the *actual* exchange with the outside (`AExch1 + AExch2`), positive into the catchment |
| `mrso` | `Prod`, the production store |
| `gw` | `Rout`, the routing store |
| `channel` | the water inside the two unit hydrographs |
| `snw`, `canopy` | identically zero; see below |

All fluxes are mm/day at every step and all stores are absolute, in mm. Every row is reported,
spinup included.

**Zero stores.** `snw` and `canopy` are reported as identically zero and named under `zero_stores`
in `run.json`: GR4J has no snow module, so all precipitation is liquid, and no interception store,
the P/E neutralisation holding no water between steps. That is what lets the budget probes score it.

**Not reported, because the model has none of them**, and named under `not_reported` in `run.json`:
`snm`, every energy flux, `stage`.

**The closure family closes because of `gwex`.** On the first gate seed of `mass/catchment-closure`,
over 11 years and 8938 mm of rain, the budget over `mrso + gw + channel` with `gwex` as a source
closes to 1e-12 mm on every step and 5e-13 mm over the record, and the exchange is −571 mm, 6.4% of
the rain: hiding it would open the budget by exactly that. The hourly kernel closes as well, to
9e-13 mm per step on the hourly case of `mass/resolution-invariance`. `run.json` carries the
whole-record residual on every run. The harness flags these `suspicious_exact`, as it must for a
budget closed to rounding. Production and routing stores are the kernel's own, and `channel` is the
account of its fluxes, held to `StateEnd` as below, so this says the kernel's bookkeeping is
consistent and that the exchange is reported.

### Reading the unit hydrographs

The kernel holds their contents in `StUH1` and `StUH2` but returns them only for the last step,
in `StateEnd`. The store at each step is therefore the account of the kernel's own fluxes,
`cumsum(PR - Q9 - Q1)` from empty, and the adapter stops with an error if that account differs
from `StateEnd` by more than 1e-8 mm (it differs by about 1e-13). The first ordinate of each array
in `StateEnd` must be left out of that comparison: it is what was released this step and it sits
in the array until the next step overwrites it. Summing `StateEnd$UH` whole overstates the store
by the last step's release (20.3 mm against 15.4 mm on airGR's example basin).

### Why the routing store is `gw`

This is a judgement and it is the choice most open to argument. #96 proposes `channel` for the
routing store and the unit hydrographs together. Here the routing store is `gw` because it is the
model's only slow store: it carries the baseflow, it is the store the groundwater exchange acts
on, and it holds water for months (it never fell below 17 mm in three years on
`momentum/routing-conservation`). AGENTS.md asks for "water in transit through a unit hydrograph"
under `channel` and "a groundwater zone" under `gw`, and `models/dhbv2` makes the same split, with
HBV's response boxes under `gw` and its unit hydrograph under `channel`. Against that, the store's
name in the literature is "routing store".

No water budget depends on this, since closure sums every store. It moves one probe's number:
as mapped here `momentum/routing-conservation` stays within 0.73 of its allowance. Under the bound
as it read before `dc9836c` (below), the store reached up to 2.9 times the allowance as mapped here,
on 9 to 11 of 1095 steps, and up to 63 times with the routing store under `channel`, on 925 to 1007
of 1095. The second mapping was not re-run under the corrected bound.

## Parameters

There is no calibration on a synthetic catchment. Each value is below.

| Parameter | Value | Where it came from |
| --- | --- | --- |
| X1, production store capacity | `soil_capacity_mm` (320 mm on most probes); 291 mm when absent | static.json, as `models/flex_lumped` and `models/sacsma_snow17` do; the fallback is the median of `Param_Sets_GR4J` |
| X2, exchange coefficient | −1.1 mm/day | median of `Param_Sets_GR4J` |
| X3, routing store capacity | 76.2 mm | median of `Param_Sets_GR4J` |
| X4, unit hydrograph time base | 1.879 days | median X4u (2.15) through `X4 = X4u / 5.995 × S^0.3` at S = 250 km2 |

`Param_Sets_GR4J` is the table of 27 generalist parameter sets airGR ships for exactly this
situation, a catchment that cannot be calibrated (Andréassian et al. 2014, *WRR* 50,
[10.1002/2013WR014761](https://doi.org/10.1002/2013WR014761)). A column-wise median is not itself
one of the 27 sets. The formula for X4 is from that data set's help page.

The set was picked before any probe ran. X2 is non-zero on purpose: with X2 = 0 the column is
identically zero and the exchange, the part of GR4J this benchmark most wants to read, is never
exercised.

S is held at 250 km2 rather than read from the case, so X4 does not follow the area told. Reading
it would be a defensible use of that published scaling, but it would change the storm response
with the stated area, which `mass/area-invariance` would then charge to the model.

**Initial state.** airGR's default levels, `Prod` 0.3 X1 and `Rout` 0.5 X3, both unit hydrographs
empty, passed explicitly so that they are the same water at every step.

## The hourly step

The issues name PT1D, but the resolution probe runs every model at its own step and the next finer
one whatever the manifest lists, so a PT1H row needs an honest answer. GR4J's own kernel cannot
give one: its hydrograph arrays hold 20 steps (40 for UH2) while X4 is 45 hours, so the kernel
would truncate them and lose water.

A PT1H row therefore runs airGR's **GR4H**, the package's own hourly twin. It has the same
structure, with its own percolation constant (21/4 against 9/4) and hydrograph exponent (1.25
against 2.5). **So at PT1H this is a closely related kernel rather than GR4J itself, and the
resolution result compares the two.** The per-day parameters are put in the units of the row:

| | Conversion | Why |
| --- | --- | --- |
| X1 | unchanged | a capacity |
| X2 | × 24^(−1/8) | exchange is X2 (R/X3)^3.5 per step, so this keeps the rate per unit time once X3 has moved |
| X3 | × 24^(1/4) | for R ≪ X3 the store releases R⁵ / (4 X3⁴) per step, so this keeps the rate per unit time |
| X4 | × 24 | a time, counted in rows |

These are derived here from the kernel's equations, not copied. Ficchì, Perrin and Andréassian
(2016, *J. Hydrol.* 538) adapt GR4J to sub-daily steps and are the place to compare them. Their
paper was not accessible when this was written, so the agreement has not been checked. No row is
resampled.

The first version of this adapter took the initial routing store from airGR's level default at
the hourly X3, which started the hourly run with 46 mm more water than the daily one. That was
the adapter's fault, and it moved the resolution result from 1.9% to 3.5% of the rain before it
was found and fixed.

## Result on the gate seeds

**20 PASS, 1 FAIL, 16 N/A**, under Docker on Apple Silicon at `3e86a79`, version `1.7.9-4bb401c.2`, on the
whole record of every probe.

| Probe | Result |
| --- | --- |
| `energy/pet-consistency` | PASS: 0.96 of potential when wet, 0.36 when dry, 0.63 over the record |
| `mass/antecedent-monotonicity` | PASS: 0.23 of the 60 mm storm |
| `mass/area-invariance` | PASS: exactly 0; `mrro` never reads the area |
| `mass/catchment-closure` | PASS: 0.0000% of rain, largest row 1.1e-12 mm; ET at most 0.683 of potential |
| `mass/causality` | PASS: 0 before the storm, +0.736 of the added rain after |
| `mass/dry-down` | PASS: 20 mm over 730 dry days |
| `mass/extreme-event-closure` | PASS: worst event residual 8.4e-13 mm against 0.001 mm |
| `mass/extreme-rain` | PASS: returns 0.99 of the 997 mm added |
| `mass/human-abstraction` | **FAIL (VIOLATION)**: all 380 mm unaccounted for |
| `mass/multi-decadal-drift` | PASS: storage changes 6e-14 mm over the final 1825 days |
| `mass/phase-counterfactual` | PASS: exactly 0; the model never reads temperature |
| `mass/precipitation-counterfactual` | PASS: of +1589 mm, 0.178 evaporates, 0.717 runs off, 0.016 stays, the rest leaves through `gwex` |
| `mass/resolution-invariance` | PASS: 1.91% of rain, GR4J against GR4H |
| `mass/response-nonnegativity` | PASS: largest dip -1.0e-05 of the storm |
| `mass/runoff-bounds` | PASS: 0.40 of the rain, inside [0.05, 0.99], declared exchange -0.07 of rain |
| `mass/spinup-cycle-invariance` | PASS: 0.00% between spin-up lengths |
| `mass/steady-state` | PASS: 0.380 mm/day under 2.5 mm/day of rain |
| `mass/time-origin-invariance` | PASS: exactly 0; the model never reads the date |
| `mass/ungauged-basin-closure` | PASS: 0.0000% of rain over 12 catchments; ET at most 0.909 of potential |
| `mass/warming-response` | PASS: mrro -0.226 and -0.282 per unit of demand |
| `momentum/routing-conservation` | PASS: within 0.73 of the allowance, the tightest of any model this probe scores; see below |

Passes that say little: `area-invariance`, `phase-counterfactual` and `time-origin-invariance`, for
the reasons in the table, and `state_bounds` on `snw` and `canopy`, which holds on a column of zeros.

N/A: 13 INCOMPLETE, the outputs it lacks: `aquifer-recharge-ordering`, `groundwater-datum-invariance`, `gw-sw-exchange-consistency`
(`gw_sw_exchange`), `froude-regime`, `stage-discharge-monotonic`,
`uniform-flow-friction-consistency` (`stage`), `snowpack-mass-closure` (`snm`), and the six energy
probes that need radiative or heat outputs. 3 INCOMPATIBLE:
`momentum/routing-lag-consistency` and `momentum/wave-celerity-bounds`, whose geometry GR4J does not
take, X4 being a fixed parameter here and not a travel time, and `mass/exchange-response`, because
GR4J's exchange follows its own routing store and no external head, so `gwh` is not declared.

### `mass/human-abstraction`

Structural, verified. GR4J has no withdrawal term and the adapter never reads `abstr`, so the
natural and irrigated runs are identical: `evspsbl`, `mrro` and storage differ by exactly 0 on all
three seeds, and the whole 380 mm is the residual. AGENTS.md says a model that never reads the
column fails this way, and it does. `gwex` does not rescue it: the exchange is X2 (R/X3)^3.5 whatever the column holds.

### `momentum/routing-conservation`

This failed until `routing_conservation` was corrected at `dc9836c`. The allowance read only runoff
already released:

> the channel holds 6.4 mm on 2001-09-28 while the largest runoff of the preceding 30 days is
> 0.14 mm/day; a 15-day hydrograph cannot hold more than 2.2 mm of that

It now also reads the fifteen days ahead, since water held now leaves within the lag, and the model
passes. What follows is the measurement from the earlier bound, kept because it says how the store
behaves.

**The mechanism, checked against airGR's fluxes on all three seeds.** GR4J's hydrographs are upstream of a store
that retains water. After a dry spell `mrro` is about 0.1 mm/day, because the routing store
releases R⁵. The first big rain then puts most of that day's effective rainfall into the
hydrographs (on the worst seed, 80 mm of rain gives PR 7.9 mm and 6.4 mm held) while the outlet has
barely moved (`Qsim` 0.137). The next day that water has left: 5.7 mm goes to the routing store,
`channel` drops to 0.37 mm, and the store is back inside the allowance.

So the excess is transient and the kernel does not leak. It lasts 9 to 11 steps of 1095 on each
seed at 2.2 to 2.9 times the allowance, with no accumulation, and a unit pulse empties in about
2 × X4 days apart from the percolation still feeding it. The excess appears because 90% of what
leaves the hydrographs goes into the routing store rather than to the outlet, so `mrro` is not the
flow feeding this channel.

## How far each choice moves the results

Each alternative was re-run on the 13 probes that scored the model before the zero-store rule, on
the same gate seeds, under `routing_conservation` as it read before `dc9836c`, and were not re-run
since. **No parameter set changed a verdict among them**: all four gave the same 12 PASS and 1 FAIL,
that failure being routing. The closure family was not part of this comparison.

| | Runoff ratio | Resolution (worst, % of rain) | Routing, channel vs allowance on the worst day | Steady runoff (mm/day of 2.5) |
| --- | --- | --- | --- | --- |
| **As packaged** (X2 −1.1, X3 76.2, X4 1.88) | 0.40 | 1.91 | 6.4 vs 2.2 mm | 0.380 |
| X2 = 0, X3 = 90, X4 = 1.7, the values usually quoted as the Perrin et al. medians | 0.47 | 3.92 | 5.9 vs 3.4 mm | 0.500 |
| airGR's help-page set (X2 +1.012, X3 88.235, X4 2.208) | 0.55 | 5.52 | 8.6 vs 3.6 mm | 0.644 |
| X1 = 291 mm, ignoring static.json | 0.41 | 2.13 | 6.7 vs 2.3 mm | 0.380 |

Those Perrin values were not checked against the paper's table: airGR's default store levels (350,
90) and airGRteaching's stated median X2 of 0 agree with them, while airGRteaching's own example
uses X4 = 1.4. The help-page set is calibrated on one basin, and its X2 is an import.

## Packaging

`rocker/r-ver:4.4.2` is pinned by digest, airGR by commit, and jsonlite 1.8.9 by SHA-256 from CRAN's
archive. Apt's `git`, `curl` and `ca-certificates` in the throwaway source stage are not pinned, and
nothing from that stage reaches the image except the checked sources. The adapter is written in R,
so the image has no Python. It uses no random numbers, so the request's seed plays no part, and no
statistic of the record: each row depends only on the rows before it.

One ten-year case of 4015 rows takes about 0.7 s through the container, most of it starting R.

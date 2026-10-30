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

All fluxes are mm/day at every step and all stores are absolute, in mm. Every row is reported,
spinup included.

**Not reported, because the model has none of them**, and named under `not_reported` in
`run.json`: `snw` (no snow module, so all precipitation is liquid), `canopy` (no interception
store; the P/E neutralisation holds no water between steps), `snm`, every energy flux, `stage`.
Nothing is filled with zeros, so the probes that require them are N/A (INCOMPLETE).

**That includes the closure family.** `mass/catchment-closure` and its relatives all require
`snw` and `canopy`, so no closure probe scores GR4J. What those probes would have read was
measured directly instead, off the harness. On the first gate seed of `mass/catchment-closure`,
over 11 years and 8938 mm of rain, the budget over `mrso + gw + channel` with `gwex` as a source
closes to 1e-12 mm on every step and 5e-13 mm over the record. The hourly kernel closes as well, to
9e-13 mm per step on the hourly case of `mass/resolution-invariance`. It closes because of `gwex`:
the exchange is −571 mm, 6.4% of the rain, and hiding it would open the budget by exactly that.
`run.json` carries the whole-record residual on every run.

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

No water budget depends on this, since closure sums every store. It moves one probe's number and
no verdict: `momentum/routing-conservation` fails either way, at up to 2.9 times its allowance as
mapped here, on 9 to 11 of 1095 steps, and at up to 63 times with the routing store under
`channel`, on 925 to 1007 of 1095.

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

12 PASS, 1 FAIL, 24 N/A: 22 INCOMPLETE (they require a store or flux above) and 2 INCOMPATIBLE
(`momentum/routing-lag-consistency` and `momentum/wave-celerity-bounds`, whose geometry GR4J does
not take: X4 is a fixed parameter here, not a travel time).

The one failure is `momentum/routing-conservation`:

> the channel holds 6.4 mm on 2001-09-28 while the largest runoff of the preceding 30 days is
> 0.14 mm/day; a 15-day hydrograph cannot hold more than 2.2 mm of that

**The mechanism, checked against airGR's fluxes on all three seeds.** The allowance is 15 days of
the largest `mrro` of the past month, including today. GR4J's hydrographs are upstream of a store
that retains water. After a dry spell `mrro` is about 0.1 mm/day, because the routing store
releases R⁵. The first big rain then puts most of that day's effective rainfall into the
hydrographs (on the worst seed, 80 mm of rain gives PR 7.9 mm and 6.4 mm held) while the outlet has
barely moved (`Qsim` 0.137). The next day that water has left: 5.7 mm goes to the routing store,
`channel` drops to 0.37 mm, and the store is back inside the allowance.

So the excess is transient and the kernel does not leak. It lasts 9 to 11 steps of 1095 on each
seed at 2.2 to 2.9 times the allowance, with no accumulation, and a unit pulse empties in about
2 × X4 days apart from the percolation still feeding it. The excess appears because 90% of what
leaves the hydrographs goes into the routing store rather than to the outlet, so `mrro` is not the
flow feeding this channel. It is reported as the probe scores it, and the probe was not changed;
whether a rising-limb transient of this kind should count is a question for the probe's owners.

## How far each choice moves the results

Every scored probe was re-run under each alternative on the same gate seeds. **No parameter set
changes a verdict**: all four give the same 12 PASS and 1 FAIL.

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

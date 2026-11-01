# mHM under HydroTuring

mHM, the mesoscale Hydrologic Model (Samaniego, Kumar and Attinger 2010,
doi:10.1029/2008WR007327; Kumar et al. 2013), with its routing module mRM, from
[mhm-ufz/mHM](https://github.com/mhm-ufz/mHM) at tag `v5.13.4`, commit
`5f2de41bf244e9bd0d4b7a27806e0c96ff8186db`, built against FORCES `v0.10.0`, commit
`1e0194777d35ddc4b0158aab16522dcaabc4c378`, the version `version_forces.txt` names. Proposed
by Valmik Nahata as [Flood-Lab/HydroTuring#58](https://github.com/Flood-Lab/HydroTuring/issues/58).
GPL-3.0-or-later (mHM), LGPL-3.0 (FORCES); the licences are in the image at `/model`, and the
image is never pushed to a registry.

## What runs

Upstream's sources, compiled unmodified with upstream's CMake as a subdirectory of
`CMakeLists.txt`, on `python:3.11-slim` pinned by digest (Debian 13). The Apptainer image is
110 MB, and a case's run needs no network.

The executable is `ht_mhm_driver.f90`, not upstream's `mhm`. It is upstream's `mhm_driver.f90`
with `mhm_eval`'s loop written out, through the same `mhm_interface_run_*` calls in the same
order, plus one record per row read from upstream's global arrays. It exists because upstream's
netCDF writer gives a day's states only as the day's mean, or else an hourly table too large for
a run's 64 MB tmpfs, and the contract needs each store at the end of its row. It computes nothing
except the area weights. It runs one domain, the first land-cover scene and no optimisation, which
is what every case here is.

### The driver against upstream's entry point

Checked on upstream's shipped case: `test_domain` as it ships (Mosel basin, 34 cells at 24 km,
1990 to 1993, 35,064 hourly steps), with its own `mhm.nml`, `mhm_parameter.nml` and output
namelists, once with `/model/mhm` and once with `/model/ht_mhm_driver`. **Verified identical:**

| Output | Compared | Result |
| --- | --- | --- |
| `daily_discharge.out` | bytes | identical |
| `discharge.nc` | 4 variables | max difference 0 |
| `mHM_Fluxes_States.nc` | 30 variables, 48 records | max difference 0 |
| `mRM_Fluxes_States.nc` | 9 variables, 1461 records | max difference 0 |

The driver's own record was then compared with upstream's executable writing hourly netCDF (the
shipped case with only the output interval changed), over 1990, the part of the record in the first
land-cover scene, since the driver reads scene 1 for its sealed weights and a case here has one
scene:

| Record | max difference | Upstream's largest value |
| --- | --- | --- |
| snowpack at the end of the row | 7.1e-15 mm | 36.8 mm |
| interception at the end of the row | 2.2e-16 mm | 1.02 mm |
| aET, sum over the row | 2.2e-15 mm | 4.78 mm |
| Q generated, sum over the row | 1.8e-15 mm | 5.67 mm |
| preEffect, sum over the row | 1.8e-14 mm | 38.2 mm |
| melt (`Qsm`), sum over the row | 1.8e-15 mm | 5.89 mm |
| PET, sum over the row | 5.3e-15 mm | 5.85 mm |
| outlet discharge, all 1461 days, against `Qsim` in `daily_discharge.out` | 5.0e-8 m3/s (the file's print precision) | 1123 m3/s |

The four soil and groundwater records (`SWC`, `sealedSTW`, `unsatSTW`, `satSTW`) were not compared
this way, because the driver multiplies them by each cell's sealed or unsealed share, which upstream's
netCDF does not store. What supports them is the closure of every case below: residual 0.0000% of rain,
with a largest single-row residual of 2.6e-13 mm on the closure case.

## What is reported

| Column | Native | Note |
| --- | --- | --- |
| `pr` | the forcing | echoed from the CSV as text |
| `evspsbl` | `aET`: soil over the unsealed share + canopy + sealed over the sealed share | as `mHM_updateDataset` composes it, mm/day |
| `mrro` | mRM's `mRM_runoff` at the gauge, as a depth over the sub-basin | what has passed the outlet |
| `dis` | `mrro` over the stated `area_km2` | scaled, not mRM's m3/s on the 286 km2 domain |
| `snm` | `L1_preEffect` | rain and melt leaving the snow module; that rain has already passed interception |
| `mrso` | soil horizons × unsealed share + `L1_sealSTW` × sealed share | end of row, mm |
| `snw` | `L1_snowPack` | end of row, mm |
| `canopy` | `L1_inter` | end of row, mm |
| `gw` | (`L1_unsatSTW` + `L1_satSTW`) × unsealed share | end of row, mm |
| `channel` | cumulative `L1_total_runoff` − cumulative `mrro` | water generated and not yet past the outlet; mRM's own store is not differenced |

Every water store mHM has is reported, so `zero_stores` is empty. `not_reported` in `run.json`: `stage`
(mRM has no water level) and the energy terms (mHM has no energy balance). There is no `gwex`: every
geological unit in the shipped table is non-karstic, `karst_loss_min` is 1.0 in the driver's record, so no
water leaves below ground. `abstr` is not read, because neither mHM nor mRM has a withdrawal process
(no reference to one in `src/mHM` or `src/mRM` at the pin).

**A limit on what closure says.** `channel` is a difference of two cumulative fluxes, so the closure
residual is the land surface's budget alone. A loss inside mRM would not appear as a residual. It would
appear as a `channel` that keeps growing, which `momentum/routing-conservation` bounds and
`mass/multi-decadal-drift` watches. Both pass.

## Domain and parameters

No calibration exists on a synthetic catchment, and these results describe one documented point in
parameter space. MPR runs as upstream ships it: upstream's global parameters, unchanged, on real 500 m
morphological grids.

`derive_domain.py` cuts the domain at build time from upstream's `test_domain`, following upstream's flow
directions from one cell of its river network (row 92, column 65 of the 500 m grid). That gives a natural
sub-basin of 285.75 km2, 31 level-1 cells at 4 km, chosen as the nearest of the candidates of that size to
the whole basin's land cover (34% forest, 3% urban cells, 63% pervious; sealed share 1.68% after MPR). No
morphological value is changed. **The stated area enters only `dis`**, as in `models/wflow_sbm` and
`models/lisflood`.

| Item | Value | Source |
| --- | --- | --- |
| Global parameters (every transfer function, 146 marked free or fixed in the file) | as shipped | `test_domain/mhm_parameter.nml` |
| `processCase` 1, 2, 3, 4, 6, 7, 9 | 1 | `test_domain/mhm.nml` |
| `processCase(5)` | 0: PET is an input, multiplied by mHM's aspect factor | upstream's option, and its test domain's |
| `processCase(8)` | 3: adaptive routing with slope-dependent celerity | `test_domain/mhm.nml` |
| Karstic loss | off | every unit in the shipped geology table |
| `resolution_Hydrology`, `resolution_Routing` | 4000 m | **choice**: upstream's 24 km gives this sub-basin one or two cells and no river |
| `timestep` | 1 h | every configuration upstream ships |
| `iFlag_soilDB` | 1 (upstream's test domain: 0) | **choice**: only option 1 takes the horizon depths from the namelist, which the capacity needs |
| `nSoilHorizons_mHM`, `tillageDepth`, `soil_Depth(1)` | 2, 200 mm, 200 mm | `test_domain/mhm.nml` |
| `soil_Depth(2)` | from `soil_capacity_mm`; 1000 mm if none is stated | **from static.json** |
| `fracSealed_cityArea` | 0.6 | `test_domain/mhm.nml` |
| `evap_coeff` | 1.0 every month (upstream's test domain: 0.72 to 1.5 by month) | **choice**: no calendar; 1.0 is upstream's value for May to November |
| `fnight_prec` | 0.5 every month | annual mean of upstream's table, flat |
| `fnight_pet` | 0.1 every month | upstream's, already constant |
| `fnight_temp` | −1.84 every month | annual mean of upstream's table, flat |
| LAI table | each class at its own annual mean, every month | **choice**, see below |
| Land cover | `lc_1991.asc`, one scene for all time | one of upstream's two scenes |

**Capacity.** mHM's capacity is saturated water content times depth over its horizons, so the depth of the
second horizon carries the stated capacity: `depth = 200 + (capacity − 0.0084 − 0.43855 × 200) / 0.34592`.
The two slopes are read from the driver's own record on this sub-basin (87.71 mm over the first 200 mm,
176.80 mm over the next 511.11 mm), and 0.0084 mm is the sealed retention, 0.5 mm over 1.68%. A capacity
below what the tillage layer already holds cannot be met, and the second horizon stops at 10 mm.

**No calendar.** Upstream's LAI and day/night ratios are monthly tables for Central Europe. Nothing in a case
says when this catchment leafs out, and keeping them would hand the model a seasonal cycle the forcing does
not have, so each is flat at its annual mean. `time-origin-invariance` passes in part because of this.

**PET.** mHM multiplies the input by an aspect factor, 0.94 on average over this sub-basin, from upstream's
`minCorrectionFactorPET` and `maxCorrectionFactorPET`. So mHM's demand is 6% below `pet`, and that is mHM.

### What the choices move

On the closure case (seed 598896396, 1095 rows, the last 730 summed; every run closes to 0.000% of rain):

| Variant | P | ET | Q | Runoff ratio | Peak `channel` (mm) |
| --- | --- | --- | --- | --- | --- |
| as packaged, 4 km | 1594.6 | 863.0 | 705.4 | 0.442 | 2.84 |
| same sub-basin at 8 km | 1594.6 | 864.2 | 704.4 | 0.442 | 2.13 |
| whole Mosel basin at 24 km, uniform forcing | 1594.6 | 844.8 | 723.4 | 0.454 | 15.68 |
| upstream's monthly LAI table | 1594.6 | 861.8 | 706.6 | 0.443 | 2.83 |
| `fnight_temp` = 0 | 1594.6 | 863.2 | 705.8 | 0.443 | 2.84 |
| routing `processCase(8)` = 1 | 1594.6 | 863.0 | 705.6 | 0.442 | 1.70 |
| middle of upstream's range for 113 of the 120 free parameters | 1594.6 | 897.9 | 662.5 | 0.415 | 3.69 |

The 4 km against 8 km row is the closest thing here to a test of what proposal #58 wanted tested: MPR moves
the same sub-basin's runoff by 0.06% of the rain when its hydrology grid doubles. That is one case, and the
adapter's run, not a probe.

## The step

**Declared `[PT1D, PT1H]`.** mHM reads daily or hourly forcing and steps hourly in both. A daily row is split
into hours by mHM's own disaggregation with the flat ratios above; the adapter does not resample. Depths per row
go in, fluxes come back as mm/day. Nothing finer is available: `request.json`'s other steps are refused.

**`timestep = 24`, the other value the namelist names, does not run correctly.** On the shipped case it exits 0,
and the daily discharge skill falls from KGE 0.746, NSE 0.767 at `timestep = 1` to KGE −3.69, NSE −25.1. Verified as
a measurement; why it does so was not traced. So the package never uses it, which answers the question in the
proposal's review. It also means **`mass/resolution-invariance` compares the same hourly mHM under two ways of
delivering the forcing**. Its 2.6% measures mHM's disaggregation, not a dependence on the model's step.

## Evaluation

Gate seeds, all 37 probes, Docker on Apple Silicon at `3e86a79`: **21 PASS, 1 FAIL, 15 N/A** (13 INCOMPLETE, 2 INCOMPATIBLE) at version
`5.13.4-ht.1`.

| Probe | Verdict | Reading | Mechanism |
| --- | --- | --- | --- |
| `mass/human-abstraction` | **FAIL** | all 380 mm of prescribed `abstr` unaccounted for; closure itself 0.0000% | **verified, property of the model**: no withdrawal process in `src/`, so the adapter does not read the column. AGENTS.md scores this as a failure |
| `energy/pet-consistency` | PASS | 0.80 of demand wet, 0.33 dry | |
| `mass/antecedent-monotonicity` | PASS | wetter catchment runs off 19.1 mm more of 60 mm | |
| `mass/area-invariance` | PASS | departure 0 | **trivial**, see below |
| `mass/catchment-closure` | PASS | residual 0.0000%, ET 0.853 of PET | land budget, see above |
| `mass/causality` | PASS | nothing moves before the storm | |
| `mass/dry-down` | PASS | 24 mm over 730 dry days | |
| `mass/extreme-event-closure` | PASS | residual 0.0000%, worst event 1.7e-13 mm | |
| `mass/extreme-rain` | PASS | returns 1.00 of the 997 mm added | |
| `mass/multi-decadal-drift` | PASS | storage changes −1.1e-13 mm over the last five years | |
| `mass/phase-counterfactual` | PASS | 3.3% of the rain | |
| `mass/precipitation-counterfactual` | PASS | added rain: ET +0.134, runoff +0.838, storage +0.028 | |
| `mass/resolution-invariance` | PASS | worst `mrro` 2.6% of rain | **partly trivial**, see "The step" |
| `mass/response-nonnegativity` | PASS | largest dip −9.8e-5 of the 120 mm | |
| `mass/runoff-bounds` | PASS | 0.54 of the rain, inside [0.11, 1.07] | |
| `mass/snowpack-mass-closure` | PASS | residual 0.96% | |
| `mass/spinup-cycle-invariance` | PASS | worst departure 0.00% | |
| `mass/steady-state` | PASS | runoff 0.684 of 2.5 mm/day, settles to 0.32% (`gw`) | |
| `mass/time-origin-invariance` | PASS | departure 0 | **partly by configuration**: the flat tables remove the calendar mHM would otherwise read |
| `mass/ungauged-basin-closure` | PASS | residual 0.0000%, ET 0.791 of PET | |
| `mass/warming-response` | PASS | runoff −0.411, ET +0.407 per unit of demand | |
| `momentum/routing-conservation` | PASS | `channel` within 0.02 of a 15-day hydrograph | |

N/A (INCOMPLETE, 13): `energy/evaporative-partition`, `latent-heat-et-consistency`, `radiation-consistency`,
`snowmelt-energy-water`, `soil-heat-storage-consistency`, `surface-energy-closure`, `mass/aquifer-recharge-ordering`,
`exchange-response`, `groundwater-datum-invariance`, `gw-sw-exchange-consistency`, `momentum/froude-regime`,
`stage-discharge-monotonic`, `uniform-flow-friction-consistency`: no energy terms, no `gwex`, no groundwater exchange,
no stage. INCOMPATIBLE (2): `momentum/routing-lag-consistency` and `wave-celerity-bounds`. The first needs channel-length attributes the adapter does not read, the second a prescribed `q_in` and a reach's geometry; the sub-basin's own network stands in for those here, and mHM has no way to take them.

**Trivial pass that matters.** Issue #58 chose mHM because `mass/area-invariance` looked like the probe that decides MPR's
scale claim. Here it does not. The stated area enters only `dis`, as it must when the domain is a fixed sub-basin
of real grids, so depths cannot depend on it. The probe's "departure 0" says nothing about MPR. The 4 km and 8 km row in
the table above is the only evidence here about MPR's flux-matching.

## An alternative parameter set

Every parameter upstream marks free (flag 1) moved to the middle of upstream's stated range, except that the script's pattern misses the seven free indexed `GeoParam(i,:)` entries, which stay at upstream's values: 113 of 120
(`midrange_parameters.py`, in the image as `/model/sensitivity/midrange.nml`; not a calibration, not upstream's). Full
evaluation on gate seeds: **15 PASS, 7 FAIL, 15 N/A**. `mass/human-abstraction` fails
as before. The other six (`catchment-closure`, `extreme-event-closure`, `multi-decadal-drift`, `precipitation-counterfactual`,
`spinup-cycle-invariance`, `ungauged-basin-closure`) fail only `state_bounds`, with `mrso` above the stated capacity, for
example 353 mm against 320 mm. **That is the package and not mHM**: the slopes that turn a capacity into a depth are read at upstream's
parameters, and the pedotransfer parameters moved, so the same depth holds more. Closure stays at 0.0000% on all of them, the other
criteria still pass, and `mass/resolution-invariance` (2.5%), `mass/phase-counterfactual` (0.9%) and
`mass/antecedent-monotonicity` (0.24 of the storm) barely move. Nothing about the verdicts as packaged depends on the capacity mapping
being more general than this.

## Limits

- One documented point in parameter space, on one sub-basin of one basin's morphology, none of it calibrated for these
  catchments. A soil capacity is met through the second horizon's depth at this sub-basin's soil mix.
- Forcing is uniform over the domain, so nothing tests spatial organisation or MPR's handling of heterogeneous forcing.
- `iFlag_soilDB` is 1 where the shipped test domain has 0. The driver check above used the shipped value; it was not rerun at 1.
- `stage` is not reported, so `momentum/stage-discharge-monotonic` is N/A; mRM has no water level to give.

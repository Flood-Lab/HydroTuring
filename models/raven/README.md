# Raven (HBV-EC emulation) under HydroTuring

Raven, the hydrological modelling framework of Craig et al. (2020, *Environ. Model. Softw.* 129,
[10.1016/j.envsoft.2020.104728](https://doi.org/10.1016/j.envsoft.2020.104728)), in its emulation of
HBV-EC, packaged from v4.12 at commit `81011f061e53ef1289ebcab36dc3622cea05c76d`. Proposed by Yan Chen as
[Flood-Lab/HydroTuring#116](https://github.com/Flood-Lab/HydroTuring/issues/116), and narrowed by the
proposer there to the HBV-EC emulation rather than Raven in general. Artistic-2.0; the image carries the licence
at `/model/LICENSE`.

**The result that most needs reading is the gap below.** The packaged configuration scores 15 FAILs. Two
changes that keep the rest of the package score 8 and 7. No configuration was chosen for its score. The
gap is one quantity, the runoff ratio, crossing one threshold, and six FAILs stay whichever
configuration is run.

## What runs

Upstream's own CMake target, built unmodified at the pinned commit, and upstream's own five files
`benchmarking/_InputFiles/Salmon_HBV/raven-hbv-salmon.rv*`, baked into the image byte for byte.
`ht_adapter.py` writes the case's forcing, applies a named list of edits to those files, runs the
executable once and reads what it wrote. It computes no hydrology. Base image `python:3.11-slim` pinned by
digest (Debian 13); the build refuses any other commit. The evaluated image used cmake 3.31.6, g++ 14.2.0
and NetCDF 4.9.3 from the Debian archive, and needs no network at run time. NetCDF is for
`:WriteNetCDFFormat`: Raven's CSV prints six significant figures, which is too coarse for a closure on a
0.02 mm drizzle and for a one-in-a-million area probe. That was an adapter fault in the first evaluation, when
`mass/area-invariance` failed on `evspsbl` at 4.9e-06, and it is why the version is `.2`.

## Every edit to the template

Each is asserted to match exactly once and is written to `run.json` under `edits_from_template`.

| File | Template | Used | Why |
| --- | --- | --- | --- |
| rvi | `:StartDate`, `:Duration`, `:TimeStep` | the case's | the case's record and step |
| rvi | `:Evaporation`, `:OW_Evaporation PET_FROMMONTHLY` | `PET_DATA` | the case prescribes `pet`; the template derives PET from Salmon River normals |
| rvi | `:WriteForcingFunctions`, `:EvaluationMetrics` | `:WriteNetCDFFormat` and three `:CustomOutput` | output only |
| rvp | soil thickness 2.036937 m (`HBV_PARA_17`, 203.4 mm) | `soil_capacity_mm` / (1000 x `POROSITY`) | capacity from the case. **Moves the verdicts; see the gap** |
| rvp | `MAX_CAPACITY`, `MAX_SNOW_CAPACITY` 10000 mm, meaning no limit | `canopy_capacity_mm` | capacity from the case. **Moves the verdicts; see the gap** |
| rvh | area 4250.6 km2, latitude 54.4848 | `area_km2`, `latitude_deg` | the case's. Elevation stays 843 m, equal to the gauge's, so no orographic correction acts |
| rvc | initial outflow 1 m3/s | 0 | Raven would release Salmon River's 1 m3/s in the first step without taking it from a store |
| rvt | `:RainCorrection` 1.141608, `:SnowCorrection` 1.024278 | 1.0, 1.0 | gauge-undercatch factors calibrated at Salmon River. At those values Raven takes in 14% more rain than the case gave it |
| rvt | monthly normals, Salmon River series | the case's `tas`, `pr`, `pet` | `pet` also fills `OW_PET` |

The hydrologic process list, the soil model, every other parameter and the routing are the template's.

## What is reported

| Column | Raven | Note |
| --- | --- | --- |
| `pr` | the forcing | echoed from the CSV as text |
| `snm` | `BETWEEN ATMOS_PRECIP and PONDED_WATER` + `BETWEEN SNOW_LIQ and PONDED_WATER` | differenced from Raven's running sum |
| `evspsbl` | `AET` custom output | canopy, canopy snow, soil and lake evaporation |
| `mrro` | `q_sim` over `area_km2` | |
| `dis` | `q_sim` | m3/s |
| `mrso` | `Soil Water[0]` | the top soil only |
| `gw` | `Soil Water[1]` + `Soil Water[2]` | the fast and slow reservoirs |
| `snw` | `Snow` + `Snow Melt (Liquid)` | |
| `canopy` | `Canopy` + `Canopy Snow` | these exist; `CANEVP_ALL` empties them in the step they fill, so they read 0 at every row's end |
| `channel` | rivulet, channel and reservoir storage, `Surface Water`, `Ponded Water` | the unit hydrograph's water is in rivulet storage |

`zero_stores` and `not_reported` are both empty: every store is real and is reported. The adapter stops if
`WatershedStorage.nc` holds a variable it has not mapped, or if the two glacier stores hold water. No flux
is filled in. `abstr` is not read, and the template has no withdrawal process to give it to. There is no
`gwex`: the template exchanges nothing across the catchment boundary.

## Parameters

No calibration exists for the synthetic catchment. Everything below is the template's own, calibrated by
upstream for Salmon River near Prince George, B.C. (file headers: Juliane Mai and James Craig), and is
unchanged. Results describe that one documented point in parameter space, with capacities from the case. Salmon
River has about 650 mm of precipitation and 430 mm of PET a year in that benchmark; the closure probe, for one,
has about 810 mm and 750 mm.

| Parameter | Value | Template's label |
| --- | --- | --- |
| `RAINSNOW_TEMP`, `RAINSNOW_DELTA` | 0.05984519 C, 2.0 C | `HBV_PARA_01`, constant |
| `MELT_FACTOR`, `MIN_MELT_FACTOR` | 4.072232, 2.2 mm/d/K | `HBV_PARA_02`, constant |
| `REFREEZE_FACTOR` | 2.001574 mm/d/K | `HBV_PARA_03` |
| `SNOW_SWI` | 0.03473693 | `HBV_PARA_04` |
| `HBV_MELT_FOR_CORR`, `HBV_MELT_ASP_CORR` | 0.4452843, 0.48 | `HBV_PARA_18`, constant |
| `POROSITY` | 0.09985144 | `HBV_PARA_05` |
| `FIELD_CAPACITY`, `SAT_WILT` | 0.5060520, 0.04505643 | `HBV_PARA_06`, `HBV_PARA_14` |
| `HBV_BETA` | 3.438486 | `HBV_PARA_07` |
| `MAX_CAP_RISE_RATE` | 18.94145 mm/d | `HBV_PARA_16` |
| fast `MAX_PERC_RATE` | 38.32455 mm/d | `HBV_PARA_08` |
| fast `BASEFLOW_COEFF`, `BASEFLOW_N` | 0.4606565 /d, 1.877607 | `HBV_PARA_09`, `1 + HBV_PARA_15` |
| slow `BASEFLOW_COEFF` | 0.06303738 /d | `HBV_PARA_10` |
| fast and slow reservoir thickness | 100 m | constant, effectively unbounded |
| `TIME_CONC`, `TIME_TO_PEAK` | 2.277781, 1.138890 d | `HBV_PARA_11` and half of it |
| `ADIABATIC_LAPSE`, `PRECIP_LAPSE` | 0.5718813, 4.873686 | `HBV_PARA_13`, `HBV_PARA_12`; inactive at gauge elevation |
| `TFRAIN`, `TFSNOW` | 0.88, 0.88 | throughfall fractions |
| `FOREST_COV`, `IMPERM` | 1, 0 | constant |
| `LAKE_PET_CORR` | 1.0 | not set by the template; Raven's default (`LandUseClass.cpp:181`) |
| initial `SOIL[2]` | 0.50657 mm | template |

From the case: `area_km2`, `soil_capacity_mm`, `canopy_capacity_mm`, `latitude_deg`, and the forcing `pr`,
`tas`, `pet`.

## The step

Raven's rates are per day and its solver takes `:TimeStep` directly, so nothing is rescaled. Forcing is
written at the case's step, as rates. The manifest lists PT1D, the issue's step and the template's, and
PT1H, which the resolution probe exercises. The adapter accepts finer steps, but nothing has tested them, so
they are not claimed. At a sub-daily step `RAINSNOW_HBV` and `POTMELT_HBV` read `temp_daily_ave`
(`UpdateForcings.cpp:1098`, `PotentialMelt.cpp:83`), which Raven takes over the day, so an hourly run's phase and
melt see later hours of the same day (read from the source, not tested). No hourly probe checks that.

## Evaluation

Gate seeds, all 37 probes, under Docker on Apple Silicon at `3e86a79`, version `4.12-81011f0.2`: **7 PASS, 15 FAIL, 15 N/A** (13 INCOMPLETE, 2 INCOMPATIBLE).
No ERROR rows. A rerun on the merged tree reproduced every row of the previous run. The whole suite takes
about a minute on a node. Closure holds everywhere it is measured, with a residual of 0.0000% of the rain
and Raven's own mass-balance error at 0, at PT1D and PT1H alike. **Raven conserves water.** The FAILs say
something else about where this configuration sends it.

### The mechanisms

Three properties of the template produce nearly all of it.

**Lake evaporation with no lake.** `:LakeEvaporation LAKE_EVAP_BASIC SLOW_RESERVOIR ATMOSPHERE` takes the
full `OW_PET` from the slow reservoir whenever it holds water (`OpenWaterEvap.cpp:286`, `LAKE_PET_CORR` 1).
No lake fraction scales it, and it does not compete with the soil's or the canopy's evaporation: unlike
them it never subtracts what has already evaporated. So total evaporation can exceed PET. That is not
something this adapter introduced. Upstream's Salmon River benchmark, run unmodified, gives ET 29472 mm on
PET 24530 mm over its 57 years, 1.20 of PET, of which the lake term is 4993 mm, or 40% of everything the slow
reservoir receives. Removing the term gives 24478 mm, 1.00 of PET.

**Capillary rise that returns the first water.** `RISE_HBV` lifts `MAX_CAP_RISE_RATE` x (1 - S/Smax) from the
fast reservoir back into the top soil (`CapillaryRise.cpp:139`), up to 18.9 mm/d. Rain generates recharge at
(S/Smax)^3.44 of itself (`Infiltration.cpp:392`), so on a soil at a third of capacity a 60 mm storm generates
1 mm, and it goes straight back. On `mass/antecedent-monotonicity`, wet variant, seed 1295520324, the storm
row takes `Soil Water[0]` from 96.9 to 151.5 mm, and the reservoirs do not rise at all. Only storms on a wet
soil, generating more than that rate, reach the reservoirs.

**No drainage from the top soil.** Only evaporation empties it, and recharge goes only through rain. So
without rain, runoff stops once the reservoirs are empty, with 89 to 200 mm still in the soil.

In series, these are two sinks of nearly fixed size on what the soil generates. Full record, 4 years,
`mass/area-invariance` control, seed 66780378, in mm:

| Run | Generated | Back by capillary rise | Percolated | Lake evaporation | Slow outflow | Q/P |
| --- | --- | --- | --- | --- | --- | --- |
| **As packaged, soil 320 mm** | 684 | 300 | 384 | 275 | 109 | 0.032 |
| Template soil and canopy, 203 mm | 848 | 289 | 559 | 374 | 185 | 0.054 |
| No lake evaporation | 684 | 300 | 384 | 0 | 385 | 0.113 |
| No capillary rise | 498 | 0 | 498 | 403 | 96 | 0.028 |
| Neither | 498 | 0 | 498 | 0 | 497 | 0.146 |

A 24% change in what is generated is a 70% change in what leaves, which is why a modest parameter moves this
number so much.

### The gap

Raven has been run at all four corners of the two capacities the adapter takes from the case, and with
lake evaporation removed. FAIL counts on the same gate seeds:

| | Canopy from the case, 2 mm | Canopy from the template, 10000 mm |
| --- | --- | --- |
| **Soil from the case, 320 mm** | **15, as packaged** | 16 |
| **Soil from the template, 203 mm** | 8 | 8 (`raven_alt`) |

The packaged capacities without `:LakeEvaporation` (`raven_nolake`): 7.

So the soil capacity drives it, and the canopy capacity trades one or two probes for others. The soil
capacities are the values most probes declare; `mass/resolution-invariance` declares 120 mm and
`mass/ungauged-basin-closure` draws its own. Verdict by probe:

| Probe | Packaged | Soil only from case | Canopy only from case | `raven_alt` | `raven_nolake` |
| --- | --- | --- | --- | --- | --- |
| `energy/pet-consistency` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `mass/antecedent-monotonicity` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `mass/dry-down` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `mass/human-abstraction` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `mass/multi-decadal-drift` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `mass/ungauged-basin-closure` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `mass/area-invariance` | FAIL | FAIL | PASS | PASS | PASS |
| `mass/causality` | FAIL | FAIL | PASS | PASS | PASS |
| `mass/extreme-rain` | FAIL | FAIL | PASS | PASS | PASS |
| `mass/response-nonnegativity` | FAIL | FAIL | PASS | PASS | PASS |
| `mass/spinup-cycle-invariance` | FAIL | FAIL | PASS | PASS | PASS |
| `mass/catchment-closure` | FAIL | FAIL | PASS | PASS | PASS |
| `mass/runoff-bounds` | FAIL | FAIL | PASS | PASS | FAIL |
| `mass/warming-response` | FAIL | FAIL | PASS | FAIL | PASS |
| `mass/resolution-invariance` | FAIL | PASS | FAIL | PASS | PASS |
| `mass/snowpack-mass-closure` | PASS | FAIL | PASS | FAIL | PASS |
| `mass/phase-counterfactual` | PASS | PASS | FAIL | PASS | PASS |
| `mass/time-origin-invariance` | PASS | FAIL | PASS | PASS | PASS |

The top six fail everywhere, and are what this emulation does, whichever capacities it has. The rest of
the count is:

- **Five probes differ only by one gate.** On `area-invariance`, `causality`, `extreme-rain`,
  `response-nonnegativity` and `spinup-cycle-invariance`, the assertion the probe exists for passes at all
  corners, and `non_degenerate` fails because the worst seed's runoff ratio is 0.0044 to 0.0106 against a
  floor of 0.02 (0.01 on spin-up) at 320 mm, and 0.056 to 0.118 at 203 mm. `runoff-bounds` is close kin: its own
  bound passes (0.19 of the rain), and the same gate fails on runoff's correlation with rain, r 0.039 against 0.05.
  This is mostly one number counted several times. A probe with a dry forcing shows it and a probe with a wet one does not:
  the same packaged configuration has a ratio of 0.18 on `catchment-closure`.
- **Two probes are close to their limits.** `catchment-closure` is ET at 1.014 of PET on one seed of five (the
  others 0.905 to 0.961), and passes at 203 mm. `warming-response` is a runoff response of -0.091 per unit of
  demand against -0.1, and `raven_alt` misses too, at -0.0997 on one seed of five; there is little runoff left for
  warmer air to take away.
- **The canopy capacity moves three.** With unlimited capacity, 12% of the precipitation evaporates at once on
  `snowpack-mass-closure`, whose case says the canopy holds nothing, so that probe's budget is open by 12.0000%.
  At 2 mm the interception cap applies per step, so the resolution probe sees more interception at PT1H than
  at PT1D; see below. `phase-counterfactual` moves between 2.8%, 4.6% and 5.3% against a 5% limit.
- **One sits exactly on the gate.** `time-origin-invariance` passes at a ratio of 0.0259 as packaged and fails at
  0.0200 with the larger canopy.

Removing the lake closes the gap by another route, because it acts on the same number: on the probes that
fail this gate the ratio at least triples, and all but `ungauged-basin-closure` clear it. The six FAILs that are always there stay, and `runoff-bounds` stays,
which suggests the slow reservoir's smoothing rather than the lake sets that correlation (inferred).

### Why the packaged configuration stands

`raven_nolake` is not the template: it deletes a process that upstream's HBV-EC emulation lists. `raven_alt` is
nearer to the template's numbers, but it ignores two inputs the case declares. Its 203 mm soil leaves the
case's stated capacity on `mass/ungauged-basin-closure` (`mrso` up to 202 against 124.358, on 2086 steps), and its
unlimited canopy opens the snowpack budget on a case that declares no canopy at all. The benchmark's rule is that
capacities come from the case where it gives them (as in `models/flex_lumped` and `models/sacsma_snow17`), and the
adapter has done so since its first full evaluation, whose runoff ratios are these. So it is kept, and neither
alternative was picked for having fewer FAILs. Strictly, what is evaluated is the template with named edits, and
the two capacity edits are the ones with consequences. A maintainer who reads the issue's "HBV-EC emulation" as
excluding them should archive the `raven_alt` column instead, knowing that those PASSes come from a soil other than
the one the case declares.

### Verdicts

Gate-seed worst values. Verified means a targeted run or the source; inferred is marked.

| Probe | Verdict | Measured | Mechanism |
| --- | --- | --- | --- |
| `energy/pet-consistency` | FAIL | ET 1.89 of demand with the soil wettest | lake evaporation on top of soil evaporation. Verified: 1.04 without it, the rest being canopy evaporation, which `CANEVP_ALL` does not limit to PET |
| `mass/antecedent-monotonicity` | FAIL | +0.0 mm from a 60 mm storm, 1.2 mm expected; ratio 0.0136 | capillary rise returns what is generated, and the lake takes what survives. Verified: removing both gives +1.42 mm, and either alone does not suffice |
| `mass/area-invariance` | FAIL | ratio 0.0106; the invariance itself holds (1.0e-13) | gate, the two sinks. Verified |
| `mass/catchment-closure` | FAIL | ET 1.014 of PET on the worst seed; closure and bounds hold | lake evaporation. Verified: passes without it |
| `mass/causality` | FAIL | ratio 0.0070; causality itself holds | gate, the two sinks. Verified |
| `mass/dry-down` | FAIL | `mrro` exactly 0 on all 730 dry days, cv 0; `dry_down` itself holds | no drainage from the top soil, and the reservoirs are empty at the start. Verified on all three gate seeds |
| `mass/extreme-rain` | FAIL | ratio 0.0103; the monotone response holds, 1.00 of 997 mm returned | gate, the two sinks. Verified |
| `mass/human-abstraction` | FAIL | all 380 mm missed | the column is not read; no withdrawal process. Verified |
| `mass/multi-decadal-drift` | FAIL | ET 1.584 of PET; closure, drift and bounds hold | lake evaporation at 56% of PET on seed 1349990263. Verified: 1.003 without it, the excess being canopy evaporation |
| `mass/resolution-invariance` | FAIL | `evspsbl` moves 10.3% of the rain, PT1H against PT1D; limit 10% | below. Verified |
| `mass/response-nonnegativity` | FAIL | ratio 0.0044; the assertion holds | gate, the two sinks. Verified |
| `mass/runoff-bounds` | FAIL | r 0.039 between runoff and 7-day rain; bounds hold | runoff is rare and leaves through a slow reservoir at 0.063/d. Inferred: r stays low without the lake (0.032) and rises at 203 mm (0.124) |
| `mass/spinup-cycle-invariance` | FAIL | ratio 0.0077 against 0.01; the invariance holds | gate, the two sinks. Verified |
| `mass/ungauged-basin-closure` | FAIL | ET 1.539 of PET and ratio 0.0044 | lake evaporation and the two sinks. Verified |
| `mass/warming-response` | FAIL | runoff -0.091 per unit of demand, -0.1 expected; ratio 0.0132 | runoff is already nearly gone. Q/P 0.0118 falls to 0.0048 on seed 732983943. Verified: without the lake it is -0.201 and passes |
| `mass/extreme-event-closure` | PASS | residual 1.0e-12 mm | Raven's bookkeeping |
| `mass/phase-counterfactual` | PASS | 2.8% of the rain | real: snow and rain take different paths |
| `mass/precipitation-counterfactual` | PASS | partitioned; 0.55 of the top rung returned | real |
| `mass/snowpack-mass-closure` | PASS | 0.0000% | depends on the case's zero canopy capacity reaching the model |
| `mass/steady-state` | PASS | runoff 0.032 mm/d under 2.5 mm/d of rain | settles, by evaporating nearly all of it |
| `mass/time-origin-invariance` | PASS | departure 0 | |
| `momentum/routing-conservation` | PASS | within 0.09 of what a 15-day hydrograph can hold | a 2.3-day triangular hydrograph holds little |

N/A: `routing-lag-consistency` and `wave-celerity-bounds` are INCOMPATIBLE, because this configuration has
`ROUTE_NONE` and a `TIME_CONC` parameter, so the case's channel geometry has nothing to enter; the energy probes,
the groundwater-exchange probes, `mass/exchange-response` and the stage-based momentum probes are INCOMPLETE,
because the variables are not computed.

### Resolution

On seed 1760592044, whole case, in mm of the 368 mm of rain:

| | Canopy | Soil | Lake | Total |
| --- | --- | --- | --- | --- |
| PT1D | 28 | 47 | 47 | 122 |
| PT1H | 42 | 65 | 51 | 158 |
| PT1D, `:SuppressCompetitiveET` | 28 | 71 | 43 | 142 |
| PT1H, `:SuppressCompetitiveET` | 42 | 70 | 51 | 163 |

That is 36 mm, 9.7% of the rain, over the case; the probe's 10.3% is over its scored window. Both parts come from
Raven's process rules. The soil term comes from competition: at PT1D the day's interception uses up the
day's whole PET budget, and at PT1H only the rainy hours'; with competition off, the two steps give the same soil
evaporation, 71 and 70 mm. The canopy term comes from `MAX_CAPACITY` applying per step: `CANEVP_ALL` empties the
canopy every step, so 2 mm a day becomes 2 mm an hour. The adapter's budget closes exactly at both steps, and neither
`pr` nor `pet` is reshaped (`PET_DATA` is not a daily PET method, `Evaporation.cpp`). So this is a property of the
model, 0.3 points over a limit. At this probe's 120 mm soil and unlimited canopy it passes (at 203 mm and 2 mm it is
11.8%).

## The alternative configuration

`raven_alt`, the template's own 203 mm soil and unlimited canopy, is the alternative parameter set the benchmark
asks for: **14 PASS, 8 FAIL, 15 N/A**. The verdicts do move, as tabulated above. `raven_nolake` was a mechanism check
and is **15 PASS, 7 FAIL, 15 N/A**. The corner runs are 6 PASS, 16 FAIL and 14 PASS, 8 FAIL. None of those
directories is kept. Their rows are not archived.

## Caveats

- The parameters are Salmon River's. No PASS or FAIL here is about Raven's skill, or about Raven as a framework: its
  other emulations were not run, and the three mechanisms above belong to this configuration of it.
- The lake term takes a depth set by PET, so it grows with PET while the slow reservoir's inflow does not. A
  climate with nearly twice Salmon River's PET makes it large. It is present at Salmon River too (1.20 of PET), so
  the probes expose it rather than cause it.
- Whether HBV-EC's own lake term is meant to act on the whole catchment is not settled here. Raven's source applies
  no fraction, and the template sets no correction; a `LAKE_PET_CORR` below 1 would be a different parameter set
  and was not tried. That is a question for the upstream authors.
- The soil capacity is met through thickness, which also moves the tension storage (46% of it) that evaporation sees.
- apt packages in the build come from the Debian archive of the pinned base image, not pinned one by one.
  Versions are in the build log.

# HBV-EC (Raven emulation, Alouette2 tutorial model) under HydroTuring

HBV-EC as emulated by Raven, the hydrological modelling framework of Craig et al. (2020, *Environ. Model.
Softw.* 129, [10.1016/j.envsoft.2020.104728](https://doi.org/10.1016/j.envsoft.2020.104728)), packaged from
v4.12 at commit `81011f061e53ef1289ebcab36dc3622cea05c76d`. The configuration is Raven's own HBV-EC tutorial
model, `Alouette2` from [RavenTutorialFiles.zip](https://raven.uwaterloo.ca/files/RavenTutorialFiles.zip).
Proposed by Yan Chen as [Flood-Lab/HydroTuring#116](https://github.com/Flood-Lab/HydroTuring/issues/116), and
renamed and moved to the tutorial model following
[Flood-Lab/HydroTuring#207](https://github.com/Flood-Lab/HydroTuring/issues/207). Artistic-2.0; the image
carries Raven's licence at `/model/LICENSE`.

## Why the name and the configuration changed

This package was `raven` until 2026-11-19. Issue #207 raised two points, and both stand.

- **Raven is a framework, not a model.** It can assemble millions of model configurations, and HBV-EC is one.
  A verdict here is about one configuration, so the package is named for it.
- **The parameters were not a calibration.** `raven` ran `benchmarking/_InputFiles/Salmon_HBV` from Raven's
  repository, and its README called those parameters calibrated for the Salmon River. They are not. Those
  files exist to catch code changes that break backward compatibility. Their parameters were sampled at
  random, and they score an NSE of about 0.1 on the Salmon River, where a calibrated model should exceed 0.8.
  The issue suggested the tutorial's Alouette2 HBV-EC model instead, and that is what runs here.

The `raven` rows, version `4.12-81011f0.2`, stay in `models/result.csv` as the record of that configuration.
With the `models/raven` directory gone, the standings no longer count them. On the same gate seeds the
tutorial model passes 18 of 23 scored probes, against 8 of 23 for the benchmark files. The section "Against
the benchmark files" below gives the probe-by-probe comparison.

## What runs

Upstream's own CMake target, built unmodified at the pinned commit. Four of the tutorial's files,
`Alouette2/Alouette2.rv{i,p,h,c}`, are fetched at build time and pinned by SHA-256, then baked into the image
byte for byte. The archive's URL carries no version, so a change to any of the four fails the build. The archive
was last modified on 2025-01-30, and its SHA-256 as fetched on 2026-11-19 was `0544b97e…57fec91bd`.

`ht_adapter.py` writes the case's forcing as the gauge file, applies a named list of edits to the other four
files, runs the executable once, and reads what it wrote. It computes no hydrology. The base image is
`python:3.11-slim`, pinned by digest (Debian 13), and the build refuses any other commit. NetCDF output is used
for its doubles, because Raven's CSV carries six significant figures (`ht_adapter.py`). The container needs no
network at run time.

## Every edit to the tutorial

Each edit is asserted to match exactly once and is written to `run.json` under `edits_from_template`. The
tutorial files have Windows line endings, and the adapter reads them without.

| File | Tutorial | Used | Why |
| --- | --- | --- | --- |
| rvi | `:StartDate`, `:Duration`, `:TimeStep` | the case's | the case's record and step |
| rvi | `:Evaporation`, `:OW_Evaporation PET_FROMMONTHLY` | `PET_DATA` | the case prescribes `pet`; the tutorial derives PET from Alouette monthly normals |
| rvi | `:EvaluationMetrics NASH_SUTCLIFFE` | `:WriteNetCDFFormat` and three `:CustomOutput` | output only |
| rvi | `:DefineHRUGroups … Glacier` | without `Glacier` | that reporting group is empty once the glacier HRUs are removed |
| rvp | topsoil thickness 1.0 m (200 mm at `POROSITY` 0.2) | `soil_capacity_mm` / (1000 x `POROSITY`) | capacity from the case |
| rvp | `MAX_CAPACITY`, `MAX_SNOW_CAPACITY` 10000 mm, meaning no limit | `canopy_capacity_mm` | capacity from the case |
| rvh | 23 HRUs, 201.8 km2, latitude 49.3 | 21 HRUs, each area scaled to `area_km2`, `latitude_deg` | the case's area and latitude. HRUs 34 and 35, 1.4% of the area, are glacier, and the case has no glacier. Without the case's glacier extent, Raven would melt ice the case never gave it |
| rvh | `Band7` and `Glacier` HRU groups | without 34 and 35 | follows from the row above |
| rvc | initial outflow 1 m3/s | 0 | Raven would release Alouette's 1 m3/s in the first step without taking it from a store |
| rvc | initial `SOIL[2]` for 23 HRUs | the 21 kept | follows from the glacier HRUs' removal; the 6.50657 mm stays |
| rvt | gauge at 582.358 m, `:RainCorrection` 1.25885, `:SnowCorrection` 1.04283, monthly normals, Alouette series | gauge at the HRUs' area-weighted mean elevation, 653.94 m; corrections 1.0; the case's `tas`, `pr`, `pet` | see below |

Every HRU keeps its elevation, land use, aspect and share of the area. The process list, the soil model,
every other parameter and the routing are the tutorial's.

**The gauge's elevation.** `OROCORR_HBV` corrects each HRU's temperature, precipitation and PET by its height
above the gauge (`OrographicCorrections.cpp`). Temperature drops at `ADIABATIC_LAPSE`, 4 C/km. Precipitation
is multiplied by 1 + 0.08 per km of height (`HBVEC_LAPSE_RATE`, Raven's default). PET is multiplied by
1 - 0.0005 per m (`HBV_PET_ELEV_CORR`). All three are linear in the height. The precipitation factor steepens
above `HBVEC_LAPSE_ELEV`, which the tutorial leaves at Raven's 5000 m default, far above the top band at
1500 m; the adapter stops if either changes. So with the gauge at the HRUs' area-weighted mean
elevation, the catchment as a whole receives exactly the case's precipitation and PET, spread over the seven
bands the way the tutorial spreads them. Neither multiplier reaches its floor of zero. At the tutorial's
582.358 m, the catchment would take in 0.6% more precipitation than the case gave it, and see 3.6% less PET.
Checked on `mass/catchment-closure`, seed 1105848510, from Raven's own `ForcingFunctions` output: PET 8381.96 mm
against the case's 8381.98, with the difference down to how the rows are offset.

**The gauge corrections** are undercatch factors for the Alouette gauge. The case's rain is exact, so at the
tutorial's values Raven would take in 26% more rain than it was given.

## What is reported

| Column | Raven | Note |
| --- | --- | --- |
| `pr` | the forcing | echoed from the CSV as text |
| `snm` | `BETWEEN ATMOS_PRECIP and PONDED_WATER` + `BETWEEN SNOW_LIQ and PONDED_WATER` | differenced from Raven's running sum |
| `evspsbl` | `AET` custom output | canopy, canopy snow, soil and lake evaporation |
| `mrro` | `q_sim` over `area_km2` | |
| `dis` | `q_sim` | m3/s |
| `mrso` | `Soil Water[0]` | the top soil only |
| `gw` | `Soil Water[1]` + `Soil Water[2]` | the fast and slow reservoirs, shared by all HRUs (`:AggregatedVariable`) |
| `snw` | `Snow` + `Snow Melt (Liquid)` | |
| `canopy` | `Canopy` + `Canopy Snow` | |
| `channel` | rivulet, channel and reservoir storage, `Surface Water`, `Ponded Water` | the unit hydrograph's water is in rivulet storage |

Every store is Raven's watershed average over the HRUs. `zero_stores` and `not_reported` are both empty, so
every store is real and is reported. The adapter stops if `WatershedStorage.nc` holds a variable it has not
mapped, or if the two glacier stores hold water. No flux is filled in. `abstr` is not read, because the
tutorial has no withdrawal process to give it to. There is no `gwex`, because the tutorial exchanges nothing
across the catchment boundary. The tutorial's two passive tracers, `SNOWMELT` and `GLACIER`, move no water and
are left as they are.

## Parameters

No calibration exists for the synthetic catchment. Everything below is the tutorial's own, for the Alouette
River in British Columbia, and is unchanged. Results describe that one parameter set, with area and capacities
from the case.

| Parameter | Value | Note |
| --- | --- | --- |
| `RAINSNOW_TEMP`, `RAINSNOW_DELTA` | 1.0 C, 2.0 C | |
| `ADIABATIC_LAPSE` | 4.0 C/km | across seven bands, 161 m to 1500 m |
| `SNOW_SWI` | 0.05 | |
| `MELT_FACTOR`, `MIN_MELT_FACTOR` | 3.76, 2.2 mm/d/K | |
| `HBV_MELT_FOR_CORR`, `HBV_MELT_ASP_CORR` | 1.0, 0.48 | aspects 0, 180 and about 126 degrees |
| `REFREEZE_FACTOR` | 2.0 mm/d/K | |
| `POROSITY`, `FIELD_CAPACITY`, `SAT_WILT` | 0.2, 0.7, 0.0 | |
| `HBV_BETA` | 1.0 | recharge in proportion to soil wetness |
| `MAX_CAP_RISE_RATE` | 0.0 mm/d | no capillary rise |
| fast `MAX_PERC_RATE` | 5.0 mm/d | |
| fast `BASEFLOW_COEFF`, `BASEFLOW_N` | 0.2 /d, 1.1 | |
| slow `BASEFLOW_COEFF` | 0.05 /d | |
| fast and slow reservoir thickness | 100 m | effectively unbounded |
| `TIME_CONC`, `TIME_TO_PEAK` | 0.8, 0.35 d | triangular unit hydrograph |
| `TFRAIN`, `TFSNOW` | 0.88, 0.88 | Raven reads them as 12% of rain and of snow intercepted on the forest-covered share (`VegetationClass.cpp:369`), each capped every step at what the canopy can still hold (`PartitionPrecip.cpp`). The case's canopy capacity sets both `MAX_CAPACITY` and `MAX_SNOW_CAPACITY`, so a step with rain and snow can hold up to twice it. `CANEVP_ALL` then evaporates all of it in the same step, whatever the PET |
| `IMPERM` | 0.3% to 33% by band | the highest band is a third impermeable |
| `FOREST_COV` | 1 on the 7 forest HRUs, 0 on the 14 open ones | 48% of the area is forest |
| `LAKE_PET_CORR` | 1.0 | not set by the tutorial; Raven's default |

From the case: `area_km2`, `soil_capacity_mm`, `canopy_capacity_mm`, `latitude_deg`, and the forcing `pr`,
`tas`, `pet`.

## The step

Raven's rates are per day and its solver takes `:TimeStep` directly, so nothing is rescaled. Forcing is
written at the case's step, as rates. The manifest lists PT1D, the tutorial's step, and PT1H, which the
resolution probe exercises. The adapter accepts finer steps, but nothing has tested them, so they are not
claimed. At a daily step Raven warns that the gauge's `TEMP_AVE` series was overwritten. It is copied from the
daily series Raven copied from that same series (`Gauge.cpp:85`, `:677`), so the values do not change. At an hourly
step `RAINSNOW_HBV` and `POTMELT_HBV` read `temp_daily_ave`, which Raven takes over the day, so an hourly run's
phase and melt see later hours of the same day. That was read from the source and has not been tested, and no
hourly probe checks it.

## Evaluation

Gate seeds, all 40 probes, under Docker on Apple Silicon at `95f99b4`, version `4.12-81011f0-alouette2.1`:
**18 PASS, 5 FAIL, 17 N/A** (INCOMPLETE or INCOMPATIBLE). There are no ERROR rows. The whole suite takes about
five minutes. Closure holds everywhere it is measured, with a residual of 0.0000% of the rain. On every case, hourly ones
included, the adapter's own budget closes to under 1e-7 mm and Raven's mass-balance error is under 1e-9 mm.

### Verdicts

Gate-seed worst values.

| Probe | Verdict | Measured | Mechanism |
| --- | --- | --- | --- |
| `energy/pet-consistency` | FAIL | ET 1.78 of demand with the soil wettest | lake evaporation on top of soil evaporation |
| `mass/dry-down` | FAIL | `mrro` exactly 0 over the 730 dry days, cv 0; `dry_down` itself holds | runoff stops 3 to 9 days after the last rain. HBV recharges the reservoirs only from new water, and lake evaporation empties the slow reservoir. The top soil keeps 93 to 157 mm, and only evaporation takes it |
| `mass/human-abstraction` | FAIL | all 380 mm missed | the column is not read, and there is no withdrawal process |
| `mass/multi-decadal-drift` | FAIL | ET 1.593 of PET; closure, drift and bounds hold | lake evaporation, on a humid case with PET about 0.56 of the rain |
| `mass/ungauged-basin-closure` | FAIL | ET 1.516 of PET on the worst seed | lake evaporation, on the seed whose PET is half the rain |
| `mass/antecedent-monotonicity` | PASS | the wetter catchment runs off 6.8 mm more from the 60 mm storm | |
| `mass/area-invariance` | PASS | worst departure 9.6e-14 | |
| `mass/catchment-closure` | PASS | ET 0.888 of PET on the worst seed; runoff 0.21 to 0.26 of the rain | |
| `mass/causality` | PASS | nothing moves before the storm, and +0.48 of the added rain after it | |
| `mass/extreme-event-closure` | PASS | worst event residual 6.0e-11 mm | |
| `mass/extreme-rain` | PASS | 0.95 of the 997 mm added at the top is returned | |
| `mass/interception-response` | PASS | capacity 2 to 10 mm: evaporation +15.4 mm, runoff -17.3 mm | |
| `mass/phase-counterfactual` | PASS | turning snow into rain moves the volumes by at most 0.6% of the rain | |
| `mass/precipitation-counterfactual` | PASS | partitioned; 0.57 of the top rung returned | |
| `mass/resolution-invariance` | PASS | `evspsbl` moves 3.66% of the rain, PT1H against PT1D | |
| `mass/response-nonnegativity` | PASS | largest dip 2.0e-21 of the added rain | |
| `mass/runoff-bounds` | PASS | runoff 0.29 of the rain | |
| `mass/snowpack-mass-closure` | PASS | 0.0000% | |
| `mass/spinup-cycle-invariance` | PASS | departure 0.00% | |
| `mass/steady-state` | PASS | runoff 0.057 mm/d under 2.5 mm/d of rain | settles, by evaporating nearly all of it |
| `mass/time-origin-invariance` | PASS | departure 0 | |
| `mass/warming-response` | PASS | runoff -0.144 and -0.204 per unit of demand | |
| `momentum/routing-conservation` | PASS | within 0.06 of what a 15-day hydrograph can hold | |

N/A: `routing-lag-consistency` and `wave-celerity-bounds` are INCOMPATIBLE, because this configuration has
`ROUTE_NONE` and a `TIME_CONC` parameter, so the case's channel geometry has nothing to enter. The energy
probes, the groundwater-exchange probes, `mass/exchange-response`, `mass/routing-network-closure` and the
stage-based momentum probes are INCOMPLETE or INCOMPATIBLE, because those variables are not computed.

### Lake evaporation with no lake

The tutorial keeps the term issue #207 expected to fail. `:LakeEvaporation LAKE_EVAP_BASIC SLOW_RESERVOIR
ATMOSPHERE` takes the full `OW_PET` from the slow reservoir whenever it holds water (`OpenWaterEvap.cpp`,
`LAKE_PET_CORR` 1). No lake fraction scales it, and it does not subtract what the soil and canopy have already
evaporated. So total evaporation can exceed PET, and it does on the humid cases. Where PET is high against the
rain, the soil leaves less to the slow reservoir (inferred), and ET stays under PET: at most 0.888 on
`mass/catchment-closure`, where PET is 0.88 to 0.99 of the rain.

Removing that one line from the process list, as a mechanism check that is not archived, leaves one FAIL in
23: `mass/human-abstraction`. Without it, ET is 0.82 of PET on `mass/multi-decadal-drift` and on
`mass/ungauged-basin-closure`. `energy/pet-consistency` sees evaporation at 0.88 of demand on a wet soil and
0.33 on a dry one. `mass/dry-down` drains 18 mm over the 730 dry days instead of stopping. So four of the five
failures come from this one term.

## The tutorial's own capacities

Capacities come from the case wherever it gives them. That is the benchmark's rule, as in `models/flex_lumped`
and `models/sacsma_snow17`, and the packaged `raven` followed it too. As a mechanism check that is not archived,
the tutorial's own 200 mm topsoil and unlimited canopy give 16 PASS, 6 FAIL and 18 N/A. The same five fail,
and so does `mass/snowpack-mass-closure`. That case declares no canopy, and the unlimited canopy takes 5.8% of
the precipitation past the snow module, which opens its budget. `mass/interception-response` becomes N/A,
because a model that ignores the case's canopy capacity cannot be asked how it responds to one. So the
capacities move one verdict, and the case's capacities are the ones that pass it.

## Against the benchmark files

The same gate seeds, on the same 40 probes, with `raven` from the archive:

| Probe | `raven`, Salmon benchmark files | `hbvec_raven`, Alouette2 |
| --- | --- | --- |
| `mass/antecedent-monotonicity` | FAIL, +0.0 mm from the 60 mm storm on the wetter catchment | PASS, +6.8 mm |
| `mass/area-invariance` | FAIL, runoff ratio 0.0106 | PASS |
| `mass/catchment-closure` | FAIL, ET 1.014 of PET | PASS, ET 0.888 of PET at most |
| `mass/causality` | FAIL, runoff ratio 0.0070 | PASS |
| `mass/extreme-rain` | FAIL, runoff ratio 0.0103 | PASS |
| `mass/resolution-invariance` | FAIL, `evspsbl` 10.3% of the rain | PASS, 3.66% |
| `mass/response-nonnegativity` | FAIL, runoff ratio 0.0044 | PASS |
| `mass/runoff-bounds` | FAIL, runoff against 7-day rain r 0.039 | PASS |
| `mass/spinup-cycle-invariance` | FAIL, runoff ratio 0.0077 | PASS |
| `mass/warming-response` | FAIL, runoff -0.091 per unit of demand | PASS, -0.144 at least |

Nothing moved the other way. The five that fail here failed there too.

Five of the ten that moved failed only the `non_degenerate` runoff-ratio floor. That
package's runoff ratio was 0.004 to 0.011 on the drier cases: capillary rise of up to 18.9 mm/d returned the
water the soil generated, and lake evaporation took most of what reached the slow reservoir. The tutorial has
no capillary rise, and its `HBV_BETA` of 1.0 against 3.44 lets a soil at a third of capacity pass a third
of the rain to the reservoirs rather than about 2%. Its runoff ratio is 0.07 to 0.13 on the probes the floor
decided before, against a floor of 0.02, so the floor no longer decides anything. The lake term and the missing withdrawal process are in both configurations, and so are the
five failures that come from them.

## Caveats

- The parameters are the Alouette River's. No PASS or FAIL here is about Raven's skill, or about Raven as a
  framework. Its other emulations were not run, and the mechanisms above belong to this configuration.
- The glacier HRUs were removed because the case declares no glacier. With a glacier, `GMELT_HBV` melts ice at
  the potential-melt rate once the snow is gone, and the case would not have given Raven that water.
- Whether HBV-EC's lake term is meant to act on the whole catchment is not settled here. Raven's source applies
  no fraction, and the tutorial sets no correction. A `LAKE_PET_CORR` below 1 would be a different parameter
  set and was not tried.
- The soil capacity is met through thickness, which also moves the tension storage that evaporation sees
  (`FIELD_CAPACITY` 0.7 of it).
- apt packages in the build come from the Debian archive of the pinned base image, and are not pinned one by
  one. Versions are in the build log.

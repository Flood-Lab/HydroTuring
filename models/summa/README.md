# SUMMA under HydroTuring

SUMMA, the Structure for Unifying Multiple Modeling Alternatives (Clark et
al. 2015, *WRR*, [10.1002/2015WR017198](https://doi.org/10.1002/2015WR017198)),
compiled from [CH-Earth/summa](https://github.com/CH-Earth/summa) release
v4.0.0 (`f787fa5`, 2026-09-04) and run as one lumped HRU. Submitted as
[Flood-Lab/HydroTuring#30](https://github.com/Flood-Lab/HydroTuring/issues/30).

SUMMA solves the coupled conservation equations for water and energy in a
vegetation canopy, a layered snowpack, a layered soil column and an aquifer
with one implicit solver, and it reports its latent, sensible and ground heat
fluxes. It is the first submission that can be asked about both budgets and
the identity between them, rather than being declared INCOMPLETE on the
energy probes.

v4.0.1 was released on the day this was packaged. It adds brackets around
`iden_ice/iden_water` in four files (a last-bit change in floating point) and
reworks the power-law conductivity profile, which the decisions used here do
not select. Nothing in v4.0.0 blocked the build or the run, so the image stays
on the release named in the submission.

## Licence

SUMMA is GPL-3.0-or-later. The image carries `COPYING` at `/model/LICENSE`
and `LICENSE.txt` beside it. The adapter runs the unmodified executable.

## Build

A two-stage Dockerfile. Debian bookworm with gfortran, netCDF-Fortran,
OpenBLAS and CMake compiles the pinned tag (the build stops if the tag no
longer points at `f787fa5`) with the release's own CMake script in Release
mode, about 15 s at `-j4`; the runtime stage adds Python 3 with numpy 1.24 and
netCDF4 1.6. No SUNDIALS: the release documents it as optional, and the
decisions used here select SUMMA's own backward-Euler solver (`num_method
homegrown`), which needs none. No NextGen, no OpenWQ. The image is 371 MB and
builds on arm64 and on amd64 (checked with `docker buildx --platform
linux/amd64`). SUMMA's `-v` banner prints no version, because CMake's `git`
calls run outside the checkout; the pin is enforced by the commit check.

## What SUMMA is given

The release ships exactly one complete setup, the CAMELS basin 01073000 case
under `test_ngen/gauge_01073000/settings/SUMMA`. Its decisions, its default
parameter tables (`localParamInfo.txt`, `basinParamInfo.txt`) and its Noah-MP
lookup tables are copied into the image unedited and used as they are:

| Decision | Option | Decision | Option |
| --- | --- | --- | --- |
| `num_method` | `homegrown` | `fDerivMeth` | `analytic` |
| `f_Richards` | `mixdform` | `bcLowrSoiH` | `drainage` |
| `groundwatr` | `bigBuckt` | `spatial_gw` | `localColumn` |
| `hc_profile` | `constant` | `infRateMax` | `GreenAmpt` |
| `surfRun_SE` | `homegrown_SE` | `subRouting` | `timeDlay` |
| `bcUpprTdyn` | `nrg_flux` | `bcLowrTdyn` | `zeroFlux` |
| `nrgConserv` | `enthalpyForm` | `stomResist` | `BallBerry` |
| `soilStress` | `NoahType` | `LAI_method` | `monTable` |
| `canopySrad` | `BeersLaw` | `canopyEmis` | `difTrans` |
| `veg_traits` | `Raupach_BLM1994` | `windPrfile` | `logBelowCanopy` |
| `astability` | `louisinv` | `alb_method` | `conDecay` |
| `snowIncept` | `lightSnow` | `snowLayers` | `CLM_2010` |
| `compaction` | `anderson` | `thCondSnow` | `jrdn1991` |
| `thCondSoil` | `funcSoilWet` | `soilCatTbl`, `vegeParTbl` | `ROSETTA`, `USGS` |

Where the probe says nothing about the catchment, the test case's attributes
stand: USGS class 15 (mixed forest), ROSETTA class 3 (loam: porosity 0.399,
residual water 0.061, saturated conductivity 1.39e-6 m/s), a 10 m measurement
height, a 16 m canopy top, flat for radiation, longitude 0 with local time.
From `static.json`, only quantities with a direct SUMMA counterpart:

| Probe attribute | SUMMA quantity | How |
| --- | --- | --- |
| `latitude_deg`, `area_km2` | `latitude`, `HRUarea` | as given |
| `snow_threshold_degC` | `tempCritRain` | + 273.15 K |
| `soil_capacity_mm` | soil column depth | capacity / porosity, so the column holds at most the capacity; the shipped layer thicknesses down to that depth, a remainder under half the layer above merged into it (320 mm: 0.802 m in 5 layers; 180 mm: 0.451 m in 4; 120 mm: 0.301 m in 3). `rootingDepth` kept inside the column |
| `canopy_capacity_mm` | `refInterceptCapRain`, `refInterceptCapSnow` | capacity divided by the largest monthly LAI + SAI of the class (5.0 for mixed forest, so 0.4 kg m-2 per unit area); zero is a bare surface, USGS class 19, which SUMMA gives no leaves or stems |

The probes give no elevation, so the HRU sits at sea level. The cold state is
the test case's (283.16 K, 0.3 volumetric water, no snow, an empty canopy)
except the aquifer, which starts empty: under the shipped
`aquiferBaseflowRate` the shipped 0.4 m would drain within hours and put
400 mm of runoff into the spinup.

## Forcing the probes do not generate

SUMMA wants shortwave and longwave down, wind, pressure and specific humidity.
They are mocked from each forcing row alone: nothing reads the calendar, the
clock or another row, and `run.json` labels every one.

* **Net radiation** of the row: the probe's `rn` where it supplies one (the
  three energy probes that need it); otherwise the net radiation
  Priestley-Taylor (alpha 1.26) needs to produce the row's `pet` at the row's
  temperature.
* **Shortwave and longwave**: split so that a reference surface at air
  temperature, with albedo 0.23 and SUMMA's soil emissivity 0.96, would have
  exactly that net radiation and no longwave deficit. Longwave down is the
  reference surface's own emission, plus the net radiation where it is
  negative; shortwave carries a positive net radiation through the albedo.
  The split is linear in the net radiation, so the same day given as 24 hourly
  rows or as one daily row delivers the same radiant energy.
* **Humidity** 70 percent relative humidity; **pressure** a standard
  atmosphere at sea level; **wind** 2 m/s at the measurement height.

SUMMA then computes its own net radiation from that forcing, with its own
albedo (a wet loam is about 0.15), emissivity and surface temperature. That is
the central packaging fact for the energy probes: SUMMA's partition is scored
against a net radiation it never received as such. The section on the energy
probes measures how far apart the two are and how much of each failure the
gap explains. Every value above moves some result; the sensitivity section
says which.

## Time

Rows go in at the step they came at, as SUMMA's data step, and nothing is
resampled. SUMMA's forcing stamps are period-ending, so an hour is stamped at
its end. A day is stamped at 23:00 rather than 24:00. SUMMA computes the day's
mean cosine of the solar zenith angle in `CLRSKY_RAD` (`sunGeomtry.f90`,
called from `derivforce.f90`) over a window that starts at the stamp's hour
minus the data step. A 24:00 stamp is hour 0 of the next day, so the window
starts at hour -24, lies wholly before the day the routine integrates, finds
no daylight and returns zero; `vegSWavRad.f90` then sets the shortwave
reaching the canopy and the ground to zero. With 24:00 stamps every day of
the ten-year closure record had a zero zenith cosine and SUMMA absorbed no
shortwave: evaporation fell from 0.67 to 0.42 of demand on the closure probe,
and the latent-heat probe's energy residual rose from 4.8 to 97 percent of the
net radiation. Stamped at 23:00 the window covers the day's daylight. This is
a limitation of SUMMA at a daily data step worth reporting upstream.

SUMMA reads the clock itself, and that is the model's physics rather than the
adapter's: its solar geometry, and its phenology, which interpolates LAI and
SAI from a monthly table by day of year.

## What is reported

Fluxes are the step means SUMMA accumulates over its sub-steps, as rates;
states are end of step.

| Column | SUMMA output | Sign and units |
| --- | --- | --- |
| `pr` | the forcing, echoed as the text it arrived as | mm/day |
| `evspsbl` | `-(scalarTotalET + scalarSnowSublimation + scalarCanopySublimation)`; SUMMA's total ET leaves sublimation out | positive upward, kg m-2 s-1 to mm/day; net deposition makes it negative |
| `sbl` | `-(scalarSnowSublimation + scalarCanopySublimation)` | positive upward; frost deposition is negative |
| `mrro` | `averageRoutedRunoff`: surface runoff plus aquifer baseflow after the time-delay histogram | m/s to mm/day |
| `channel` | cumulative `averageInstantRunoff` minus `averageRoutedRunoff` | mm; SUMMA normalises the histogram to sum to one |
| `hfls` | `-scalarLatHeatTotal` | W m-2, positive away from the surface (SUMMA's are positive downward) |
| `hfss` | `-scalarSenHeatTotal` | W m-2, positive away from the surface |
| `hfg` | `scalarGroundNetNrgFlux` | W m-2, positive into the ground |
| `mrso` | `scalarTotalSoilWat` (liquid plus ice) plus the cumulative `scalarSoilCompress` since the cold state | mm |
| `snw` | `scalarSWE`, ice plus liquid in the pack, snow without a layer included | mm |
| `canopy` | `scalarCanopyLiq + scalarCanopyIce` | mm |
| `gw` | `scalarAquiferStorage` | m to mm |

`mrso` carries a term SUMMA's soil balance keeps outside the volumetric water
content: water stored by compressing the soil matrix under `specificStorage`
(1e-6 per metre of head). It is small but not nothing. A dry column that
rewets takes up to 3 mm into it in a day, and two rainless years release
3.9 mm from it; left out, it shows as exactly that residual in the water
budget.

`hfg` is the net energy flux into the snow-soil column through its top: the
soil surface when there is no snow, the snow surface when there is, so the
snowpack's heat storage and melt energy are inside it, as the correction rule
for storage above the soil surface requires. The canopy's net energy flux is
not ground and is not in it; it is 0.03 to 0.07 percent of the net radiation
on the energy probes. `scalarSfcMeltPond` is not added to `snw`: it records
melt SUMMA has already passed to the soil in the same step. There is no ponded
surface store in these decisions.

## The budgets as SUMMA keeps them

Water. Recomputed from the reported columns, the residual is at most 1.1e-6
mm on any step and 2.6e-6 mm over any record of the 108 archived cases, with
one exception: on one day of one latent-heat seed a trace of snowfall
(4.3e-5 mm/day at -3.8 C) never reaches the pack, and the record is short by
4e-5 mm. The `closure` criterion reports 0.0000 percent on every seed of every
probe that asks. SUMMA's own per-step balance diagnostics stay below
3.3e-9 for soil mass, 1e-9 for aquifer mass and 2.8e-4 for the energy of any
domain, in SUMMA's units.

Energy. SUMMA's own net radiation minus `hfls + hfss + hfg` is the heat
carried by precipitation at the wet-bulb temperature (SUMMA's advective heat
flux) plus the canopy's net energy flux. With no canopy nothing else remains.
With a canopy a term of up to 1.1 W/m2 on a single day remains, 0.2 W/m2 on
average in absolute value, which averages to under 0.001 W/m2 over the
record, as the heat stored in the canopy air space (a SUMMA state) would.
Against its own net radiation the surface budget closes to 0.09 to 0.24
percent on every seed of every energy probe.

## The energy probes: SUMMA's net radiation is not the probe's

Rerunning every gate seed with SUMMA's own series kept (`SUMMA_HT_DIAG=1`)
gives each energy criterion twice: as scored, against the probe's `rn`, and
against the net radiation SUMMA computed from the mocked forcing.

| Probe (seeds) | Mean `rn` | SUMMA's net radiation | Mean abs gap | As scored, against `rn` | Against SUMMA's own net radiation |
| --- | --- | --- | --- | --- | --- |
| latent-heat-et-consistency (5) | 89.2 to 91.1 W/m2 | 85.2 to 86.9 | 5.3 to 5.5 | `energy_closure` 4.68 to 4.84 %, passes | 0.18 to 0.19 % |
| evaporative-partition (3) | 75.7 to 80.4 | 72.3 to 76.9 | 5.0 to 5.2 | `energy_closure` 4.41 to 4.70 %, passes; `partition_shift` sums to -174, -262, -341 W m-2 day against limits of 141, 204, 267, fails | SUMMA's own net radiation fell by 183, 272, 352 W m-2 day over the drought window, all of it longwave; what remains is the precipitation heat term, +9 to +11 |
| surface-energy-closure (5) | 77.2 to 81.3 | 84.1 to 89.2 | 7.2 to 8.3 | `energy_closure_by_phase`: 16 to 22 of 28 blocks fail | 0 of 28 blocks fail on every seed; the worst block's residual is 0.7 to 1.6 W/m2 against a 2 W/m2 allowance |

The daily closure passes against `rn` by 0.16 to 0.32 percentage points, a
margin that depends on the mock (see the sensitivity table). The hourly phase
test fails on the gap alone. By day the bare loam's albedo is about 0.15
against the reference 0.23, so SUMMA absorbs more shortwave than the reference
surface, and its skin averages 1.3 K above the air (up to 2.7 K over a
daylight block) and emits more; by night it averages 0.35 K below (up to
1.2 K) and emits less. The night allowance is 2 W/m2, and a surface 0.4 K off
the air temperature is already 2 W/m2 of longwave. A model that computes its
own surface temperature cannot meet that against a net radiation prescribed
for a surface at air temperature, and a row-only mock cannot know SUMMA's
temperature in advance.

`partition_shift` is an identity only while net radiation is fixed. SUMMA's
is not: in the drought run the drying surface is 0.2 to 0.45 K warmer over
the four months, emits more longwave, and its net radiation falls by what the
three fluxes fail to cancel, to within the precipitation heat term. The
energy is conserved; it is the probe's premise, a net radiation independent
of the surface, that SUMMA does not share.

`flux_identity` fails on the model's own constants. SUMMA converts evaporated
water at `LH_vap` = 2.501e6 J/kg, the latent heat of vaporisation at 0 C,
whatever the temperature, and sublimated water at `LH_sub` = 2.8347e6, which
is exactly the probe's value. The reported latent heat equals
`LH_vap * E_liquid + LH_sub * E_ice` to 8e-4 W/m2. The probe asks for
2.501e6 - 2361 T, so every liquid-evaporation step with the air more than
about 5 C from freezing leaves the 0.5 percent tolerance: 970 to 1018 of 3650
days on the latent-heat probe and 151 to 183 of 1095 on the partition probe,
every one of them above 5.3 C in absolute temperature. That is what
`reference_constant_lambda` is built to show. The harness message leads with
two further findings about `sbl`: SUMMA's net sublimation is negative on 34 to
230 days, when frost deposits on the pack, and it is non-zero on 3 to 33 days
when canopy ice sublimates after the ground snow has gone, which the
criterion's snow test (ground `snw` only) cannot see.

## Verdict

```
### HydroTuring `summa` v4.0.0-f787fa5.1

FAIL (VIOLATION) · 8/18 probes passed · suite 0.1.0
```

It passes `pet-consistency`, `area-invariance`, `causality`, `extreme-rain`,
`response-nonnegativity`, `time-origin-invariance` (bit for bit),
`warming-response` and `routing-conservation`. For each failure, the
mechanism, and whether it is the model or a choice the packaging had to make:

| Probe | Failing criterion (worst seed) | Mechanism | Model or packaging |
| --- | --- | --- | --- |
| `energy/latent-heat-et-consistency` | `flux_identity`; `state_bounds` canopy up to 2.09 mm on 4 days | a latent heat of vaporisation held at its 0 C value; deposition and canopy-ice sublimation as above. The canopy: liquid above capacity drains at `canopyDrainageCoeff` (0.005 s-1), so a wet day ends a few hundredths of a mm above it | model (constants, drainage law); the canopy capacity is mapped from `static.json` |
| `energy/evaporative-partition` | `partition_shift`; `flux_identity`; `state_bounds` canopy (2 of 3 seeds) | SUMMA's net radiation responds to the drought through its surface temperature; constant latent heat; finite canopy drainage | model; the size of the shift residual depends on the wind mock |
| `energy/surface-energy-closure` | `energy_closure_by_phase`, 19 of 28 blocks | the gap between SUMMA's net radiation and `rn`: albedo by day, surface temperature by day and night; SUMMA's own budget passes every block | the mock cannot deliver `rn`, meeting the model's own surface temperature |
| `mass/catchment-closure` | `state_bounds` canopy up to 2.10 mm, 7 days | finite canopy drainage above capacity | model |
| `mass/steady-state` | `steady_state`: soil water 14.3 %, canopy 10.4 %, runoff 1.0 %; -0.038 mm/day unplaced | under constant weather the phenology still cycles LAI and SAI through the year, so transpiration, soil water and interception cycle with it; the forest evaporates 2.54 of the 2.5 mm/day rain, runoff is 0.0006 mm/day and the soil is still drying 14 mm a year in the third year (1.8 mm of it out of the soil's elastic storage), which is the unplaced residual | model (its phenology keys on the day of year); the bare surface passes |
| `mass/resolution-invariance` | runoff differs by 22.0 % of the rain (18.8 to 22.0 across seeds) | at the hourly step Green-Ampt infiltration excess makes 53.5 mm of surface runoff from bursts up to 43 mm/h, none at the daily step; and the daily step, fed a daily-mean temperature and humidity, draws 23 W/m2 of sensible heat from the air into evaporation where the hourly step returns 11 W/m2 to it, so ET is 42 mm lower hourly. Both steps receive nearly the same radiation (104 and 96 W/m2 of shortwave on average) | model (intensity-dependent infiltration, step-dependent turbulent exchange); how large depends on the humidity mock |
| `mass/antecedent-monotonicity` | +0.0002 to +0.0014 of the storm (0.02 needed) | the wetter month's extra 120 mm is evaporated before the storm (135 against 17 mm of ET in those 30 days, demand 157 mm), so both runs meet the storm with 121 and 123 mm in an 802 mm-deep column and neither drains within the month | model at this demand; passes with 90 percent humidity or a bare surface |
| `mass/dry-down` | runoff rises 7.2 % between weeks 1 and 2 (1 of 3 seeds) | a delayed drainage pulse: after the last rains the bottom layer keeps wetting (volumetric water 0.154 to 0.158 over four weeks) while the top dries, so free drainage out of the column rises from 0.0021 to 0.0028 mm/day, on a runoff of a few thousandths of a mm | model (Richards redistribution) |
| `mass/runoff-bounds` | `non_degenerate`: runoff/rain correlation 0.018 (1 of 5 seeds; 0.08 to 0.21 on the others) | 87 to 94 percent of the runoff leaves in March to May as melt drains through the loam column and the aquifer; summer and autumn rain is held and evaporated | model and configuration (soil class, vegetation); the bounds themselves pass |
| `mass/phase-counterfactual` | `non_degenerate`: correlation 0.038 (1 of 3 seeds) | the same melt-season release; `phase_invariance` itself passes, 3.4 to 3.7 percent | model and configuration |

## Sensitivity

The first gate seed of each probe, one choice changed at a time from the
evaluated configuration, set with the `SUMMA_HT_*` variables the adapter
documents. P is pass and F fail for the whole probe; the number is the
quantity that moves. The table was computed before the soil's elastic storage
was added to `mrso`. That moves `mrso` by at most 4 mm; of the quantities
below it moves only the steady-state variation, by about two points (12.3 to
14.3 percent on this seed), and no verdict.

| Probe, quantity | evaluated | clear-sky split | RH 0.5 | RH 0.9 | wind 1 | wind 4 | albedo 0.15 | albedo 0.30 | bare surface | 24:00 stamp |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| latent-heat, energy residual % | F 4.8 | F 9.5 | F 4.6 | F 5.4 | F 6.1 | F 3.7 | F 13.4 | F 4.4 | F 7.4 | F 97.0 |
| surface-energy-closure, blocks failing of 28 | F 19 | F 24 | F 23 | F 22 | F 13 | F 18 | F 15 | F 25 | F 19 | F 19 |
| evaporative-partition, shift residual (limit 0.05) | F 0.064 | F 0.065 | F 0.063 | F 0.065 | F 0.104 | F 0.043 | F 0.064 | F 0.064 | F 0.326 | F 0.064 |
| pet-consistency, wet-soil ET / demand | P 0.93 | F 0.87 | F 1.03 | F 0.59 | P 0.81 | F 1.03 | P 0.90 | P 0.95 | F 0.59 | F 0.62 |
| catchment-closure, canopy overshoot mm | F 0.023 | F 0.023 | F 0.012 | F 1.94 | F 0.029 | F 0.015 | F 0.023 | F 0.022 | F 0 | F 0.027 |
| steady-state, largest variation | F 12.3 % | F 10.7 % | F 9.2 % | F 34 % | F 15.4 % | F 7.3 % | F 11.5 % | F 10.9 % | P 0.0 % | F 14.7 % |
| resolution-invariance, % of rain | F 22.0 | F 21.7 | F 31.9 | P 7.5 | F 13.2 | F 31.6 | F 21.7 | F 22.3 | P 3.4 | P 9.8 |
| antecedent-monotonicity, share of storm | F 0.000 | F 0.000 | F 0.000 | P 0.051 | F 0.000 | F 0.000 | F 0.000 | F 0.000 | P 0.173 | F 0.017 |
| phase-counterfactual, % of rain | P 3.7 | P 3.5 | P 0.5 | F 5.3 | P 3.8 | P 3.7 | P 3.6 | P 3.7 | P 3.4 | P 3.3 |
| warming-response | P | P | P | P | P | P | P | P | P | P |

What that says. `flux_identity` fails in every column (900 to 1160 days
wherever shortwave reaches SUMMA, 424 under the midnight stamp), and so does
the hourly phase test: those are the model. The daily energy residual passes
or fails on the mock, between 3.7 and 13.4 percent, and its pass in the
evaluated configuration is not robust. `pet-consistency` passes in four of ten
columns and fails on either side: with dry air or strong wind the forest
evaporates more than the probe's demand, with moist air or no leaves too
little. The step dependence and the missing antecedent response shrink below
their limits only with humid air or a bare surface, because both remove the
evaporation that erases the difference. The bare surface passes steady state
because it has no phenology, and in exchange evaporates too little, condenses
more than it evaporates on some days (negative ET) and leaves a third of the
drought shift uncancelled. The evaluated albedo and wind are the FAO-56
reference surface's and the humidity is the value this repository's other
radiation mock uses; all three were fixed before this table was made. The
split was changed once, for the reason below.

## Faults found in the adapter before the verdict

Each was found by a failing criterion or a budget that did not close, traced,
fixed and re-run.

* **A midnight daily stamp** discarded all shortwave, as described under
  Time.
* **The first radiation split** took Brutsaert's clear-sky longwave and gave
  shortwave whatever net radiation was left, on any row, including night hours
  whose shortwave SUMMA throws away because its sun is down: the hourly
  surface probe failed all 28 blocks with 42 W/m2 residuals at night.
  Restricting shortwave to rows with positive net radiation cut that to 24
  blocks, but the clear-sky split is not linear in the row, so the same day
  given hourly received 150 W/m2 of shortwave against 215 given daily. The
  linear split replaced it; the resolution failure did not move with it (22.2
  against 22.4 percent in the same configuration), so it is not the split's.
* **Snow on the canopy.** With only the rain capacity mapped, SUMMA's default
  snow capacity held 8 to 9 mm of intercepted snow against the stated 2 mm.
  Both capacities now come from the probe's.
* **A thin bottom layer merged upward in the wrong place.** `layers[-2] +=
  layers.pop()` reads the target before the pop and writes it after, so the
  120 mm column became 0.376 m deep and held up to 145 mm. Only the 120 mm
  columns were affected.
* **The melt pond counted twice.** Adding `scalarSfcMeltPond` to `snw` left
  step residuals of 2.8 to 4 mm on melt days, cancelling the next day.
* **The soil's elastic storage left out.** With `mrso` as liquid plus ice
  only, rain on a dry column left a residual of up to 3.3 mm that day, and the
  two rainless years of the dry-down a drift of -3.9 mm; both equal
  `scalarSoilCompress` over the step to rounding. It is now inside `mrso`.

## Running it

```bash
ht verify-adapter --model summa
ht run --model summa --gate-seeds --csv models/result.csv
```

One ten-year daily case (4015 rows) takes about 9 s in the container, 8.5 s
of it in SUMMA; the hourly surface case (384 rows) 0.4 s. Every probe is
evaluated on its full record. `run.json` carries the mapping, the mocks, the
water residual, the soil's elastic storage change, SUMMA's balance
diagnostics and the net-radiation gap for every case.

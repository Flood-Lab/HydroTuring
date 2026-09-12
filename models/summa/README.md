# SUMMA under HydroTuring

SUMMA, the Structure for Unifying Multiple Modeling Alternatives (Clark et
al. 2015, *WRR*, [10.1002/2015WR017198](https://doi.org/10.1002/2015WR017198)),
compiled from [CH-Earth/summa](https://github.com/CH-Earth/summa) release
v4.0.0 (`f787fa5`, 2026-09-04) and run as one lumped HRU. Proposed by Yuanhang
Liu (Independent Researcher) as
[Flood-Lab/HydroTuring#30](https://github.com/Flood-Lab/HydroTuring/issues/30).

SUMMA solves the coupled conservation equations for water and energy in a
vegetation canopy, a layered snowpack, a layered soil column and an aquifer
with one implicit solver, and it reports its latent, sensible and ground heat
fluxes. It is the first submission that can be asked about both budgets and
the identity between them, rather than being N/A (INCOMPLETE) on the
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

A two-stage Dockerfile. Both stages start from the same `debian:bookworm-slim`
index, pinned by digest (it lists linux/amd64 and linux/arm64 images). The
build stage installs gfortran, netCDF-Fortran, OpenBLAS and CMake. It compiles
the pinned tag with the release's own CMake script in Release mode, which takes
about 13 s at `-j4`. The build stops if the tag no longer points at `f787fa5`.
CMake runs inside the checkout, so its `git` calls find the fetched tag. The
`-v` banner therefore reads `Version: v4.0.0` and
`Git Branch: tags/v4.0.0-0-gf787fa5` with the full hash. The first version ran
CMake outside the checkout, and its banner printed no version.

The runtime stage adds Python 3 with numpy and netCDF4. Its apt packages still
resolve from the bookworm archive at build time, so a Debian point release can
move their patch versions. The archived evaluation ran with:
- libgfortran 12.2.0
- netCDF-C 4.9.0 and netCDF-Fortran 4.5.4
- OpenBLAS 0.3.21
- Python 3.11.2, numpy 1.24.2 and netCDF4 1.6.2

The build uses no SUNDIALS. The release documents it as optional, and the
decisions used here select SUMMA's own backward-Euler solver (`num_method
homegrown`), which needs none. It also uses no NextGen and no OpenWQ.

The image is 371 MB. It was built and evaluated on arm64. The first version's
Dockerfile was also built for linux/amd64 with `docker buildx`.

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
stand. That means:
- USGS class 15 (mixed forest);
- ROSETTA class 3 (loam: porosity 0.399, residual water 0.061, saturated
  conductivity 1.39e-6 m/s);
- flat for radiation;
- longitude 0 with local time.

**Heights.** Canopy top and bottom are the heights SUMMA's own vegetation
table gives the class. They are `HVT` and `HVB` in `TBL_MPTABLE.TBL`, which
`pOverwrite` makes the defaults:
- mixed forest: 16 m and 10 m;
- barren class 19: 0 m and 0 m.

The measurement height written to the attributes is the one SUMMA applies.
`derivforce.f90` raises it to the canopy top plus `minMeasHeight` (1 m)
whenever `mHeight` is below the canopy top. The adapter therefore writes the
shipped 10 m, raised the same way: 17 m over the forest and 10 m over bare
ground. `run.json` records SUMMA's own `scalarAdjMeasHeight`, and it was 17 m
or 10 m, as written, in every archived case.

The first version wrote a 10 m height with a 16 m canopy top for every class.
SUMMA therefore used 17 m everywhere, including the bare hourly surface probe,
while the attributes and this page said 10 m.

From `static.json`, only quantities with a SUMMA counterpart are used:

| Probe attribute | SUMMA quantity | How |
| --- | --- | --- |
| `latitude_deg` | `latitude` | as given; required, and the adapter stops if it is missing rather than invent one (it drives SUMMA's solar geometry) |
| `area_km2` | `HRUarea` | as given |
| `snow_threshold_degC` | `tempCritRain` | translated from an air temperature to the wet-bulb temperature SUMMA compares it with; see below |
| `soil_capacity_mm` | soil column depth | See the soil note below the table. |
| `canopy_capacity_mm` | `refInterceptCapRain`, `refInterceptCapSnow` | capacity divided by the largest monthly LAI + SAI of the class (5.0 for mixed forest, so 0.4 kg m-2 per unit area); zero is a bare surface, USGS class 19, which SUMMA gives no leaves or stems |

Soil column depth is the capacity divided by the porosity, so the column holds
at most the capacity. The shipped layer thicknesses are kept down to that
depth, and a remainder under half the layer above is merged into it.
`rootingDepth` is kept inside the column.
- 320 mm: 0.802 m in 5 layers.
- 180 mm: 0.451 m in 4 layers.
- 120 mm: 0.301 m in 3 layers.

The probes give no elevation, so the HRU sits at sea level. The cold state is
the test case's (283.16 K soil at a matric head of -1 m, no snow, an empty
canopy) except the aquifer, which starts empty. Under the shipped
`aquiferBaseflowRate`, the shipped 0.4 m would drain within hours and put
400 mm of runoff into the spinup.

SUMMA does not start from the water content in the file. `check_icond.f90`
recomputes each soil layer's liquid water from its matric head through the van
Genuchten curve (`volFracLiq` in `soil_utils.f90`, with the ROSETTA loam's
parameters). The loam therefore starts at 0.3247 volumetric water, not the
shipped file's 0.3. The adapter writes that value, so the cold state it hands
SUMMA says what SUMMA starts from; SUMMA's output is the same bit for bit
either way. The water diagnostics in `run.json` check the first step against
that state.

### Rain and snow

The probes define their snow threshold as an air temperature: precipitation
on a day below it is snow. SUMMA's `tempCritRain` is not an air temperature.
`derivforce.f90` compares it with each row's wet-bulb temperature and splits
linearly over a ramp of `tempRangeTimestep` (2 K, shipped and kept) centred
there.

The adapter therefore translates the definition. `tempCritRain` is the
wet-bulb temperature SUMMA itself computes for air at the probe's threshold,
under this adapter's own humidity mock (70 percent) and standard pressure.
The adapter ports the two routines SUMMA uses:
- `SPHM2RELHM`, capped at saturation as `derivforce.f90` caps it;
- `WETBULBTMP` from `convert_funcs.f90`, a Newton iteration with SUMMA's own
  starting value, step, tolerance and iteration limit, and SUMMA's constants
  (610.8 Pa, 0.622, 273.16 K).

On every row of the 16 cases run with SUMMA's own series kept, the port
reproduces SUMMA's `scalarTwetbulb` to 4.0e-7 K. For the 0 C threshold every probe uses,
`tempCritRain` is **271.4656 K**, a wet bulb of -1.69 C.

The value is a constant for a case, so the mock stays row-only and
calendar-free. Every row carries the same relative humidity, so a row's wet
bulb falls below `tempCritRain` when its air falls below the threshold. What
remains is SUMMA's own 2 K ramp around it: a row at the threshold gets half
rain and half snow. Changing the humidity mock moves the translation with it,
as the definition requires: 270.297 K at 50 percent, 272.602 K at 90 percent.
`SUMMA_HT_THRESHOLD=air` restores the first two versions' direct mapping,
threshold + 273.15 K, for the sensitivity table.

What the translation changes, on the evaluated runs:
- **Snow share.** On the probes with a snow season, snow is 0.25 to 0.36 of
  SUMMA's precipitation, against 0.25 to 0.37 under the probe's own rule:
  0.93 to 1.03 times the rule's share. The direct mapping gave 1.22 to 1.41
  times, with snow at air temperatures up to 2.9 C.
- **Phase-counterfactual `warm` variant.** This variant lifts every day below
  the threshold to 1 C above it. SUMMA's snow falls from 0.27 to 0.36 of the
  control's precipitation to 0.021 to 0.027, so 90 to 93 percent of the
  control snow becomes rain. Under the direct mapping, `warm` kept 0.32 to
  0.41 of its precipitation as snow, 0.91 of the control's, so that probe did
  not test SUMMA's response to the phase change.
- **No snow where the probe's rule gives none.** The second version said SUMMA
  makes only a little snow there. That was wrong: under the direct mapping,
  `warm` kept 0.32 to 0.41 as snow, and hourly resolution-invariance rows made
  0.04 to 0.25. With the translation, the ramp still gives some snow to rows
  just above the threshold: 0.021 to 0.027 in `warm`, and 0.008 to 0.13 at the
  hourly step of resolution-invariance, where the daily step makes none.

### Rain on a freezing canopy

Centring the ramp on the probe's threshold also puts rain on days whose air is
a little below freezing, which the direct mapping had kept as snow. With the
translation, SUMMA rains at air temperatures down to about -1.1 C and snows up
to about +1.1 C. On those days SUMMA, run with the shipped setup's decisions
and default parameters, keeps the rain on its canopy as ice without limit.

**How.** Five steps in SUMMA's code, each set by the shipped setup:
1. **All rain is intercepted.** The shipped `modelDecisions.txt` does not set
   `cIntercept`, so SUMMA keeps its backwards-compatible default, `unDefined`
   (`mDecisions.f90`). `vegLiqFlux.f90` then sets rain throughfall to zero,
   whatever `refInterceptCapRain` is. The storage-proportional option,
   `storageFunc`, is not active.
2. **The water freezes on a cold canopy.** `tempAdjust.f90` splits canopy water
   into ice and liquid by canopy temperature, so a canopy just below 0 C holds
   its water almost entirely as ice.
3. **Only liquid drains.** `vegLiqFlux.f90` drains liquid above
   `scalarCanopyLiqMax`; ice does not drain.
4. **The ice capacity bounds only snowfall.** `scalarCanopyIceMax` (exposed
   LAI + SAI times `refInterceptCapSnow`, set in `coupled_em.f90`; 0.96 mm in
   February) limits the interception of falling snow in `canopySnow.f90`. It
   does not limit ice made from intercepted rain.
5. **Nothing unloads the ice.** The shipped decisions leave `snowUnload` at its
   default, `meltDripUnload`. That option unloads `snowUnloadingCoeff` times the
   ice, and the shipped `snowUnloadingCoeff` is 0 (`localParamInfo.txt`; SUMMA's
   own range runs to 1.5e-6 s-1). The ice leaves only by sublimation and, once
   it melts, with the melt drip.

A positive `snowUnloadingCoeff`, or the `windUnload` decision, would remove ice
from the canopy. Both are SUMMA options. Neither is changed here, because the
evaluation runs the shipped setup untuned.

**Traced** on latent-heat seed 788749541 with SUMMA's canopy fluxes written
out:
- On 2002-02-11 the air is -0.01 C and the canopy -0.20 C. SUMMA intercepts all
  30.7 mm of rain and freezes 30.4 mm of it; of 31.2 mm of snow it intercepts
  0.53 mm.
- Ice gain, sublimation and the liquid left (26.8, 4.2 and 0.25 mm) add up to
  31.2 mm: all of the rain plus about 0.5 mm of snow.
- Canopy ice reaches 26.8 mm, 27.9 times `scalarCanopyIceMax`.
- The ice sublimates at 2.3 to 4.2 mm/day for four days. When the air warms it
  melts and drains, and some ice goes with the melt drip.
- Over that record, freezing of intercepted rain adds 234 mm of canopy ice on
  the days ice grows, and intercepted snow adds 111 mm.

**On the five probes that score canopy bounds**, the canopy peaks at:
- 14.3 to 27.0 mm on latent-heat, on 23 to 36 days above capacity per seed;
- 13.4 to 20.3 mm on catchment-closure, on 30 to 42 days;
- 14.5 to 23.0 mm on the precipitation counterfactual, on 20 to 63 days;
- 5.3 to 18.9 mm on evaporative-partition, on 2 to 5 days;
- 13.5 to 20.1 mm on human-abstraction, on 22 to 51 days, the same in the
  natural and irrigated runs.

Every peak falls on a day within 0.7 C of freezing with 25 to 71 mm of
precipitation, and most days above capacity are within 1 C of freezing.
Under the direct mapping, the largest excess was 0.10 mm of liquid.

Cases of probes that do not score canopy bounds reach as high, and no verdict
moves with them. The archive maximum is 28.5 mm, on pet-consistency.
Warming-response, runoff-bounds and phase-counterfactual reach 21.6 to
26.3 mm.

**Attribution.** The failure is SUMMA with the shipped setup's decisions and
default parameters, reached through the translated threshold. It does not
come from the adapter's capacity mapping. Total interception ignores
`refInterceptCapRain`, and `refInterceptCapSnow` caps only snowfall, so SUMMA's
own default capacities would still accrete ice from frozen rain. The mocks set
how large the ice gets (see the sensitivity table). The threshold
translation, centred on the probe's threshold, is what brings rain to sub-zero
days. No placement of SUMMA's 2 K ramp avoids both rain below freezing and
snow above it.

## Forcing the probes do not generate

SUMMA wants shortwave and longwave down, wind, pressure and specific humidity.
They are mocked from each forcing row alone: nothing reads the calendar, the
clock or another row, and `run.json` labels every one.

**Net radiation of the row.** Where the probe supplies `rn` (the three energy
probes that need it), that is used. Otherwise the row's `pet` is multiplied by
the Priestley-Taylor (alpha 1.26) conversion from potential evaporation to
net radiation. The conversion is evaluated at one fixed reference temperature,
20 C, where it is 33.0 W/m2 per mm/day. 20 C is FAO-56's standard temperature
for its constant latent heat of 2.45 MJ/kg. It was fixed before any run with
it, and the sensitivity table sets 10 C, 30 C and the old per-row factor
beside it. The target is then linear in `pet`, so a day gives the same net
radiation whether it arrives as 24 hourly rows or as one daily row.

**Correction to the first version (`4.0.0-f787fa5.1`, commit `cd7a70c`).**
That version evaluated the conversion at each row's own temperature, and this
page and its commit message said the mock was step-consistent. It was not.
The factor falls as the air warms: 48.2 W/m2 per mm/day at 5 C, 41.4 at 10 C,
28.5 at 30 C. Hourly `pet` peaks in the warm hours. On resolution-invariance
the daily step therefore received 7 to 8 percent more net radiation and
shortwave than the same days given hourly. Over the first gate seed's full
record that was 79.6 against 74.0 W/m2 of target net radiation. With the fixed
reference the two are identical: 64.2 W/m2 over the full record, 64.7 over the
scored window.

**Shortwave and longwave.** They are split so that a reference surface at
air temperature, with albedo 0.23 and SUMMA's soil emissivity 0.96, would
have exactly the target net radiation and no longwave deficit. Longwave down
is the reference surface's own emission, plus the net radiation where it is
negative. Shortwave carries a positive net radiation through the albedo.

The shortwave is linear in the target. The longwave is not: it carries the
air's emission, sigma T^4, which a daily-mean temperature slightly
underestimates. On resolution-invariance the hourly longwave averages about
0.4 W/m2 above the daily. That is the step dependence the mock still has.

**Humidity** is 70 percent relative humidity at the row's temperature.
**Pressure** is a standard atmosphere at sea level.

**Wind** is 2 m/s at the height SUMMA applies: 17 m over the forest, 10 m over
bare soil. It is a round value, not a height-corrected reference: FAO-56's
2 m/s is measured 2 m above short grass.

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
its end. A day is stamped at 23:00 rather than 24:00.

SUMMA computes the day's mean cosine of the solar zenith angle in `CLRSKY_RAD`
(`sunGeomtry.f90`, called from `derivforce.f90`). It integrates over a window
that starts at the stamp's hour minus the data step. A 24:00 stamp is hour 0
of the next day, so the window starts at hour -24. It then lies wholly before
the day the routine integrates, finds no daylight and returns zero, and
`vegSWavRad.f90` sets the shortwave reaching the canopy and the ground to zero.

With 24:00 stamps every day of the ten-year closure record had a zero zenith
cosine, and SUMMA absorbed no shortwave:
- on the closure probe, evaporation fell from 0.68 to 0.45 of demand;
- on the latent-heat probe, the energy residual rose from 4.5 to 97 percent of
  the net radiation.

Stamped at 23:00 the window covers the day's daylight. This is a limitation of
SUMMA at a daily data step worth reporting upstream.

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
| `mrro` | `averageRoutedRunoff`: surface runoff plus aquifer baseflow after the time-delay histogram | m/s to mm/day |
| `channel` | cumulative `averageInstantRunoff` minus `averageRoutedRunoff` | mm; SUMMA normalises the histogram to sum to one |
| `hfls` | `-scalarLatHeatTotal` | W m-2, positive away from the surface (SUMMA's are positive downward) |
| `hfss` | `-scalarSenHeatTotal` | W m-2, positive away from the surface |
| `hfg` | `scalarGroundNetNrgFlux` | W m-2, positive into the ground |
| `mrso` | `scalarTotalSoilWat` (liquid plus ice) plus the cumulative `scalarSoilCompress` since the cold state | mm |
| `snw` | `scalarSWE`, ice plus liquid in the pack, snow without a layer included | mm |
| `canopy` | `scalarCanopyLiq + scalarCanopyIce` | mm |
| `gw` | `scalarAquiferStorage` | m to mm |

`sbl` is not reported. SUMMA's snow and canopy sublimation are net fluxes over
each step: on a step when more frost deposits than sublimates they are
negative. They are carried inside `evspsbl`, so deposition lowers it. The
contract's `sbl` is the non-negative share of `evspsbl` that left as ice, and
no non-negative share describes a step of net deposition. Clipping at zero
would misstate those steps, so the column is left out, as the contract
allows. The third version reported the signed flux as `sbl`, and its
`flux_identity` numbers are kept under What changed.

`mrso` carries a term SUMMA's soil balance keeps outside the volumetric water
content: water stored by compressing the soil matrix under `specificStorage`
(1e-6 per metre of head). It is small but not nothing. A dry column that
rewets takes up to about 3 mm into it in a day, and two rainless years release
about 4 mm from it. Left out, it shows as exactly that residual in the water
budget.

`hfg` is SUMMA's `scalarGroundNetNrgFlux`. That is absorbed shortwave, net
longwave, turbulent heat and precipitation heat at the ground surface
(`vegNrgFlux.f90:981`). SUMMA passes it as the upper boundary flux of its
snow-soil column (`iLayerNrgFlux(0)`, `snowSoilNrgFlux.f90:218`).

The ground surface is the top layer of that column (`mLayerTempTrial_1`,
`vegNrgFlux.f90:259`): the soil when there is no snow, the top snow layer
when there is. With snow, SUMMA also switches the ground's latent heat to
sublimation (`vegNrgFlux.f90:548-550`). So when snow lies, the surface
exchanging turbulent and radiative energy with the atmosphere is the snow
surface.

The contract names no separate snowpack heat term, so the flux is taken
there, with the pack's heat storage and melt energy inside it. Mapping it to
the soil-snow interface instead would leave the pack's storage and melt as
an unexplained surface-budget residual on every snow day. The only snow-free
energy probe, surface-energy-closure, is unaffected: its ground is always
soil.

The canopy is treated differently, and plainly so. Its heat storage and the
phase change of water intercepted on it are also energy held above the soil
surface, but they are in no column. SUMMA's `scalarCanopyNetNrgFlux` is not
added to `hfg`, which stays the flux SUMMA itself calls ground heat.

With the canopy freezing rain, that term is no longer small on single days:
- On the traced latent-heat record it reaches -117 W/m2 on the day 30 mm of
  rain froze there, as the canopy gives up the fusion heat.
- It exceeds 5 W/m2 in size on 49 days of that record, and averages
  -0.17 W/m2.
- Over the energy probes' records it is -0.27 to -0.07 percent of the net
  radiation.

No verdict moves with it: the daily closure criteria pass by 0.39 to 0.77
points. `run.json` reports the term.

`scalarSfcMeltPond` is not added to `snw`: it records melt SUMMA has already
passed to the soil in the same step. There is no ponded surface store in these
decisions.

## The budgets as SUMMA keeps them

**Water.** The archived run has 126 cases. Recomputed from the reported
columns, starting from the cold state SUMMA itself starts from, every step
of every case closes to 7.3e-6 mm and every record to 8.6e-6 mm, with four
exceptions. The first step, which the third version did not check, closes
to 2.7e-8 mm. Each exception is a single step, off by
2.6e-5 to 7.9e-5 mm, on a day with only a trace of precipitation (4.3e-5 to
0.0023 mm/day):
- one latent-heat seed;
- two `warm` phase-counterfactual cases;
- one runoff-bounds seed.

The `closure` criterion reports 0.0000 percent on every seed of every probe
that asks. SUMMA's own per-step balance diagnostics stay below 3.5e-9 for soil
mass, 1.0e-9 for aquifer mass and 2.8e-4 for the energy of any domain, in
SUMMA's units.

**Energy.** SUMMA's own net radiation minus `hfls + hfss + hfg` is,
algebraically, the canopy's net energy flux minus the heat carried by
precipitation at the wet-bulb temperature (SUMMA's advective heat flux). With
no canopy only the precipitation term remains. Against its own net radiation,
the surface budget closes to 0.04 to 0.15 percent on every seed of every
energy probe.

## The energy probes: SUMMA's net radiation is not the probe's

Every gate seed was rerun with SUMMA's own series kept (`SUMMA_HT_DIAG=1`).
Each energy criterion is then given twice: as scored, against the probe's
`rn`, and against the net radiation SUMMA computed from the mocked forcing.

| Probe (seeds) | Mean `rn` | SUMMA's net radiation | Mean abs gap | As scored, against `rn` | Against SUMMA's own net radiation |
| --- | --- | --- | --- | --- | --- |
| latent-heat-et-consistency (5) | 89.2 to 91.1 W/m2 | 85.2 to 86.9 | 5.3 to 5.5 | `energy_closure` 4.37 to 4.61 %, passes | 0.04 to 0.13 % |
| evaporative-partition (3) | 75.7 to 80.4 | 72.4 to 76.9 | 4.9 to 5.2 | `energy_closure` 4.23 to 4.52 %, passes; `partition_shift` sums to -174, -262, -341 W m-2 day against limits of 141, 204, 267, fails | 0.05 to 0.10 %; SUMMA's own net radiation fell by 183, 272, 352 W m-2 day over the drought window, all of it longwave, and what remains is the precipitation heat term, +9 to +11 |
| surface-energy-closure (5) | 77.2 to 81.3 | 83.9 to 89.0 | 7.0 to 8.2 | `energy_closure_by_phase`: 17 to 22 of 28 blocks fail | 0 of 28 blocks fail on every seed; the worst block's residual is 0.7 to 1.6 W/m2 against a 2 W/m2 allowance |

The daily closure passes against `rn` by 0.39 to 0.77 percentage points. That
margin depends on the mock (see the sensitivity table).

The hourly phase test fails on the gap alone, at the 10 m measurement height
SUMMA applies over bare soil. It failed 16 to 22 blocks in the first version,
when SUMMA applied 17 m. The probe is warm, so the threshold translation does
not touch it.
- By day the bare loam's albedo is 0.145 to 0.148 against the reference 0.23,
  so SUMMA absorbs more shortwave than the reference surface. Its skin averages
  1.1 to 1.3 K above the air over the sunlit hours, and emits more.
- By night the skin averages 0.6 to 0.7 K below the air, and emits less.
- The night allowance is 2 W/m2, and a surface 0.4 K off the air temperature is
  already 2 W/m2 of longwave.

A model that computes its own surface temperature cannot meet that against a
net radiation prescribed for a surface at air temperature. A row-only mock
cannot know SUMMA's temperature in advance.

`partition_shift` is an identity only while net radiation is fixed. SUMMA's
is not. In the drought run the drying surface is 0.2 to 0.45 K warmer over
the four months and emits more longwave. Its net radiation falls by what the
three fluxes fail to cancel, to within the precipitation heat term. The
energy is conserved; it is the probe's premise, a net radiation independent
of the surface, that SUMMA does not share.

`flux_identity` fails on the model's own constants. SUMMA converts evaporated
water at `LH_vap` = 2.501e6 J/kg, the latent heat of vaporisation at 0 C,
whatever the temperature. It converts sublimated water at `LH_sub` = 2.8347e6,
which is exactly the probe's value. The reported latent heat equals
`LH_vap * E_liquid + LH_sub * E_ice` to 9.1e-4 W/m2. The probe asks for 2.501e6 - 2361 T.

The adapter reports no `sbl`, so the criterion does two things:
- On snow-free steps, where no pack lies and none could fall, it holds latent
  heat to that equality.
- Where a pack is or could be present, it asks only that latent heat lie
  between the liquid and ice conversions of the reported evaporation.

SUMMA meets the interval on every such step. It fails the equality:
- on 931 to 962 of 3650 days on the latent-heat probe, 928 to 960 of them with
  the air more than 5.3 C from freezing;
- on 141 to 171 of 1095 on the partition probe, all but one of them above
  5.3 C.

The other 2 to 9 latent-heat days a seed are days when the canopy's ice
sublimates, about 1 mm/day, with no snow on the ground and none falling.
SUMMA converts that water at `LH_sub`, while the criterion's equality asks for
the liquid value. They are the worst steps, 7.3 to 9.4 times the tolerance,
and the ice comes from the rain frozen on the canopy (Rain on a freezing
canopy).

That is what `reference_constant_lambda` is built to show. The third version,
which reported the signed sublimation as `sbl`, failed the same criterion
through the split instead; its numbers are under What changed.

## Verdict

```
### HydroTuring `summa` v4.0.0-f787fa5.4

FAIL (VIOLATION) · 8/21 probes passed · suite 0.1.0
```

It passes eight probes:
- `pet-consistency`, with wet-soil evaporation 0.93 to 0.99 of demand against
  0.7;
- `area-invariance`, `causality`, `extreme-rain` and `response-nonnegativity`;
- `time-origin-invariance`, bit for bit;
- `warming-response` and `routing-conservation`.

Against `4.0.0-f787fa5.3` no verdict changed. The fourth version's changes
leave SUMMA's output the same bit for bit, and `flux_identity` still fails
without `sbl`. The count was out of 20 because the suite gained
`mass/human-abstraction`, which SUMMA fails because it has no human water use.
It is out of 21 since `mass/extreme-event-closure` merged. SUMMA fails that
probe on the canopy's `state_bounds` alone, holding up to 25 mm of ice against
a 2 mm capacity on 68 steps, while every wet event's water budget closes to
within 3e-4 of its allowance.

In `4.0.0-f787fa5.3` the threshold translation flipped no probe verdict
against `4.0.0-f787fa5.2`. It moved one failure inside a probe and changed the
size and mechanism of another:
- **Phase-counterfactual.** `phase_invariance` began failing on seed 25625370,
  and the `non_degenerate` failure on seed 1039529598 passed.
- **Canopy `state_bounds`.** On the four probes that scored canopy bounds then,
  the excess grew from hundredths of a millimetre of liquid to 4 to 25 mm of
  ice.

For each failure: the mechanism, and whether it is the model or a choice the
packaging had to make.

| Probe | Failing criterion (worst seed) | Mechanism | Model or packaging |
| --- | --- | --- | --- |
| `energy/latent-heat-et-consistency` | `flux_identity` (5 of 5 seeds); `state_bounds` canopy 12.3 to 25.0 mm above the 2 mm capacity (5 of 5) | a latent heat of vaporisation held at its 0 C value, and on 2 to 9 days a seed canopy ice sublimating at `LH_sub` with no snow on the ground (above). The canopy: rain frozen on it near 0 C, up to 27.0 mm of ice (section on the freezing canopy); a few warm days also end a few hundredths of a mm above capacity as liquid drains at 0.005 s-1 | model (constants; drainage law). The canopy ice is SUMMA with the shipped setup's decisions and default parameters (all rain intercepted, `snowUnloadingCoeff` 0), reached through the translated threshold |
| `energy/evaporative-partition` | `partition_shift` (3 of 3); `flux_identity` (3 of 3); `state_bounds` canopy 4.1 and 16.9 mm (2 of 3) | SUMMA's net radiation responds to the drought through its surface temperature; constant latent heat; canopy ice | model; the size of the shift residual depends on the wind mock. The canopy ice is attributed as for latent-heat |
| `energy/surface-energy-closure` | `energy_closure_by_phase`, 20 of 28 blocks (17 to 22 across seeds) | the gap between SUMMA's net radiation and `rn`: albedo by day, surface temperature by day and night; SUMMA's own budget passes every block | the mock cannot deliver `rn`, meeting the model's own surface temperature |
| `mass/catchment-closure` | `state_bounds` canopy 11.4 to 18.3 mm above capacity (5 of 5) | rain frozen on a sub-zero canopy: peaks of 13.4 to 20.3 mm on days of 31 to 45 mm at 0.1 to 0.5 C, on 30 to 42 days a seed | SUMMA with the shipped setup's decisions and default parameters, reached through the translated threshold; 0.02 mm under the direct mapping |
| `mass/precipitation-counterfactual` | `state_bounds` canopy 15.7 to 16.7 mm above capacity (3 of 3 seeds; up to 21 mm across the variants) | canopy ice as above, more in the wetter variants. The partition the probe is about passes: on the worst seed evaporation takes 0.24 to 0.30 of the added or removed rain, runoff 0.67 to 0.73 and storage 0.03, summing to 1.000 in each variant; on every seed runoff rises along the ladder, returning 0.70 to 0.73 of the rain added at the top | SUMMA with the shipped setup's decisions and default parameters, reached through the translated threshold; the bare surface passes |
| `mass/human-abstraction` | `human_abstraction`: all 380 mm of the prescribed withdrawal are missing from the budget difference, 100 % against a 5 % limit (3 of 3 seeds); `state_bounds` canopy 11.6 to 18.1 mm above capacity (3 of 3) | SUMMA has no abstraction or water-use process, so nothing reads the `abstr` column. The natural and irrigated runs are identical, and evaporation, runoff and storage each change by 0.0 mm. The canopy ice is the same in both runs, as on the other probes that score canopy bounds | model: SUMMA has no human water use, and the adapter does not invent one, so this is not a packaging gap. The canopy ice is attributed as for catchment-closure |
| `mass/steady-state` | `steady_state`: soil water 10.2 %, canopy 10.6 %, runoff 1.1 %; -0.029 mm/day unplaced (3 of 3) | under constant weather the phenology still cycles LAI and SAI through the year, so transpiration, soil water and interception cycle with it; the forest evaporates 2.53 of the 2.5 mm/day rain, runoff is 0.0008 mm/day and the soil is still drying 10.7 mm a year in the third year, which is the unplaced residual | model (its phenology keys on the day of year); the bare surface passes |
| `mass/resolution-invariance` | runoff differs by 20.9 % of the rain (17.1 to 20.9 across seeds) | at the hourly step Green-Ampt infiltration excess makes 54.2 mm of surface runoff from bursts up to 43 mm/h, against 0.5 mm at the daily step. The daily step, fed a daily-mean temperature and humidity, draws 33 W/m2 of sensible heat from the air into evaporation, where the hourly step returns 5 W/m2 to it, so ET is 41 mm lower hourly. Both steps receive the same net radiation (64.7 W/m2 target, 66.0 in SUMMA). The hourly step also makes some snow the daily step does not: 0.008, 0.13 and 0.029 of its precipitation on the three seeds against none daily (0.04 to 0.25 against 0 to 0.009 under the direct mapping), from sub-zero April hours and SUMMA's ramp. Retained snow lowers hourly runoff, against the direction of the deviation, so it does not explain the failure | model (intensity-dependent infiltration, step-dependent turbulent exchange); how large depends on the humidity and wind mocks |
| `mass/antecedent-monotonicity` | +0.0002 to +0.0025 of the 60 mm storm (0.02 needed) | the wetter month's extra 120 mm is evaporated before the storm: 137 to 156 mm of ET in those 30 days against 19 to 42 mm in the drier run, on 131 to 157 mm of demand. Both runs meet the storm with soil water within 2 to 5 mm of each other (120 to 132 mm in an 802 mm-deep column), and neither drains within the month | model at this demand; passes with 90 percent humidity or a bare surface |
| `mass/dry-down` | runoff rises 8.2 % between weeks 1 and 2 (1 of 3 seeds) | a delayed drainage pulse: after the last rains the column's free drainage keeps rising for four weeks, from 0.0022 to 0.0028 mm/day, on a runoff of a few thousandths of a mm. At the first version SUMMA's layer water showed the bottom layer wetting while the top dried | model (Richards redistribution) |
| `mass/runoff-bounds` | `non_degenerate`: runoff/rain correlation 0.014 (1 of 5 seeds; 0.12 to 0.21 on the others) | 90 percent of that seed's runoff leaves in March to May as melt drains through the loam column and the aquifer. Its snow is now the probe's own: 0.340 of the precipitation against the rule's 0.343 | the humidity mock decides it, now through evaporation alone: at 90 percent humidity the seed passes (0.071), at 50 percent it fails further (0.007); under the direct mapping it failed at 0.019; the bounds themselves pass |
| `mass/phase-counterfactual` | `phase_invariance`: runoff changes by -5.2 % of the rain when snow falls as rain (1 of 3 seeds; 3.5 and 4.4 % on the others; limit 5 %) | the `warm` variant now turns 93 percent of that seed's control snow into rain (snow 0.286 to 0.021 of precipitation), and SUMMA's runoff responds by just over the limit. `non_degenerate` passes on every seed (0.085 to 0.206) | model, close to the limit: the seed passes at 90 percent humidity (4.9 %) and 1 m/s wind (5.0 %) and fails at 50 percent humidity (10.0 %) or 4 m/s (7.1 %). Under the direct mapping `warm` kept 0.91 of the snow and the criterion passed at 2.4 to 3.3 %, testing nothing |

## Sensitivity

The first gate seed of each probe was rerun, changing one choice at a time
from the evaluated configuration, with the `SUMMA_HT_*` variables the adapter
documents. P is pass and F fail for the whole probe; the number is the
quantity that moves. The Priestley-Taylor columns change only the probes
without `rn` and are left blank for the three that supply it. "Direct
threshold" is the first two versions' mapping, `tempCritRain` = threshold +
273.15 K.

| Probe, quantity | evaluated | direct threshold | PT at row temperature | PT at 10 C | PT at 30 C | clear-sky split | RH 0.5 | RH 0.9 | wind 1 | wind 4 | albedo 0.15 | albedo 0.30 | bare surface | 24:00 stamp |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| latent-heat, energy residual % | F 4.5 | F 4.8 | | | | F 9.3 | F 4.1 | F 5.3 | F 6.0 | F 3.4 | F 13.1 | F 4.6 | F 7.4 | F 96.7 |
| surface-energy-closure, blocks failing of 28 | F 20 | F 20 | | | | F 26 | F 23 | F 21 | F 12 | F 18 | F 15 | F 28 | F 20 | F 20 |
| evaporative-partition, shift residual (limit 0.05) | F 0.064 | F 0.064 | | | | F 0.065 | F 0.063 | F 0.065 | F 0.104 | F 0.043 | F 0.064 | F 0.064 | F 0.347 | F 0.064 |
| pet-consistency, wet-soil ET / demand; runoff-rain r (min 0.05) | P 0.94; 0.090 | P 0.87; 0.062 | P 0.99; 0.089 | P 0.99; 0.087 | P 0.91; 0.092 | F 0.86; 0.038 | F 1.13; 0.045 | F 0.56; 0.169 | P 0.79; 0.114 | F 1.05; 0.064 | P 0.92; 0.091 | P 0.96; 0.088 | F 0.49; 0.170 | F 0.68; 0.125 |
| catchment-closure, canopy excess mm | F 17.4 | F 0.024 | F 17.2 | F 17.3 | F 17.7 | F 17.3 | F 23.8 | F 5.9 | F 14.3 | F 19.7 | F 17.5 | F 17.3 | F 0 | F 20.8 |
| precipitation-counterfactual, canopy excess mm | F 15.7 | F 0.071 | F 13.3 | F 15.4 | F 15.8 | F 15.8 | F 17.2 | F 8.7 | F 11.5 | F 18.5 | F 15.7 | F 15.7 | P 0 | F 17.5 |
| steady-state, largest variation | F 10.6 % | F 10.6 % | F 14.3 % | F 16.1 % | F 10.7 % | F 10.9 % | F 14.3 % | F 36.0 % | F 16.0 % | F 7.3 % | F 10.7 % | F 12.2 % | P 0.0 % | F 14.7 % |
| resolution-invariance, % of rain | F 20.9 | F 21.0 | F 21.9 | F 21.3 | F 20.7 | F 20.7 | F 32.6 | P 6.4 | F 11.9 | F 30.7 | F 20.8 | F 21.0 | P 1.2 | F 10.6 |
| antecedent-monotonicity, share of storm | F 0.000 | F 0.000 | F 0.000 | F 0.000 | F 0.000 | F 0.000 | F 0.000 | P 0.051 | F 0.000 | F 0.000 | F 0.000 | F 0.000 | P 0.160 | F 0.019 |
| runoff-bounds, runoff / rain; r | P 0.41; 0.21 | P 0.43; 0.20 | P 0.40; 0.21 | P 0.40; 0.21 | P 0.42; 0.21 | P 0.41; 0.13 | P 0.33; 0.19 | P 0.59; 0.27 | P 0.46; 0.23 | P 0.36; 0.20 | P 0.42; 0.21 | P 0.40; 0.21 | P 0.59; 0.24 | P 0.61; 0.25 |
| phase-counterfactual, % of rain; r | F 5.2; 0.17 | P 3.3; 0.12 | F 5.6; 0.16 | F 5.2; 0.17 | F 5.2; 0.17 | F 5.9; 0.05 | F 10.0; 0.13 | P 4.9; 0.20 | P 5.0; 0.18 | F 7.1; 0.16 | F 5.1; 0.17 | F 5.3; 0.17 | F 6.6; 0.18 | F 5.9; 0.18 |
| warming-response | P | P | P | P | P | P | P | P | P | P | P | P | P | P |

Because the evaluated failures of runoff-bounds and phase-counterfactual sit
on seeds other than the first, the threshold and humidity columns were also
run on every gate seed of both probes. The seeds that fail in either mapping:

| Gate seed, quantity | evaluated (RH 0.7) | direct threshold | RH 0.5 | RH 0.9 |
| --- | --- | --- | --- | --- |
| runoff-bounds 607076768, runoff-rain r (min 0.05) | F 0.014 | F 0.019 | F 0.007 | P 0.071 |
| snow share of precipitation (probe's rule 0.343) | 0.340 | 0.432 | 0.340 | 0.340 |
| share of runoff in March to May | 0.90 | 0.92 | 0.94 | 0.77 |
| runoff / rain | 0.41 | 0.44 | 0.33 | 0.58 |
| phase-counterfactual 25625370, runoff change % of rain (limit 5) | F 5.2 | P 3.3 | F 10.0 | P 4.9 |
| control snow share (probe's rule 0.307) | 0.286 | 0.388 | 0.285 | 0.286 |
| phase-counterfactual 1039529598, runoff-rain r (min 0.05) | P 0.085 | F 0.039 | F 0.040 | P 0.112 |
| runoff change % of rain | 4.4 | 2.8 | F 13.4 | 3.7 |

Every other gate seed of the two probes passes in all four columns.

With the translated threshold, the snow share no longer moves with humidity,
so the humidity columns now act through evaporation alone. Runoff-bounds'
failing seed still fails with the probe's own snowpack, and passes only when
moist air leaves more rain to run off.

What the table says:

- **Model findings.** `flux_identity` fails in every column, on 853 to 1175
  days wherever shortwave reaches SUMMA and 387 under the midnight stamp. The
  hourly phase test fails in every column too. Both are the model.
- **The threshold mapping.** It moves the canopy excess by three orders of
  magnitude on the closure probes, from 0.02 to 0.07 mm to 16 to 17 mm. It
  decides phase-counterfactual on the first seed: 5.2 percent evaluated
  against 3.3 percent direct. It changes no other verdict in the table.
- **The canopy excess** under the translated threshold also depends on the
  mocks. Moist air cuts it to a third (5.9 mm at 90 percent humidity), dry air or
  strong wind raises it (23.8 mm, 19.7 mm), and a bare surface removes it.
- **The daily energy residual** passes or fails on the mock, between 3.4 and
  13.1 percent, so its pass in the evaluated configuration is not robust.
- **The Priestley-Taylor reference** flips no verdict at 10 C, 30 C or the
  row's own temperature. It moves resolution-invariance between 20.7 and 21.9
  percent, steady state between 10.6 and 16.1 percent, pet-consistency's
  wet-soil ratio between 0.91 and 0.99, and the phase response between 5.2
  and 5.6 percent.
- **`pet-consistency`** passes in 8 of 14 columns and fails on either side:
  - with a clear-sky split or dry air, runoff barely follows the rain (r under
    0.05);
  - with dry air or strong wind, the forest evaporates more than the demand;
  - with moist air, no leaves or the midnight stamp, it evaporates too little.
- **The phase response** sits at the 5 percent limit on the first seed in most
  columns. It passes with moist air or light wind and fails with dry air,
  strong wind, a bare surface, a clear-sky split or a midnight stamp.
- **The step dependence and the missing antecedent response** fall below their
  limits only with humid air or a bare surface, because both remove the
  evaporation that erases the difference.
- **The bare surface** passes steady state (it has no phenology) and the
  precipitation counterfactual (it has no canopy). In exchange it evaporates
  too little, condenses more than it evaporates on some days (negative ET),
  and leaves a third of the drought shift uncancelled.
- **How the evaluated values were chosen.** Each was fixed before any table was
  made.
  - The albedo is FAO-56's reference-grass value.
  - The wind is a round 2 m/s, not a height-corrected reference.
  - The humidity is the value this repository's other radiation mock uses.
  - The 20 C reference is FAO-56's standard temperature.
  - The threshold translation follows from the probes' definition and SUMMA's
    own wet-bulb routine.
  - The split was changed once, for the reason below.

## What changed in 4.0.0-f787fa5.4

The pull request's review raised three points about the adapter, and main
gained a twentieth probe. Each change and its effect:

- **`sbl` is no longer reported.**
  - Why: SUMMA's snow and canopy sublimation are net fluxes, negative when
    frost deposits. The contract's `sbl` is a non-negative share of
    `evspsbl`, and no such share describes a step of net deposition. Clipping
    at zero would misstate those steps.
  - Effect: `evspsbl` is unchanged, with deposition inside it.
    `flux_identity` still fails on both energy probes that score it, now on
    the snow-free equality: 931 to 962 latent-heat days a seed and 141 to 171
    partition days.
  - History: with the signed flux reported as `sbl`, the third version failed
    the same criterion at 4.95 to 5.64 times the tolerance on latent-heat and
    4.69 to 5.15 on partition. Its message led with `sbl` negative on 174
    steps and non-zero where the criterion saw no snow on 36 (latent-heat's
    worst seed); on partition's worst seed the counts were 25 and 3.
- **`hfg` stays at the top of the snow-soil column.** The reason is now given
  with SUMMA's Fortran under What is reported.
- **The water diagnostics start from the cold state.**
  - Before: the first step was never checked.
  - Finding: checking it showed that SUMMA does not start from the file's
    water content. `check_icond.f90` recomputes it from the matric head, which
    gives 0.3247 for loam, not 0.3.
  - Change: the adapter now writes that value. SUMMA's output is the same bit
    for bit, and the first step closes (numbers under The budgets as SUMMA
    keeps them).
- **The proposer is named**: Yuanhang Liu (Independent Researcher), here and
  in CONTRIBUTORS.
- **The branch merged main again**, which brings `mass/human-abstraction`.
- **No verdict changed.** SUMMA fails the new probe, so the count is 8 of 20.

## What changed in 4.0.0-f787fa5.3

A second review found that the adapter mistranslated the rain-snow threshold,
which the second version had documented but not corrected. It also found two
smaller faults. Each fault, the change made and its effect:

- **The threshold was an air temperature written into a wet-bulb parameter.**
  - Change: `tempCritRain` is now the wet-bulb temperature SUMMA computes for
    air at the threshold under the adapter's humidity mock (271.4656 K), as
    described under Rain and snow. The direct mapping is kept as a sensitivity
    column.
  - Effect on snow: SUMMA's snow share went from 1.22–1.41 to 0.93–1.03 times
    the probe's rule. The `warm` phase-counterfactual variant now removes 90
    to 93 percent of SUMMA's control snow instead of 9 percent.
  - Effect on phase-counterfactual: the verdict is unchanged, but
    `phase_invariance` now fails on one seed (5.2 percent) and
    `non_degenerate` now passes on every seed.
  - Effect on runoff-bounds: the failing seed's correlation went from 0.019 to
    0.014.
  - Effect on the canopy: rain now falls on sub-zero days. SUMMA, with the
    shipped setup's decisions and default parameters, intercepts all of it and
    keeps the ice that forms. Canopy `state_bounds` on four probes fails by 4
    to 25 mm where it failed by 0.06 to 0.10 mm.
  - Effect on energy: daily closure against `rn` fell to 4.23–4.61 percent,
    from 4.41–4.84.
- **This page said SUMMA made only a little snow where the probe's rule makes
  none.** The claim is withdrawn, with the numbers under Rain and snow.
- **The steady-state figure in the second version's change log was wrong**
  (12.4 percent). It is corrected below to 14.3 percent.
- **The resolution-invariance row now gives the step-dependent snow share.**
- **The canopy-ice mechanism was first described through an inactive branch of
  SUMMA's rain interception.** It is now traced through the options the shipped
  setup selects, under Rain on a freezing canopy.
- **The branch merged main twice more.** The first merge brought wflow_sbm.
  The second brought the fix that opens the antecedent-monotonicity window on
  the storm, and LF line endings.
  - Summa's antecedent-monotonicity row is re-scored under the new window. It
    still fails: 0.0002 to 0.0025 of a 60 mm storm, where the old window gave
    0.0002 to 0.0013 of 116 mm.
  - `models/result.csv` is main's file byte for byte, with summa's rows
    appended.

No probe verdict changed. The count is 8 of 19, as it was.

## What changed in 4.0.0-f787fa5.2

A review of the first version found three documentation faults that each
bore on how a verdict is attributed, and three small ones. Each fault, the
change made and its effect:

- **The Priestley-Taylor inversion was not step-consistent**, as described
  under Forcing.
  - Change: the conversion is evaluated at a fixed 20 C.
  - Effect: no verdict moved. The resolution-invariance deviation went from
    22.0 to 21.0 percent. The steady-state variation went from 14.3 to 10.6
    percent. pet-consistency's wet-soil ratios went from 0.87-0.99 to
    0.82-0.94 of demand. The dry-down rise went from 7.2 to 8.2 percent on its
    one failing seed.
- **SUMMA applied 17 m where the attributes said 10 m**, bare soil included.
  - Change: the heights now come from SUMMA's vegetation table, the attributes
    carry the height SUMMA applies, and `run.json` records
    `scalarAdjMeasHeight`.
  - Effect: the bare hourly probe now runs at 10 m. It fails 17 to 22 blocks
    against `rn` where it failed 16 to 22, and still none against SUMMA's own
    net radiation.
- **The rain-snow split is on the wet bulb.**
  - Change: this was documented, and the runoff-bounds and
    phase-counterfactual attributions named the humidity mock. The threshold
    itself was translated only in the third version.
- **Canopy storage** is stated plainly to be outside `hfg`; the mapping is
  unchanged.
- **`latitude_deg` is required** instead of defaulting to 40 N. Every probe
  supplies it.
- **Both Dockerfile stages are pinned by digest.** CMake runs inside the
  checkout, so the banner prints the release, and the Dockerfile comment about
  the version string is corrected.

## Faults found in the adapter before the verdict

Each was found by a failing criterion or a budget that did not close, traced,
fixed and re-run.

* **A midnight daily stamp** discarded all shortwave, as described under
  Time.
* **The first radiation split** took Brutsaert's clear-sky longwave and gave
  shortwave whatever net radiation was left, on any row. That included night
  hours, whose shortwave SUMMA throws away because its sun is down, so the
  hourly surface probe failed all 28 blocks with 42 W/m2 residuals at night.
  Restricting shortwave to rows with positive net radiation cut that to 24
  blocks. But the clear-sky split is not linear in the row, so the same day
  given hourly received 150 W/m2 of shortwave against 215 given daily. The
  linear split replaced it. The resolution failure did not move with it (22.2
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
  only, rain on a dry column left a residual of up to 3.3 mm that day. The
  two rainless years of the dry-down left a drift of -3.9 mm. Both equal
  `scalarSoilCompress` over the step to rounding, and it is now inside `mrso`.
* **The snow threshold as an air temperature**, as described under What
  changed in 4.0.0-f787fa5.3.

## Running it

```bash
ht verify-adapter --model summa
ht run --model summa --gate-seeds --csv models/result.csv
```

One ten-year daily case (4015 rows) takes about 10 s in the container, most
of it in SUMMA; the hourly surface case (384 rows) about 0.5 s. Every probe
is evaluated on its full record. `run.json` carries, for every case:
- the mapping and the mocks;
- the heights from the table and the height SUMMA applied;
- the translated `tempCritRain`;
- the cold-state storage, the water residual from the first step on, and the
  soil's elastic storage change;
- SUMMA's snow share against the probe rule's, and the canopy's net energy
  flux;
- SUMMA's balance diagnostics and the net-radiation gap.

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

From `static.json`, only quantities with a direct SUMMA counterpart are used:

| Probe attribute | SUMMA quantity | How |
| --- | --- | --- |
| `latitude_deg` | `latitude` | as given; required, and the adapter stops if it is missing rather than invent one (it drives SUMMA's solar geometry) |
| `area_km2` | `HRUarea` | as given |
| `snow_threshold_degC` | `tempCritRain` | + 273.15 K; a wet-bulb threshold in SUMMA, see below |
| `soil_capacity_mm` | soil column depth | See the soil note below the table. |
| `canopy_capacity_mm` | `refInterceptCapRain`, `refInterceptCapSnow` | capacity divided by the largest monthly LAI + SAI of the class (5.0 for mixed forest, so 0.4 kg m-2 per unit area); zero is a bare surface, USGS class 19, which SUMMA gives no leaves or stems |

Soil column depth is the capacity divided by the porosity, so the column holds
at most the capacity. The shipped layer thicknesses are kept down to that
depth, and a remainder under half the layer above is merged into it.
`rootingDepth` is kept inside the column.
- 320 mm: 0.802 m in 5 layers.
- 180 mm: 0.451 m in 4 layers.
- 120 mm: 0.301 m in 3 layers.

**Rain and snow.** SUMMA does not compare `tempCritRain` with the air
temperature. `derivforce.f90` compares it with the wet-bulb temperature and
splits linearly over a ramp of `tempRangeTimestep`. That ramp is 2 K in the
shipped `localParamInfo.txt`, and it is kept as shipped. The wet bulb comes
from the humidity mock below (70 percent). So the probe's air-temperature
threshold becomes a wet-bulb threshold with a 2 K ramp, and SUMMA makes snow
at air temperatures up to 2.9 C.

On the probes with a snow season, snow is 0.32 to 0.46 of SUMMA's
precipitation. The probe's own rule (air below the threshold) gives 0.25 to
0.37 on the same forcing, so SUMMA gets 1.22 to 1.41 times as much snow. On
the five latent-heat seeds, snow falls on 126 to 155 days with the air at or
above 0 C, 678 to 860 mm in all. Humidity moves the share. On the five
runoff-bounds seeds it is:
- 0.38 to 0.43 at the evaluated 70 percent;
- 0.45 to 0.49 at 50 percent;
- 0.33 to 0.37 at 90 percent;
- 0.30 to 0.34 under the probe's rule.

Where the probe's rule gives no snow at all, SUMMA still makes a little: up to
4.4 percent of the precipitation on dry-down. `run.json` reports the share for
every case.

The probes give no elevation, so the HRU sits at sea level. The cold state is
the test case's (283.16 K, 0.3 volumetric water, no snow, an empty canopy)
except the aquifer, which starts empty. Under the shipped
`aquiferBaseflowRate`, the shipped 0.4 m would drain within hours and put
400 mm of runoff into the spinup.

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
shortwave than the same days given hourly: 79.6 against 74.1 W/m2 of target
net radiation on the first gate seed. With the fixed reference the two are
identical, 64.2 and 64.2 W/m2.

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
- on the closure probe, evaporation fell from 0.65 to 0.42 of demand;
- on the latent-heat probe, the energy residual rose from 4.8 to 97 percent of
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
rewets takes up to about 3 mm into it in a day, and two rainless years release
about 4 mm from it. Left out, it shows as exactly that residual in the water
budget.

`hfg` is SUMMA's ground heat flux: the net energy flux into the snow-soil
column through its top. That top is the soil surface when there is no snow and
the snow surface when there is, so the snowpack's heat storage and melt energy
are inside `hfg`.

The canopy is treated differently, and plainly so. Its heat storage and the
phase change of water intercepted on it are also energy held above the soil
surface, but they are in no column. SUMMA's `scalarCanopyNetNrgFlux` is not
added to `hfg`, which stays the flux SUMMA itself calls ground heat. That term
is 0.03 to 0.07 percent of the net radiation over the energy probes' records,
and about 1 W/m2 at most on a single day. It moves no verdict, and `run.json`
reports it.

`scalarSfcMeltPond` is not added to `snw`: it records melt SUMMA has already
passed to the soil in the same step. There is no ponded surface store in these
decisions.

## The budgets as SUMMA keeps them

**Water.** The archived run has 120 cases. Recomputed from the reported
columns, their residual is at most 1.1e-6 mm on any step and 2.7e-6 mm over
any record, with two exceptions:
- On one day of one latent-heat seed, a trace of snowfall (4.3e-5 mm/day at
  -3.8 C) never reaches the pack. The record is short by 4.1e-5 mm.
- On one melt day of one variant of the precipitation counterfactual, with
  406 mm of snow on the ground, the step is off by 8e-6 mm.

The `closure` criterion reports 0.0000 percent on every seed of every probe
that asks. SUMMA's own per-step balance diagnostics stay below 3.2e-9 for soil
mass, 9.7e-10 for aquifer mass and 3.5e-4 for the energy of any domain, in
SUMMA's units.

**Energy.** SUMMA's own net radiation minus `hfls + hfss + hfg` is,
algebraically, the canopy's net energy flux minus the heat carried by
precipitation at the wet-bulb temperature (SUMMA's advective heat flux). With
no canopy only the precipitation term remains. Against its own net radiation,
the surface budget closes to 0.09 to 0.24 percent on every seed of every
energy probe.

## The energy probes: SUMMA's net radiation is not the probe's

Every gate seed was rerun with SUMMA's own series kept (`SUMMA_HT_DIAG=1`).
Each energy criterion is then given twice: as scored, against the probe's
`rn`, and against the net radiation SUMMA computed from the mocked forcing.

| Probe (seeds) | Mean `rn` | SUMMA's net radiation | Mean abs gap | As scored, against `rn` | Against SUMMA's own net radiation |
| --- | --- | --- | --- | --- | --- |
| latent-heat-et-consistency (5) | 89.2 to 91.1 W/m2 | 85.2 to 86.9 | 5.3 to 5.5 | `energy_closure` 4.68 to 4.84 %, passes | 0.18 to 0.19 % |
| evaporative-partition (3) | 75.7 to 80.4 | 72.3 to 76.9 | 5.0 to 5.2 | `energy_closure` 4.41 to 4.70 %, passes; `partition_shift` sums to -174, -262, -341 W m-2 day against limits of 141, 204, 267, fails | SUMMA's own net radiation fell by 183, 272, 352 W m-2 day over the drought window, all of it longwave; what remains is the precipitation heat term, +9 to +11 |
| surface-energy-closure (5) | 77.2 to 81.3 | 83.9 to 89.0 | 7.0 to 8.2 | `energy_closure_by_phase`: 17 to 22 of 28 blocks fail | 0 of 28 blocks fail on every seed; the worst block's residual is 0.7 to 1.6 W/m2 against a 2 W/m2 allowance |

The daily closure passes against `rn` by 0.16 to 0.59 percentage points. That
margin depends on the mock (see the sensitivity table).

The hourly phase test fails on the gap alone, now at the 10 m measurement
height SUMMA applies over bare soil. It failed 16 to 22 blocks when SUMMA
applied 17 m.
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
`LH_vap * E_liquid + LH_sub * E_ice` to 8.4e-4 W/m2. The probe asks for
2.501e6 - 2361 T, so every liquid-evaporation step with the air more than
about 5 C from freezing leaves the 0.5 percent tolerance:
- 970 to 1018 of 3650 days on the latent-heat probe;
- 151 to 183 of 1095 on the partition probe;
- every one of them above 5.3 C in absolute temperature.

That is what `reference_constant_lambda` is built to show. The harness
message leads with two further findings about `sbl`:
- SUMMA's net sublimation is negative on 34 to 230 days, when frost deposits
  on the pack.
- It is non-zero on 3 to 33 days when the criterion's snow test (ground `snw`
  only) sees no snow. Some of those days are canopy ice sublimating after the
  ground snow has gone, and a few follow snow that the wet-bulb split let fall
  above 0 C.

## Verdict

```
### HydroTuring `summa` v4.0.0-f787fa5.2

FAIL (VIOLATION) · 8/19 probes passed · suite 0.1.0
```

It passes eight probes:
- `pet-consistency`, with wet-soil evaporation 0.82 to 0.94 of demand against
  0.7;
- `area-invariance`, `causality`, `extreme-rain` and `response-nonnegativity`;
- `time-origin-invariance`, bit for bit;
- `warming-response` and `routing-conservation`.

Against `4.0.0-f787fa5.1` no verdict changed. The count is now out of 19
because the suite gained `mass/precipitation-counterfactual`, which SUMMA
fails on its canopy alone.

For each failure: the mechanism, and whether it is the model or a choice the
packaging had to make.

| Probe | Failing criterion (worst seed) | Mechanism | Model or packaging |
| --- | --- | --- | --- |
| `energy/latent-heat-et-consistency` | `flux_identity` (5 of 5 seeds); `state_bounds` canopy up to 2.09 mm on 4 days (5 of 5) | a latent heat of vaporisation held at its 0 C value; deposition and canopy-ice sublimation as above. The canopy: liquid above capacity drains at `canopyDrainageCoeff` (0.005 s-1), so a wet day ends a few hundredths of a mm above it | model (constants, drainage law); the canopy capacity is mapped from `static.json` |
| `energy/evaporative-partition` | `partition_shift` (3 of 3); `flux_identity` (3 of 3); `state_bounds` canopy (2 of 3) | SUMMA's net radiation responds to the drought through its surface temperature; constant latent heat; finite canopy drainage | model; the size of the shift residual depends on the wind mock |
| `energy/surface-energy-closure` | `energy_closure_by_phase`, 20 of 28 blocks (17 to 22 across seeds) | the gap between SUMMA's net radiation and `rn`: albedo by day, surface temperature by day and night; SUMMA's own budget passes every block | the mock cannot deliver `rn`, meeting the model's own surface temperature |
| `mass/catchment-closure` | `state_bounds` canopy up to 2.10 mm, 7 days (5 of 5) | finite canopy drainage above capacity | model |
| `mass/precipitation-counterfactual` | `state_bounds` canopy up to 2.07 mm on 3 days (3 of 3 seeds; up to 2.11 mm across the variants) | the same drainage law. Every day above capacity is liquid, on heavy summer or June rain (13 to 110 mm/day at 9.6 to 26.5 C), and the wetter variants overshoot more. The partition the probe is about passes: on the worst seed evaporation takes 0.23 to 0.28 of the added or removed rain, runoff 0.69 to 0.74 and storage about 0.03, summing to 1.000 in each variant, and on every seed runoff rises along the ladder, returning 0.71 to 0.75 of the rain added at the top | model (drainage law); the bare surface passes |
| `mass/steady-state` | `steady_state`: soil water 10.2 %, canopy 10.6 %, runoff 1.1 %; -0.029 mm/day unplaced (3 of 3) | under constant weather the phenology still cycles LAI and SAI through the year, so transpiration, soil water and interception cycle with it; the forest evaporates 2.53 of the 2.5 mm/day rain, runoff is 0.0008 mm/day and the soil is still drying 10.7 mm a year in the third year, which is the unplaced residual | model (its phenology keys on the day of year); the bare surface passes |
| `mass/resolution-invariance` | runoff differs by 21.0 % of the rain (18.0 to 21.0 across seeds) | at the hourly step Green-Ampt infiltration excess makes 54.1 mm of surface runoff from bursts up to 43 mm/h, against 0.5 mm at the daily step. The daily step, fed a daily-mean temperature and humidity, draws 33 W/m2 of sensible heat from the air into evaporation, where the hourly step returns 4 W/m2 to it, so ET is 41 mm lower hourly. Both steps now receive the same net radiation (64.7 W/m2 target, 66.1 in SUMMA); under the first version's per-row factor the daily step received 7 to 8 percent more and the deviation was 22.0 % | model (intensity-dependent infiltration, step-dependent turbulent exchange); how large depends on the humidity and wind mocks |
| `mass/antecedent-monotonicity` | +0.0002 to +0.0013 of the storm (0.02 needed) | the wetter month's extra 120 mm is evaporated before the storm: 137 to 156 mm of ET in those 30 days against 19 to 42 mm in the drier run, on 131 to 157 mm of demand. Both runs meet the storm with soil water within 2 to 5 mm of each other (120 to 132 mm in an 802 mm-deep column), and neither drains within the month | model at this demand; passes with 90 percent humidity or a bare surface |
| `mass/dry-down` | runoff rises 8.2 % between weeks 1 and 2 (1 of 3 seeds) | a delayed drainage pulse: after the last rains the column's free drainage keeps rising for four weeks, from 0.0022 to 0.0028 mm/day, on a runoff of a few thousandths of a mm. At the first version SUMMA's layer water showed the bottom layer wetting while the top dried | model (Richards redistribution) |
| `mass/runoff-bounds` | `non_degenerate`: runoff/rain correlation 0.019 (1 of 5 seeds; 0.08 to 0.20 on the others) | 92 percent of that seed's runoff leaves in March to May as melt drains through the loam column and the aquifer. The wet-bulb split gives that seed 0.43 of its precipitation as snow, where the probe's rule gives 0.34 | the humidity mock decides it, through the wet-bulb snow split and evaporation: at 90 percent humidity the seed passes (0.065) and at 50 percent it fails further (0.014); the bounds themselves pass |
| `mass/phase-counterfactual` | `non_degenerate`: correlation 0.039 (1 of 3 seeds) | the same melt-season release: 95 percent of that seed's runoff in March to May, snow 0.46 of its precipitation against 0.37 under the probe's rule; `phase_invariance` itself passes, 2.4 to 3.3 percent | the humidity mock decides it, the same way: 0.093 at 90 percent humidity, 0.041 at 50 percent |

The first version attributed the last two failures to "model and
configuration (soil class, vegetation)". The per-seed runs below show that
both turn on the humidity mock.

## Sensitivity

The first gate seed of each probe was rerun, changing one choice at a time
from the evaluated configuration with the `SUMMA_HT_*` variables the adapter
documents. P is pass and F fail for the whole probe; the number is the
quantity that moves. The three Priestley-Taylor columns change only the
probes without `rn` and are left blank for the three that supply it.

| Probe, quantity | evaluated (20 C) | PT at row temperature | PT at 10 C | PT at 30 C | clear-sky split | RH 0.5 | RH 0.9 | wind 1 | wind 4 | albedo 0.15 | albedo 0.30 | bare surface | 24:00 stamp |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| latent-heat, energy residual % | F 4.8 | | | | F 9.5 | F 4.6 | F 5.4 | F 6.1 | F 3.7 | F 13.4 | F 4.4 | F 8.3 | F 97.0 |
| surface-energy-closure, blocks failing of 28 | F 20 | | | | F 26 | F 23 | F 21 | F 12 | F 18 | F 15 | F 28 | F 20 | F 20 |
| evaporative-partition, shift residual (limit 0.05) | F 0.064 | | | | F 0.065 | F 0.063 | F 0.065 | F 0.104 | F 0.043 | F 0.064 | F 0.064 | F 0.347 | F 0.064 |
| pet-consistency, wet-soil ET / demand; runoff-rain r (min 0.05) | P 0.87; 0.062 | P 0.93; 0.060 | P 0.92; 0.059 | P 0.85; 0.063 | F 0.84; 0.019 | F 0.98; 0.014 | F 0.54; 0.153 | P 0.74; 0.081 | F 0.98; 0.043 | P 0.86; 0.063 | P 0.89; 0.061 | F 0.48; 0.121 | F 0.62; 0.087 |
| catchment-closure, canopy overshoot mm | F 0.024 | F 0.023 | F 0.023 | F 0.024 | F 0.024 | F 0.013 | F 2.15 | F 0.030 | F 0.016 | F 0.024 | F 0.023 | F 0 | F 0.027 |
| precipitation-counterfactual, canopy overshoot mm | F 0.071 | F 0.070 | F 0.070 | F 0.071 | F 0.071 | F 0.061 | F 3.65 | F 0.076 | F 0.064 | F 0.071 | F 0.070 | P 0 | F 0.074 |
| steady-state, largest variation | F 10.6 % | F 14.3 % | F 16.1 % | F 10.7 % | F 10.9 % | F 14.3 % | F 36.0 % | F 16.0 % | F 7.3 % | F 10.7 % | F 12.2 % | P 0.0 % | F 14.7 % |
| resolution-invariance, % of rain | F 21.0 | F 22.0 | F 21.4 | F 20.7 | F 20.8 | F 31.2 | P 6.4 | F 11.9 | F 30.8 | F 20.8 | F 21.1 | P 1.8 | F 10.6 |
| antecedent-monotonicity, share of storm | F 0.000 | F 0.000 | F 0.000 | F 0.000 | F 0.000 | F 0.000 | P 0.043 | F 0.000 | F 0.000 | F 0.000 | F 0.000 | P 0.155 | F 0.017 |
| runoff-bounds, runoff / rain; r | P 0.43; 0.20 | P 0.42; 0.21 | P 0.42; 0.20 | P 0.44; 0.20 | P 0.43; 0.11 | P 0.38; 0.16 | P 0.59; 0.26 | P 0.48; 0.22 | P 0.38; 0.19 | P 0.44; 0.21 | P 0.43; 0.20 | P 0.60; 0.24 | P 0.63; 0.23 |
| phase-counterfactual, % of rain; r | P 3.3; 0.12 | P 3.7; 0.12 | P 3.3; 0.12 | P 3.3; 0.12 | P 3.3; 0.07 | P 0.8; 0.08 | P 4.9; 0.19 | P 3.4; 0.13 | P 3.5; 0.11 | P 3.2; 0.12 | P 3.2; 0.12 | P 2.8; 0.17 | P 3.3; 0.12 |
| warming-response | P | P | P | P | P | P | P | P | P | P | P | P | P |

The first seed of runoff-bounds and of phase-counterfactual passes in every
column, so the humidity columns were also run on every gate seed of both. The
seeds that fail in the evaluated configuration:

| Failing gate seed, quantity | RH 0.7 (evaluated) | RH 0.5 | RH 0.9 |
| --- | --- | --- | --- |
| runoff-bounds 607076768, runoff-rain r (min 0.05) | F 0.019 | F 0.014 | P 0.065 |
| snow share of precipitation (probe's rule 0.343) | 0.432 | 0.488 | 0.369 |
| share of runoff in March to May | 0.92 | 0.96 | 0.79 |
| runoff / rain | 0.44 | 0.40 | 0.58 |
| phase-counterfactual 1039529598, runoff-rain r | F 0.039 | F 0.041 | P 0.093 |
| snow share of precipitation (probe's rule 0.372) | 0.456 | 0.516 | 0.392 |
| share of runoff in March to May | 0.95 | 0.97 | 0.83 |

Every other gate seed of the two probes passes in all three columns.

Moister air gives both failing seeds less snow, closer to the probe's own
rule, and less of their runoff in spring. It also evaporates less, so more of
the rain runs off; these runs cannot separate the two channels. Either way the
two verdicts are set by the humidity mock, not by the soil class or the
vegetation.

What the table says:

- **Model findings.** `flux_identity` fails in every column, on 901 to 1169
  days wherever shortwave reaches SUMMA and 424 under the midnight stamp. The
  hourly phase test fails in every column too. Both are the model.
- **The daily energy residual** passes or fails on the mock, between 3.7 and
  13.4 percent, so its pass in the evaluated configuration is not robust.
- **The Priestley-Taylor reference** flips no verdict at 10 C, 30 C or the
  row's own temperature. It moves resolution-invariance between 20.7 and 22.0
  percent, steady state between 10.6 and 16.1 percent, and pet-consistency's
  wet-soil ratio between 0.85 and 0.93.
- **`pet-consistency`** passes in 7 of 13 columns and fails on either side:
  - with a clear-sky split, dry air or strong wind, runoff barely follows the
    rain (r under 0.05);
  - with moist air, no leaves or the midnight stamp, the forest evaporates too
    little.
- **The step dependence and the missing antecedent response** fall below their
  limits only with humid air or a bare surface, because both remove the
  evaporation that erases the difference.
- **The bare surface** passes steady state (it has no phenology) and the
  precipitation counterfactual (it has no canopy). In exchange it evaporates
  too little, condenses more than it evaporates on some days (negative ET),
  and leaves a third of the drought shift uncancelled.
- **At 90 percent humidity the canopy** holds up to 3.65 mm on the
  precipitation counterfactual and 2.15 mm on the closure probe. Those days
  have heavy precipitation within a few tenths of a degree of freezing, with
  snow on the ground; this was not traced further. In the evaluated
  configuration every overshoot is liquid, on a warm day.
- **How the evaluated values were chosen.** Each was fixed before any table was
  made.
  - The albedo is FAO-56's reference-grass value.
  - The wind is a round 2 m/s, not a height-corrected reference.
  - The humidity is the value this repository's other radiation mock uses.
  - The 20 C reference is FAO-56's standard temperature.
  - The split was changed once, for the reason below.

## What changed in 4.0.0-f787fa5.2

A review of the first version found three documentation faults that each
bore on how a verdict is attributed, and three small ones. Each fault, the
change made and its effect:

- **The Priestley-Taylor inversion was not step-consistent**, as described
  under Forcing.
  - Change: the conversion is evaluated at a fixed 20 C.
  - Effect: no verdict moved. The resolution-invariance deviation went from
    22.0 to 21.0 percent. The steady-state variation went from 12.4 to 10.6
    percent. pet-consistency's wet-soil ratios went from 0.87-0.99 to 0.82-0.94 of demand. The
    dry-down rise went from 7.2 to 8.2 percent on its one failing seed.
- **SUMMA applied 17 m where the attributes said 10 m**, bare soil included.
  - Change: the heights now come from SUMMA's vegetation table, the attributes
    carry the height SUMMA applies, and `run.json` records
    `scalarAdjMeasHeight`.
  - Effect: the bare hourly probe now runs at 10 m. It fails 17 to 22 blocks
    against `rn` where it failed 16 to 22, and still none against SUMMA's own
    net radiation.
- **The rain-snow split is on the wet bulb.**
  - Change: this is now documented with its numbers. The runoff-bounds and
    phase-counterfactual attributions name the humidity mock, and the
    sensitivity section has the per-seed humidity runs above. The shipped
    `tempRangeTimestep` is unchanged.
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

## Running it

```bash
ht verify-adapter --model summa
ht run --model summa --gate-seeds --csv models/result.csv
```

One ten-year daily case (4015 rows) takes about 9 s in the container, most
of it in SUMMA; the hourly surface case (384 rows) about 0.4 s. Every probe
is evaluated on its full record. `run.json` carries, for every case:
- the mapping and the mocks;
- the heights from the table and the height SUMMA applied;
- the water residual and the soil's elastic storage change;
- the snow share of precipitation and the canopy's net energy flux;
- SUMMA's balance diagnostics and the net-radiation gap.

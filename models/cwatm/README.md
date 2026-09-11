# CWatM under HydroTuring

The Community Water Model of IIASA's Water Security group (Burek et al.
2020, *GMD*, [10.5194/gmd-13-3267-2020](https://doi.org/10.5194/gmd-13-3267-2020)),
the ISIMIP global hydrological model, packaged from
[iiasa/CWatM](https://github.com/iiasa/CWatM) at release 1.11 (`5baaadd`).
Proposed by [@kawh1111](https://github.com/kawh1111) in
[Flood-Lab/HydroTuring#28](https://github.com/Flood-Lab/HydroTuring/issues/28).

CWatM is a daily grid model driven by one settings file and netCDF maps:
degree-day snow, interception, a three-layer van Genuchten soil under an
Arno infiltration curve with preferential flow and capillary rise, a linear
groundwater reservoir, a triangular runoff-concentration lag inside each
cell, kinematic-wave routing between cells, and optional lakes, reservoirs,
water demand and MODFLOW. It is the first process-based global model in the
benchmark, and the first whose every store is a declared variable the
adapter can read rather than reconstruct.

## Verdict

**FAIL (INCOMPLETE)**, 12 of 18 probes passed, on the gate seeds and the full
record of every probe (`ht run --model cwatm --gate-seeds`). INCOMPLETE is
the worst reason: CWatM reports no energy fluxes. Three probes it can be
asked fail as VIOLATION, and they are not alike.

- `resolution-invariance` is the model: CWatM has no step other than a day,
  and no choice made here moves it.
- `pet-consistency` and `response-nonnegativity` are packaging choices as
  much as model results. Both were run with `preferentialFlow = True`, which
  this package takes from the settings template of
  [iiasa/CWatM-Earth-30min](https://github.com/iiasa/CWatM-Earth-30min) at
  `e9dfd99`, the repository its parameter maps come from. The pinned model
  repository ships the switch off in its own 30′ templates
  (`Tutorials/General/06_Watercycle/settings_CWatM_template_30min.ini`,
  `Tutorials/General/09_Calibration_renovation/settings_templates_CWatM/settings_CWatM_template_30min.ini`,
  `Toolkit/Calibration/templates_CWatM/settings_CWatM_template.ini`) and in
  its global 30′ test setup
  (`pytest/settings/30min/global_30min/settings_global_30min.ini`); its
  regional 1 km and 1′ setups mostly switch it on. That one switch alone
  turns both probes into PASS on the same seeds: evaporation on wet soil
  rises to 0.708–0.712 of demand, and the only dip left, 0.048 mm/day on
  seed 1713476937 (4.0e-4 of the added storm), is within the probe's
  tolerance of 1e-3 of the storm. The crop coefficients and the land-cover split
  move `pet-consistency` as well, and the land cover and the
  runoff-concentration lag move `response-nonnegativity` (sensitivity
  table). The switch stays as it was run: picking the setting that passes
  after seeing the result would be tuning to the probe, and which of the two
  upstream settings to follow is the maintainer's call.

| Probe | Result | Mechanism |
| --- | --- | --- |
| `energy/evaporative-partition`, `energy/latent-heat-et-consistency`, `energy/surface-energy-closure` | INCOMPLETE | CWatM reports no latent, sensible or ground heat flux (and the last probe is hourly) |
| `energy/pet-consistency` | VIOLATION: evaporation on the wettest fifth of soil days is 0.695 of demand on the worst seed (0.695–0.705; at least 0.7) | a land cover's transpiration, bare-soil and interception evaporation together cannot exceed `crop_correct × cropKC × ETRef`, and the fraction-weighted crop coefficient is 0.689; snow evaporation is added on top of that cap, which is why the ratio sits just above 0.689 (Sensitivity). **A packaging choice**: with `preferentialFlow = False` it passes at 0.708–0.712, and higher crop coefficients or more forest pass too |
| `mass/resolution-invariance` | VIOLATION: runoff differs by 52.1 % of the rain between PT1H and PT1D (38.6–52.1 %), evaporation by 0.6 % | CWatM has no dt (below); its groundwater reservoir releases `recessionCoeff × storage` per step, so at PT1H it drains 24 times too fast |
| `mass/response-nonnegativity` | VIOLATION: runoff 0.29 mm/day below the control on 2002-02-12, eight days after 120 mm was added (seed 1713476937; the other two never dip) | on wetter soil preferential flow takes a larger share of a later storm, and its interflow part replaces surface runoff; runoff concentration releases interflow through a slower kernel than surface runoff, so the next day carries less (below). **A packaging choice**: with `preferentialFlow = False` it passes (largest dip 0.048 mm/day, within tolerance), and the verdict also moves with land cover and with the slope that sets the lag |
| `mass/catchment-closure` | PASS | residual 0.029 % of the rain; `mrso` within 320 mm; ET 0.43–0.48 of PET |
| `mass/time-origin-invariance` | PASS | identical to floating point under a 28-year shift; residual 0.050 %. The shift keeps the day of year, so it cannot see CWatM's day-of-year snow terms (Parameters) |
| `mass/area-invariance`, `mass/causality` | PASS | identical to floating point |
| `mass/dry-down` | PASS | 83 mm drains in two rainless years against a 406 mm bound |
| `mass/extreme-rain` | PASS | returns 1.00 of the 997 mm added at the top of the ladder |
| `mass/phase-counterfactual` | PASS | turning snow into rain moves runoff by at most 1.8 % of the rain |
| `mass/steady-state` | PASS | settles to within 0.8 % |
| `mass/runoff-bounds` | PASS | runoff 0.61 of the rain, inside [0.11, 1.07] |
| `mass/antecedent-monotonicity` | PASS | the wetter catchment runs off 0.15 of the storm more |
| `mass/warming-response` | PASS | runoff −0.25 and −0.28, evaporation +0.27 and +0.31 per unit of demand, warmer and cooler |
| `momentum/routing-conservation` | PASS | the runoff-concentration store stays within 0.50 of the 15-day bound |

### Preferential flow, runoff concentration and the dip

`soil.py` computes preferential flow, water that bypasses the soil matrix,
as `availWaterInfiltration × relSat ^ preferentialFlowConstant`: the relative
saturation of the upper two layers raised to the fourth power, taken before
infiltration and surface runoff. Of what leaves the soil column that way or
by percolation, the fraction `percolationImp` (0.166) becomes interflow and
the rest recharges groundwater. The Arno curve then limits infiltration on
what is left, and the remainder is surface runoff.

Stepping the worst seed's two runs and reading CWatM's variables on
2002-02-11, a 46 mm day a week after the added storm, the pulse run (upper
two layers at 173 mm before the rain, against 163) differs from the control
by:

| CWatM variable, pulse − control | mm |
| --- | --- |
| preferential flow | +6.71 |
| infiltration | −5.58 |
| surface runoff (`directRunoff`) | −1.13 |
| interflow | +1.12 |
| runoff generated that day, surface plus interflow | −0.02 |
| groundwater recharge | +5.60 |
| baseflow | +0.37 |

The extra recharge comes out of infiltration, not out of runoff, and raises
baseflow; it cannot lower the flow. What changes is the lag the storm's
runoff enters. `runoff_concentration.py` releases each component through a
triangular kernel whose peak is `0.5 + 0.6 × 50000 / (86400 × √tanslope)`
days (3.65 days at the global-median slope) times a per-component factor,
clamped: grassland surface runoff peaks at 1.82 days and releases 0.44 of
itself on the next day, forest surface runoff at 3.0 days (0.17), interflow
at 3.65 days (0.11). Turning 1.13 mm of surface runoff into 1.12 mm of
interflow changes the next day by −0.38 mm if the surface runoff was
grassland's and −0.06 mm if it was forest's. On 2002-02-12 the pulse run
runs off 4.63 mm against the control's 4.92, a dip of 0.29 mm/day despite
0.37 mm/day more baseflow. More water arrived, less left that day.

The kernels explain the variants in the sensitivity table. With preferential
flow off almost nothing moves between kernels (one seed still dips
0.048 mm/day, within tolerance), and with runoff concentration off both
components leave on the day they are generated; both pass. At the
10th-percentile slope the surface kernels are clamped to 3.0 days and
interflow's to 4.0, the same shift costs the next day only 0.08 mm, the
extra baseflow outweighs it, and the probe passes. At the 90th-percentile
slope grassland surface runoff leaves on the day it is generated while
interflow peaks at 1.70 days, and the dip moves to the 46 mm day itself and
grows to 1.40 mm/day. Land cover works through the same kernels: all forest, whose
surface kernel is close to interflow's, never dips; all grassland dips
1.33 mm/day.

## Sensitivity

Two of the three VIOLATIONs move with choices this package had to make. Each
row is the packaged adapter with one choice changed, scored on the same gate
seeds by the harness (an image layered on the packaged one; not archived).
Ranges are across seeds; a verdict is shown where it differs from PASS, and
a dip is reported only where the perturbed run falls below the control.

| Variant | `pet-consistency`: ET on wet soil / demand | `response-nonnegativity`: worst dip | `warming-response`: runoff per unit demand, warmer / cooler | `catchment-closure`: residual; ET / PET | `phase-counterfactual` |
| --- | --- | --- | --- | --- | --- |
| **as packaged**: forest share 0.51, crop coefficients 0.86 / 0.52, template defaults, one snow zone, preferential flow on | **0.695–0.705, FAIL** | **0.29 mm/day, FAIL** | −0.25 / −0.28 | 0.029 %; 0.43–0.48 | 1.8 % |
| `preferentialFlow = False`, as in the pinned model repository's 30′ templates | 0.708–0.712 | 0.048 mm/day on one seed, within tolerance | −0.30 / −0.34 | 0.019 %; 0.45–0.50 | 3.2 % |
| crop coefficients 1.0 for forest and grassland | 0.99–1.00 | 0.16 mm/day, FAIL | −0.19 / −0.24 | 0.055 %; 0.52–0.58 | 1.6 % |
| the template's example calibration: `SnowMeltCoef` 0.0027, `crop_correct` 1.11, `preferentialFlowConstant` 4.5, `arnoBeta_add` 0.19, `factor_interflow` 2.8, `recessionCoeff_factor` 5.278, `runoffConc_factor` 0.1 | 0.78–0.79 | 0.22 mm/day, FAIL | −0.21 / −0.25 | 0.050 %; 0.46–0.51 | 2.7 % |
| all grassland | 0.53, FAIL | 1.33 mm/day, FAIL | −0.26 / −0.29 | 0.028 %; 0.37–0.42 | 1.9 % |
| all forest | 0.86–0.87 | no dip | −0.23 / −0.28 | 0.032 %; 0.50–0.53 | 1.5 % |
| seven snow elevation zones, as in the templates | 0.685–0.694, FAIL | 0.29 mm/day, FAIL | −0.25 / −0.28 | 0.027 %; 0.43–0.47 | 1.7 % |
| `includeRunoffConcentration = False` | 0.695–0.705, FAIL | no dip | −0.25 / −0.28 | 0.029 %; 0.43–0.48 | 1.9 % |
| `tanslope` at its 10th percentile, 0.0018 | 0.695–0.705, FAIL | no dip | −0.25 / −0.28 | 0.029 %; 0.43–0.48 | 1.8 % |
| `tanslope` at its 90th percentile, 0.084 | 0.695–0.705, FAIL | 1.40 mm/day, on the 46 mm day, FAIL | −0.25 / −0.28 | 0.029 %; 0.43–0.48 | 1.8 % |

(The template calibration's `soildepth_factor` of 1.28 is left at 1 because
the adapter sizes the soil column from `soil_capacity_mm`. With preferential
flow off, `extreme-rain`, `dry-down` and `antecedent-monotonicity` were also
rerun and still pass.)

`preferentialFlow` is the one switch that moves both packaging-dependent
verdicts at once. It is on because the CWatM-Earth-30min template that
supplies the other settings and every parameter map has it on; the pinned
model repository's 30′ templates and its global 30′ test setup have it off.
Choosing between the two by which one passes would be tuning to the probe,
so the evaluation is recorded as run and this table says what the other
setting gives.

`pet-consistency` otherwise follows the crop coefficients. A land cover's
transpiration, bare-soil and interception evaporation are capped at
`crop_correct × cropKC × ETRef`, 0.689 of `pet` with the global medians.
Snow evaporation, up to `snowEvapFactor` (0.4) × `minCropKC` (0.2) ×
`ETRef`, is added to `totalET` for the whole cell on top of that cap
(`landcoverType.py:1017`). On the three gate seeds, on the criterion's own wet-soil days and counting
a day as snowy when `snw` is above zero at its end, evaporation is
0.684–0.686 of demand on snow-free days and 0.718–0.763 on days with snow,
and the largest daily ratio is 0.769, which is 0.689 + 0.08: the excess over
the cap is snow evaporation. An upstream slip sits next to it:
`snow_frost.py:789` subtracts `self.var.snowEvap`, still zero inside the
elevation-zone loop, instead of the zone's own `snowEvap`, so bare-soil
potential is never reduced by the snow evaporation drawn from the same
demand. Because potential transpiration is computed as the land cover's
potential minus bare-soil potential, fixing the slip would move that demand
to transpiration rather than lower the cap, so it is documented here, not
counted as the cause. The probe's floor of 0.7 sits just above the cap;
coefficients of 1, the example calibration's `crop_correct` or more forest
lift the cap and pass, more grassland lowers it. With preferential flow off
the probe also passes, 0.008–0.012 above its floor; that run was not taken
apart.

`response-nonnegativity` is the interflow share of preferential flow meeting
a slower runoff-concentration kernel. The kernel peaks scale with
`1/√tanslope`, and the slope is a global median, so the verdict moves with
it: with the lag off, or at the 10th-percentile slope, the probe passes; at
the 90th-percentile slope the dip is almost five times larger. It moves with
land cover through the same kernels, from 1.33 mm/day for all grassland to
no dip for all forest.

`resolution-invariance` is not among the choices; see the next sections.
Closure, the warming response and the phase counterfactual hold in every
variant.

## Licence

CWatM is GPL-3.0. The image is built from the pinned source for this
evaluation, carries the licence at `/opt/cwatm/LICENSE`, and is not pushed
to a registry.

## What the adapter reports

All values are CWatM's own variables, read after every step. CWatM works in
metres of water over the cell; land-cover variables are its own
fraction-weighted sums.

| Column | CWatM variable |
| --- | --- |
| `pr` | the forcing, echoed |
| `evspsbl` | `totalET`: transpiration, bare-soil, open-water, interception and snow evaporation |
| `mrro` | `runoff`: surface runoff, interflow and baseflow after the runoff-concentration lag, i.e. what leaves the cell |
| `dis` | the same over the catchment area, m3/s |
| `mrso` | `sum_soil` = `sum_w1 + sum_w2 + sum_w3` (+ `sum_topwater`, zero without paddy fields) |
| `snw` | `SnowCover`; the degree-day scheme holds no liquid water in the pack |
| `canopy` | `sum_interceptStor` |
| `gw` | `storGroundwater` |
| `channel` | `gridcell_storage`: runoff generated but still inside the runoff-concentration lag |

With abstraction, inflow and MODFLOW off, CWatM exchanges no water with the
outside, so there is no `gwex`. Every step the adapter checks
P − ET − Q against the change in the reported stores, and that CWatM's own
total water storage `tws` equals their sum; both numbers go to `run.json`.

## One grid cell

A lumped catchment is one active CWatM cell (a 1 × 1 mask on a 0.5° grid,
local drain direction 5) whose `CellArea` is the catchment's area. The
forcing is written as netCDF stacks with one record per forcing row, and the
settings file is generated per case. CWatM is then run through its own entry
points (`parse_configuration`, `checkifDate`, `CWATModel`, `ModelFrame`),
with the frame stepped by the adapter so that the stores can be read after
each step. No CWatM code is changed.

| Switch | Setting | Why |
| --- | --- | --- |
| `calc_evaporation` | False | CWatM reads reference ET (`ETMaps`) and open-water evaporation (`E0Maps`) directly instead of computing Penman–Monteith; both are the probe's `pet` |
| `includeWaterDemand`, `includeIrrigation` | False | no probe prescribes abstraction; paddy and irrigated fractions are zero |
| `includeWaterBodies` | False | no lakes or reservoirs in a lumped case |
| `modflow_coupling` | False | CWatM's own linear groundwater reservoir instead |
| `includeRouting` | False | one cell has no river network to route along; CWatM still initialises the routing module, so the channel geometry and `lakeEvaFactor` are written as placeholders no flux reads |
| `includeRunoffConcentration`, `CapillarRise` | True | the within-cell lag and capillary rise from groundwater, as in every 30′ template of CWatM and of CWatM-Earth-30min |
| `preferentialFlow` | True | as in the CWatM-Earth-30min template at `e9dfd99`, where the parameter maps come from. The pinned model repository's own 30′ templates and its global 30′ test setup ship it False, and that switch alone turns `pet-consistency` and `response-nonnegativity` into PASS (Verdict, Sensitivity) |
| `inflow`, `calc_environflow`, `waterquality`, `includeGlaciers`, `usepySnowClim` | False | not part of a lumped water balance |
| `NumberSnowLayers` | 1 | the template's seven elevation zones need a relative-elevation distribution the probe does not give; a lumped case has one elevation |

Sealed and open-water fractions are zero, so `E0Maps` is read but used by no
flux.

## Parameters

Four attributes of `static.json` have an unambiguous counterpart in CWatM:

- `area_km2` is `CellArea`.
- `snow_threshold_degC` is `TempSnow`, the temperature below which
  precipitation falls as snow.
- `canopy_capacity_mm` is the interception capacity of forest and grassland,
  held constant through the year (and `minInterceptCap` is lowered to it if
  the capacity is below CWatM's 1 mm floor).
- `soil_capacity_mm` sizes the soil column. CWatM's soil store is bounded by
  saturation, not by field capacity, so the column is made to hold exactly
  the catchment's capacity at saturation: `StorDepth1` and `StorDepth2` keep
  the ratio of their global medians and are scaled together until
  Σ θs × depth over the three layers, weighted over forest and grassland,
  equals `soil_capacity_mm`. On the 320 mm probes that is 0.14 m and 0.56 m.
  CWatM's smallest column (a 5 cm top layer and the 5 cm minimum it allows
  the other two) holds 67.9 mm at saturation with these medians, so a
  catchment with `soil_capacity_mm` below about 68 mm stops the adapter with
  an error, which the harness records as ERROR. The merged probes use 120,
  180 and 320 mm, and the code is left as it is.

Everything else has no counterpart in the probe. Five calibration factors
take the neutral values the [CALIBRATION] comments of the 30′ templates
document beside the example calibration (`preferentialFlowConstant` 4,
`arnoBeta_add` 0.1, `factor_interflow` 1, `recessionCoeff_factor` 1,
`runoffConc_factor` 1), and `crop_correct` and `soildepth_factor` are 1, no
scaling. `SnowMeltCoef` has no documented neutral value, and 0.004 m/°C/day
is a choice: the CWatM-Earth-30min template sets 0.0027 and the model
repository's 30′ test setups 0.0034. 0.004 is the value the templates carry
commented out under [SNOW] (`#SnowMeltCoef = 0.004`, beside "Snow melt
coefficient: default: 4.0") and the one the model repository's Bhima 1 km
test setup sets (`pytest/settings/1km/Bhima/settings_Bhima.ini:171`). Snow,
frost, runoff concentration and Arno constants are the template's. Every
quantity CWatM reads from a map is the median over the 67,130 land cells
(cells with a valid drain direction) of CWatM's own 30′ input set,
[iiasa/CWatM-Earth-30min](https://github.com/iiasa/CWatM-Earth-30min) at
`e9dfd99`, written back as a one-cell map:

| Quantity | CWatM map | Median | 10th–90th percentile |
| --- | --- | --- | --- |
| saturated conductivity, layers 1–2 / 3 (cm/day) | `ksat1..3` | 3.06 / 2.92 | 2.00–4.43 / 2.41–3.44 |
| forest saturated conductivity, layers 1–2 | `forest_ksat1..2` | 3.44 | 2.03–4.80 |
| θs, layers 1–2 / 3 | `thetas1..3` | 0.447 / 0.423 | 0.40–0.50 / 0.38–0.47 |
| forest θs, layers 1–2 | `forest_thetas1..2` | 0.488 | 0.45–0.53 |
| θr, layers 1–2 / 3 | `thetar1..3` | 0.097 / 0.099 | 0.08–0.11 |
| van Genuchten α, layers 1–2 / 3 | `alpha1..3` | 0.038 / 0.042 | 0.013–0.050 / 0.030–0.050 |
| van Genuchten λ, layers 1–2 / 3 | `lambda1..3` | 0.158 / 0.154 | 0.13–0.23 / 0.11–0.22 |
| storage depths (m), ratio only | `storageDepth1`, `storageDepth2` | 0.267, 1.068 | 0.12–0.30, 0.48–1.20 |
| impeded-percolation fraction | `percolationImp` | 0.166 | 0–0.95 |
| crop group | `cropgrp` | 3 | 2–3.7 |
| crop coefficient, forest / grassland (annual mean) | `cropCoefficient*_10days` | 0.857 / 0.516 | 0.57–1.09 / 0.20–0.89 |
| root depth, forest / grassland (m) | `maxRootDepth` | 1.0 / 0.5 | constant |
| root fraction in layers 1–2, forest / grassland | `rootFraction1` | 0.640 / 0.679 | 0.28–0.93 / 0.44–0.99 |
| groundwater recession (1/day) | `recessionCoeff` | 0.0069 | 0.0002–0.061 |
| specific yield | `specificYield` | 0.05 | 0.01–0.23 |
| elevation standard deviation (m) | `elvstd` | 57 | 11–298 |
| tangent of slope | `tanslope` | 0.0122 | 0.0018–0.084 |
| relative elevation, 12 quantiles (m) | `dzRel_hydro1k` | −41.9 … 115.2 | |
| forest share of forest + grassland | `fractionLandcover` (area-weighted) | 0.508 | |

The crop coefficient and interception capacity maps are 10-day cycles; the
adapter uses the annual mean, constant through the year, because a global
median of a seasonal cycle mixes the hemispheres and describes no catchment.
`derive_parameters.py` in this directory recomputes the table from the
CWatM-Earth-30min maps.

CWatM's snow scheme also reads the calendar, and this package does not
flatten it. The melt coefficient varies as `SnowMeltCoef + 0.0005 ×
sin((doy − 81) × 0.9856°)` m/°C/day, a northern-hemisphere phase whatever
`latitude_deg` says (`snow_frost.py:680`; `SeasonalSnowMeltSin`, which would
shift it, is not set). Between days 166 and 259 an extra "ice melt" of
`IceMeltCoef` (0.007 m/°C/day) × `Tavg` × a half-sine acts on any snow cover
(`snow_frost.py:681-684` and `:767-772`). Both follow the day of year of
CWatM's own calendar, which advances one day per forcing row (next section).
On the daily probes that is the probe's calendar; the time-origin probe's
28-year shift keeps the day of year, so it cannot detect either term, and a
southern-hemisphere catchment would melt out of phase.

## The timestep

CWatM 1.11 has no sub-daily step. `miscInitial.py` sets `DtSec = 86400.0` as
a literal, and although `DtDay` multiplies the meteorological inputs and the
degree-day melt, the soil module's conductivities and sub-steps, the
groundwater recession, the interception and evaporation, the frost index
and the runoff-concentration peak times are all per day, with nothing that
rescales them. The adapter therefore gives CWatM one step per forcing row at
every timestep, as the depth that row carries, and CWatM's calendar
advances one day per row. At PT1H that is a day of drainage, percolation and
evaporation applied to each hour of forcing, and a calendar running 24 times
too fast: the hourly resolution case, 960 rows from 1 April 2000, takes
CWatM's day-of-year snow terms through 960 days, three of their
day-166-to-259 ice-melt windows included. No parameter is rescaled: the
resolution probe is there to measure a model without a dt, and this is one.

Which of those per-day rates carries the dependence was checked by hand,
outside the package. On the resolution probe's gate seed 1253639930, over
the whole 40-day record with its spinup (284 mm of rain), the packaged
model runs off 165 mm at PT1H and 78 mm at PT1D, a difference of 30.5 % of
the rain, with evaporation within 0.2 %; on the 30 days the probe scores,
the same seed differs by 38.6 %. Putting only the groundwater recession
coefficient in per-hour units,
1 − (1 − k)^(1/24), brings the hourly runoff to 72 mm (2.3 %); putting the
soil conductivities in per-hour units as well gives 69 mm (3.4 %). The step
dependence is CWatM's linear groundwater reservoir releasing 0.69 % of its
storage per step instead of per day: at PT1H it drains 24 times too fast and
the water it should have held for months leaves inside the window. CWatM has
no setting that would do this rescaling, so the adapter does not.

## The water balance

CWatM's budget closes to a few millimetres over ten years, not to floating
point, and the residual is the model's. In `soil.py` the capillary rise from
groundwater is added to the third soil layer before recharge is computed;
when it exceeds the percolation that would have recharged the groundwater,
recharge is set to zero and `capRiseFromGW` is reduced, but the water already
added to the soil is not taken back, so the excess is created. On the
closure probe's first gate seed the adapter's per-step check finds at most
0.008 mm in a step and −2.6 mm (0.03 % of the rain) over the record; with
`CapillarRise = False` the same record closes to 2e-13 mm. It is far inside
the closure probe's 5 % and is reported here because it is there.

## Running it

```bash
ht verify-adapter --model cwatm
ht run --model cwatm --gate-seeds --csv models/result.csv
```

The image fetches only the `cwatm` package and the licence at the pinned
commit, compiles CWatM's routing and runoff-concentration kernel from the
`t5.cpp` it ships (the repository carries prebuilt binaries for x86-64 only,
so this is what lets the image run on arm64 as well), and installs numpy,
scipy, netCDF4, pandas and rasterio. The base image is pinned by digest
(`python:3.11-slim@sha256:9534e5a8…`, Python 3.11.16 on Debian trixie), the
one the archived evaluation was built on: rebuilt with the pin, every
filesystem layer but the adapter copy (whose file changed only in comments)
is identical, and the output is byte-identical. g++ and the
transitive Python dependencies are not pinned. It is 847 MB. In the container
at one CPU a three-year record (1,095 rows) takes about 2.5 s and a ten-year
record
about 9 s, most of it CWatM opening its netCDF stacks once per step, well
inside the 60 s budget of the shortest probe.

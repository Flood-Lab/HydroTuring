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
asked fail as VIOLATION: one on the model's missing timestep, one on its
preferential-flow term and one on its crop coefficients. The verdicts of
the last two move with the land-cover and crop-coefficient choices this
package had to make, and the sensitivity table below says how far.

| Probe | Result | Mechanism |
| --- | --- | --- |
| `energy/evaporative-partition`, `energy/latent-heat-et-consistency`, `energy/surface-energy-closure` | INCOMPLETE | CWatM reports no latent, sensible or ground heat flux (and the last probe is hourly) |
| `energy/pet-consistency` | VIOLATION: evaporation on the wettest fifth of soil days is 0.695 of demand on the worst seed (0.695–0.705; at least 0.7) | a land cover's evaporation cannot exceed `crop_correct × cropKC × ETRef`, and the fraction-weighted crop coefficient is 0.689; the probe's floor sits on the ceiling the crop coefficients set. **A packaging choice**: see the sensitivity table |
| `mass/resolution-invariance` | VIOLATION: runoff differs by 52.1 % of the rain between PT1H and PT1D (38.6–52.1 %), evaporation by 0.6 % | CWatM has no dt (below); its groundwater reservoir releases `recessionCoeff × storage` per step, so at PT1H it drains 24 times too fast |
| `mass/response-nonnegativity` | VIOLATION: runoff 0.29 mm/day below the control on 2002-02-12, eight days after 120 mm was added (seed 1713476937; the other two never dip) | preferential flow (below) sends a storm on wetter soil to groundwater faster than saturation excess adds surface runoff. The mechanism is the model's; **whether this storm exposes it moves with the land cover**: see the sensitivity table |
| `mass/catchment-closure` | PASS | residual 0.029 % of the rain; `mrso` within 320 mm; ET 0.43–0.48 of PET |
| `mass/time-origin-invariance` | PASS | identical to floating point under a 28-year shift; residual 0.050 % |
| `mass/area-invariance`, `mass/causality` | PASS | identical to floating point |
| `mass/dry-down` | PASS | 83 mm drains in two rainless years against a 406 mm bound |
| `mass/extreme-rain` | PASS | returns 1.00 of the 997 mm added at the top of the ladder |
| `mass/phase-counterfactual` | PASS | turning snow into rain moves runoff by at most 1.8 % of the rain |
| `mass/steady-state` | PASS | settles to within 0.8 % |
| `mass/runoff-bounds` | PASS | runoff 0.61 of the rain, inside [0.11, 1.07] |
| `mass/antecedent-monotonicity` | PASS | the wetter catchment runs off 0.15 of the storm more |
| `mass/warming-response` | PASS | runoff −0.25 and −0.28, evaporation +0.27 and +0.31 per unit of demand, warmer and cooler |
| `momentum/routing-conservation` | PASS | the runoff-concentration store stays within 0.50 of the 15-day bound |

### Preferential flow and the dip

`soil.py` computes preferential flow, water that bypasses the soil matrix
straight to groundwater, as `availWaterInfiltration × relSat ^
preferentialFlowConstant`, with the relative saturation of the upper two
layers raised to the fourth power. It is taken before infiltration and
before surface runoff, and the Arno curve then limits infiltration on what is
left. On a wetter soil both terms grow, and the fourth power grows faster.
Stepping the worst seed's two runs and reading CWatM's variables on
2002-02-11, a 46 mm day a week after the added storm: the pulse run's upper
two soil layers hold 173 mm before the rain against 163, its preferential
flow is 27.6 mm against 20.9, its infiltration capacity has fallen by
5.6 mm, and so its surface runoff is 8.1 mm against 9.3. Of the extra
6.7 mm of preferential flow, 5.6 mm recharges the groundwater reservoir,
which releases 0.69 % of its storage a day, and the rest becomes interflow;
after the runoff-concentration lag the next day's runoff is 4.63 mm
against 4.92. More water arrived, less left that day. With
`preferentialFlow = False` the same seeds pass the probe; they also pass
with the cell all forest, while with it all grassland the dip is
1.33 mm/day.

## Sensitivity

Two of the three VIOLATIONs move with choices this package had to make. Each
row is the packaged adapter with one choice changed, scored on the same gate
seeds by the harness (an image layered on the packaged one; not archived).
Ranges are across seeds; a verdict is shown where it differs from PASS.

| Variant | `pet-consistency`: ET on wet soil / demand | `response-nonnegativity`: worst dip | `warming-response`: runoff per unit demand, warmer / cooler | `catchment-closure`: residual; ET / PET | `phase-counterfactual` |
| --- | --- | --- | --- | --- | --- |
| **as packaged**: forest share 0.51, crop coefficients 0.86 / 0.52, template defaults, one snow zone, preferential flow on | **0.695–0.705, FAIL** | **0.29 mm/day, FAIL** | −0.25 / −0.28 | 0.029 %; 0.43–0.48 | 1.8 % |
| crop coefficients 1.0 for forest and grassland | 0.99–1.00 | 0.16 mm/day, FAIL | −0.19 / −0.24 | 0.055 %; 0.52–0.58 | 1.6 % |
| the template's example calibration: `SnowMeltCoef` 0.0027, `crop_correct` 1.11, `preferentialFlowConstant` 4.5, `arnoBeta_add` 0.19, `factor_interflow` 2.8, `recessionCoeff_factor` 5.278, `runoffConc_factor` 0.1 | 0.78–0.79 | 0.22 mm/day, FAIL | −0.21 / −0.25 | 0.050 %; 0.46–0.51 | 2.7 % |
| all grassland | 0.53, FAIL | 1.33 mm/day, FAIL | −0.26 / −0.29 | 0.028 %; 0.37–0.42 | 1.9 % |
| all forest | 0.86–0.87 | no dip | −0.23 / −0.28 | 0.032 %; 0.50–0.53 | 1.5 % |
| seven snow elevation zones, as in the template | 0.685–0.694, FAIL | 0.29 mm/day, FAIL | −0.25 / −0.28 | 0.027 %; 0.43–0.47 | 1.7 % |
| `preferentialFlow = False` | 0.71 | no dip | −0.30 / −0.34 | 0.019 %; 0.45–0.50 | 3.2 % |

(The template calibration's `soildepth_factor` of 1.28 is left at 1 because
the adapter sizes the soil column from `soil_capacity_mm`. With preferential
flow off, `extreme-rain`, `dry-down` and `antecedent-monotonicity` were also
rerun and still pass.)

`pet-consistency` is decided by how much of the probe's demand CWatM's crop
coefficients let a land cover use. With the global medians the ceiling is
0.689 of `pet` and the probe asks for 0.7; any choice that raises the
fraction-weighted coefficient (coefficients of 1, the example calibration's
`crop_correct`, more forest) passes, and more grassland fails by more. The
verdict here is a statement about the global-median crop coefficients as much
as about CWatM. `response-nonnegativity` is CWatM's preferential flow: with
the term off the dip is gone whatever else holds. Whether this storm exposes
it depends on the land cover, from 1.33 mm/day for grassland to none for
forest on these seeds. `resolution-invariance` is not among the choices; see
the next sections. Closure, the warming response and the phase counterfactual
hold in every variant.

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
| `includeRunoffConcentration` | True | the within-cell lag, as in CWatM's shipped template |
| `preferentialFlow`, `CapillarRise` | True | as in CWatM's shipped 30′ template |
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

Everything else has no counterpart in the probe. Calibration factors take
the neutral defaults CWatM's settings template documents beside its example
calibration (`SnowMeltCoef` 0.004, `crop_correct` 1, `soildepth_factor` 1,
`preferentialFlowConstant` 4, `arnoBeta_add` 0.1, `factor_interflow` 1,
`recessionCoeff_factor` 1, `runoffConc_factor` 1). Snow, frost, runoff
concentration and Arno constants are the template's. Every quantity CWatM
reads from a map is the median over the 67,130 land cells (cells with a
valid drain direction) of CWatM's own 30′ input set,
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

## The timestep

CWatM 1.11 has no sub-daily step. `miscInitial.py` sets `DtSec = 86400.0` as
a literal, and although `DtDay` multiplies the meteorological inputs and the
degree-day melt, the soil module's conductivities and sub-steps, the
groundwater recession, the interception and evaporation, the frost index
and the runoff-concentration peak times are all per day, with nothing that
rescales them. The adapter therefore gives CWatM one step per forcing row at
every timestep, as the depth that row carries, and CWatM's calendar
advances one day per row. At PT1H that is a day of drainage, percolation and
evaporation applied to each hour of forcing. No parameter is rescaled: the
resolution probe is there to measure a model without a dt, and this is one.

Which of those per-day rates carries the dependence was checked by hand,
outside the package. On the resolution probe's gate seed 1253639930, over
the whole 40-day record with its spinup (284 mm of rain), the packaged
model runs off 165 mm at PT1H and 78 mm at PT1D, a difference of 30.5 % of
the rain, with evaporation within 0.2 %; on the 30 days the probe scores,
the same seed differs by 38.6 %.
Putting only the groundwater recession coefficient in per-hour units,
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
scipy, netCDF4, pandas and rasterio. It is 847 MB. In the container at one
CPU a three-year record (1,095 rows) takes about 2.5 s and a ten-year record
about 9 s, most of it CWatM opening its netCDF stacks once per step, well
inside the 60 s budget of the shortest probe.

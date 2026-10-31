# HYPE under HydroTuring

HYPE, the Hydrological Predictions for the Environment model that SMHI runs operationally
(Lindström et al. 2010, *Hydrol. Res.* 41, [10.2166/nh.2010.007](https://doi.org/10.2166/nh.2010.007)),
packaged from release 5.36.0 of 2026-06-24. Proposed by Yan Chen as
[Flood-Lab/HydroTuring#113](https://github.com/Flood-Lab/HydroTuring/issues/113). LGPL-3.0; the
image carries both licence texts at `/model/`.

## What runs

Upstream's `hype` executable, built unmodified from the versioned release archive with upstream's
`make hype`. `ht_adapter.py` writes HYPE's input files, runs it once, and reads one basin output
file. It computes no hydrology.

Pins: the source archive `hype_5_36_0_src.tgz` (sha256 `d90bc367…f43a`) and `demo2.zip` (sha256
`117fdd72…9024`) are checked at build time, and the base image `python:3.11-slim` is pinned by
digest (Debian 13). The evaluated image used gfortran 14.2.0 at build time, from the Debian
archive rather than a pin. No network at run time.

The configuration is upstream's own `demo2`: its `par.txt`, its `GeoClass.txt` and the `GeoData.txt`
row of subbasin 3672, the largest of its subbasins with neither a lake class nor a wetland (12.9
km2). The issue names no configuration. demo2 was in the first image, built before any probe had run
(13:53 against first results at 14:04 on 2026-10-30), so it was not picked for its score. Two reasons
keep it: it is upstream's discharge demo, with `Qobs.txt` and a performance criterion on `cout`, where
`demo` is a nitrogen and phosphorus demo with no observed flow; and a lake-free row leaves no store
that HYPE holds and this adapter does not read back. `demo` scores two more FAILs, which is not an
argument either way; "The alternative configuration" gives the reasons for each.

Every difference from those files is in `configure()` and is listed in `run.json`:

- `MAINDOWN` 0, because one subbasin is the whole catchment. `AREA` from `area_km2`.
- `RIVLEN` from `main_channel_length_km` and a new `LOC_RIVLEN` column from
  `centroid_channel_length_km`, when the case supplies them. Otherwise demo2's 5700 m main river, and
  HYPE's own default local river, the square root of the land area (`sw_proc.f90:237`).
- Every soil depth in `GeoClass.txt` (the layer bottoms, the tile depth and the stream depth)
  multiplied by one factor so that the area-weighted pore volume equals `soil_capacity_mm`. demo2
  holds 483.14 mm in this row, so the probes' 320 mm gives a factor of 0.662.
- `Pobs.txt`, `Tobs.txt` and an `Xobs.txt` series of `repo`, the documented replacement for HYPE's
  temperature formula under `petmodel 0` (`atm_proc.f90:2091`). `epot` then equals the case's `pet`
  to the output's precision.
- `info.txt` is written for the case: dates, step, and the basin output below. demo2's lake-stage
  criteria and `readobsid` go with it.

## What is reported

| Column | HYPE output | Note |
| --- | --- | --- |
| `pr` | the forcing | echoed from the CSV as text |
| `snm` | `cprf` + `melt` | rain reaches the soil directly, and HYPE holds no liquid in the pack |
| `evspsbl` | `evap` | mm per step over the step length |
| `mrro` | `cout` over `area_km2` | |
| `dis` | `cout` | m3/s |
| `mrso` | `soim` | every soil layer, including water below wilting point |
| `snw` | `snow` | |
| `canopy` | 0 | named under `zero_stores` |
| `channel` | `clrv` + `cmrv` over `area_km2` | local and main river, including their translation queues and dead volume |

## What is not reported

`zero_stores`: `canopy`. HYPE has no interception store. `icpe` is a flux set by a precipitation
correction (`pcluse`), which demo2 leaves at zero, so nothing is intercepted at all.

`not_reported`: `gw`. demo2 uses no aquifer, and the three soil layers are all under `mrso`.

There is no `gwex`. `abstr` is not read: HYPE has irrigation and water transfers, in demo3's
configuration rather than demo2's, and the adapter does not map the column onto either.

## Parameters

No calibration exists for the synthetic catchment. Everything is demo2's `par.txt`, calibrated by
upstream for its own Swedish demo basin, unchanged at PT1D. Results describe that one documented point
in parameter space. Landuse-indexed values run Lake, Mire, Forest, Grass, Other, Agriculture;
soil-indexed values run Organic, Clay, Coarse, Moraine, Thin, Water.

| Parameter | Value | Meaning |
| --- | --- | --- |
| `ttpi` | 1 | half-width of the rain/snow interval, deg C |
| `lp` | 1 | share of field capacity below which evaporation is limited |
| `ttmp` | 0.2, 0.2, 0.5, 0.2, 0.2, 0.2 | threshold temperature for melt and evaporation, by landuse |
| `cmlt` | 2.7, 2.7, 2.2, 2.7, 2.7, 2.7 | melt factor, mm/C per step, by landuse |
| `cevp`, `cevpam`, `cevpph` | 0.1-0.18, 0.5, 45 | HYPE's own PET formula; unused, since `repo` replaces it |
| `rrcs1`, `rrcs2`, `rrcs3` | 0.3-0.7, 0.001-0.04, 0.0002 | soil runoff recession, upper, lower, slope term |
| `trrcs`, `srrcs` | 0.15, 0.15 | tile and surface runoff recession, per step |
| `mperc1`, `mperc2` | 100, 100 | percolation caps, mm per step |
| `damp`, `rivvel` | 0.5, 0.5 | river damping share and velocity, m/s |
| `tcalt` | 0.6 | temperature lapse with class elevation; `GeoData.txt` gives no deviation |
| `frost`, `sfrost` | 2, 1 | frost depth parameters |
| `macrate`, `mactrinf`, `mactrsm`, `srrate` | 0.3, 15, 0.5, 0.1 | macropore and surface share of infiltration above the threshold |
| `wcwp` | 0.2, 0.4, 0.04, 0.15, 0.15, 0.1 | wilting point, by soil |
| `wcfc` | 0.2, 0.11, 0.19, 0.23, 0.1, 0.1 | field capacity above wilting point, by soil |
| `wcep` | 0.04, 0.01, 0.1, 0.08, 0.03, 0.04 | effective porosity, by soil |
| `epotdist` | 3 | evaporation depth distribution |
| `Qmean`, `rivvel1-3`, `rivwidth1-3`, `limqprod`, `gldepi`, `gratp`, `gratk`, `sdnsnew`, `snowdensdt` | demo2's | lake, river-width and snow-density terms that do not change the water budget of this row |

The class mix of subbasin 3672 is 3% class 4 (forest, coarse), 55% class 5 (forest, moraine), 29% class
6 (forest, thin), 5% class 9 (grass, moraine), 1% class 11 (agriculture, clay) and 7% class 12
(agriculture, coarse). `SLOPE_MEAN` 5.4, `ELEV_MEAN` 211.4, `parreg` 1 with no `cevpcorr`.

From the case: `area_km2`, `soil_capacity_mm`, `main_channel_length_km`,
`centroid_channel_length_km`, and the forcing `pr`, `tas`, `pet`.

## The step

PT1D is the only declared step: it is the issue's, and the only one demo2's parameters describe. HYPE's
rates are per timestep and it does not rescale them. The resolution probe asks for the next finer
step of every model, so at PT1H the adapter rescales: `cmlt`, `cevp`, `mperc1`, `mperc2`, `mactrinf`,
`rrcs3` by dt, and `rrcs1`, `rrcs2`, `trrcs`, `srrcs` as 1 - (1 - k)^dt. Nothing else is resampled.
Hourly-only probes come back N/A. `steplength` takes sub-daily steps; PT1H is the only one exercised.

## Evaluation

Gate seeds, all 37 probes, on the merged tree: **15 PASS, 8 FAIL, 14 N/A** (13 INCOMPLETE, 1
INCOMPATIBLE). No ERROR rows; a rerun on the merged tree reproduced every row of the previous run. Closure
holds on every probe that measures it, with a budget residual of at most 0.0001% of the rain. The evaluation
is on the full record (`window_days: full`), and the whole suite takes under a minute on a node.

| Probe | Verdict | What was measured | Mechanism |
| --- | --- | --- | --- |
| `mass/catchment-closure` | FAIL | `mrso` leaves [0, 320] on 18 steps, up to 363 | supersaturation, below. Verified. Closure, ET and fidelity pass |
| `mass/extreme-event-closure` | FAIL | `mrso` up to 380 on 58 steps | supersaturation. Verified. Event closure passes (5.3e-05 mm) |
| `mass/precipitation-counterfactual` | FAIL | `mrso` up to 345 on 16 steps | supersaturation. Verified |
| `mass/spinup-cycle-invariance` | FAIL | `mrso` up to 324 on 10 steps | supersaturation. Verified. The invariance itself passes |
| `mass/ungauged-basin-closure` | FAIL | `mrso` up to 135 against 124.358 on 76 steps | supersaturation. Verified |
| `mass/human-abstraction` | FAIL | all 380 mm of `abstr` missed, and `mrso` up to 344 | no withdrawal in this configuration, column not read; then supersaturation. Verified |
| `mass/area-invariance` | FAIL | `channel` moves by 57.75 relative | local river length is sqrt(area). Verified |
| `momentum/routing-lag-consistency` | FAIL | xlarge lag 9.00 d against [1.17, 7.17] | all water crosses both reaches in series. Verified; depends on a mapping choice |

### Supersaturation

The state bound is the only thing failing on five of these probes, and it is a property of HYPE.
HYPE does not cap a soil layer at its pore volume. Infiltration is not limited by the room left in the
top layer (`soil_proc.f90`, `calculate_infiltration_flow_diversion`). Water above pore volume in that
layer leaves as saturated overland flow at only `srrcs` per step (`soil_proc.f90:2039`), which demo2
sets at 0.15. And macropore flow that finds the lower layers full is added to the top layer after that
overflow has been taken (`soilmodel0.f90:397`, then `:431`; `soil_proc.f90`, `add_macropore_flow`), so
the end-of-step store can stand above capacity whatever `srrcs` is. `soim` reports what HYPE holds.

The check, on `mass/catchment-closure` with seed 598896396 (the probe's worst is another seed):

| Run | Max `mrso` (mm) | Rows above capacity | ET (mm) | Q (mm) |
| --- | --- | --- | --- | --- |
| As packaged, capacity 320 | 341.70 | 9 | 4766 | 4146 |
| `srrcs` = 1, so overflow leaves in its own step | 328.53 | 6 | 4764 | 4148 |
| demo2's depths, capacity 483.14 | 494.14 | 2 | 4957 | 3957 |

So the overflow rate is part of it and not all of it, and the excess also appears at upstream's own
depths. The adapter's scaling is not what causes it.

### Area

HYPE's local river defaults to the square root of the land area (`sw_proc.f90:237`), and the
area probe supplies no channel length. Ten times the area gives a local stream 3.16 times longer, so
water takes longer to cross it and the depth held in `channel` grows. On seed 66780378 the mean `channel`
is 0.334 mm at 250 km2 and 0.864 mm at 2500 km2, and `cout` per unit area moves too. With `LOC_RIVLEN`
fixed at 15811 m in both runs `channel` agrees to 6e-07 relative. A larger catchment having longer
streams is a sensible default and it is upstream's, so this is left as a model property. The probe's
57.75 is the relative departure on a row where the control's store is near zero.

### Lag

The case gives a main-channel length L and a centroid length Lc nested inside it. HYPE has two reaches in
series, and every drop of local runoff crosses both, so `RIVLEN` = L and `LOC_RIVLEN` = Lc gives a
1.5 L path at `rivvel` 0.5 m/s. On the 10,000 km2 catchment that is 480 km, and the peak comes 9 days
after the storm. demo2's own subbasin is built the same way, with a 5.7 km main river and a 3.6 km
default local one on 12.9 km2, which is 1.63 times its main-river length, so the literal mapping is not
out of proportion with upstream's practice. What is out of proportion is a single 10,000 km2 subbasin: a
real HYPE set-up would divide it.

The verdict does depend on this mapping, which was in the first full evaluation and has not been changed
since, so no alternative here was tried in order to move the verdict. Lag in days, seed 169671131, against the probe's bounds:

| Mapping | small | medium | large | xlarge |
| --- | --- | --- | --- | --- |
| **As packaged, `RIVLEN` = L, `LOC_RIVLEN` = Lc** | 0 | 1 | 5 | 9 |
| `RIVLEN` = Lc, `LOC_RIVLEN` = L - Lc | 0 | 1 | 3 | 6 |
| `RIVLEN` = L, `LOC_RIVLEN` = 1 m | 0 | 1 | 3 | 5 |
| `RIVLEN` = Lc, `LOC_RIVLEN` = 1 m | 0 | 0 | 2 | 3 |
| Lengths ignored, 5700 m and sqrt(A) | 0 | 1 | 1 | 2 |

The maintainer should read this FAIL as a statement about HYPE's serial reaches under one geometry
mapping, and not as one about its routing physics.

### Passes

| Probe | Measured | Note |
| --- | --- | --- |
| `energy/pet-consistency` | 0.91 of demand when wet, 0.30 when dry | real: `evap` is limited by `epot` |
| `mass/antecedent-monotonicity` | +10.9 mm from the 60 mm storm | real |
| `mass/causality` | nothing moves before the storm | a forward model at a daily step |
| `mass/dry-down` | 10 mm over 730 dry days against a 322 mm bound | real |
| `mass/extreme-rain` | 1.00 of the top rung returned | real |
| `mass/multi-decadal-drift` | every store settles | real |
| `mass/phase-counterfactual` | 3.0% of the rain | real: rain and melt take different paths |
| `mass/resolution-invariance` | worst `mrso_end` 2.23% of the rain, PT1H against PT1D | depends on the adapter's rescaling, above |
| `mass/response-nonnegativity` | largest dip -7.6e-07 of the added rain | rounding |
| `mass/runoff-bounds` | 0.52 of the rain | real |
| `mass/snowpack-mass-closure` | residual 0.0000% | HYPE holds no liquid in snow |
| `mass/steady-state` | runoff 0.500 mm/d under 2.5 mm/d of rain | real |
| `mass/time-origin-invariance` | departure 0 | close to trivial: the calendar enters through the PET formula, which `repo` replaces. Inferred |
| `mass/warming-response` | runoff -0.298, ET +0.301 per unit of demand | real |
| `momentum/routing-conservation` | within 0.05 of what a 15-day hydrograph can hold | real |

N/A: the energy probes (`hfls`, `hfss`, `hfg`, `rlus`, `ts` and the rest are not computed), the
groundwater-exchange probes (no aquifer here), `mass/exchange-response` (no `gwex`), and the
stage-based momentum probes.

## The alternative configuration

The other parameter set is upstream's `demo`, subbasin 609 of the same release, through the same adapter:
**13 PASS, 10 FAIL, 14 N/A**, against 15, 8 and 14. It was evaluated once, as the alternative the benchmark
asks for, and its directory is not kept. Its score does not make it better or worse as a choice. Most of
the difference comes from things the demo contains and the packaged row does not:

| Probe | `hype` | `demo`, 609 | Why it moves |
| --- | --- | --- | --- |
| `mass/multi-decadal-drift` | PASS | FAIL, ET 1.100 of PET | `cevpcorr` 0.1 for `parreg` 2 raises PET by 10% (`atm_proc.f90:2128`). Verified: `epot` 22087 mm on a `pet` of 20080, and with `cevpcorr` 0 ET is 20080 mm |
| `energy/pet-consistency` | PASS | FAIL, 1.02 of demand | the same 10%. Inferred from the above |
| `mass/ungauged-basin-closure` | FAIL, bounds | FAIL, bounds and ET 1.087 | the same 10%. Inferred from the above |
| `momentum/routing-conservation` | PASS | FAIL, 0.1 mm held at no runoff | dead river volume, `deadl` and `deadm` 0.003. Verified: with both at 0 the minimum `channel` goes from 0.052 mm to 0 |
| `mass/extreme-event-closure` | FAIL, bounds | FAIL, bounds and 3.88 mm event residual on 669 of 1316 events | internal lake, 11.7% of the subbasin, whose store this adapter does not read. Verified: moving that class to forest takes the worst row residual from 20.9 mm to 7.7e-05 mm. A package of `demo` would have to report the lake, so this row says nothing about HYPE |
| `mass/spinup-cycle-invariance` | FAIL, 324 of 320 | PASS | a margin of 4 mm on demo2; inferred to be demo's different soils and classes. Supersaturation still shows in demo on five other probes |

Seven of the eight FAILs are FAILs under `demo` too, for the same reason: supersaturation on five, the
area result, and the lag result, with `mass/human-abstraction` counted among the five. That agreement is
the useful part of the comparison. Of `demo`'s three new FAILs, two come from its 10% PET correction and
one from its dead river volume. Its lake, which this adapter was never written to read, adds a second
reason on a probe that already failed.

## How much each choice moves the result

| Choice | Effect |
| --- | --- |
| Soil capacity from the case rather than demo2's 483 mm | catchment-closure ET 4766 vs 4957 mm, Q 4146 vs 3957 mm of 8938 mm rain. Exceedance appears at either depth |
| Channel lengths from the case | decides the lag FAIL, table above |
| Local river left at HYPE's default on the area probe | decides the area FAIL |
| `repo` for PET | required by the contract: the case prescribes `pet` |
| Subbasin 3672 of demo2 rather than 609 of demo | table above |
| demo2's other lake-free subbasins | not run. They share `par.txt` and differ in class mix |

## Caveats

- The parameters are a Swedish demo basin's. No PASS or FAIL here is about HYPE's skill.
- A single subbasin carries the whole catchment, at every area a probe asks for. HYPE is normally a
  network of them, and upstream inflow, lake routing and water transfers are not exercised at all. demo2's
  own row 3672 has an upstream neighbour, which `MAINDOWN` 0 and a one-row `GeoData.txt` remove.
- The soil capacity is met through depth, which also moves the tile and stream depths and the
  layer thicknesses that percolation and evaporation see.
- The apt packages in the build are taken from the Debian archive of the pinned base image, not pinned
  one by one; the versions are in the build log. The source and demo are checked by checksum.
- `basinoutput decimals 6` and eight significant figures bound the budget at about 1e-5 mm a row, which is
  why the area run with a fixed local stream agrees to 2e-05 on `mrro` rather than to floating point.

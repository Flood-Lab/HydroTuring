# SHUD under HydroTuring

SHUD, the Simulator for Hydrologic Unstructured Domains, by Lele Shu and contributors, from
[SHUD-System/SHUD](https://github.com/SHUD-System/SHUD) at tag `cpu-accel-v1.1.1`, commit
`c7404df0ebaffd67ea267bce32fff36ba36f6cfe`, the revision named in
[Flood-Lab/HydroTuring#106](https://github.com/Flood-Lab/HydroTuring/issues/106), proposed by
Yuanhang Liu. MIT; the image carries the licence at `/model/LICENSE`.

**Read this first.** At this tag SHUD loses 30 to 60% of the rain in these cases, and that accounts for most
of the verdict. The loss is one line of upstream's `ET()`, it appears on upstream's own example with no
HydroTuring code involved, and **upstream has since fixed it** in commit `0958682c78da` (2026-10-13, upstream
issue #15), which postdates the tag and every tagged release so far. See "The headline" and "The tag or the
fix". The minimal mesh does not cause it, and neither does the adapter.

## What runs

Upstream's `shud` executable, compiled unmodified with upstream's Makefile against CVODE 6.0.0, the release its
`./configure` asks for (tarball checked by sha256), on `python:3.11-slim` pinned by digest (Debian 13). Upstream's
Makefile compiles its optional BoomerAMG wrapper into every binary, so Hypre, MPI and OpenBLAS are installed though
the default SPGMR path used here never calls them. The Apptainer image is 369 MB. `ht_adapter.py` writes the text
inputs `shud` reads and reads the binary outputs it writes.

Upstream's writer averages whatever it prints, states included, over the output interval, so a store's end-of-row
value exists only if printed at the solver's own step: 144 records per day per variable, too much for a run's tmpfs.
The output files are therefore named pipes, each reduced as it is written to one value per row: the last sample for a
store, the mean for a flux.

## What is reported

| Column | Native | Note |
| --- | --- | --- |
| `pr` | the forcing | echoed from the CSV as text |
| `evspsbl` | `eleveta`, mean over the mesh | canopy + soil + ground + pond evaporation and transpiration, mm/day |
| `mrro` | `rivqdown` of the reach, over the mesh's area | m3/day to mm/day |
| `dis` | `mrro` over the stated `area_km2` | scaled, not the mesh's reach |
| `mrso` | `eleyunsat` × specific yield, mean over the mesh | end of row, mm of water |
| `gw` | `eleygw` × specific yield, mean over the mesh | end of row, mm of water |
| `snw` | `eleysnow` | end of row, mm |
| `canopy` | `eleyic` | end of row, mm |
| `channel` | `eleysurf` (pond) + `rivystage` × width × length over the mesh | end of row, mm |
| `stage` | `bed_elevation_m`, or the mesh's own bed, + `rivystage` | diagnostic; the mesh's reach |

Every store SHUD has is reported, so `zero_stores` is empty. `not_reported` in `run.json`: `snm`, because
`qEleNetPrep` is net of canopy capture and the snow module's own output is not printed, and the energy terms, because
SHUD computes no surface energy budget. Without `snm`, `mass/snowpack-mass-closure` is N/A (INCOMPLETE). There is no
`gwex`: the mesh is closed and its only outlet is the reach. `abstr` is not read: SHUD has a prescribed-flux boundary
for an element's groundwater, but no withdrawal process, and making one of it would be the adapter's invention.

## The mesh

Upstream runs a catchment on a mesh `rSHUD` builds from rasters. There is no raster here, so the mesh is a representative
stretch of river with a hillslope on each bank: **four triangles and one reach**, sized from upstream's example catchment,
North Fork Cache Creek (`input/ccw`: 196 km2, 1147 elements, 103 reaches). Forcing is uniform. **The stated area enters only
`dis`**, as in `models/wflow_sbm` and `models/lisflood`.

| Item | Value | Source |
| --- | --- | --- |
| Reach length | 1167 m | median length of a type 1 reach in `ccw.sp.riv` |
| Hillslope length | 776 m | catchment area / (2 × total river length) |
| Element area | 0.45 km2 (mean in the example 0.17 km2) | from the two above |
| Hillslope slope | 0.25 | median element slope |
| River bed slope | 0.076 | median reach slope |
| Outlet elevation | 450 m | the example's |

## Parameters

No calibration exists on a synthetic catchment, and these results describe one documented point in parameter space. Everything
is upstream's example: the median over its 1147 elements for the soil and geology tables, its dominant land cover, and its
`ccw.cfg.calib` applied as written, except where the last column says otherwise.

| Table | Value | Source |
| --- | --- | --- |
| Soil `KsatV`, `ThetaS`, `ThetaR`, `InfD`, `Alpha`, `Beta`, `hAreaF`, `macKsatV` | 0.2615 m/d, 0.4056, 0.01, 0.1 m, 3.277 /m, 1.189, 0.01, 26.15 m/d | element medians, `ccw.para.soil` |
| Geology `KsatH`, `KsatV`, `ThetaS`, `ThetaR`, `vAreaF`, `macKsatH`, `Dmac` | 2.215 m/d, 0.2215 m/d, 0.4117, 0.01, 0.01, 221.5 m/d, 1 m | element medians, `ccw.para.geol` |
| Land cover `ALBEDO`, `VEGFRAC`, `ROUGH`, `RZD`, `SOILDGRD`, `IMPAF` | 0.288, **0.7255**, 0.04, 0.4 m, 0.1, 0 | class 11 of 17 in `ccw.para.lc`, that of 1020 of 1147 elements; `RZD` capped below |
| LAI, melt factor | 1.489, 0.0022 | means of class 11's column in `ccw.tsd.lai` and of `ccw.tsd.mf`, flat over the case so no calendar enters |
| `ccw.cfg.calib`, all keys | as written, including `SOIL_KINF` 0.01 | **except `RIV_WDTH+`**: 0 instead of 50, which belongs to a network draining 196 km2 |
| River type: depth, bank slope, width, sinuosity, Manning, `Cwr`, `KsatH`, `BedThick` | see below, 0, 2 m, 1, 0.04, 0.6, 0.1 m/d, 0.1 m | upstream's river type 1 |
| `ABSTOL`, `RELTOL`, `MAX_SOLVER_STEP`, `LSM_STEP`, `INIT_SOLVER_STEP` | 1e-4, 1e-4, 10 min, 60, 1 | `ccw.cfg.para` |
| `INIT_MODE` | 2 (upstream's own guess: unsaturated 0.3 and saturated 0.4 of the aquifer, river 0.2 of its depth) | the example reads a file for a catchment that does not exist here |

**Capacity.** SHUD's whole subsurface is one aquifer, its storage specific yield × depth, so `depth = soil_capacity_mm / Sy`,
with `Sy` = 0.4017. A stated 320 mm gives 0.80 m where the example has 30 m. The river's depth keeps the example's proportion,
5.5 m in 30 m. **This is a real departure from what upstream runs**, and "What the choices move" measures it.

**Roots.** A capacity below about 160 mm would put the roots below the aquifer, which upstream refuses
(`RootReachLevel` must be positive, `Element.cpp`). So root depth is the smaller of the example's 0.4 m and 90% of the aquifer.
This was an adapter fault: `ht.1` passed the example's value through and two probes came back as ERROR. See "Adapter faults".

## Mocked forcing

SHUD computes potential evaporation by Penman-Monteith from temperature, humidity, wind and net radiation, and has no PET input.
Wind is constant at the example's mean, 1.82 m/s. Humidity and net radiation are solved row by row from that row's `pet` and
`tas` alone, no calendar, so that SHUD's open-water and canopy demands, each linear in available energy and in vapour pressure
deficit, both equal `pet`. Labelled in `run.json`. On the closure case SHUD's printed PET over the forcing's `pet` is 0.99999988,
and printed precipitation over `pr` is 1.0000000.

## The step

**Declared `[PT1D, PT1H]`.** SHUD's output interval is a setting and its solver is adaptive under a 10 minute maximum, so a row can
be any multiple of that; these two are what the package has been checked at. A row's forcing is a step function in SHUD's own table.
Nothing is resampled. Rates go in as mm/day at every step, as SHUD reads them, and every store is read at the end of its row.

## The headline: one line of `ET()`

`src/ModelData/MD_ET.cpp` at the tag, once per solver step, so 144 times a day:

```cpp
double icStg = yEleIS[i];                       // read as a depth over the vegetated area
...
icAcc  = min(prcp - snAcc, max(0., (icMax - icStg) / DT_min));
icEvap = min(max(0., icStg / DT_min), qPotEvap[i]);
icStg += (icAcc - icEvap) * DT_min;
yEleIS[i]      = icStg * vgFrac;                // written back as a depth over the whole element
qEleE_IC[i]    = icEvap * vgFrac;
qEleNetPrep[i] = (1. - snFrac) * prcp + snMelt - icAcc * vgFrac;
```

The fluxes carry `vgFrac × (icAcc − icEvap)`, but the store is multiplied by `vgFrac` as a whole. Each call therefore destroys
`(1 − vgFrac) × yEleIS` with no flux to account for it. The canopy refills from the next step's rain, so the less it rains, the larger
the share that disappears, and **the rate depends on how often `ET()` is called**.

**Verified, four independent ways:**

1. **By the equation.** On the closure case (seed 598896396, 1095 rows, the last 730 summed), `(1 − vgFrac) × Σ yEleIS` over every call
   predicts **539.08 mm**. The canopy and snow module's actual residual, `P − netP − E_ic − ΔSnow − ΔCanopy`, is **539.08 mm**. It matches to the hundredth at every other step tried: 718.70 mm at 5 min, 814.71 mm at 2 min. Total residual
   547.1 mm, 34.3% of the rain. Everything else in the budget, pond, soil, aquifer and river, leaves 8.0 mm.
2. **With no HydroTuring code at all.** Upstream's `input/ccw`, its 1147-element mesh and its parameters and forcing, run for a year: residual
   28.2% of rain, of which the same module holds 170.3 of 171.2 mm. Upstream's mesh, on this case's weather: residual 31.0%.
3. **By removing it.** The tag with that line alone changed to `yEleIS[i] = icStg`: the case above goes to 1.19% and runoff ratio 0.09 to 0.31.
   On upstream's mesh, 31.0% to 2.4%. `VEGFRAC = 1`, which changes no code and makes the factor 1: 1.03%.
4. **Upstream agrees.** Commit `0958682c78da`, "keep interception storage on the element-area basis", says the same thing in its message:
   "Every macro step lost (1 - vgFrac) * yEleIS without booking it: 21.7 % of the land precipitation on ccw over 60 days at the default 10 min
   step, 62.3 % at 1 min, so the macro step acted as a physical parameter."

The step dependence shows here too. On the same case, `MAX_SOLVER_STEP` of 2, 5, 10 and 60 minutes gives residuals of 52.0%, 45.8%, 34.3% and 8.8%,
and the example's own value is 10.

## Is it the mesh, or the adapter?

Asked hard, because 3 PASS of 22 scored is far worse than any comparable physical model.

| Suspect | Test | Result |
| --- | --- | --- |
| A solver on a mesh it was not designed for | upstream's 1147-element mesh, same weather, same 730 days | **Not the cause of the loss**: 31.0% there, 34.3% here. Runoff ratio 0.04 there, 0.09 here |
| A thin aquifer | 30 m, the example's, on the minimal mesh | **Not the cause**: same canopy residual to the hundredth, 539.08 mm. Runoff ratio 0.093 against 0.090 |
| Units | printed P and PET over the forcing; river and element fluxes by upstream's `tau` | ratios 1.0000000 and 0.99999988; the rest of the budget closes to 1% once the canopy is out of it |
| Initial conditions | budgets start after spin-up; canopy, snow and pond start at zero | not involved: the residual is per step and does not decay. Runs at 6, 9 and 10 cycles show the same 41% |
| Spin-up | `mass/spinup-cycle-invariance`, 3652 rows | same; see below for what does not converge |
| Stores missing from the budget | pond, soil, aquifer, snow, canopy and river are all in it; `VEGFRAC = 1` leaves 1.0% | none missing |
| Reading the pipes | the adapter against an ordinary run writing the same outputs to files | identical: residual 547.10 mm both ways |
| Forcing mock | `pet` honoured; not involved in a store's arithmetic | ET over PET depends on it only as intended |

Where the minimal mesh **does** move things, it is the runoff, not the loss: it is a single short reach with little storage, so a flood
arrives and leaves within a day. That is what the rating probes read, and "Limits" says so.

## Adapter faults fixed

- **`ht.1` to `ht.2`: root depth.** `mass/resolution-invariance` (capacity 120 mm, aquifer 0.30 m) and `mass/ungauged-basin-closure`
  (capacities down to 124 mm) were ERROR because the adapter passed the example's 0.4 m roots to an aquifer shallower than that, and upstream
  exits. Fixed as above. Both are now scored, and both fail for the reasons below.

Nothing else in the adapter moved a verdict. It stays at upstream's tolerances and step: see the section on tolerances.

## Evaluation

Gate seeds, all 37 probes, Apptainer: **3 PASS, 19 FAIL, 15 N/A** (11 INCOMPLETE, 4 INCOMPATIBLE) at version `1.1.1-ht.3`, tag
`cpu-accel-v1.1.1`. About 1 s per simulated year.

Four mechanisms, used in the table:

- **Leak**: the canopy line above. Verified.
- **ET double count**: `f_etFlux` takes transpiration as `(qPotTran − qEleE_IC) × va` and ground evaporation as `(qPotEvap − Es) × vb`, where
  `qEleE_IC` and `Es` are already element-area rates, so the demand a wet canopy and a pond have satisfied is subtracted at the wrong weight. Total ET can
  reach `1 + 2v(1 − v)` of PET, 1.40 at this `VEGFRAC`. The line is unchanged on upstream's current master. Verified: at `VEGFRAC = 1` the formula gives 1 and
  `energy/pet-consistency` passes, with `et_plausible` at 1.000.
- **Solver**: `ABSTOL` is 1e-4 m, a tenth of a millimetre, and the right-hand side is discontinuous (`Flux_Infiltration` switches on the sign of
  `Ysurf + netprcp`). `uYsf` is not constrained to be positive, so the pond sits as low as −0.17 mm, and row budgets carry tenths of a millimetre of noise. Upstream's own
  example shows the same (pond −0.05 mm on its last day). Verified by tolerance: with `VEGFRAC = 0` the residual falls from 3.47% at 1e-4 to 2.47% at 1e-5 and 0.11% at 1e-6, and
  the lowest pond from −0.114 mm to −0.005 mm. With the fix and `TOL` 1e-6 together, a diagnostic that overruns the time budget on four probes, `mass/extreme-rain` passes and the pond's minimum is −0.0013 mm, which still fails a limit of −1e-6 mm.
- **Wet soil infiltrates faster**: `u_effkInfi` grows with saturation through `u_satKr` and `macKsatV × u_satn`, while the head term `1 + av / infD` has no suction in it. So a wetter soil
  takes rain in faster, and runs off less. With `hAreaF` 0.01 and `macKsatV` 26 m/d, infiltration capacity is hundreds of mm/day, so storms below saturation produce almost no runoff.
  Mechanism read from the equations and inferred; the numbers below are measured, and they hold at `TOL` 1e-6.

"Fixed" is the same evaluation with upstream's current master, commit `851fb8c`, and "`VEGFRAC` 1" the tag with that one parameter changed.

| Probe | Tag | Reading at the tag | Mechanism at the tag | Fixed | `VEGFRAC` 1 |
| --- | --- | --- | --- | --- | --- |
| `energy/pet-consistency` | **FAIL** | wet soil evaporates 1.09 of demand | ET double count, verified | FAIL 1.12 | PASS |
| `mass/antecedent-monotonicity` | **FAIL** | wetter catchment runs off −0.4 mm of a 60 mm storm; runoff ratio 0.015 | wet soil infiltrates faster, numbers measured (0.38 mm dry, 0 wet; 0.47 and 0 at `TOL` 1e-6 and `VEGFRAC` 1), mechanism inferred; a property of the model at these parameters. The leak adds the ratio | FAIL | FAIL |
| `mass/area-invariance` | **FAIL** | invariance holds, trivially; runoff ratio 0.0062 | leak, verified: only `non_degenerate` fails | PASS | PASS |
| `mass/catchment-closure` | **FAIL** | residual 52.8%; ET 1.25 of PET | leak, verified; ET double count, verified | FAIL, ET 1.27 only | FAIL, ET 1.000 only |
| `mass/causality` | **FAIL** | nothing moves before the storm, but runoff answers +0.000 of a 40 mm storm; ratio 0.017 | leak, verified: with so much rain lost the soil is too dry for the storm to run off; both alternatives answer it | PASS | PASS |
| `mass/extreme-event-closure` | **FAIL** | residual 30.6%; 1122 of 1316 events open | leak, verified; what remains is solver noise: 66 to 130 events open by 0.16 to 0.18 mm | FAIL, events only | FAIL, events only |
| `mass/extreme-rain` | **FAIL** | runoff rose by 1.06 of the rain added; ratio 0.0002 | ET falls when rain is added, measured: 111 mm more rain, soil evaporation −4.2 mm, transpiration −2.2 mm, the same at `TOL` 1e-6 (1.062). Inferred: a raised water table drains the soil faster and leaves less for ET. Absent at `VEGFRAC` 1 and `TOL` 1e-6 (0.995). Ratio: leak, verified | FAIL 1.13 | FAIL 1.09 |
| `mass/human-abstraction` | **FAIL** | all 380 mm unaccounted for; residual 31.3% | no withdrawal process, and the adapter does not read the column, verified; leak, verified | FAIL | FAIL |
| `mass/multi-decadal-drift` | **FAIL** | residual 59.5%; pond below zero on all 18250 rows; runoff ratio 0 | leak, verified: under this drizzle it takes everything ET leaves. Pond: solver, verified | FAIL | FAIL |
| `mass/phase-counterfactual` | **FAIL** | −8.3% of rain | leak, verified: snow bypasses the canopy, rain does not, so changing phase changes the loss | PASS | PASS 2.9% |
| `mass/precipitation-counterfactual` | **FAIL** | residual 32.3%; responses account for 0.78 of added rain | leak, verified | PASS | PASS |
| `mass/resolution-invariance` | **FAIL** | runoff differs 33.0% of rain between hour and day; it was 40.9% before hourly weather stopped arriving up to one 10-minute step late | leak, verified: hourly rain falls in fewer solver steps, so less is lost (seed 746687816, on that earlier image: `mrro` 73.7 mm daily, 160.7 mm hourly, of 212.5 mm). What remains is the split below ground: at `VEGFRAC` 1, `mrso_end` 24.9 and 48.4 mm, `gw_end` 78.3 and 58.8 mm, and runoff within 8.9 mm | FAIL 10.5% | FAIL 14.7% |
| `mass/response-nonnegativity` | **FAIL** | never falls below control; runoff ratio 0.015 | leak, verified. Without the leak the wet-soil mechanism appears: runoff 1 to 4 mm/day below control a month after the pulse, also at `TOL` 1e-6 | FAIL | FAIL |
| `mass/spinup-cycle-invariance` | **FAIL** | residual 41.2%; pond negative; ratio 0.009; runoff differs 294% of the 0.05 mm/day floor | leak and solver, verified. The 294% is 0.147 mm/day on one storm day of the same date, with 0.14 to 0.39 mm differences in the stores. With the fix it is larger, and at `TOL` 1e-6 it remains: the runs do not reach one periodic solution in 6 to 10 cycles. Not traced, a finding as it stands | FAIL | FAIL |
| `mass/steady-state` | **FAIL** | pond varies 6.4% and does not settle; 0.554 of 2.5 mm/day unaccounted for; runoff zero | leak, verified; pond noise: solver, verified: the pond is constant at `TOL` 1e-6. On master the store trend remains, inferred from the ET double count, because `VEGFRAC` 1 settles | FAIL | PASS |
| `mass/time-origin-invariance` | **FAIL** | invariance holds; residual 39.6%; ratio 0.011 | leak, verified | PASS | PASS |
| `mass/ungauged-basin-closure` | **FAIL** | residual 50.9%; pond negative; ratio 0.010 | leak, verified; solver, verified | FAIL | FAIL |
| `momentum/routing-conservation` | **FAIL** | `channel` reaches −0.106 mm | solver, verified: it is the pond. The limit is −1e-6 mm, out of reach at any tolerance tried | FAIL | FAIL |
| `momentum/stage-discharge-monotonic` | **FAIL** | monotone holds; rising limb too short to read a loop | leak, verified: the river is almost always dry. Without it the loop reads inverted by 2 to 9 mm, because `channel` holds the pond, which leads the river; inferred | FAIL | FAIL |
| `mass/dry-down` | PASS | 0 mm over 730 dry days | **trivial**: the reach is already dry | PASS | PASS |
| `mass/runoff-bounds` | PASS | 0.24 of rain, inside [0.10, 1.07] | | PASS | PASS |
| `mass/warming-response` | PASS | runoff −0.336, ET +0.381 per unit of demand | | PASS | PASS |

**Other trivial passes.** `invariance` in `mass/area-invariance` and `mass/time-origin-invariance` passes by construction: area enters only `dis`, and the tables are flat.
That answers the mesh-scale claim in the proposal only in the sense that it was not tested. SHUD's claim concerns mesh and solver flexibility, which a fixed four-element mesh does not exercise.

N/A (INCOMPLETE, 11): the six energy probes that need heat fluxes, `mass/aquifer-recharge-ordering`, `exchange-response`, `groundwater-datum-invariance`, `gw-sw-exchange-consistency` and
`snowpack-mass-closure`. N/A (INCOMPATIBLE, 4): `momentum/froude-regime`, `routing-lag-consistency`, `uniform-flow-friction-consistency`, `wave-celerity-bounds`, whose reach geometry and `q_in` the package does not consume.

## The tag or the fix

Evaluated at upstream's current master, `851fb8c7e146a533a9545f786e8f68fd8a8d817f` (2026-10-15), same adapter, same parameters: **8 PASS, 14 FAIL**. Closure on the closure case goes from 52.8% to 0.33%,
and the leak's failures flip as the table says. The tag's one-line patch used for diagnosis gives the same 8 and 14, which also shows upstream's fix and that line are the same repair. The fix is not in any tag: `cpu-accel-v1.1.2` and `v2.2` predate it.

This package stays on the tag the proposal named. **Whether the archived verdict should describe the tag, or a commit carrying the fix, is the maintainer's decision.** What it reports either way: the leak was a real
defect of the released model, upstream found it, and what remains after it, in both, is the ET double count and the solver noise.

## What the choices move

Closure case, seed 598896396, 1095 rows, the last 730 summed; P = 1594.6 mm. Residual is what the budget leaves, as a share of rain.

| Variant | ET | Q | ΔS | Residual (mm) | Residual | Runoff ratio | Lowest pond (mm) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **as packaged** | 887.5 | 144.3 | 15.8 | 547.1 | 34.31% | 0.090 | −0.112 |
| `VEGFRAC` 0 | 828.5 | 721.5 | −10.8 | 55.4 | 3.47% | 0.452 | −0.114 |
| `VEGFRAC` 0.5 | 902.6 | 179.8 | 8.1 | 504.1 | 31.62% | 0.113 | −0.110 |
| `VEGFRAC` 1 | 1024.3 | 535.8 | 18.0 | 16.5 | 1.03% | 0.336 | −0.137 |
| `MAX_SOLVER_STEP` 2 min | 733.9 | 30.1 | 1.2 | 829.4 | 52.02% | 0.019 | −0.027 |
| `MAX_SOLVER_STEP` 5 min | 782.8 | 69.4 | 12.3 | 730.1 | 45.79% | 0.043 | −0.085 |
| `MAX_SOLVER_STEP` 60 min | 1046.7 | 393.7 | 13.7 | 140.5 | 8.81% | 0.247 | −0.116 |
| `ABSTOL`, `RELTOL` 1e-5 | 891.7 | 142.1 | 16.2 | 544.6 | 34.15% | 0.089 | −0.010 |
| `ABSTOL`, `RELTOL` 1e-6 | 899.3 | 140.0 | 16.3 | 539.0 | 33.80% | 0.088 | −0.002 |
| aquifer 30 m, not 0.80 m | 1048.8 | 148.2 | −150.8 | 548.3 | 34.38% | 0.093 | −0.179 |
| river depth half the aquifer's | 870.0 | 160.8 | 12.0 | 551.8 | 34.60% | 0.101 | −0.110 |
| hillslope slope 0.10, not 0.25 | 887.0 | 156.2 | −8.4 | 559.7 | 35.10% | 0.098 | −0.107 |
| river 52 m wide, not 2 m | 890.7 | 141.9 | 15.6 | 546.4 | 34.27% | 0.089 | −0.098 |
| one-line patch of the tag | 1071.2 | 487.1 | 17.2 | 19.0 | 1.19% | 0.305 | −0.109 |
| upstream's mesh, parameters and 30 m aquifer, tag | 712.5 | 63.8 | 323.4 | 494.9 | 31.0% | 0.040 | |
| upstream's mesh, parameters and 30 m aquifer, patched | 850.8 | 63.4 | 641.9 | 38.4 | 2.4% | 0.040 | |

The canopy residual is 539.08 mm in all of the mesh, aquifer, river and tolerance rows. What sets the runoff at this point is the leak and the solver step.
The geometry and the aquifer hardly move it. A third of the budget going missing makes the rest of the table small by comparison, so read `VEGFRAC` 1 before reading anything else into these runoff ratios.

## An alternative parameter set

`VEGFRAC` = 1 changes one parameter, inside what upstream's table allows, and turns the leak off without any code change: **10 PASS, 12 FAIL**. That flips `energy/pet-consistency`, `mass/area-invariance`, `causality`, `phase-counterfactual`,
`precipitation-counterfactual`, `steady-state` and `time-origin-invariance`. The other failures, solver noise and the wet-soil mechanism, stay. It is reported to show what is left. It is not a calibration, and it is not offered as the model's configuration: choosing it because it passes would be tuning.

## Limits

- One point in parameter space, one mesh, uniform forcing. Nothing here tests SHUD's unstructured-mesh coupling across scales, which is the claim the proposal makes.
- The aquifer carries the stated capacity, so the model runs 0.3 to 1.7 m deep where its example runs 30 m. This changes `mrso` and `gw` and their state bounds, not the loss.
- Upstream's `ABSTOL` of 1e-4 m is a tenth of a millimetre, and a number of criteria here are at that scale. Stating those failures as "solver" does not excuse them: it is the solver setting upstream ships. `TOL` 1e-6 costs 4 to 6× the run time and overran the time budgets of four probes in the diagnostic run.
- `stage` is the mesh's reach and `dis` is scaled to the stated area, so a rating reads the mesh's short reach, not a 250 km2 river. The reach's response is within a day, shorter than a row.
- `abstr` could be run through SHUD's prescribed groundwater boundary; that was not done, because it is not a withdrawal process upstream provides for this purpose.

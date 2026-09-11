# wflow_sbm under HydroTuring

The SBM land-hydrology model of Deltares' [Wflow.jl](https://github.com/Deltares/Wflow.jl),
at the v1.0.4 release (commit `82df720`), requested in
[Flood-Lab/HydroTuring#29](https://github.com/Flood-Lab/HydroTuring/issues/29).
The model paper (van Verseveld et al. 2024, *GMD*,
[10.5194/gmd-17-3199-2024](https://doi.org/10.5194/gmd-17-3199-2024)) evaluates v0.7.3;
this is v1.0.4, as the request asked, and the TOML format and the model code have both
changed since.

wflow_sbm couples the SBM soil column (an unsaturated store over soil layers above a
saturated store with a pseudo water table) to Gash or Rutter interception, HBV snow, and
kinematic waves for lateral subsurface, overland and river flow. It is the first physical,
distributed model submitted, and the first adapter written in Julia rather than Python.

## Licence

MIT. The image carries Wflow's licence at `/opt/wflow/LICENSE`.

## The catchment Wflow is given

Wflow is distributed: a TOML configuration, a netCDF of static maps and a netCDF of
forcing. The adapter gives a lumped catchment to Wflow as one representative cell and runs
it through Wflow's own time loop (`Wflow.Model`, `Wflow.run_timestep!`), reading the
model's variables after every step.

- **One representative cell at the resolution of Wflow's Moselle model.** Hillslope length
  has no lumped definition, so it comes from the model Deltares ships to test Wflow, as
  every other parameter the catchment does not name does. The cell is the square with the
  mean cell area of that schematisation as Wflow computes it: a 0.00833° grid, cell lengths
  through Wflow's `lattometres`, 607 × 927 m on average, 562 446 m² and a side of
  749.96 m. Every flux and store is a depth over the cell, so **the catchment's area
  enters only `dis`**.
- The cell sits in a 2 × 2 grid in metres (`cell_length_in_meter__flag`) with the other three
  cells inactive. Wflow reads the cell size from the coordinate spacing and drops dimensions
  of length one, so one active cell needs a second coordinate on each axis.
- The cell is **land, river and outlet** at once: local drain direction 5 (a pit), river
  mask 1. Wflow gives a pit its diagonal as flow length, 1060.6 m, and the matching flow
  width, so the land and subsurface kinematic waves drain a hillslope of that length.
- **The river** runs the same length and is 6.27 m wide, which makes it cover 1.183 % of the
  cell, the share the Moselle model's rivers cover of its catchment: river length × width
  summed over its 5 809 river cells (333.1 km²), over its 28 158 km² of active cells. Its
  median river width, 30 m, would make the river 5.7 % of a cell this size.
- The outlet has no cell downstream of it, so Wflow passes the overland and lateral
  subsurface flow that leave the outlet cell out of the map, not into the river. **`mrro` is
  everything that leaves the cell**: river discharge at the outlet plus the overland and
  subsurface outflow of the outlet cell. On the closure probe's first seed that is 80.5 %
  overland flow, 17.4 % lateral subsurface flow and 2.0 % river.

**Why not one cell with the catchment's area.** Adapter version `1.0.4-ht.1` did that. A
250 km² catchment became one 15.8 km cell draining a 22 km hillslope, and the tenfold area
of `mass/area-invariance` a 71 km one. Wflow computes lateral flow per unit length, so the
runoff depth followed the area it was told, and the overland wave held water for months at
low flow. Two of that version's four violations (area-invariance, routing-conservation) were
that geometry and pass with this one. It stays in the sensitivity tables below as the "ht.1
geometry" row.

## Parameters

What the catchment names is used; nothing is derived from the forcing.

| Wflow parameter | Value | Source |
| --- | --- | --- |
| soil thickness | `soil_capacity_mm` / (θs − θr) | Wflow defines the soil water capacity as thickness × (θs − θr) (`soil.jl`), so the column holds exactly `soil_capacity_mm` |
| maximum canopy storage `Cmax` | `canopy_capacity_mm` | static.json |
| snowfall threshold `TT`, melt threshold `TTM` | `snow_threshold_degC` | static.json |
| degree-day factor `Cfmax` | `degree_day_factor_mm_per_C_day` | static.json |
| cell side; flow length | 749.96 m; 1060.6 m | Moselle: the square with its mean cell area; the pit cell's diagonal |
| river width; length | 6.27 m; 1060.6 m | Moselle river-area share of the cell, 1.183 %; the pit cell's flow length |
| θs, θr | 0.435808, 0.165628 | Moselle test model, median |
| `KsatVer`, `f` | 256.983 mm/day, 0.00335494 /mm | Moselle, median |
| Brooks–Corey `c` per layer | 9.38966, 9.70355, 10.0999, 10.0956 | Moselle, median per layer |
| soil layers | 100, 300, 800 mm and the rest | Wflow default |
| rooting depth | 387 mm | Moselle, median (Wflow caps it at 0.99 × thickness) |
| `KsatHorFrac` | 100 | Moselle, median |
| land slope, Manning n | 0.0741422, 0.4768 | Moselle, median |
| canopy gap fraction | 0.305933 | Moselle: median over cells of Wflow's exp(−Kext × LAI) at each cell's annual-mean LAI, held constant |
| Gash `EoverR` | 0.11 | Moselle, median |
| `TTI`, `WHC` | 2.0 °C, 0.1 | Moselle, median |
| `PathFrac`, `InfiltCapPath`, `rootdistpar`, `cf_soil` | 0, 5 mm/day, −500, 0.038 | Moselle, median |
| `MaxLeakage` | 0 | Moselle, median |
| river slope, Manning n, bankfull depth | 0.0017048, 0.03, 1 m | Moselle, median over river cells |
| Feddes heads, air-entry pressure, capillary rise | −100 … −16000 cm, −10 cm, 2000 mm | Wflow defaults |

"Moselle" is the schematisation Deltares ships to test Wflow (`staticmaps-moselle.nc`,
wflow-artifacts v1.0.0): medians over its 50 063 active cells, or its 5 809 river cells.
`baseflow_coefficient` and `latitude_deg` have no counterpart and are unused (the probe gives
potential evaporation, so latitude plays no part). Every value used is written to `run.json`.

**Switched off**, because the probe's catchment has none of it: reservoirs and lakes,
glaciers, paddies, irrigation and all water demand, floodplains, lateral snow transport
(which at a pit would carry snow out of the map), the frozen-soil infiltration reduction,
open water outside the river (`WaterFrac` 0), and leakage out of the saturated store
(`MaxLeakage` 0, as in the test model). Leaf area index is not cyclic, so nothing in the
model follows the calendar.

## What the adapter reports

| Column | What it is |
| --- | --- |
| `pr` | the forcing, echoed as given |
| `evspsbl` | Wflow's total actual evapotranspiration (`actevap`): interception, soil evaporation, transpiration, open water |
| `mrro` | river `q_av` at the outlet + overland `q_av` + lateral subsurface `ssf` out of the outlet cell, as a depth rate over the cell |
| `dis` | `mrro` over the catchment's area: `mrro × area_km2 / 86.4`, m3/s |
| `gwex` | minus the leakage out of the saturated store; identically zero with `MaxLeakage` 0 |
| `mrso` | the SBM soil column: unsaturated store over all layers plus the saturated store |
| `gw` | zero: the `sbm` model type has no store below the soil column |
| `snw` | dry snow plus the liquid water held in the pack |
| `canopy` | canopy storage |
| `channel` | water in the land and river kinematic waves, as a depth over the cell |

**Why the saturated store is `mrso` and not `gw`.** SBM's saturated zone is the lower part of
the same soil column, below a pseudo water table `zi`, bounded by the same thickness and
porosity; the unsaturated and saturated stores together can never hold more than
thickness × (θs − θr), which is the capacity `soil_capacity_mm` names and the state bound the
probes apply to `mrso`. The contract's `gw` is groundwater *below* the soil column, and the
`sbm` type has none (Wflow's `sbm_gwf` type adds an aquifer there). Reporting the saturated
store as `gw` would count the soil capacity twice in the dry-down and runoff bounds, which
add the initial `gw` to the capacity.

## The step

Wflow takes the step as `timestepsecs` and rescales its per-day parameters itself, so rows go
in at the case's step as depths and fluxes come back as rates; the manifest lists PT1D and
PT1H. Wflow labels forcing at the end of its interval, so a row stamped with the start of an
interval is written at start + step. Two things change with the step inside Wflow itself:

- **Interception** is Gash at steps of 23 hours or more and the modified Rutter model below
  that (`sbm.jl`). Gash evaporates the day's interception within the day and carries no
  canopy store, so `canopy` is identically zero at the daily step; Rutter carries it.
- The kinematic waves sub-step at a fixed 3600 s (land) and 900 s (river) whatever the model
  step, while the lateral subsurface wave is solved once per model step.

## Packaging

Both stages of the Dockerfile use `julia:1.12.6` pinned by the digest of its multi-architecture
index (`sha256:3688355d…4a2489`), the Julia version Wflow's committed Manifest was resolved
with. Wflow is fetched at the pinned commit and installed with exactly the versions of Wflow's
own `Manifest.toml` at that commit, plus JSON 1.8.0, Parsers 3.0.0 and StructUtils 2.8.5 for the
adapter; the committed `Manifest.toml` here records that set. The adapter is a small package
(`src/WflowSbmAdapter.jl`) with a PrecompileTools workload that runs a daily and an hourly
synthetic case at build time, so the Wflow code a case needs is compiled into the image:
without it a container spent 25 to 30 seconds compiling before its first step. The depot is
read-only at run time (`JULIA_DEPOT_PATH=/tmp/julia-depot:/opt/julia-depot`) and compiled for a
generic CPU of the architecture; Julia writes nothing to `/tmp` under the harness's isolation.
The adapter writes only its scratch schematisation there (static maps, forcing and TOML,
under a megabyte for a ten-year daily case), in a temporary directory it removes after the run.

Built on aarch64 (Docker Desktop, Apple silicon); image 2.05 GB. One ten-year daily case takes
4.9 s in the container (3.4 s inside Julia, peak RSS 548 MB), a forty-day hourly case 2.1 s.

## Running it

```bash
ht verify-adapter --model wflow_sbm
ht run --model wflow_sbm --gate-seeds
```

## Result

**FAIL (INCOMPLETE), 13 of 18 probes passed**, adapter `1.0.4-ht.2`, on the gate seeds, every
case on the full record (`window_days: full`).

| Probe | Verdict | Reason | Detail |
| --- | --- | --- | --- |
| `energy/evaporative-partition` | FAIL | INCOMPLETE | does not report `hfls`, `hfss`, `hfg` |
| `energy/latent-heat-et-consistency` | FAIL | INCOMPLETE | does not report `hfls`, `hfss`, `hfg` |
| `energy/pet-consistency` | PASS | OK | evaporation 0.96 of demand when the soil is wettest, 0.15 when driest |
| `energy/surface-energy-closure` | FAIL | INCOMPLETE | does not report `hfls`, `hfss`, `hfg` |
| `mass/antecedent-monotonicity` | FAIL | VIOLATION | a wet month before the storm adds no runoff (−0.0011 and −0.0010 of the storm on 2 of 3 seeds) |
| `mass/area-invariance` | PASS | OK | identical to floating point at ten times the area |
| `mass/catchment-closure` | PASS | OK | residual 1.8e-4 to 2.1e-4 of the rain; runoff ratio 0.45 to 0.52; ET 0.53 to 0.63 of demand |
| `mass/causality` | PASS | OK | identical to floating point before the storm; 1.00 of it runs off after |
| `mass/dry-down` | FAIL | VIOLATION | on one seed runoff rises 3.6 % between weekly blocks with no rain |
| `mass/extreme-rain` | PASS | OK | returns 1.00 of the rain added at every rung |
| `mass/phase-counterfactual` | PASS | OK | snow as rain moves the volumes by 0.6 to 3.8 % of the rain |
| `mass/resolution-invariance` | PASS | OK | PT1D against PT1H: 0.4 to 1.1 % of the rain |
| `mass/response-nonnegativity` | PASS | OK | largest dip 7.5e-4 of the added rain (limit 1e-3) |
| `mass/runoff-bounds` | PASS | OK | runoff 0.49 to 0.52 of the rain, inside the bounds |
| `mass/steady-state` | PASS | OK | nothing varies; the budget balances to 1e-14 mm/day |
| `mass/time-origin-invariance` | PASS | OK | identical to floating point in 1972 and 2000 |
| `mass/warming-response` | PASS | OK | runoff falls by 0.27 to 0.31 per unit of added demand |
| `momentum/routing-conservation` | PASS | OK | the channel holds at most 0.17 of what a 15-day hydrograph of recent runoff allows |

The three energy-flux probes are INCOMPLETE because wflow_sbm computes no latent, sensible
or ground heat flux; that is the model declining to be asked, not a failure.

**What the geometry change moved**, from `1.0.4-ht.1` (one cell with the catchment's area) to
`1.0.4-ht.2` (this representative cell), on the same gate seeds:

| Probe | ht.1 | ht.2 |
| --- | --- | --- |
| `mass/area-invariance` | VIOLATION: `mrro` departs by 9.2 to 13.9 relative | PASS: identical |
| `momentum/routing-conservation` | VIOLATION: 3.5 to 4.9 times the bound | PASS: 0.09 to 0.17 of it |
| `mass/dry-down` | VIOLATION on seed 1200831778: runoff too steady (cv 0.090) | VIOLATION on the same seed: runoff rises 3.6 % |
| `mass/antecedent-monotonicity` | VIOLATION on two seeds | VIOLATION on the same two seeds |
| `mass/catchment-closure` | PASS, residual 3e-6 of the rain | PASS, residual 2e-4 of the rain |

Everything else passes in both.

**The budget.** Wflow's own land, overland and subsurface mass balances close to 2e-13 mm,
3e-16 m3/s and 9e-11 m3 per step. The adapter's budget over the reported stores closes to
1.8e-4 to 2.1e-4 of the rain on the closure probe and 2.7e-4 to 3.2e-4 on the time-origin
control, and all of that is the river kinematic wave. Its lateral inflow is rain on the river
surface minus open-water evaporation from it, which is negative on 2 877 of the 4 015 days of
the closure probe's first seed; the wave floors its discharge at 1e-30 m3/s rather than drying
the channel, and the water that creates, 1.71 mm over the record, is exactly the adapter's
residual. On days with non-negative inflow the error is below 2e-17 m3/s. It is larger than in
ht.1 (3e-6) because the river covers 1.18 % of a 750 m cell rather than 0.27 % of a 15.8 km one,
and it stays more than two orders of magnitude inside the 5 % limit.

**Steps.** Between PT1D and PT1H the volumes move by 0.4 to 1.1 % of the rain (runoff 0.7 %,
end-of-record soil storage 1.1 % on the worst seed) although interception changes scheme with
the step: Wflow puts its per-day parameters in the units of the step itself, and its kinematic
waves sub-step at the same fixed lengths at both.

### The two violations

Each was diagnosed with the adapter's own `simulate()`, a per-step trace of the soil column and
the harness's own criteria on every gate seed. Neither depends on the cell geometry. Both move
with the mapping of the stated soil capacity onto Wflow's soil thickness and the Moselle
rooting depth, which the sensitivity tables below vary, so both are statements about SBM under
that mapping rather than about SBM alone.

**`mass/antecedent-monotonicity`.** The generator puts a 60 mm storm on 19 July. In the wet
variant 120 mm falls over the 20 days before 9 July, then both variants have ten dry days
before the storm. The criterion counts its month from 9 July, the start of the dry days, so the
"storm" it reports is the 60 mm plus the rain of the following weeks (71 and 89 mm).

On the two failing seeds, 1295520324 and, in brackets, 1802472438:
- On 18 July the wet column holds 220.5 mm against the dry column's 218.2 mm (226.8 against
  220.7): only 2.3 mm (6.1 mm) wetter. It transpired 83 mm (82 mm) of the 120 mm while the rain
  fell and 36 mm (35 mm) more in the ten dry days, against 2 to 3 mm and 0.4 mm in the dry run.
- Both columns sit just above 215 mm, which is (soil thickness − rooting depth) × (θs − θr) =
  (1184 − 387) × 0.270. SBM transpires from the unsaturated store only in rooted layers, and
  from the saturated store only where the water table is above the rooting depth (a sigmoid
  with `rootdistpar` −500). Water below the Moselle roots is out of evaporation's reach, and
  both runs are drained down to it.
- The storm lifts both to 273 to 282 mm, far from the 320 mm at which SBM makes
  saturation-excess runoff. Over the 30 days from the storm they run off 1.70 against 1.63 mm
  (1.44 against 1.41 mm). Over the criterion's month from 9 July the wet run returns 1.57
  against 1.64 mm (1.46 against 1.55 mm).

The third seed passes because its wet column arrives 57 mm wetter (276.5 against 219.1 mm),
and 17.4 and 25.3 mm more rain fall within the next nine days. That saturates it (318.8 mm) and
returns 40 mm more.

The verdict holds under every mapping tried, but whether a storm crosses the threshold depends
on where the mapping puts that evaporation floor against the stated capacity:
- With the Moselle model's own 2000 mm soil thickness and θs lowered to keep the 320 mm
  capacity, the floor rises to (2000 − 387) × 0.16 = 258 mm. The wet column on the second seed
  then reaches 319.1 mm, returns 0.026 of the storm, and that seed passes.
- With roots through 99 % of the column there is no floor. The wet column arrives 100 mm wetter
  but still far from saturation, and the two seeds still fail (0.0009 and 0.016).

**`mass/dry-down`.** Two of the three seeds drain 38 and 45 mm against a 322 mm bound and
never rise. On seed 1200831778:
- The rain stops with the water table at 467 mm, below the 387 mm roots, and 56 mm in the
  unsaturated store.
- Transpiration draws on the rooted part of that store, which drains from 48 to 20 mm. As it
  empties, transpiration falls from 0.85 to 0 mm/day in February, while potential evaporation
  rises from 0.9 to 2.3 mm/day.
- Gravity drainage from the unsaturated store into the saturated store continues (0.37 falling
  to 0.04 mm/day), so the saturated store refills from 190.0 to 192.4 mm and the water table
  rises 9 mm.
- Lateral subsurface flow follows the water table, not the total storage, and here it is all of
  the runoff. It rises 11 %, and runoff rises 3.6 % between weekly blocks 6 and 7, with 17.9 mm
  drained over the two years.

In ht.1 the same seed failed the other way: its 22 km hillslope turned the same drainage into a
trickle too steady to count as varying (cv 0.090). With the Moselle thickness the water table
starts at the roots (389 mm), and runoff rises 7.6 % by the same route. With roots through 99 %
of the column, transpiration reaches the saturated store and empties it (0.001 to 0.004 mm
drains in two years), and all three seeds pass. So this is a statement about SBM with a water
table below its roots.

### Sensitivity

On the gate seeds, through the adapter's `simulate()` with one setting each (`apply_overrides!`),
scored with the harness's criteria:

- **ht.1 geometry** is one cell with the catchment's own area, its river 30 m wide.
- **roots 0.99 × thickness** sets the rooting depth to Wflow's own cap.
- **Moselle thickness** is 2000 mm with θs lowered so the column still holds the stated
  capacity.

The land roughness and horizontal-conductivity settings tried on ht.1 moved neither remaining
violation there and were not repeated.

`momentum/routing-conservation`, channel over its bound (limit 1), seeds 417694852,
924646966, 1431599080:

| Setting | | | | largest channel store |
| --- | --- | --- | --- | --- |
| as shipped (ht.2) | 0.17 | 0.12 | 0.09 | 20 to 24 mm |
| ht.1 geometry | 3.49 | 4.79 | 4.91 | 95 to 101 mm |
| roots 0.99 × thickness | 0.16 | 8.56 | 0.14 | 20 to 24 mm |
| Moselle thickness | 0.13 | 0.18 | 0.12 | 20 to 24 mm |

The 8.56 is one step, the first scored day. Deep roots have dried the column to 68 mm, runoff
over the preceding 30 days is 0.00001 mm/day, and the channel still holds a 0.002 mm trace.

`mass/area-invariance`, largest relative departure of `mrro` (limit 1e-6), seeds 66780378,
573732492, 1080684606:

| Setting | | | |
| --- | --- | --- | --- |
| as shipped (ht.2) | 0 | 0 | 0 |
| ht.1 geometry | 10.1 | 13.9 | 9.2 |
| roots 0.99 × thickness | 0 | 0 | 0 |
| Moselle thickness | 0 | 0 | 0 |

`mass/antecedent-monotonicity`, extra runoff as a share of the storm (at least 0.02), seeds
1295520324, 1802472438, 161940905:

| Setting | | | |
| --- | --- | --- | --- |
| as shipped (ht.2) | −0.0011 | −0.0010 | 0.337 |
| ht.1 geometry | −0.00005 | 0.0007 | 0.354 |
| roots 0.99 × thickness | 0.0009 | 0.016 | 0.022 |
| Moselle thickness | 0.0017 | 0.026 (pass) | 0.205 |

`mass/dry-down`, seeds 186927550, 693879664, 1200831778 (a weekly rise may not exceed 2 %,
runoff cv must reach 0.1):

| Setting | | | 1200831778 |
| --- | --- | --- | --- |
| as shipped (ht.2) | pass | pass | fail: runoff rises 3.6 % (cv 0.23) |
| ht.1 geometry | pass | pass | fail: cv 0.090 |
| roots 0.99 × thickness | pass | pass | pass (0.002 mm drains) |
| Moselle thickness | pass | pass | fail: runoff rises 7.6 % |

The mapping of the saturated store to `mrso` rather than `gw` moves no result: `mrso` stays
inside its bound on every closure seed, and the dry-down and runoff bounds, which would add an
initial `gw` to the capacity, pass or fail here for reasons that do not involve it.

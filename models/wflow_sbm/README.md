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
forcing. The adapter builds the smallest valid schematisation of a lumped catchment and
runs it through Wflow's own time loop (`Wflow.Model`, `Wflow.run_timestep!`), reading the
model's variables after every step.

- **One active cell**, a square whose area is `area_km2`, in a 2 × 2 grid in metres
  (`cell_length_in_meter__flag`) with the other three cells inactive. Wflow reads the cell
  size from the coordinate spacing and drops dimensions of length one, so one active cell
  needs a second coordinate on each axis.
- The cell is **land, river and outlet** at once: local drain direction 5 (a pit), river
  mask 1. Wflow gives a pit the diagonal as its flow length (√2 × the side) and the matching
  flow width, so the land and subsurface kinematic waves drain a hillslope of that length;
  the river is given the same length.
- The outlet has no cell downstream of it, so Wflow passes the overland and lateral
  subsurface flow that leave the outlet cell out of the map, not into the river. In a
  one-cell catchment that is most of the water, so **`mrro` is everything that leaves:
  river discharge at the outlet plus the overland and subsurface outflow of the outlet
  cell**. The river itself only receives rain falling on its own surface.

## Parameters

What the catchment names is used; nothing is derived from the forcing.

| Wflow parameter | Value | Source |
| --- | --- | --- |
| soil thickness | `soil_capacity_mm` / (θs − θr) | Wflow defines the soil water capacity as thickness × (θs − θr) (`soil.jl`), so the column holds exactly `soil_capacity_mm` |
| maximum canopy storage `Cmax` | `canopy_capacity_mm` | static.json |
| snowfall threshold `TT`, melt threshold `TTM` | `snow_threshold_degC` | static.json |
| degree-day factor `Cfmax` | `degree_day_factor_mm_per_C_day` | static.json |
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
| river slope, Manning n, width, bankfull depth | 0.0017048, 0.03, 30 m, 1 m | Moselle, median over river cells |
| river length | √2 × cell side | the pit cell's own flow length |
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
| `mrro` | river `q_av` at the outlet + overland `q_av` + lateral subsurface `ssf` out of the outlet cell, as a depth rate over the catchment |
| `dis` | the same outflow in m3/s |
| `gwex` | minus the leakage out of the saturated store; identically zero with `MaxLeakage` 0 |
| `mrso` | the SBM soil column: unsaturated store over all layers plus the saturated store |
| `gw` | zero: the `sbm` model type has no store below the soil column |
| `snw` | dry snow plus the liquid water held in the pack |
| `canopy` | canopy storage |
| `channel` | water in the land and river kinematic waves, over the catchment area |

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

`julia:1.12.6` (the Julia version Wflow's committed Manifest was resolved with; multi-arch).
Wflow is fetched at the pinned commit and installed with exactly the versions of Wflow's own
`Manifest.toml` at that commit, plus JSON 1.8.0, Parsers 3.0.0 and StructUtils 2.8.5 for the
adapter; the committed `Manifest.toml` here records that set. The adapter is a small package
(`src/WflowSbmAdapter.jl`) with a PrecompileTools workload that runs a daily and an hourly
synthetic case at build time, so the Wflow code a case needs is compiled into the image:
without it a container spent 25 to 30 seconds compiling before its first step. The depot is
read-only at run time (`JULIA_DEPOT_PATH=/tmp/julia-depot:/opt/julia-depot`), compiled for a
generic CPU of the architecture, and a run under the harness's isolation writes nothing to
`/tmp`.

Built on aarch64 (Docker Desktop, Apple silicon); image 2.05 GB. One ten-year daily case takes
3.9 s in the container (2.6 s inside Julia, peak RSS 550 MB), a forty-day hourly case 2.1 s.

## Running it

```bash
ht verify-adapter --model wflow_sbm
ht run --model wflow_sbm --gate-seeds
```

## Result

**FAIL (INCOMPLETE), 11 of 18 probes passed**, on the gate seeds, every case on the full
record (`window_days: full`).

| Probe | Verdict | Reason | Detail |
| --- | --- | --- | --- |
| `energy/evaporative-partition` | FAIL | INCOMPLETE | does not report `hfls`, `hfss`, `hfg` |
| `energy/latent-heat-et-consistency` | FAIL | INCOMPLETE | does not report `hfls`, `hfss`, `hfg` |
| `energy/pet-consistency` | PASS | OK | evaporation 0.97 of demand when the soil is wettest, 0.16 when driest |
| `energy/surface-energy-closure` | FAIL | INCOMPLETE | does not report `hfls`, `hfss`, `hfg` |
| `mass/antecedent-monotonicity` | FAIL | VIOLATION | the wetter catchment runs off −0.003 and +0.06 mm more from 71 and 89 mm storms (2 of 3 seeds) |
| `mass/area-invariance` | FAIL | VIOLATION | `mrro` departs by 9.2 to 13.9 relative when the area is ten times larger |
| `mass/catchment-closure` | PASS | OK | residual 2e-6 to 3e-6 of the rain; runoff ratio 0.43 to 0.51; ET 0.55 to 0.65 of demand |
| `mass/causality` | PASS | OK | identical to floating point before the storm |
| `mass/dry-down` | FAIL | VIOLATION | on one seed runoff is nearly constant (cv 0.090 < 0.1) |
| `mass/extreme-rain` | PASS | OK | returns 1.00 of the rain added at the top rung |
| `mass/phase-counterfactual` | PASS | OK | snow as rain moves the volumes by 0.9 to 3.8 % of the rain |
| `mass/resolution-invariance` | PASS | OK | PT1D against PT1H: 0.6 to 2.0 % of the rain |
| `mass/response-nonnegativity` | PASS | OK | never dips below the control |
| `mass/runoff-bounds` | PASS | OK | runoff 0.48 to 0.51 of the rain, inside the bounds |
| `mass/steady-state` | PASS | OK | nothing varies; the budget balances to 1e-13 mm/day |
| `mass/time-origin-invariance` | PASS | OK | identical to floating point in 1972 and 2000 |
| `mass/warming-response` | PASS | OK | runoff falls by 0.27 to 0.31 per unit of added demand |
| `momentum/routing-conservation` | FAIL | VIOLATION | the channel holds up to 4.9 times what a 15-day hydrograph of recent runoff allows |

The three energy-flux probes are INCOMPLETE because wflow_sbm computes no latent, sensible
or ground heat flux; that is the model declining to be asked, not a failure.

**The budget.** Wflow's own land, overland and subsurface mass balances close to 1e-13 mm,
1e-12 m3/s and 3e-8 m3 per step, and the adapter's budget over the reported stores to 3e-6 of
the rain. What little remains is the river kinematic wave. Its lateral inflow includes rain on
the river surface minus open-water evaporation from it, which is negative on 882 of 1095 scored days of the routing probe's first seed;
the wave floors its discharge at 1e-30 m3/s rather than drying the channel, and the water that
creates is exactly Wflow's reported river mass-balance error. On days with non-negative inflow
that error is below 3e-15 m3/s.

**Steps.** Between PT1D and PT1H runoff moves by 0.6 to 2.0 % of the rain although interception
changes scheme with the step (Gash daily, Rutter hourly): Wflow puts its per-day parameters in
the units of the step itself, and its kinematic waves sub-step at the same fixed lengths at
both.

### The four violations

Each was diagnosed with the adapter's own `simulate()` and the harness's own criteria on every
gate seed; the sensitivity tables are below. Two are about the domain the adapter had to
choose, two are the model.

**`momentum/routing-conservation`: the one-cell domain.** The channel store is the overland
kinematic wave (98 % of it on average, the river holding the rest), and 98.6 % of the runoff
leaves as overland flow out of the outlet cell (lateral subsurface flow 0.9 %, river 0.5 %).
A kinematic wave stores A = αQ^0.6 per unit length, so its residence time, L·α·Q^−0.4, grows
without bound as the flow falls: in a dry September it held 0.25 mm while releasing
0.004 mm/day. That is Wflow's routing. That it crosses the bound is the domain: one cell makes
the flow path the diagonal of the whole catchment, 22 km for 250 km2. A 1 km cell, the
resolution of Wflow's Moselle model, passes on all three seeds; Wflow's default land roughness
(0.072 against the Moselle median 0.4768) and a tenfold horizontal conductivity do not.

**`mass/area-invariance`: the one-cell domain.** `mrro` departs by 9.2 to 13.9 relative,
`channel` by 5.2 to 7.8, `evspsbl` by 0.2 to 2.1, `mrso` by about 0.01. The area changes
nothing but the cell: its side is √area, so at 2500 km2 the land and subsurface waves drain a
hillslope 3.2 times longer. Running the tenfold case on the control's cell gives output
identical to floating point on all three seeds. Wflow computes lateral flow per unit length, so
a catchment schematised as one cell cannot be area-invariant.

**`mass/antecedent-monotonicity`: the model.** On two seeds the storm falls on 9 July into a
column 99 to 101 mm below its 320 mm capacity (after a rainless month) or 47 to 60 mm below it
(after 120 mm), and neither run saturates within the scored month (peaks 274 to 291 mm). Below
saturation SBM makes runoff only by infiltration excess (never, at 257 mm/day), by exfiltration
and by lateral subsurface flow, and the wetter column spends its extra water on transpiration
and soil evaporation instead: on seed 1295520324, 85 mm of evaporation during the wet month
against 4 mm, 107 mm against 69 mm in the month after the storm, and both runs end that month
with the same storage (224 against 223 mm) and the same 0.22 mm of runoff. The third seed
passes because its wet column starts 1.6 mm below capacity and saturates, returning 41 mm
more. The same two seeds fail and the same one passes with a 1 km cell, Wflow's default
roughness and a tenfold horizontal conductivity, so this is the SBM column under the capacity
the catchment states, not the domain.

**`mass/dry-down`: the model, in two forms.** On seed 1200831778 the rain stops with the column
64 mm below capacity and 0.18 mm on the land, and for two years the only outflow is lateral
flow from the saturated store at an almost constant 0.0029 mm/day, so runoff has a cv of 0.090
against the 0.1 the probe asks for. The drain is otherwise sound: 2.3 mm against a 322 mm
bound, and nothing rises. The other two seeds stop within 8 to 12 mm of capacity with 6 to 8 mm
on the land and recede from their overland wave (cv 3.2 and 3.8). Where lateral flow is made to
carry the runoff (a 1 km cell, a tenfold horizontal conductivity) the same seed fails the other
way: in mid-February evaporation collapses from 0.59 to 0.04 mm/day, the unsaturated store keeps
draining into the saturated store, the water table rises 7 mm, and lateral flow, which follows
the water table rather than the total storage, rises 7 %, so runoff rises 3.3 % between weekly
blocks with no rain.

### Sensitivity

On the gate seeds, through the adapter's `simulate()` with one override each; "1 km cell"
changes the cell side and nothing else, so it no longer has the catchment's area.

`momentum/routing-conservation`, channel over its bound (limit 1), seeds 417694852,
924646966, 1431599080:

| Setting | | | | largest channel store |
| --- | --- | --- | --- | --- |
| as shipped | 3.49 | 4.79 | 4.91 | 95 to 101 mm |
| 1 km cell | 0.13 | 0.13 | 0.11 | 24 to 28 mm |
| land Manning n 0.072 (Wflow default) | 1.66 | 1.91 | 1.74 | 40 to 47 mm |
| KsatHorFrac × 10 | 2.33 | 2.78 | 2.64 | 95 to 101 mm |

`mass/antecedent-monotonicity`, extra runoff as a share of the storm (at least 0.02), seeds
1295520324, 1802472438, 161940905:

| Setting | | | |
| --- | --- | --- | --- |
| as shipped | −0.00005 | 0.0007 | 0.354 |
| 1 km cell | −0.0008 | −0.0007 | 0.318 |
| land Manning n 0.072 | −0.00005 | 0.0007 | 0.371 |
| KsatHorFrac × 10 | −0.0002 | −0.0006 | 0.328 |

`mass/dry-down`, seeds 186927550, 693879664, 1200831778 (a weekly rise may not exceed 2 %,
runoff cv must reach 0.1):

| Setting | | | 1200831778 |
| --- | --- | --- | --- |
| as shipped | pass | pass | fail: cv 0.090 (largest weekly rise −0.05 %) |
| 1 km cell | pass | pass | fail: runoff rises 3.3 % (cv 0.18) |
| land Manning n 0.072 | pass | pass | fail: cv 0.018 (rise 1.5 %) |
| KsatHorFrac × 10 | pass | pass | fail: runoff rises 2.4 % (cv 0.14) |

`mass/area-invariance`, largest relative departure of `mrro` (limit 1e-6), seeds 66780378,
573732492, 1080684606:

| Setting | | | |
| --- | --- | --- | --- |
| as shipped | 10.1 | 13.9 | 9.2 |
| tenfold area on the control's cell | 0 | 0 | 0 |

The mapping of the saturated store to `mrso` rather than `gw` moves no result: `mrso` stays
inside its bound on every closure seed, and the dry-down and runoff bounds, which would add an
initial `gw` to the capacity, pass or fail here for reasons that do not involve it.

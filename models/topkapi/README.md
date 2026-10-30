# TOPKAPI

TOPKAPI (Liu and Todini, 2002, doi:10.5194/hess-6-859-2002) integrates the kinematic wave over
each grid cell into three non-linear reservoirs, soil, overland and channel, each solving
`dV/dt = a - b V^alpha`. This is [sahg/PyTOPKAPI](https://github.com/sahg/PyTOPKAPI) at tag
`v0.4.0`, commit `f7b953eebcb201403f80970846b429bb1081739c`, BSD-3-Clause, by Theo Vischel and
Scott Sinclair. Submission: issue #123.

The image vendors `pytopkapi/` from that commit, unmodified, and calls `pytopkapi.run`. Nothing
is reimplemented.

## Packaging

- **No GDAL.** `parameter_utils/create_file.py`, which turns rasters into a parameter table, is
  the only code needing it, and `import pytopkapi` does not reach it. The adapter writes the two
  parameter files that `pretreatment.read_global_parameters` and `read_cell_parameters` read,
  and the HDF5 forcing `model.run` reads, which is what `create_file.py` would have produced.
- **Python 3.8 and networkx 1.11 are forced.** `pretreatment.sort_cell` takes `np.array()` of
  `nx.topological_sort()`, which is a list only before networkx 2.0, and networkx 1.11 does not
  import on Python 3.9 or later. Python 3.8 is end-of-life; the image's digest is pinned.
- **matplotlib** is installed because `import pytopkapi` imports its plotting tools. It is not used.
- `pytopkapi.run(ini, quiet=True, parallel_exec=False)`: serial, so nothing depends on the
  host's core count.

## The grid

Upstream's smallest, the four cells of its continuity test `tests/test_continuity/4cells`
(Parak, 2006): one hillslope cell draining into three channel cells in a line, 1 km square, 4 km²
in all. It is a representative hillslope and reach, not a map of the stated catchment.

**The stated area enters only `dis`**, `dis = mrro * area_km2 / 86.4`, as `models/wflow_sbm` and
`models/lisflood` do and for their reason: both found that building a cell from the stated area
routes the same depth differently. TOPKAPI would do the same, since `b_s`, `b_o` and `b_c` all
contain `X`. So the depths do not depend on the area, and `mass/area-invariance` passes by that
construction; it is not evidence about TOPKAPI's own geometry. For the same reason `stage` is not
reported: the channel's depth `V_c / (W X_c)` belongs to the 4 km² grid, and would not correspond
to a `dis` scaled to the stated area.

Forcing is uniform over the grid. The adapter reports the catchment as a whole, so nothing here
tests the cell-level conservation the issue mentions.

## Parameters

No calibration exists on a synthetic catchment. Every value is from `4cells/` unless stated.

| Parameter | Value | Source |
| --- | --- | --- |
| `X`, cell side | 1000 m | `global_parameters.dat` |
| `Dt` | the case's step | the request |
| `alpha_s` | 3 | `global_parameters.dat` |
| `alpha_o`, `alpha_c` | 1.66666667 | `global_parameters.dat`; Manning |
| `A_thres` | 1e6 m² | `global_parameters.dat` |
| `W_min`, `W_max` | 1 m, 10 m | `global_parameters.dat` |
| `tan_beta`, by cell | 0.09, 0.08, 0.07, 0.05 | `cell_parameters.dat` |
| `tan_beta_channel` | -, 0.08, 0.07, 0.05 | `cell_parameters.dat` |
| `L`, by cell | 0.5, 0.75, 1, 1.25 m, times `fac_L` | `cell_parameters.dat`; `fac_L` below |
| `Ks` | 0.001 | `cell_parameters.dat`; see "Upstream's `Ks`" |
| `theta_r` | 0.04, 0.05, 0.06, 0.07 | `cell_parameters.dat` |
| `theta_s` | 0.45, 0.5, 0.55, 0.6 | `cell_parameters.dat` |
| `n_o` | 0.1, 0.2, 0.3, 0.4 | `cell_parameters.dat` |
| `n_c` | -, 0.04, 0.05, 0.06 | `cell_parameters.dat` |
| `kc` | 1 | `cell_parameters.dat` |
| `psi_b`, `lambda` | 146.6 mm, 0.322 | `cell_parameters.dat` |
| initial saturation, overland, channel flow | 2.4 %, 0, 0.001 m³/s | `cell_parameters.dat`; forgotten over the spinup |
| `fac_L` | `soil_capacity_mm` / 423.75 | the catchment's capacity; upstream's own multiplier on `L` |
| `fac_Ks`, `fac_n_o`, `fac_n_c` | 1 | `4cells.ini` |
| `solve_s`, `solve_o`, `solve_c` | 1 | upstream's worked example, `model-simulation.ini` |

`fac_L` makes the grid's mean `(theta_s - theta_r) L` equal the catchment's capacity, keeping the
differences between cells: at 320 mm the cells hold 155, 255, 370 and 500 mm.

`ETr`, the crop demand, and `ETo`, the open-water demand on channels, are both `pet`. The probe
gives one demand.

Nothing is rescaled to the step. `Dt` is in seconds in the ODE solvers and in the fluxes, and the
model carries no constant per day.

## Choices, and what each moves

Each was rerun on the thirteen scored probes on the gate seeds, with no container. The first
row is the verdict below.

| Choice | Alternative | What moves |
| --- | --- | --- |
| `4cells` parameters | - | 12 PASS, `momentum/routing-conservation` FAIL |
| `fac_Ks` = 1, `fac_n_c` = 1 | 60 and 1.7, the worked example's calibration for the Liebenbergsvlei | Routing-conservation PASS. `energy/pet-consistency` FAIL (0.36 of demand when wet, limit 0.7), `mass/extreme-rain` FAIL (1.27 of added rain), `mass/response-nonnegativity` FAIL (-0.279 mm/day) |
| `4cells` soils | Liebenbergsvlei soils from `d8/` on these four cells: `Ks` 6.06e-5, `alpha_s` 2.5, gentler slopes | Routing-conservation PASS. `mass/antecedent-monotonicity` FAIL (0.0005 of the storm, limit 0.02), `mass/response-nonnegativity` FAIL (-0.474 mm/day) |
| `fac_L` from the capacity | `fac_L` = 1, 424 mm | Routing-conservation still FAIL (9.1 mm against 6.6 mm); `mass/area-invariance` FAIL on `non_degenerate` (correlation 0.044, limit 0.05) |
| Overland store in `channel` | In `mrso` | Routing-conservation PASS, nothing else moves. See the FAIL below |
| `mrro` is everything crossing the outlet | `/Channel/Qc_out` alone, which is what upstream's plotting tool reads | 11 of 13 FAIL, runoff ratio about 0.005 |
| PT1H native | PT1D native | Resolution pairs PT1D with PT1H rather than PT1H with PT1M: worst 0.76% of rain against 0.04% |

The alternative parameter sets' FAILs were not diagnosed, since they are not in the verdict. They
are recorded because they show that this verdict is a statement about one parameter set: three
sets drawn from upstream's own files give three different patterns of failure.

**The solver.** `solve = 0`, the full Runge-Kutta-Fehlberg of `4cells.ini`, cannot run these
records. Its tolerances are relative to the store, so on a near-empty channel under real inflow
(measured: `V0` = 6e-13 m³, inflow 0.067 m³/s) it cannot take steps of useful length, its
minimum being `10e-10` s: a ten-year daily case had not passed day 744 after 280 s, where `solve = 1` finishes all
4015 rows in 4 s. Todini's quasi-analytical solution falls back on the same solver when it
fails, so the stall is still reachable. It did not happen on the gate seeds; a fresh seed that
hits it would be an ERROR at the probe's time budget.

## Reported, and not reported

| Column | From |
| --- | --- |
| `pr` | echoed |
| `evspsbl` | `/ET_out`, cell mean, plus `/Channel/Ec_out` over the grid |
| `mrro` | `/Channel/Qc_out` plus `/Q_down` at the outlet cell, the two outflow terms of upstream's own continuity test |
| `dis` | `mrro * area_km2 / 86.4` |
| `mrso` | `/Soil/V_s` over the grid |
| `channel` | `/Overland/V_o` plus `/Channel/V_c` over the grid: water generated as runoff and not yet out, as in `models/wflow_sbm` and `models/lisflood` |

Upstream writes the initial state to row 0 of each result and row t's outcome to row t+1, so
result row t+1 is the end of forcing row t.

Not reported, named in `run.json`, and the reason 24 probes are N/A:

- `snw`: no snow store; precipitation is rain whatever the temperature.
- `canopy`: no interception store.
- `gw`: no store below the soil store.
- `stage`: see the grid. `q_in` is not read, though upstream's external flow could carry it.
- `gwex` and every energy term: the model has none. Withdrawals (`abstr`) are not read.

Every budget probe requires `snw` and `canopy`, so `mass/catchment-closure` and its kin are
N/A (INCOMPLETE). Nothing here reports a column of zeros to reach them.

## Results

Evaluated under Apptainer on gate seeds: 12 PASS, 1 FAIL, 24 N/A (22 INCOMPLETE, 2 INCOMPATIBLE).

| Probe | Result | Measured |
| --- | --- | --- |
| `energy/pet-consistency` | PASS | 0.98 of potential when wet, 0.36 when dry, 0.59 over the record |
| `mass/antecedent-monotonicity` | PASS | 11.0 mm more from the 60 mm storm, 0.18 of it |
| `mass/area-invariance` | PASS | identical, by construction: see the grid |
| `mass/causality` | PASS | nothing moves before the storm; `mrro` returns 0.948 of the added rain |
| `mass/dry-down` | PASS | 45 mm over 730 dry days against a 322 mm bound |
| `mass/extreme-rain` | PASS | returns 1.00 of 526 mm added |
| `mass/phase-counterfactual` | PASS | 0.0%; the model never reads temperature |
| `mass/resolution-invariance` | PASS | worst 0.043% of rain, `mrro`, PT1H against PT1M |
| `mass/response-nonnegativity` | PASS | largest dip 9.2e-4 of the added 120 mm, limit 1e-3 |
| `mass/runoff-bounds` | PASS | 0.50 of rain, inside [0.12, 1.05] |
| `mass/steady-state` | PASS | 0.843 mm/day under 2.5 mm/day of rain, nothing varying |
| `mass/warming-response` | PASS | per unit of demand, warmer: `mrro` -0.493, `evspsbl` +0.501 |
| `momentum/routing-conservation` | **FAIL** | 6.2 mm in `channel` against a 4.3 mm bound, one seed of three; the others at 0.93 and 0.96 of it |

### The failure

`momentum/routing-conservation` allows `channel` at most fifteen days of the largest runoff of
the preceding month. On seed 1431599080 a seven-day 48 mm/day storm lands on a catchment dry for
weeks, and on 2003-02-03 `channel` holds 6.2 mm while the outlet has shown at most 0.27 mm/day.

Verified from upstream's own results, day by day:

| Date | Rain | Soil, mm, cells 0 to 3 (capacity 155, 255, 370, 500) | Overland, mm, cells 0 to 3 | Outlet runoff, mm/day |
| --- | --- | --- | --- | --- |
| 02-01 | 64.6 | 154, 177, 200, 244 | 0, 0, 0, 0 | 0.08 |
| 02-02 | 53.9 | 155, 255, 269, 298 | 4.4, 3.6, 0, 0 | 0.08 |
| 02-03 | 48.0 | 155, 255, 370, 365 | 4.3, 14.0, 6.4, 0 | 0.27 |
| 02-04 | 48.0 | 155, 255, 370, 500 | 4.6, 11.9, 25.1, 10.7 | 6.16 |
| 02-05 | 52.2 | 155, 255, 370, 500 | 4.7, 12.9, 21.9, 41.3 | 43.66 |

The soils saturate from the shallowest, upslope, downwards. Overland water forms on a saturated
cell and enters the next cell's soil store, which takes all of it until it too is full. So
water sits in the overland stores for a day before the outlet cell, the deepest, fills and any of
it leaves. This is not a stuck store: the overland stores turn over in under a day (43.6 mm/day
out of 4.3 mm on cell 0, 78.7 out of 14.0 on cell 1), and they hold 0.03 mm three days after the
rain stops. The channel proper holds 0.03 mm that day and at most 0.06 mm in the record.

So the violation is a real feature of the model, surface water ahead of the outlet because
run-on infiltrates downstream, measured on a bound built for a store that the outlet's runoff
feeds. Whether that is what the probe means to catch is for the maintainers; the probe was left
alone. Moving the overland store into `mrso` passes it and is rejected, because that water is
not in the soil and the other distributed models report it under `channel`.

## Findings beyond the verdict

**The budget, which no scored probe tests.** Run by hand on the closure probe's five gate seeds
at PT1D, the budget closes to -1.2e-5 to -7.2e-5 of rain: `Q/P` 0.43 to 0.49, `E/PET` 0.56 to 0.64,
`mrso` 7.6 to 320 mm, `channel` 0 to 29 mm. At PT1H on the resolution probe's case, -8e-9. On the
first seed three of 16060 cell-steps account for nearly all of it, the largest 0.87 mm of water
made on a day the outlet cell's soil reached capacity. Upstream sets the soil's outflow to `Q_max` when the store
overflows, then takes `a_o = max(0, a_s - dV/dt - Qs_out)`; where `Q_max` exceeds what is left,
the clip makes up the difference. The row was traced; the reading of the clip is from the code.

**Evaporation can pass the demand slightly.** Channel cells evaporate `ETo` from the channel on
top of `ks * kc * ETr` from the soil, and the catchment's largest `E/PET` is 1.006. It stays inside
`energy/pet-consistency`'s tolerance on this grid because `W/X` is about 0.01.

**The channel barely carries anything on this grid.** A channel cell passes only `W/X` of its
hillslope flow into the channel and the rest to the next cell, so about 1% of the water reaches
the channel and the outlet's non-channel drainage carries the other 99%. The routing-network case
the issue argues for is not really exercised here; it would need a grid with many more cells, at
a cost in run time.

**Upstream's `Ks`.** The same number is used in two units. `model.py` compares it with a rain
rate in mm/s in the Green-Ampt infiltration, as its docstring says, while
`pretreatment.compute_cell_param` multiplies it by `L` in metres to give a discharge in m³/s, which
needs m/s. So `Ks` = 0.001 is 3.6 mm/h at the surface and 3.6 m/h along the slope. This is read
from the code, not measured. It is left as upstream has it.

**The time step in the issue.** The issue lists PT1H and PT6H. `Dt` is free, and upstream's
fixtures use 1 h and 6 h. The contract has no PT6H.

## Reproducing

```bash
sbatch --chdir="$PWD" --job-name=ht-build-topkapi scripts/hpc/build_model.sbatch topkapi
sbatch --chdir="$PWD" --job-name=ht-topkapi scripts/hpc/run_model.sbatch topkapi
```

A ten-year daily case takes about 4 s, and a minute-step case of 17 days about 20 s, outside the
container. The image is
about 155 MB. The default window applies; most scored probes widen it to their whole record.

# mizuRoute under HydroTuring

mizuRoute v3.1.1 (ESCOMP, tag commit `28514dea3fd28a2f820c7ecaa9014f820c1bb220`,
Apache-2.0), with ParallelIO at `pio2_6_6` (`1d516b82768ef9d2cade524724c38ab63fbe08e4`),
the tag upstream's `.gitmodules` names. Upstream's stand-alone executable, built with its
own Makefile, unmodified. Issue #132.

## Verdict

| Probe | Verdict | Note |
| --- | --- | --- |
| `momentum/wave-celerity-bounds` | PASS | Diffusive wave, 0.9 to 1.6% below the Manning `dQ/dA` at all three states on all three gate seeds, and rising with the state. Mechanism verified, below |
| The other 36 probes | N/A (INCOMPLETE) | Each needs a flux, store or diagnostic this model does not have: `mrro`, `pr`, `evspsbl`, `stage`, a soil or snow store, a heat flux. Most are also daily. No probe was stopped by anything else |

**One pass, on one probe, in one of five schemes, and the scheme decides it.** The table
in "The five schemes" is as much the result as the verdict is. These results describe one
documented point in parameter space.

### Why it passes

`diffusive_wave` in `dfw_route.f90` recomputes the wave speed at every step: `ck` is
`celerity(Qbar, ...)` at a three-point mean of the reach's current discharge, with the width,
slope and roughness in the network file. That is the state dependence the probe asks for,
and `celerity` is `5/3 S^0.3 Q^0.4 B^-0.4 n^-0.6`, Manning with the hydraulic radius taken as
the depth. Evaluated at the three states on the three seeds it is 0.704, 0.929, 1.226 m/s,
then 0.716, 0.945, 1.247, then 0.741, 0.977, 1.289, which is 0.4 to 1.1% above the probe's
exact-radius `dQ/dA`. The observed speeds sit about 2% below `ck` at every state on every seed.
So the state dependence is read from the source and confirmed by the numbers; the 2% offset
is measured and not traced to a term in `solve_ade`. It is not close to the 5% line.

The pass is not trivial, since a fixed speed fails, as `models/rvic` does on this probe.
It is narrow: the probe reads the first moment only, so it says nothing about the numerical
diffusivity described under "The five schemes", and it has no junction or lateral inflow.

The reach-indexed output the issue argues for, `emits.routing`, is not declared: no merged
probe requests `routing.csv` or supplies `routing_network`, so the adapter would have no
way to be tested on it. mizuRoute's native output carries all three variables, and
wiring them up is the natural next step when `mass/routing-network-closure` merges.

## What the adapter does

It writes a control file, a river network and a runoff file, runs `route_runoff`, and
reads the history file. It computes no routing.

| Item | Value | Source |
| --- | --- | --- |
| `route_opt` | 5, diffusive wave | `route/settings/SAMPLE.control`. Chosen from that file before any scheme was run. The documented default is 1, IRF, which reads neither width nor roughness and so could not declare them |
| Network | Two reaches: a headwater draining into the declared reach | Adapter. Under the default `hw_drain_point` of 2 a headwater does no routing: its `REACH_Q` is its local runoff (`dfw_route.f90`), so the inflow arrives whole at the top of the second reach |
| Headwater length | 100 m | Enters only KWT, and is the same in both reach lengths. Not read by the evaluated scheme |
| `length`, `slope`, `width`, `man_n` | `reach_length_m`, `slope`, `width_m`, `manning_n` | Case, read with `hydGeometryOption 0` |
| `sideSlope` | 0 | The declared section is rectangular |
| `doesBasinRoute` | 0 | `q_in` is a river inflow, so no hillslope hydrograph |
| `q_in` | Over the HRU `area` of `area_km2`, in mm/s, with `units_qsim mm/s` | mizuRoute multiplies a runoff depth rate by the HRU area to get a discharge, so the area cancels and has no effect on `dis` |
| `param.nml` | `route/ancillary_data/param.nml.default` unchanged | Upstream. `fshape`, `tscale`, `velo`, `diff` and the namelist `mann_n` and `wscale` do not reach the evaluated scheme, whose geometry comes from the file |
| `floodplain` | F | Upstream's default: unlimited bank depth, which is the declared section |
| Molecules | 20 for KW and DW, 2 for MC | Fixed in `init_model_data.f90` |
| `pio_netcdf_type` | `netcdf` | Upstream documents netCDF-only ParallelIO as sufficient for a small domain |

`dis` is `DWroutedRunoff` of the declared reach. `channel` is `DWvolume` summed over both
reaches, over the same area, in mm. It is upstream's own volume, and cumulative `q_in - dis`
matches it to 1.7e-9 of the inflow on every gate-seed case.

The history file holds a row per step, `outputFrequency 1`, so nothing is aggregated.

### Stores

`channel` is reported and is the only store declared. `snw`, `canopy`, `mrso` and `gw` are
left out and named under `not_reported`.

Version `3.1.1-ht.1` reported `snw` and `canopy` as zero under `zero_stores`. That was
withdrawn. The zero-store rule in `AGENTS.md` is about a model whose control volume includes
the land surface and which has no such compartment: rain falls into it, so a zero `snw` says
all of that water went elsewhere, and the budget tests it. mizuRoute's control volume starts at
its inflow. It never sees weather, and any snow or canopy water is upstream of its boundary,
so its value is unknown rather than zero, which is the distinction the rule itself draws for
a discharge-only network. A zero here would also not be able to cost the model anything,
because no budget covers it, so it would only ever be able to help on a probe about a store
the model does not have. No verdict moved: every probe that reads `snw` or `canopy` also
needs a flux or store this model lacks, so all 36 were N/A (INCOMPLETE) before and after.

### Determinism and the future

The runs are deterministic and read a row at a time. Each reach starts empty, which is
what the 20-day spinup is for.

## The five schemes

The adapter reads `HT_MIZU_ROUTE_OPT`, which the harness never sets, so the comparison is
the same executable on the same network and three gate seeds, scored with the same criterion.
Residual is `(c_obs - c_kin) / c_kin` at low, medium, high.

| `route_opt` | Scheme | Criterion, per seed | Residuals | Observed diffusivity (m2/s) |
| --- | --- | --- | --- | --- |
| 5 | Diffusive wave, **evaluated** | PASS, PASS, PASS | -1.6 to -0.9% | 890 to 3000 |
| 3 | Euler kinematic wave | PASS, PASS, PASS | +1.1 to +1.9% | up to 1300 |
| 2 | Lagrangian kinematic wave | PASS, PASS, PASS | +1.3 to +2.0% | under 2 |
| 4 | Muskingum-Cunge | FAIL, FAIL, PASS | +1.0, +1.3, +5.2% on the failing seeds; +4.75% high on the third | 62 to 121 |
| 1 | Impulse response | FAIL, FAIL, FAIL | No delay at any length | n/a |

So the verdict moves with the scheme. Each row below is a statement about what was run:

- **Muskingum-Cunge** crosses the 5% line by 0.22% at the high state on two seeds and stays
  inside it on the third. Its delay at the high state is below the kinematic schemes' by
  0.12 h. That is measured and not traced to a term in `mc_route.f90`.
- **Impulse response.** The reach's output is the same at 4 km and 20 km. Upstream's own
  `<desireId>` debug print shows why: the reach's hydrograph is `1, 0, 0, ...` at both
  lengths, where `make_uh` with `velo 1.5` would spread it over hours. The only line that sets
  that value is the override for lakes in `process_ntopo.f90`, guarded by
  `NETOPO_in(iSeg)%islake`, and with `is_lake_sim F` nothing assigns `islake`, which has
  no default in `dataTypes.f90`. So this looks like an uninitialised flag in upstream,
  with `MODE=fast` under gfortran 14. The hydrograph is observed; the cause is read from the
  source and not confirmed by a patched build. Either way, this row is not a measurement of
  a fixed-celerity router, and IRF would read neither width nor roughness in any case.

The observed diffusivity is a diagnostic the criterion reports and does not score. The
physical `Q/(2BS)` on the first seed is 38, 76 and 153 m2/s, and diffusive wave reads about
20 times that at this step, 898, 1602 and 2802 m2/s. Most of it is numerical. A verdict
on celerity says nothing about that: the first moment is what the probe measures because
it survives it.

## Step

`dt_qsim` is a configuration key. Serving the same `q_in` at finer steps on the first gate
seed, by hand since no probe does, with the same centroid estimate:

| Step | `c_obs` low, medium, high (m/s) | Residual |
| --- | --- | --- |
| PT1H | 0.690, 0.911, 1.202 | -1.6, -1.3, -1.0% |
| PT15M | 0.691, 0.912, 1.204 | -1.5, -1.2, -0.9% |
| PT5M | stops | `basinUH: cannot identify the maximum number of bins for the tdh` |
| PT1M | stops | Same |

PT5M and PT1M stop because `basinUH` builds the hillslope hydrograph with at most 1000 bins
even under `doesBasinRoute 0`, and the shipped `tscale` of 86400 s needs more at those
steps. A smaller `tscale` would probably clear it without touching routing, but that has not
been tried, so they are not declared. Upstream also accepts `dt_qsim 86400`. PT1D is not
declared because no daily probe supplies `q_in`, so it could not be run or checked here.

## Limits

- This is one reach and one upstream pour point. The network machinery, junctions, lateral
  inflow and domain decomposition, which is what the issue is mainly interested in, are not
  exercised by any merged probe.
- About 1.5 to 5 s per case.

## References

Mizukami et al. (2016), GMD 9(6), 2223-2238, https://doi.org/10.5194/gmd-9-2223-2016.
Cortés-Salazar et al. (2023), HESS 27, 3505-3524, https://doi.org/10.5194/hess-27-3505-2023.

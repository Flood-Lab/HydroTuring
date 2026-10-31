# LISFLOOD-FP under HydroTuring

**Licence: the archive carries no licence file for LISFLOOD-FP itself.** It holds only
`LICENCE-3RD-PARTY.txt`, for its dependencies, and `lisflood.cpp` says "copyright Bristol
University Hydrology Research Group 2008". The Zenodo record lists GPL-2.0. The issue
quotes the project page as GPL-3.0 for non-commercial use, which was not checked here. They
disagree, so the maintainers should settle this before anything is published. The image is
built locally and never pushed.

LISFLOOD-FP v8.1, CPU FV1 solver: the two-dimensional shallow water equations by first-order
finite volumes. The Zenodo archive `LISFLOOD-FP v8.1.zip`, record 6912932, MD5
`77fa78a9c074d7baf9fbb2157fac4090`, unmodified, built with upstream's CMake. This is not the
JRC LISFLOOD rainfall-runoff model, which is `models/lisflood`. Issue #77.

## Verdict

| Probe | Verdict | Mechanism |
| --- | --- | --- |
| `momentum/wave-celerity-bounds` | FAIL (VIOLATION) | At a cell as wide as the reach the bed is a staircase whose steps are comparable to the water depth, so the low state is not a Manning reach. Verified, below |
| The other 36 probes | N/A (INCOMPLETE) | Each needs a flux, store or diagnostic this configuration does not report: `mrro`, `pr`, `evspsbl`, `stage`, a soil or snow store, a heat flux. Most are also daily. No probe was stopped by anything else |

**The verdict is a statement about this grid and not about the solver in the limit.** The
same executable passes at a 25 m cell on all three gate seeds, and it converges on its own
Manning rating as the cell shrinks. The cell size was set before any scoring and has been
kept, because moving it after seeing the score would be tuning to pass. The table under "The
grid" is as much the result as the verdict is. These results describe one documented point
in parameter space.

The proposal's own target, the spatial probe of #75, is not merged, so what the issue
actually asks for (spatial water balance, walls and openings) is not tested here. This
scores the one merged probe that can feed a hydraulic solver.

### The failure

Residual is `(c_obs - c_kin) / c_kin`, on the evaluated grid.

| Gate seed | `c_obs` low, medium, high (m/s) | `c_kin` low, medium, high (m/s) | Residual | Bed step, cell depth at low (m) |
| --- | --- | --- | --- | --- |
| 397273707 | 0.620, 0.879, 1.193 | 0.701, 0.923, 1.214 | -11.5%, -4.8%, -1.8% | 0.105, 0.149 |
| 904225821 | 0.625, 0.893, 1.212 | 0.713, 0.939, 1.234 | -12.3%, -4.9%, -1.8% | 0.112, 0.151 |
| 1411177935 | 1.217, 0.843, 1.217 | 0.738, 0.972, 1.278 | +64.9%, -13.3%, -4.8% | 0.166, 0.116 |

A raster bed is piecewise constant, so the declared slope becomes one step of `S dx` at
every face. FV1 evaluates each face at `max(z_left, z_right)` (`swe/fv1/modifiedvars.cpp`),
so the depth it sees on the downstream side of a face is the cell depth minus the step.

- **Verified: the settled flow is shallower than normal depth, in proportion to the step.**
  Separating the fixed boundary excess using the two reach lengths, the depth in the uniform
  part of the reach is below the wide-channel normal depth by 0.46 to 0.72 of the step on
  the evaluated grid, and by a similar fraction at 25 m (0.42 to 0.67) and 10 m (0.47 to
  0.63). On the first seed at the low state that is 0.149 m against 0.224 m. So the error in
  depth is first order in the cell size.
- **Verified: on the third seed the low state has disconnected.** That seed has the steepest
  slope, and its depth rasters, dumped on the settled plateau with upstream's `saveint` at `ht.2`, show
  cells 0.116 m deep below a 0.166 m step. Every face in the reach is dry on its downstream
  side, so no cell can feel the one below it. At the medium and high states on that seed, and
  at all three on the first seed, no face is dry.
- **Inferred: the disconnected cells behave like a cascade of weirs.** That is why the third
  seed's low state goes fast where the other two seeds go slow, and why it breaks the ordering.
  A weir's `Q ~ h^(3/2)` gives `1.5 q/h` = 1.17 m/s against the 1.217 m/s observed. That
  estimate is the only evidence for the weir reading.
- **Measured, not traced to a term:** where the faces stay connected the wave is slow, and
  increasingly so as the step grows against the depth.

It is the numerics of a stepped bed, not the friction law. The probe's own reference runs a
first-order scheme on about the same 80 m mesh and passes, and it has flat cells with the
slope as an explicit source term (`models/reference_saint_venant/ht_adapter.py`).

So this is a property of the configuration, on top of a property of the solver: FV1's bed is
piecewise constant, and this grid does not resolve flow this shallow over it.

### The grid

The final image, with the cell size set by `HT_LFP_DX`, which the harness never sets, and
scored with the probe's own criterion. The evaluated rows reproduce the harness's report.

| Cell | Seed | Residual, low, medium, high | Criterion | One 20 km run, alone |
| --- | --- | --- | --- | --- |
| Reach width, 80 to 89 m (**evaluated**) | 397273707 | -11.5, -4.8, -1.8% | FAIL | 18 s |
| | 904225821 | -12.3, -4.9, -1.8% | FAIL | 18 s |
| | 1411177935 | +64.9, -13.3, -4.8% | FAIL | 23 s |
| 50 m | 397273707 | -4.2, -2.5, -0.5% | PASS | 32 s |
| | 904225821 | -4.6, -2.8, -0.7% | PASS | 36 s |
| | 1411177935 | -8.6, -5.05, -2.2% | FAIL | 31 s |
| 25 m | 397273707 | -2.5, -0.1, +0.9% | PASS | 104 s |
| | 904225821 | -2.4, -0.6, +0.9% | PASS | 105 s |
| | 1411177935 | -4.0, -1.6, +0.2% | PASS | 101 s |
| 10 m | 397273707 | +0.3, +1.0, +1.5% | PASS | 670 s |

At 10 m, run on the first seed only and at `ht.1` of this adapter, which differs by the
millisecond ramp described below, `c_obs` is 0.703, 0.932, 1.233 m/s. The solver's friction takes the hydraulic radius
as the depth, which gives a rating of `5/3 S^0.3 q^0.4 n^-0.6` = 0.704, 0.929, 1.226 m/s. So
the solver converges on its own Manning rating to within 0.6%, and the positive residual at
10 m is that wide-channel rating against the probe's exact-radius `dQ/dA`.

10 m takes more than twice the probe's `max_runtime_s` of 300 s, so it could not be archived
as it stands. 25 m would fit.

### Why the cell was not moved

Reach width was chosen before anything was scored: it makes each raster cell the section it
stands for, and it is about the probe's own 80 m reference mesh. Nothing in the upstream
manual gives a rule for this. A cell small enough to pass could be justified by the
convergence above, but 25 m would be a value picked after seeing that it passes, and 50 m
does not pass. The maintainers may reasonably prefer a resolved grid, and the rule would
want to come from the flow rather than from the score: the bed step small against the
normal depth at the first row's discharge. Setting `HT_LFP_DX`'s default is a one-line
change followed by a re-run.

## What the adapter does

It writes the raster, boundary files and parameter file upstream reads, runs the
executable, and reads the mass file back. It computes no hydraulics. Standard library only.

| Item | Value | Source |
| --- | --- | --- |
| Solver | `fv1` | Issue #77 |
| Precision | `_NUMERIC_MODE=1`, double | Upstream's documented switch, the one departure from `config.default.cmake`. A 3.8e6 s record in single precision resolves time to about 0.25 s |
| Grid | One row of square cells of about the declared width, fitting `reach_length_m` exactly: 45 to 50 cells on 4 km, 226 to 249 on 20 km | Adapter, see "The grid" |
| Bed | `S (L - x)` at cell centres, falling to zero at the outfall | `slope`, `reach_length_m`. Nothing reads an absolute level, so `bed_elevation_m` is not consumed or declared |
| `fpfric` | `manning_n` | Case |
| Inflow | `QVAR` on the west face, `q_in / width_m`, held across each row | `q_in` is already a discharge, so it is divided by the width and nothing else |
| Outflow | `FREE` on the east face, no slope given | Upstream's uniform-flow boundary. In FV1 it copies the inside state to the outside |
| Banks | North and south faces closed, the default | They carry no friction, so the section is wide-channel whatever the cell size |
| `initial_tstep` | 10 s, the largest solver step | Upstream's default. The step sits at this cap on the settled plateaus. `cfl` is FV1's default of 0.5 |
| `massint` | 60 s | Adapter. Report sampling only, it does not affect the solution |
| `DepthThresh` | 1e-3 m | Upstream's default |

`dis` is the mean of the `Qout` samples in a row, each carried to the declared width by
`width_m / dx`. `channel` is the last `Vol` in the row, carried the same way and spread over
`area_km2`. Both come from the solver's own mass file. Cumulative `q_in - dis` matches the
store to 4.7e-5 of the inflow on the gate seeds. That is good enough to trust the store, and
it is limited by the 60 s sampling. No probe scores it here.

`area_km2` enters `channel` only, never `dis`.

### Stores

`channel` is reported and is the only store declared. `snw`, `canopy`, `mrso` and `gw` are
left out and named under `not_reported`.

The zero-store rule in `AGENTS.md` is about a model whose control volume includes the land
surface and which has no such compartment: rain falls into it, so a zero `snw` says all of
that water went elsewhere, and the budget tests it. Here the control volume starts at the
inflow face. It never sees weather, and any snow or canopy water is upstream of its boundary,
so its value is unknown rather than zero. Version `8.1-ht.1` reported both as zero and was
never archived; no verdict depends on the difference, since every probe that reads them also
needs something this configuration lacks.

`stage` is available from this solver and is not declared. The probes that read it are
daily and supply no `q_in`, so declaring it would only move them from INCOMPLETE to
INCOMPATIBLE.

### Step

PT1H, the one step a merged probe can feed this configuration at. The solver's own step is
adaptive and in seconds, and the row step is only an output setting, so nothing is
resampled. The issue's PT1M reporting belongs to the unmerged #75, and `STEP_SECONDS` holds
only PT1H for that reason.

### Determinism and the future

The solver is deterministic, reads boundary values at the current time, and ignores the seed.
Nothing depends on later rows. Upstream interpolates linearly between the points it is given,
so each row is written as two points, the ramp between rows falling in the first millisecond
of the later row. No row reads a later row's value. Version `8.1-ht.2` had that ramp in the
last millisecond of the earlier row, which was a leak of one millisecond and is why this is
`ht.3`.

### Trivial passes

None: nothing passed. The failure is not trivial either, since a fixed-speed router fails
this probe differently (see `models/rvic`).

## Limits

- One row of cells is the smallest raster that can carry the reach. The two-dimensional
  terms, wetting and drying of a floodplain, and everything #75 asks about are not exercised.
- The last cell is two to four times as deep as its neighbours at the settled states (0.43 m
  against 0.12 m at the third seed's low state) and the first is slightly deep. Both are the
  same in the two reach lengths, so the paired probe cancels them, but they do sit in `channel`.
  Measured and not traced further.
- DG2, upstream's second-order solver, has sloped cells and might not have this staircase.
  The issue names FV1, and DG2 was not run.
- About 7 to 23 s per case on the evaluated grid.

## References

Sharifian, Kesserwani, Chowdhury, Neal and Bates (2023), GMD 16, 2391-2413,
https://doi.org/10.5194/gmd-16-2391-2023. LISFLOOD-FP 8.1 archive,
https://doi.org/10.5281/zenodo.6912932.

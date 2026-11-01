# RVIC under HydroTuring

RVIC 1.1.1 (UW-Hydro, tag `rvic.1.1.1`, commit `afde3841d9b8de243c88a9e539cb3bf13b0fe313`,
GPL-3.0), the routing model of Lohmann et al. (1996, 1998): a linearised Saint-Venant
impulse response, built once before the run and convolved with the inflow. Issue #119.

## Verdict

Gate seeds, Docker on Apple Silicon at `3e86a79`, version `1.1.1-ht.4`: **0 of 1 scored**. One
probe can score this model, and it fails it, by design of the model rather than by any choice
made here.

| Probe | Verdict | Mechanism |
| --- | --- | --- |
| `momentum/wave-celerity-bounds` | FAIL (VIOLATION) | The impulse response has one velocity for the whole record, so the wave travels at the same speed at 8, 16 and 32 m3/s: 0.734, 0.746 and 0.783 m/s across all three states on the three gate seeds, against a Manning `dQ/dA` that rises by 73%. Verified, below |
| The other 36 probes | N/A (INCOMPLETE) | Each needs a flux, store or diagnostic this model does not have: `mrro`, `pr`, `evspsbl`, `stage`, a soil or snow store, a heat flux. Most are also daily. No probe was stopped by anything else |

These results describe one documented point in parameter space. The verdict does not
depend on the point, which is the finding: see "Alternative parameters".

### The failure

`rvic.core.make_uh` builds one Green's function per cell from one `velocity` and one
`diffusion`, and `Rvar.convolve` applies it unchanged for ever. That is what makes RVIC
fast and is stated in its own description. It is the same construction as this probe's
must-fail control, `reference_fixed_celerity`.

Verified two ways on the kept gate-seed runs of this version. The response to the pulse,
in units of the pulse, is the same at all three states to the last bit in all six runs: the
largest difference is exactly 0. And the delay can be computed without running anything,
by sampling `make_uh`'s Green's function from the velocity and diffusivity in `run.json`:
6.057, 5.954 and 5.677 h on the three seeds, each equal to the probe's measured delay to
three decimals.

| Gate seed | `c_obs` low, medium, high (m/s) | `c_kin` low, medium, high (m/s) | residual |
| --- | --- | --- | --- |
| 397273707 | 0.734, 0.734, 0.734 | 0.701, 0.923, 1.214 | +4.6%, -20.5%, -39.6% |
| 904225821 | 0.746, 0.746, 0.746 | 0.713, 0.939, 1.234 | +4.7%, -20.5%, -39.5% |
| 1411177935 | 0.783, 0.783, 0.783 | 0.738, 0.972, 1.278 | +6.1%, -19.4%, -38.8% |

The velocity set was `dQ/dA` at the low state, so the low state should be the one that
fits, and it nearly does. The offset that remains there is the hourly sampling of a
narrow kernel. On the first seed the Green's function puts weights of 0.127 and 0.872 on
hours 1 and 2 of the short reach, a centroid of 1.873 h where `L/v` is 1.585 h, and 7.930 h
where `L/v` is 7.923 h on the long reach. This moves the low state over the 5% line on the
third seed; the two higher states and the ordering fail on every seed whatever it does.

## What the adapter does

Upstream's own `rvic.parameters.parameters` and `rvic.convolution.convolution_init`,
`convolution_run` and `convolution_final`, from config files, as `rvic parameters` and
`rvic convolution` would. No routing is reimplemented. Three small things sit around them:

- **A recorder on `Rvar.get_rof`,** which `convolution_run` already calls every step,
  to read the outflow and the ring. RVIC's history tape holds no storage, and the ring
  is otherwise visible only in a restart file.
- **`WholeSecondStep`,** a subclass of upstream's `DataModel` that rounds the step to the
  whole second. Upstream takes it from two times in days since year 1, gets 3599.9999966 s
  at this date, and then refuses to run because that does not divide the impulse
  response's 3600 s. This is a floating-point problem in upstream at a sub-daily step.
- **Forcing times written on RVIC's own clock,** as the same float sum of days since year
  1 that `Dtime` accumulates, because `DataModel.read` stops if a forcing time differs
  from that clock by a microsecond. They are the row times.

### The mapping

| Item | Value | Source |
| --- | --- | --- |
| Grid | Two cells on one row, source cell draining east to the outlet | Adapter, the smallest grid RVIC accepts. Lat and lon are placeholders: distance comes from `Flow_Distance` |
| `Flow_Distance` | `reach_length_m` | Case. `make_grid_uh` gives a source cell `uh_box` convolved with the river response of the cell it drains to, which at the outlet is the Green's function at the outlet's own `Flow_Distance` |
| `VELOCITY` | Rectangular Manning `dQ/dA` at the first row's `q_in` | Case geometry. RVIC takes a constant here and its samples set it by hand. The first row is the only discharge known without reading ahead |
| `DIFFUSION` | `Q/(2BS)` at the same discharge | Low-inertia hydraulic diffusivity of the declared section, the expression the probe's README cites |
| `uh_box` | An impulse in the first hour | `q_in` is a river inflow, so no hillslope delay |
| `q_in` | Over the outlet-cell area `area_km2`, in kg m-2 s-1, on the source cell | RVIC routes a flux and multiplies back by that same area on output, so the area cancels: it has no effect on `dis` |
| `BASIN_FLOWDAYS`, `CELL_FLOWDAYS` | 50, 2 | `samples/configs/rvic.parameters.pnw.cfg` |
| `SUBSET_DAYS` | 50, not the sample's 10 | The subset path passes a float as an array shape, which every NumPy with a wheel for this Python refuses. At 50 nothing is truncated |
| `REMAP`, `AGGREGATE` | False | The grid is the routing grid, so CDO is not needed |
| `RUN_TYPE` | `drystart` | The 20-day spinup fills the ring, which holds a few hours of water here |

`channel` is the ring from its second entry on, `ring[1:].sum() * dt`: the first entry
leaves during the step, and the rest is still in transit at its end. The ring is a real
store: cumulative `q_in - dis` matches it to 2.5e-16 of the inflow on every gate-seed
case. `Rvar.get_storage` is not used, because it sums the ring over outlets and not over
time, so with one outlet it returns the outflow. That is an upstream bug, harmless to
upstream because no history field reads it.

Upstream's convention releases a row's inflow from that same row, reading the impulse
response's first ordinate, at `t = dt`, as this step's. So the water leaves up to one
step early. It is common to both reach lengths and the paired probe cancels it.

### Stores

`channel` is reported and is the only store declared. `snw`, `canopy`, `mrso` and `gw` are
left out and named under `not_reported`.

Version `1.1.1-ht.4` sets `USER` in the Dockerfile, because RVIC calls `getpass` on import and
crashed in Docker under a uid with no passwd entry. The verdict is unchanged.

Version `1.1.1-ht.1` reported `snw` and `canopy` as zero under `zero_stores`. That was
withdrawn. The zero-store rule in `AGENTS.md` is about a model whose control volume includes
the land surface and which has no such compartment: rain falls into it, so a zero `snw` says
all of that water went elsewhere, and the budget tests it. RVIC's control volume starts at
its inflow. It never sees weather, and any snow or canopy water is upstream of its boundary,
so its value is unknown rather than zero, which is the distinction the rule itself draws for
a discharge-only network. A zero here would also not be able to cost the model anything,
because no budget covers it, so it would only ever be able to help on a probe about a store
the model does not have. No verdict moved: every probe that reads `snw` or `canopy` also
needs a flux or store this model lacks, so all 36 were N/A (INCOMPLETE) before and after.

### Step

PT1H, the step in the issue. The impulse response is built at the step of the `uh_box`
file, so another step is a different parameter file and has not been checked. A daily
RVIC is a normal configuration upstream, but no daily probe supplies `q_in`, so it would
not change what can be scored.

### Determinism and the future

The parameter file is built from static attributes and the first row. The convolution
reads a row at a time. Nothing depends on the seed.

## Alternative parameters

Upstream's `pnw` sample values, `VELOCITY 1`, `DIFFUSION 2000`, on the same three seeds.
The adapter reads them from `HT_RVIC_VELOCITY` and `HT_RVIC_DIFFUSION`, which the harness
never sets.

| Parameters | `c_obs`, all states (m/s) | Residual, low to high | Verdict |
| --- | --- | --- | --- |
| Case geometry at first row (evaluated) | 0.734, 0.746, 0.783 on the three seeds | +5%, -20%, -39% | FAIL |
| Upstream sample, v = 1, D = 2000 | 1.126 on all three seeds | +53 to +61%, +16 to +22%, -7 to -12% | FAIL |

The verdict does not move. No choice of the one velocity can pass: the three targets span
a factor of 1.73 and the criterion wants them ordered.

Declaring the width, slope and roughness rests on the first row of that table. With
upstream's sample values the adapter would read none of them, and could not declare
them, so this model would be N/A (INCOMPATIBLE) on the probe.

## Limits

- This is a single reach, one source cell and one outlet, so the convolution is the whole
  model and RVIC's remapping and multi-source accounting are not exercised at all. The
  issue's conservation identity over a real grid is not tested by any merged probe.
- Python 3.7, NumPy 1.16, netCDF4 1.4, because RVIC relies on `np.float` and on
  `num2date` returning real datetimes. matplotlib is left out: RVIC imports it optionally,
  for diagnostic plots only.
- About 2 to 5 s per case.

## References

Lohmann, Nolte-Holube and Raschke (1996), Tellus A 48(5), 708-721,
https://doi.org/10.3402/tellusa.v48i5.12200. Lohmann, Raschke, Nijssen and Lettenmaier
(1998), Hydrol. Sci. J. 43(1), 131-141, https://doi.org/10.1080/02626669809492107.

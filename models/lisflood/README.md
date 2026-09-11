# LISFLOOD under HydroTuring

LISFLOOD, the distributed rainfall-runoff model the European Commission's
Joint Research Centre runs behind the European and Global Flood Awareness
Systems ([documentation](https://ec-jrc.github.io/lisflood/)), packaged from
[ec-jrc/lisflood-code](https://github.com/ec-jrc/lisflood-code) at release
5.0.0 (`025cff0`), installed as the PyPI release `lisflood-model==5.0.0`.
Requested in [Flood-Lab/HydroTuring#20](https://github.com/Flood-Lab/HydroTuring/issues/20)
by [@kawh1111](https://github.com/kawh1111).

It is the first operational distributed model in the pool, and the first
whose water balance module checks its own budget every step. The adapter
reads the same stores and fluxes that module sums, so what HydroTuring scores
is LISFLOOD's own ledger.

## Licence

EUPL-1.2. The image carries it at `/model/LICENSE`, copied from the pinned
commit.

## The image

conda-forge's `pcraster`, which LISFLOOD 5.0.0 imports throughout (its
dynamic framework, drain-direction and routing operators), is built for
linux-64, osx-64, osx-arm64 and win-64 and not for linux-aarch64, so a
multi-arch image is not possible without building PCRaster from source. The
image is pinned to `linux/amd64`: native on the CI runner, emulated on an
Apple-silicon host, buildable on both. It is `mambaorg/micromamba:2.9.0-debian13`
with the binary stack at the versions of the upstream `environment.yml` for
this tag (python 3.12, pcraster 4.4.2, libgdal-core 3.12.4, numpy 2.2.6,
numba 0.65.1, netCDF4 1.7.4), then `pip install --no-deps lisflood-model==5.0.0`.
The PyPI sdist (`lisflood_model-5.0.0.tar.gz`, sha256 `a1f7bd46…3572`) was
checked against the tag: its `src/lisflood` tree, `LICENSE` and `VERSION` are
identical to commit `025cff0`, so the model code in the image is the pinned
commit's.
The official `jrce1/lisflood` image was not used as a base: it is amd64-only
too, unpinned (`latest`), 8.6 GB uncompressed, and it runs the upstream test
suite at build time.

Two build details, both for emulated builds: micromamba's parallel package
extraction deadlocked under emulation and runs single-threaded, and the
package's `setup.py` calls `gdal-config` for a version string, so the
environment's `bin` is on `PATH` before the pip step.

## The domain

A lumped case is given as the smallest valid LISFLOOD domain: one square cell
whose area is the catchment's (`area_km2`; cell length the square root),
whose local drain direction is a pit, and which carries a channel of the
cell's length, so everything the cell generates leaves through LISFLOOD's own
overland-flow and channel kinematic wave. The clone map is a one-cell lat/lon
grid centred on `latitude_deg`, which is all LISFLOOD reads from coordinates
(the sign of its seasonal snowmelt coefficient); cell length and area are
given explicitly (`gridSizeUserDefined`). The whole cell is the rainfed
"other" land-use fraction.

Switched off, because no probe prescribes them and each needs inputs a
synthetic catchment does not have: `wateruse` (and its demand, region and
smoothing options), `riceIrrigation`, `drainedIrrigation`, `simulateLakes`,
`simulateReservoirs`, `simulatePolders`, `TransLoss`, `openwaterevapo`,
`varfractionwater`, `SplitRouting`, `MCTRouting`, `dynamicWave`, `inflow`,
`indicator`, `TransientLandUseChange`, `simulatePF`, `cropsEPIC`, and every
map and timeseries report. On: `gridSizeUserDefined`. `repMBTs`, which makes LISFLOOD compute and write
its own mass-balance error every step, was on through development (its error
never exceeded 2.1e-13 mm on any case) and is off in the packaged adapter,
which recomputes the same budget from the same terms after the run and
reports the largest step residual in `run.json` (`budget_residual_mm_max_abs`).

Initial state: soil at field capacity, groundwater zones, snow, interception
and overland flow empty, channel at half bankfull (LISFLOOD's default,
`TotalCrossSectionAreaInitValue = -9999`). No prerun: at the lower zone's
default time constant of 100 days, the probe's spinup year covers three and a
half of them.

## Parameters and where they come from

| Parameter | Value | Source |
| --- | --- | --- |
| cell area, cell and channel length | `area_km2`, its square root | static.json |
| `TempSnow` | `snow_threshold_degC` | static.json |
| `SnowMeltCoef` | `degree_day_factor_mm_per_C_day` | static.json |
| soil depths, layers 1b and 2 | test catchment means (748, 1709 mm) scaled together so that the saturated water of the three layers equals `soil_capacity_mm`; the 50 mm top layer is kept | static.json and test catchment |
| LAI | constant, the value at which LISFLOOD's interception capacity `SMax = 0.935 + 0.498 LAI - 0.00575 LAI^2` equals `canopy_capacity_mm` (2.19 for 2 mm) | static.json through the model's equation |
| `UpperZoneTimeConstant`, `LowerZoneTimeConstant` | 10, 100 days | LISFLOOD reference default |
| `GwPercValue`, `LZThreshold` | 0.5 mm/day, 10 mm | LISFLOOD reference default |
| `GwLoss` | 0 ("a closed lower boundary is recommended as a starting value") | LISFLOOD reference default |
| `b_Xinanjiang`, `PowerPrefFlow`, `CalChanMan` | 0.7, 3.5, 2.0 | LISFLOOD reference default |
| `TempMelt`, `SnowSeasonAdj`, `SnowFactor`, `TemperatureLapseRate` | 1.0 degC, 1.0, 1.0, 0.0065 | LISFLOOD reference default |
| `LeafDrainageTimeConstant`, `kdf`, `AvWaterRateThreshold`, frost constants, `beta`, `OFDepRef`, `GradMin`, `ChanGradMin`, `CourantCrit` | as in the reference settings | LISFLOOD reference default |
| `DtSecChannel`, channel routing sub-step | 21600 s (reference 3600 s) | packaging choice, to fit 60 s probe budgets under emulation; see Speed |
| soil hydraulics (theta_s, theta_r, lambda, Van Genuchten alpha, Ksat, three layers) | catchment means | test catchment |
| crop coefficient, crop group, overland Manning's n | 0.9994, 2.692, 0.0933 | test catchment means |
| hillslope gradient, elevation standard deviation | 0.236, 160 m | test catchment means |
| channel Manning's n, bottom width, bankfull depth, side slope, gradient | 0.0468, 8.15 m, 0.331 m, 1, 0.0097 | test catchment medians |

"Test catchment" is `tests/data/LF_ETRS89_UseCase` at the pinned commit,
averaged over its 2847-cell mask. Channel geometry takes medians because
channel dimensions grow with upstream area and the mean is set by the few
large-river cells. `run.json` records every value with its source.

## What the adapter reports

| Column | What it is |
| --- | --- |
| `pr` | the forcing, echoed |
| `evspsbl` | transpiration + evaporation of intercepted water + soil evaporation (`TaWB + TaInterceptionWB + ESActWB`) |
| `mrro` | channel outflow at the outlet over the step (`ChanQAvg * DtSec`), as a depth over the cell |
| `dis` | the same outflow, m3/s |
| `gwex` | minus the lower zone's loss to deep groundwater (`GwLossWB`); zero at the default `GwLoss = 0` |
| `snw` | `SnowCover`, mean of the three elevation zones (the degree-day pack holds no liquid water) |
| `canopy` | interception storage `CumInterception` (plus sealed-surface depression storage, zero here) |
| `mrso` | the three soil layers, `W1a + W1b + W2` |
| `gw` | upper and lower groundwater zones, `UZ + LZ` |
| `channel` | overland flow storage (`WaterDepth`) plus channel water (`ChanM3` over the cell) |

These are exactly the terms of `waterbalance.py`: stored water is channel
plus hillslope (`WaterDepth + SnowCover + LZ + fraction-weighted
(CumInterception + W1 + W2 + UZ) + sealed CumInterSealed`), outgoing water is
outlet flow plus `TaWB + TaInterceptionWB + ESActWB + GwLossWB`. On the
closure probe's 30-day window the budget reconstructed from the reported
columns closes to 1e-13 mm per step, and LISFLOOD's own `MBErrorMM` (checked
with `repMBTs` on) agrees.
With the daily leaf-drainage time constant of one day, interception water
that is not evaporated drains within the step, so `canopy` is zero at the
end of every daily step.

The harness flags this as `suspicious_exact` on the two probes that score
closure: a budget closed to machine precision on every step can mean a store
solved as the residual. It does not here. The stores are LISFLOOD's own
state variables, read after each step, not derived from the fluxes; they stay
inside their physical bounds (`state_bounds` passes, the soil store never
leaves `[0, soil_capacity_mm]`); and LISFLOOD's own water balance module,
which sums the same terms independently, reports the same 1e-13 mm error.
LISFLOOD removes every flux from the store it came from in float64, so it is
exact bookkeeping, as the physical reference models are.

## Inputs

The meteorological reader (`readmeteo.dynamic`) is replaced by one that sets
the same five variables in the same units: `Precipitation = pr * DtDay`,
`Tavg = tas`, and `ETRef = ESRef = EWRef = pet * DtDay`. LISFLOOD wants three
potential evaporations (reference crop, bare soil, open water, normally from
LISVAP); the probe gives one, and all three are set to it. With
`openwaterevapo` off, `EWRef` drives only evaporation of intercepted water.

## Timestep

`DtSec` is the case's step (86400 at PT1D, 3600 at PT1H); rows are fed as
given and nothing is resampled. LISFLOOD converts its per-day parameters with
`DtDay` itself; channel routing sub-steps at `DtSecChannel = 21600 s`, so
four sub-steps a day and one an hour (LISFLOOD uses the smaller of the two).

## Speed

Four probes score the ten-year daily record (4015 rows) with a 60 s budget
per container run, and on this Apple-silicon host the amd64 image runs
emulated. At LISFLOOD's reference settings a ten-year record took 102–111 s.
Four settings bring it to about 50 s; three of them leave every output
bit-identical, one changes the routed numbers.

| Setting | Effect on outputs | 395-row case, stepping time |
| --- | --- | --- |
| reference settings, numba JIT on with a warm cache | reference | 10.0 s (plus about 40 s compiling when the cache is cold, which a read-only container always is) |
| `NUMBA_DISABLE_JIT=1`: the jitted loops run as Python | bit-identical | 10.1 s |
| `repMBTs` off: no per-step mass-balance timeseries through PCRaster | bit-identical | 7.5 s |
| numexpr on one thread (60 expression calls a step on one-cell arrays) | bit-identical | 7.2 s |
| `DtSecChannel` 21600 s instead of 3600 s: 4 routing sub-steps a day instead of 24 | changes routed flow | 5.2 s |

Start-up (imports and initialisation) adds about 5 s to every run. With all
four, a ten-year record takes about 50 s in the container, a 395-row case
about 10 s, the hourly resolution month about 15 s. The routing sub-step is
the one choice here that is about this host rather than the model, and the
result below is compared probe by probe with the reference sub-step.

The probe that forced it is `mass/precipitation-counterfactual`: four
variants of a ten-year record per seed, 60 s each. Run at the reference
sub-step outside the harness, its twelve containers took 101–119 s each, so
the harness would have stopped every one; scored offline with the probe's own
criteria, they pass all five (closure exact; in the +20% variant the added
rain goes 0.13 to evaporation, 0.85 to runoff and 0.02 to storage; runoff
returns 0.85 of the rain added at the top of the ladder). At 21600 s the same
containers take 52–61 s, the upper end when other containers are competing
for the CPU, so on this host the margin is not enough (Result, below).

## Result

`ht run --model lisflood --gate-seeds` on the 19 merged probes, adapter
`5.0.0-onecell.2`: **FAIL (ERROR), 13 of 19 probes passed.**

| Probe | Result | What decides it |
| --- | --- | --- |
| `mass/resolution-invariance` | FAIL (VIOLATION) | runoff moves 13.1% of the rain between PT1H and PT1D on the first gate seed, 8.9% and 9.9% on the others (limit 10%): potential infiltration is a storage multiplied by the step length (below) |
| `mass/area-invariance` | FAIL (VIOLATION) | the water in transit departs by 7.8, 14.0 and 7.2 times its mean across the three seeds (limit 1e-6), and the routed runoff's timing with it: the kinematic waves are nonlinear in volume (below) |
| `mass/precipitation-counterfactual` | FAIL (ERROR) | a container exceeded the 60 s budget; the variants take 52–61 s under amd64 emulation, and every criterion passes when they are allowed to finish (below) |
| `energy/evaporative-partition`, `energy/latent-heat-et-consistency`, `energy/surface-energy-closure` | FAIL (INCOMPLETE) | LISFLOOD reports no latent, sensible or ground heat flux |
| the other 13 | PASS | `energy/pet-consistency` (evaporation 0.90 of demand with the soil wettest, 0.09 driest), `mass/antecedent-monotonicity`, `mass/catchment-closure` (residual 2e-15 of the rain), `mass/causality`, `mass/dry-down`, `mass/extreme-rain` (returns 1.00 of the rain added at the top), `mass/phase-counterfactual` (2.2% of the rain), `mass/response-nonnegativity`, `mass/runoff-bounds` (runoff 0.60–0.64 of the rain), `mass/steady-state` (nothing varies by more than 0.23%), `mass/time-origin-invariance` (bit-identical), `mass/warming-response`, `momentum/routing-conservation` |

The harness flags `suspicious_exact` on the two closure-scoring probes; the
section on what the adapter reports explains why exact closure is LISFLOOD's
bookkeeping and not a store solved as the residual.

### The time-budget ERROR

`mass/precipitation-counterfactual` runs four variants of a ten-year record
per seed with 60 s for each. Under amd64 emulation on this host the packaged
adapter takes 52–61 s per variant, depending on what else the machine is
doing; in the full evaluation and again in a re-run of that probe alone, a
container went just over, and the harness stops a probe at its first timeout.
At LISFLOOD's reference routing sub-step the same variants take 101–119 s.
Scored offline with the probe's own criteria and without the time limit, all
five criteria pass at both sub-steps. At the packaged one (the harness's
outputs for the two seeds it finished, the third seed run outside it): the
budget closes exactly, on the worst seed +20% rain goes 0.12 to evaporation,
0.88 to runoff and 0.01 to storage, and runoff returns 0.85 of the rain added
at the top of the ladder. At the reference sub-step: +20% rain goes 0.13,
0.85 and 0.02, and runoff again returns 0.85. The ERROR is a statement about
emulated amd64 on this machine, not about LISFLOOD's water; a native amd64
runner should fit the budget, which has not been measured here.

### The routing sub-step changes no verdict

Every probe at LISFLOOD's reference routing sub-step (adapter `.1`,
`DtSecChannel` 3600 s, scored by the same harness; `mass/warming-response`
re-scored after `criteria/response.py` changed) against the packaged one:

| Probe | 3600 s (reference) | 21600 s (packaged) |
| --- | --- | --- |
| `energy/evaporative-partition`, `energy/latent-heat-et-consistency`, `energy/surface-energy-closure` | INCOMPLETE | INCOMPLETE |
| `energy/pet-consistency` | PASS, wet ratio 0.904 | PASS, 0.904 |
| `mass/antecedent-monotonicity` | PASS | PASS |
| `mass/area-invariance` | VIOLATION, 7.8 / 14.0 / 7.2 | VIOLATION, 7.8 / 14.0 / 7.2 |
| `mass/catchment-closure` | PASS, residual 2.1e-15 | PASS, 2.4e-15 |
| `mass/causality` | PASS | PASS |
| `mass/dry-down` | PASS | PASS |
| `mass/extreme-rain` | PASS, returns 1.00 | PASS, 1.00 |
| `mass/phase-counterfactual` | PASS, 2.15% | PASS, 2.15% |
| `mass/precipitation-counterfactual` | PASS offline (containers 101–119 s) | ERROR in the harness (containers 52–61 s); PASS offline |
| `mass/resolution-invariance` | VIOLATION, 13.1 / 8.9 / 9.9% | VIOLATION, 13.1 / 8.9 / 9.9% |
| `mass/response-nonnegativity` | PASS | PASS |
| `mass/runoff-bounds` | PASS, 0.638 | PASS, 0.638 |
| `mass/steady-state` | PASS, 0.23% | PASS, 0.23% |
| `mass/time-origin-invariance` | PASS, bit-identical | PASS, bit-identical |
| `mass/warming-response` | PASS | PASS |
| `momentum/routing-conservation` | PASS | PASS |

The routed water per day is the same to the digits the probes report; what
the sub-step moves is within-day timing, which no probe here scores. The
reference run used the harness before `mass/precipitation-counterfactual`
merged; nothing else it scores changed in that merge.

## Mechanisms checked with targeted runs

Each of these was run on the probe's own first gate seed, through the
adapter's `simulate()`, one simulation per process, in the amd64 LISFLOOD
environment, at LISFLOOD's reference routing sub-step (`DtSecChannel` 3600 s)
with its own mass-balance reporting on.

### Step dependence: potential infiltration is a storage scaled by the step

`mass/resolution-invariance` serves one month at the day and the hour. On
the first gate seed, runoff at the hourly step is 13.1% of the rain higher
than at the daily step and evaporation 7.8% lower; on the other two gate
seeds runoff moves 8.9% and 9.9%, just inside the 10% limit, so the probe
fails on its worst seed and the size depends on how the month's rain is
structured inside its days. The first seed is the one analysed below. Term by term, surface
runoff gains 30% of the rain while preferential flow loses 22%, upper-zone
outflow 17% and infiltration 7%; snow and interception barely move. The
cause is one line of `soilloop.py`:

    InfiltrationPot = StoreMaxPervious * (1 - SatFraction) ** PowerInfPot * DtDay

`StoreMaxPervious` is a storage, the Xinanjiang pore space of the top soil,
and multiplying it by `DtDay` makes it a per-step capacity: rain that falls
in a few hours meets one twenty-fourth of the day's capacity at the hourly
step and runs off, the topsoil stays drier, and preferential flow
(`RelSat1 ** PowerPrefFlow`) falls with it.

| Hourly run of the same month | Runoff vs daily, % of rain | Evaporation vs daily, % of rain |
| --- | --- | --- |
| as packaged | 13.1 | 7.8 |
| that line without `* DtDay` | 0.9 | 1.0 |
| as packaged, each day's rain spread evenly over its hours | 0.06 | 0.00 |

Without the factor the two steps agree to within what the physical models
move; with uniform rain inside the day they agree exactly.

The size of the effect depends a little on the choices this package had to
make and its sign not at all. The same month, daily against hourly, with one
choice changed at a time:

| Configuration | Runoff, daily / hourly (mm) | Shift, % of rain | Evaporation, daily / hourly (mm) | Shift, % of rain |
| --- | --- | --- | --- | --- |
| as packaged (LISFLOOD reference defaults) | 158.9 / 186.7 | 13.1 | 51.0 / 34.4 | 7.8 |
| test catchment's calibrated means for the nine calibration parameters | 153.5 / 179.6 | 12.3 | 48.8 / 28.1 | 9.8 |
| soil depths not scaled to `soil_capacity_mm` (about 1080 mm of pore space) | 133.5 / 168.8 | 16.6 | 51.0 / 48.8 | 1.0 |
| LAI at the test catchment's mean (2.84) instead of from `canopy_capacity_mm` | 157.0 / 183.2 | 12.4 | 53.3 / 37.3 | 7.5 |

The calibrated means are `UpperZoneTimeConstant` 17.3 d,
`LowerZoneTimeConstant` 5210 d, `GwPercValue` 1.17, `GwLoss` 0.27,
`LZThreshold` 12.8, `b_Xinanjiang` 1.34, `PowerPrefFlow` 4.10, `CalChanMan`
1.41, `SnowMeltCoef` 4.27, averaged over the test catchment's cells. The dependence is
LISFLOOD's infiltration formulation meeting sub-daily rain intensity, not the
adapter's units: every rate LISFLOOD takes per day it converts with `DtDay`
itself, and the adapter feeds rows at the step it is given.

### Area dependence: the kinematic waves are nonlinear in volume

`mass/area-invariance` tells the same catchment as 250 and as 2500 km2.
Evaporation, soil water, groundwater and snow are bit-identical between the
two, and so is the runoff generated before routing. The routed runoff is
not: on the first gate seed its daily values move by up to 4.6 times the
control's mean runoff and the water in transit by 7.8 times its mean (14.0
and 7.2 on the other two seeds), while the year's total agrees to 0.01%. The
same depth over ten times the area is ten times the volume per metre of
channel and of overland flow, and LISFLOOD routes both with a kinematic wave
storing `A = alpha * Q ** beta` with `beta = 0.6`, so the larger volume
travels on a different schedule. Holding the cell and channel length at the
250 km2 value for both tellings (only the area changes) barely changes this:

| Cell and channel length | Runoff, largest daily difference / mean | Channel store, same | Year's runoff |
| --- | --- | --- | --- |
| square root of the stated area (as packaged): 15.8 and 50 km | 4.6 | 7.8 | 471.80 vs 471.74 mm |
| 15.8 km for both | 4.4 | 7.5 | 471.80 vs 471.83 mm |

So it is the routing's nonlinearity in volume, not the one-cell domain's
length. A lumped water budget written per unit area is invariant to the
area; a distributed model's routing is not, and LISFLOOD's own comment on its
overland flow width says as much ("results will depend on cell size").

## Running it

```bash
ht verify-adapter --model lisflood
ht run --model lisflood --gate-seeds --markdown
```

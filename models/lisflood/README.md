# LISFLOOD under HydroTuring

LISFLOOD, the distributed rainfall-runoff model the European Commission's
Joint Research Centre runs behind the European and Global Flood Awareness
Systems ([documentation](https://ec-jrc.github.io/lisflood/)), packaged from
[ec-jrc/lisflood-code](https://github.com/ec-jrc/lisflood-code) at release
5.0.0 (`025cff0`), installed from the PyPI sdist `lisflood-model==5.0.0`.
Requested in [Flood-Lab/HydroTuring#20](https://github.com/Flood-Lab/HydroTuring/issues/20)
by Yuanhang Liu ([@kawh1111](https://github.com/kawh1111)).

The adapter steps LISFLOOD's own Python framework and reads the stores and
fluxes its water balance module sums, so what HydroTuring scores is
LISFLOOD's own ledger. When a probe prescribes a withdrawal, LISFLOOD's own
water-use module makes it.

## Licence

EUPL-1.2. The image carries it at `/model/LICENSE`, copied from the pinned
commit.

## The image

conda-forge's `pcraster`, which LISFLOOD 5.0.0 imports throughout (its
dynamic framework, drain-direction and routing operators), is built for
linux-64, osx-64, osx-arm64 and win-64 and not for linux-aarch64. The image
is therefore pinned to `linux/amd64`: native on an x86-64 host, emulated on
an Apple-silicon one, buildable on both. HydroTuring's CI runs submitted
models only on manual dispatch, so no pull-request check builds it.

- Base images `mambaorg/micromamba:2.9.0-debian13` and `debian:13-slim`,
  pinned by digest.
- The binary stack at the versions of the upstream `environment.yml` for this
  tag: python 3.12, pcraster 4.4.2, libgdal-core 3.12.4, numpy 2.2.6, numba
  0.65.1, netCDF4 1.7.4. `xarray` 2026.4.0, `nine` 1.2.0 and `future` 1.0.0
  at that file's pip versions.
- LISFLOOD from the PyPI sdist `lisflood_model-5.0.0.tar.gz`, installed with
  `--require-hashes` against sha256 `a1f7bd46ac43573d264210da0bb609c9823055ab1c552b3257dcba16ba103572`.
  Its `src/lisflood` tree, `LICENSE` and `VERSION` were checked identical to
  commit `025cff0`. `--no-deps`, because its requirement list pins test
  tooling and GDAL's Python bindings, neither of which the model imports.
- Two build details, both for emulated builds: micromamba's parallel package
  extraction deadlocked under emulation and runs single-threaded, and the
  package's `setup.py` calls `gdal-config` for a version string, so the
  environment's `bin` is on `PATH` before the pip step.

The official `jrce1/lisflood` image was not used as a base: it is amd64-only
too, unpinned (`latest`), 8.6 GB uncompressed, and runs the upstream test
suite at build time.

## The domain

A lumped case is given as **one representative cell at the resolution its
parameters come from**: the shipped test catchment's 5 km grid, so a cell
5000 m square (25 km²) with a channel 5000 m long (that catchment's median
`chanlength`). Its drain direction is a pit and it carries a channel, so
everything the cell generates leaves through LISFLOOD's own overland-flow and
channel kinematic waves. The clone map is a one-cell lat/lon grid centred on
`latitude_deg`, which is all LISFLOOD reads from coordinates (the sign of its
seasonal snowmelt coefficient); cell length and area are given explicitly
(`gridSizeUserDefined`). The whole cell is the rainfed "other" land-use
fraction.

**Every flux and store is a depth over the cell, and the catchment's stated
area enters only the discharge**, `dis = mrro * area_km2 / 86.4` (mm/day
times km² to m³/s). This follows the `wflow_sbm` package, which moved to one
representative cell for the same reason.

**Why not one cell of the catchment's area.** Adapters `.1` and `.2` did
that: a 250 km² catchment became one 15.8 km cell with a 15.8 km channel. The
stated area then reached LISFLOOD through the geometry built from it:
- overland sideflow is a depth times `PixelArea` over `PixelLength`;
- the overland kinematic wave uses `PixelLength` as its flow width;
- channel storage is a cross-section times `ChanLength`, including the
  half-bankfull water it starts with.

`mass/area-invariance` failed, and targeted runs show where.
- At ten times the area the hillslope did not move: evaporation, soil water,
  groundwater and the runoff generated were identical.
- Only the routed outflow and the water in transit changed, because
  LISFLOOD's kinematic waves are nonlinear in the volume per unit length they
  carry.
- Holding the cell and channel length at 15.8 km while only the area grew
  left the departure about as large as letting both grow. No single geometric
  term is to blame: any cell built from the stated area routes the same depth
  differently.
- The channel medians used were also from a 5 km grid, laid along a 15.8 km or
  50 km channel.

The representative cell takes neither its length nor its area from the stated
area, so every depth is the same whatever the area. "Mechanisms checked with
targeted runs" below has the numbers.

Switched off, because no probe prescribes them and each needs inputs a
synthetic catchment does not have: `riceIrrigation`, `drainedIrrigation`,
`simulateLakes`, `simulateReservoirs`, `simulatePolders`, `TransLoss`,
`openwaterevapo`, `varfractionwater`, `SplitRouting`, `MCTRouting`,
`dynamicWave`, `inflow`, `indicator`, `TransientLandUseChange`,
`TransientWaterDemandChange`, `wateruseRegion`, `groundwaterSmooth`,
`simulatePF`, `cropsEPIC`, `repMBTs`, and every map and timeseries report.
On: `gridSizeUserDefined`, and `wateruse` only when the forcing carries an
`abstr` column (see "A prescribed withdrawal").

Initial state: soil at field capacity, groundwater zones, snow, interception
and overland flow empty, channel at half bankfull (LISFLOOD's default,
`TotalCrossSectionAreaInitValue = -9999`). No prerun: at the lower zone's
default time constant of 100 days, the probe's spinup year covers three and a
half of them.

## Parameters and where they come from

| Parameter | Value | Source |
| --- | --- | --- |
| cell length and area; channel length | 5000 m, 25 km²; 5000 m | shipped test catchment grid; `chanlength` median |
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
| `DtSecChannel`, channel routing sub-step | 3600 s | LISFLOOD reference default |
| soil hydraulics (theta_s, theta_r, lambda, Van Genuchten alpha, Ksat, three layers) | catchment means | test catchment |
| crop coefficient, crop group, overland Manning's n | 0.9994, 2.692, 0.0933 | test catchment means |
| hillslope gradient, elevation standard deviation | 0.236, 160 m | test catchment means |
| channel Manning's n, side slope | 0.0468, 1 | test catchment medians over all cells |
| channel bottom width, bankfull depth, gradient | 4.34 m, 0.165 m, 0.066 | test catchment medians over its 961 headwater cells |

"Test catchment" is `tests/data/LF_ETRS89_UseCase` at the pinned commit,
averaged over its 2847-cell mask. Channel dimensions and low flows grow with
upstream area. The representative cell is a pit that drains only its own
25 km², so it takes them from the catchment's 961 headwater cells, whose
upstream area is one cell: the channel's bottom width, bankfull depth and
gradient, and the environmental-flow threshold. The threshold had two
sourced definitions:
- the headwater median of `ad_dis_nat_10`, each cell's 10th-percentile natural
  discharge: 0.0700 m3/s, with quartiles 0.035 and 0.152;
- the median over all cells of `ad_dis_nat_10` per unit upstream area, on the
  cell's 25 km²: 0.0813 m3/s.

The adapter uses the first. It takes the reserve from the same cells as the
channel it is applied to, and does not assume that low flow scales with area.
The second is a sensitivity row under "A prescribed withdrawal". Channel
Manning's n stays at its all-cell median of 0.0468 (headwater median 0.0613),
and channel length at 5000 m, the grid's cell length (headwater median
3653 m). `static.json`'s `baseflow_coefficient` has no LISFLOOD
counterpart (LISFLOOD's recession is two linear zones with their own time
constants) and is not used; `run.json` lists it. `run.json` records every
value with its source.

**Required inputs.** `latitude_deg`, `area_km2` and `canopy_capacity_mm` are
required.
- The latitude centres the clone cell and sets the hemisphere of LISFLOOD's
  seasonal snowmelt coefficient.
- The area scales `dis`.
- The canopy capacity sets the LAI.

If one is missing, the adapter stops with an error before it stages anything,
rather than inventing a value. `run.json` records the latitude it used, with
its source. `snow_threshold_degC`, `degree_day_factor_mm_per_C_day` and
`soil_capacity_mm` fall back to LISFLOOD's reference defaults and the test
catchment's soil depths when absent, and `run.json` names the source used.
Every probe supplies all six.

## What the adapter reports

| Column | What it is |
| --- | --- |
| `pr` | the forcing, echoed |
| `evspsbl` | transpiration + evaporation of intercepted water + soil evaporation (`TaWB + TaInterceptionWB + ESActWB`) |
| `mrro` | channel outflow at the outlet over the step (`ChanQAvg * DtSec`), as a depth over the cell |
| `dis` | `mrro * area_km2 / 86.4`, m3/s |
| `gwex` | what leaves the reported stores to the outside, negative: the lower zone's loss to deep groundwater (`GwLossWB`, zero at the default `GwLoss = 0`) plus, when `abstr` is prescribed, what the water-use module actually withdrew (`abstraction_GW_actual_M3 + withdrawal_CH_actual_M3`, plus lake and reservoir abstraction, less return flow to the channel; the last three are zero here), over the cell |
| `snw` | `SnowCover`, mean of the three elevation zones (the degree-day pack holds no liquid water) |
| `canopy` | interception storage `CumInterception` (plus sealed-surface depression storage, zero here) |
| `mrso` | the three soil layers, `W1a + W1b + W2` |
| `gw` | upper and lower groundwater zones, `UZ + LZ` |
| `channel` | overland flow storage (`WaterDepth`) plus channel water (`ChanM3` over the cell) |

These are exactly the terms of `waterbalance.py`:
- **stored water** is channel plus hillslope: `WaterDepth + SnowCover + LZ`,
  fraction-weighted `CumInterception + W1 + W2 + UZ`, and sealed
  `CumInterSealed`;
- **outgoing water** is outlet flow plus `TaWB + TaInterceptionWB + ESActWB + GwLossWB`,
  and, with water use on, the withdrawals it counts in `wateruseCum` and
  `IrriLossCUM`.

The adapter recomputes that budget after every run from the same terms and
reports the largest step residual in `run.json` (`budget_residual_mm_max_abs`).
LISFLOOD's own reporting of it (`repMBTs`) is off. During development, with it
on, LISFLOOD's error never exceeded 2.1e-13 mm on any case.

With the daily leaf-drainage time constant of one day, interception water
that is not evaporated drains within the step, so `canopy` is zero at the end
of every daily step and non-zero at PT1H.

**Exact closure.** The harness flags `suspicious_exact` where a budget
closes to machine precision on every step, which can mean a store solved as
the residual. No reported store is solved as the budget's residual here:
- **Soil, snow, interception and groundwater** are LISFLOOD's state
  variables, each updated by removing a flux from the store it leaves.
- **Overland and channel water** work the other way round. LISFLOOD sets the
  storage from discharge through Manning's relation (`routing.py`,
  `surface_routing.py`). It then takes the step's average outflow as the
  continuity residual of that storage (`kinematic_wave_parallel_tools.py`),
  clamping a negative average to zero.

That residual is LISFLOOD's routing, local to the channel and exact by
construction. A clamp that fired would show as a budget residual, and none
has. The soil store stays inside its physical bounds (`state_bounds` passes).

## Inputs

The meteorological reader (`readmeteo.dynamic`) is replaced by one that sets
the same five variables in the same units: `Precipitation = pr * DtDay`,
`Tavg = tas`, and `ETRef = ESRef = EWRef = pet * DtDay`.
- **Three potential evaporations from one.** LISFLOOD wants three potential
  evaporations (reference crop, bare soil, open water, normally from LISVAP);
  the probe gives one, and all three are set to it. With `openwaterevapo` off,
  `EWRef` drives only evaporation of intercepted water. Total demand is about
  `pet`: potential transpiration is the crop coefficient times ET0 times one
  minus the canopy term, less interception evaporation, and bare-soil
  evaporation is ES0 times the canopy term.
- **Hourly temperature.** LISFLOOD documents `Tavg` as the daily mean even at
  sub-daily steps. At PT1H the probe's hourly temperature is fed as given,
  because the contract forbids resampling.

**Snow is more than the degree-day factor.** On top of `SnowMeltCoef`,
LISFLOOD applies terms keyed to the forcing's day of year:
- a seasonal melt coefficient of ±0.5 mm/degC/day at `SnowSeasonAdj = 1`;
- a summer ice-melt term on any pack left between days 165 and 257;
- melt enhanced by the rain depth of the step, which is itself step-dependent.

They are the model's own physics, not adapter inputs.
`mass/time-origin-invariance`'s 28-year shift preserves the day of year, so it
cannot see them.

## A prescribed withdrawal (`abstr`)

When the forcing carries an `abstr` column, LISFLOOD's water-use option
(`wateruse`) is on for that run, whatever the column holds. A run whose forcing
has no such column is configured exactly as before. The natural and irrigated
variants of `mass/human-abstraction` therefore run one configuration and
differ only in the demand.

- **The demand.** Every step, before the water-use module runs, the adapter
  sets LISFLOOD's industrial demand (`IndustrialDemandMM`) to `abstr * DtDay`.
  With `TransientWaterDemandChange` off, LISFLOOD reads its demand maps once,
  so this is where a demand that changes through the year enters. Domestic,
  livestock and energy demand are zero, there is no irrigated land use, and
  lakes, reservoirs and non-conventional sources are off.
- **Net, not gross.** `abstr` is already net of return flow, so the industrial
  consumptive-use fraction is 1 (LISFLOOD's reference value is 0.15). LISFLOOD
  then withdraws the consumptive use only, and its return flow from
  groundwater users to the channel is zero.
- **Where the water comes from.** LISFLOOD splits the demand by
  `FractionGroundwaterUsed`, 0.168: the shipped test catchment's mean of
  `fracgwusedNew`, the map the reference settings name.
  - That share is subtracted from the lower groundwater zone
    (`LZ -= abstraction_GW_actual_M3`) with no availability check, so it is
    always taken in full. LZ may fall below zero. In the packaged runs it
    never does: its lowest value, 0.0016 mm, is its first step from empty,
    and on withdrawal days the groundwater store (UZ + LZ) stays above
    10.3 mm.
  - The rest is asked of the channel, limited to the channel water above an
    environmental-flow reserve, `ChanM3 - EFlowThreshold * DtSec`.
    `EFlowThreshold` is 0.0700 m3/s, the median of `ad_dis_nat_10` over the
    test catchment's headwater cells (see the parameter notes). What the channel cannot give is recorded as a shortage
    (`areatotal_shortage_SW_M3`) and is not taken from anywhere else. The
    channel withdrawal is removed inside the kinematic-wave routing, one
    sub-step at a time.
- **What is declared.** `gwex` is minus what the module actually removed,
  `abstraction_GW_actual_M3 + withdrawal_CH_actual_M3`, never the
  prescription. Lake and reservoir abstraction are added and the return flow
  subtracted, both zero here. LISFLOOD's own water balance counts the same
  terms (`IrriLossCUM`, `wateruseCum`), and with the withdrawal on the
  adapter's budget still closes to 1.8e-13 mm per step. `run.json` carries a
  `water_use` block: the option's inputs with their sources, the prescribed
  total, what came from groundwater and from the channel, the channel
  shortage, the share withdrawn and the lowest LZ over the record.

**Zero withdrawal changes nothing.** Two other probes' cases were run as staged
and again with a column of zeros added, which switches the option on:
`mass/catchment-closure` (395 rows) and `mass/steady-state` (1460 rows). Both
result tables are byte-identical to the originals (sha256 99589bb6de08 and
5e67467799ce). The second run of each shows the option on and nothing
withdrawn.

**What LISFLOOD withdraws on `mass/human-abstraction`.** LISFLOOD's rule on
this domain is plain:
- it takes the groundwater share of the demand in full;
- it takes the rest only from channel water in the water region, above the
  environmental-flow reserve;
- it records what the channel cannot give as shortage, and takes that from
  nowhere else.

The water region here is one headwater cell, so most of the prescription is
left untaken. These figures come from the probe's own cases, windows and
criteria on its three gate seeds, run outside the harness's 60 s limit. They
run from the sourced reserve to no reserve at all:

| Gate seed | Withdrawn, sourced reserve | Residual | Withdrawn, no reserve | Residual |
| --- | --- | --- | --- | --- |
| 1129545695 | 17.0% | 82.9% | 31.6% | 68.4% |
| 1636497809 | 17.2% | 82.7% | 33.1% | 66.7% |
| 2143449923 | 17.4% | 82.7% | 33.1% | 67.1% |

- **The bracket.** LISFLOOD withdraws between 17.0% and 33.1% of the 418 mm
  prescribed over the record. That leaves a residual of 82.9% down to 66.7% of
  the 380 mm scored, against a 5% limit.
- **Where the withdrawal comes from.** The groundwater share is 70.2 mm on
  every seed. With the sourced reserve the channel gives 1.0 to 2.3 mm and
  leaves 345 to 347 mm as shortage. With no reserve it gives 62 to 68 mm.
- **Other criteria.** `closure` and `state_bounds` pass on every seed.
- **Where it shows.** Water taken from the lower zone returns as less
  baseflow, so runoff carries the withdrawal and storage barely moves
  (-0.03 mm).

The channel gives little because LISFLOOD's channel abstraction draws on the
water held in the region's channels at the start of the step, above the
reserve, not on what flows through them during the step.
- The reserve, 0.0700 m3/s held for a day, is 0.24 mm over the cell.
- On this one cell the channel and overland store averages 0.41 mm over the
  record and 0.22 mm on withdrawal days. It exceeds the reserve on 1267 of
  4015 days.
- The cell's outflow is below the threshold on 814 of 4015 days (20%), where
  a 10th-percentile flow would be undercut on about 10%. The old all-cell
  median, 0.2604 m3/s, was undercut on 2772 of 4015 days with the `.3` cell.
- In the test catchment a water region's median size is 121 cells; here it is
  one.

On seed 1129545695, runs that change one water-use input at a time show which
input sets the share:

| Run | From groundwater | From channel | Channel shortage | Withdrawn | Residual | Lowest LZ |
| --- | --- | --- | --- | --- | --- | --- |
| packaged: `EFlowThreshold` 0.0700 m3/s, headwater median | 70.2 mm | 1.0 mm | 346.8 mm | 17.0% | 82.9% | 0.0016 mm |
| `EFlowThreshold` 0.0813 m3/s, per unit upstream area | 70.2 mm | 0.7 mm | 347.0 mm | 17.0% | 83.0% | 0.0016 mm |
| no reserve: `EFlowThreshold` 0 | 70.2 mm | 61.7 mm | 286.1 mm | 31.6% | 68.4% | 0.0016 mm |
| all from groundwater: `FractionGroundwaterUsed` 1 | 418.0 mm | 0 | 0 | 100% | 0.0% | -9.2 mm |

- **LZ.** LZ is recorded at the end of every step. In the first three runs its
  lowest value is its first step, filling from empty. On withdrawal days in
  the packaged run the groundwater store (UZ + LZ) stays above 10.3 mm.
- **All from groundwater.** Asked for everything, LISFLOOD takes all of it and
  lets LZ fall to -9.2 mm, and the probe passes.
- **Budget.** Every run's budget closes to 1.8e-13 mm per step.
- **Geometry.** With the all-cell channel geometry of `.3`, no reserve gave
  43.6%. A headwater channel holds less water, so it gives less.

`FractionGroundwaterUsed` stays the test catchment's 0.168, and every other
input stays the test catchment's value; none is set to pass. What LISFLOOD
removes is what `gwex` declares.

## Timestep

`DtSec` is the case's step (86400 at PT1D, 3600 at PT1H); rows are fed as
given and nothing is resampled. LISFLOOD converts its per-day parameters with
`DtDay` itself. Channel routing sub-steps at LISFLOOD's reference
`DtSecChannel = 3600 s`: 24 sub-steps a day, one an hour.

## Speed

Three settings, each checked to leave every output byte unchanged:

| Setting | Why |
| --- | --- |
| `NUMBA_DISABLE_JIT=1` | LISFLOOD's soil and interception loops are numba-jitted with parallel loops over pixels; on one cell they gain nothing, and a read-only container cannot keep the compilation cache, so every run would spend about 40 s compiling. As Python, with a warm cache, they run at the same speed and every output column is bit-identical |
| `repMBTs` off | LISFLOOD's per-step mass-balance timeseries go through PCRaster, about a quarter of the run time; the terms the adapter reads are set regardless of the option, and every output is bit-identical |
| numexpr on one thread | its expressions run on one-cell arrays, where a thread pool has nothing to share; outputs bit-identical |

The adapter's own per-step store sums run on the numpy arrays under
LISFLOOD's vegetation-fraction wrappers.

These timings were measured with the `.3` image under amd64 emulation, on an
Apple-silicon host that was also evaluating other models (load average 15 to
16). `.4` does the same work per step: it ran the same 395- and 1460-row cases
in 14.1 and 34.7 s of wall time, beside the archive run.

| Case | Rows | Initialise | Run | Per step |
| --- | --- | --- | --- | --- |
| `mass/catchment-closure`, scored window | 395 | 5.2 s | 8.4 s | 21 ms |
| `mass/steady-state` | 1460 | 5.2 s | 30.7 s | 21 ms |
| `mass/catchment-closure`, whole ten-year record | 4015 | 5.2 s | 91.7 s | 23 ms |

A ten-year daily record takes 97 s there. That is over the 60 s budget of
`mass/precipitation-counterfactual` and `mass/human-abstraction`. Left to
finish, the first ten-year cases of a `.5` run took 87 and 92 s; that run's log
records a load average of 15.6 at its start. The harness now kills a
container at its time budget, so the archived ERROR rows record no wall time
of their own.

## Result

**FAIL (ERROR), 15 of 18 probes passed, 3 N/A (INCOMPLETE).** These are the
rows of the full gate-seed run of `5.0.0-onecell.5`, made on the emulated host
described under "Native re-run". The verdict is ERROR because two probes ran
out of time on that host. Any ERROR among the scored probes makes the verdict
FAIL (ERROR), whatever the other probes score.

The `mass/extreme-event-closure` row was added when that probe merged, from
`5.0.0-onecell.5` on the same host. Its five twenty-year cases (7,665 rows
each) finished in 134 s together, case generation included, against a 300 s
budget for each case, and every wet event closes to 1e-13 mm.

**The two ERROR rows are provisional.** They come from the emulated host and
are to be replaced by a native x86-64 evaluation; "Native re-run" gives the
commands and the row replacement.

- **N/A (INCOMPLETE), 3, not scored:** `energy/evaporative-partition`,
  `energy/latent-heat-et-consistency` and `energy/surface-energy-closure`.
  LISFLOOD reports no heat fluxes, so these probes cannot ask it anything.
  They are neither a pass nor a fail, and do not decide the verdict.
- **ERROR, 2:** `mass/precipitation-counterfactual` and
  `mass/human-abstraction`. The container exceeded the 60 s budget, because a
  ten-year record takes about 90 s under emulation (87 and 92 s when left to
  finish).
  These rows come from the emulated host.
  - Run outside the limit on all three gate seeds,
    `mass/precipitation-counterfactual` passes every criterion.
  - `mass/human-abstraction` fails `human_abstraction` there. LISFLOOD
    withdraws 17.0 to 33.1% of the prescription, from the sourced reserve to
    none. That leaves a residual of 82.9 to 66.7% against a 5% limit;
    `closure` and `state_bounds` pass.
  - On a host fast enough for the budget, the first should PASS and the second
    be VIOLATION. The model's verdict would then be FAIL (VIOLATION), with 16
    of 18 probes passed and 3 N/A.
- **VIOLATION, 1:** `mass/resolution-invariance`. Rain that falls within an
  hour runs off, so `mrro` differs by 13.0% of `pr` between PT1H and PT1D,
  against a 10% limit.
- **PASS, 15:**
  - `energy/pet-consistency`;
  - `mass/antecedent-monotonicity`, `mass/area-invariance`,
    `mass/catchment-closure`, `mass/causality`, `mass/dry-down`,
    `mass/extreme-event-closure` and `mass/extreme-rain`;
  - `mass/phase-counterfactual`, `mass/response-nonnegativity`,
    `mass/runoff-bounds`, `mass/steady-state`, `mass/time-origin-invariance`
    and `mass/warming-response`;
  - `momentum/routing-conservation`.

  The budget closes to 1e-13 mm per step. The harness flags `suspicious_exact`
  on `mass/catchment-closure` and `mass/time-origin-invariance`; "What the
  adapter reports" says why the closure is exact.

Against the `.2` rows:
- `mass/area-invariance` moved from VIOLATION to PASS with the representative
  cell.
- `mass/resolution-invariance` is VIOLATION in both (13.1% then, 13.0% now).
- `mass/precipitation-counterfactual` is ERROR on this host in both.
- `mass/human-abstraction` is new since `.2`.
- The three energy-flux probes were FAIL (INCOMPLETE) under the earlier
  roll-up and are N/A (INCOMPLETE) under main's.
- No other probe's verdict moved.

The standing counts passes out of the 18 probes that could score LISFLOOD; the
three N/A energy-flux probes are in neither number.

Against the `.3` rows, `.4` changes the environmental-flow reserve and the
channel's bottom width, bankfull depth and gradient to the headwater values.
No probe's verdict or reason moved:
- the routing-sensitive probes still pass, and `mass/area-invariance` still
  departs by exactly 0;
- `mass/resolution-invariance` is 13.0% in both;
- the two ten-year probes are ERROR on this host in both.

`.5` requires `latitude_deg`, `area_km2` and `canopy_capacity_mm` and records
the latitude it used in `run.json`. Every probe supplies them, and re-run on
all 20 probes with the same harness, no row moved against `.4` except its
version and the contract row's timing. Re-run again with main's harness from
`89f2f14`, which kills a timed-out container and counts passes out of the
scored probes, no row moved either; these are the archived rows.

## Mechanisms checked with targeted runs

**`mass/resolution-invariance`: rain within the hour runs off.** LISFLOOD's
potential infiltration is a pore-space storage multiplied by the step length,
`InfiltrationPot = StoreMaxPervious * (1 - SatFraction)^PowerInfPot * DtDay`
(`soilloop.py`). Summed over a day's 24 hourly steps that capacity is the daily
one, but a storm arrives in a few of those hours. What does not infiltrate in
its hour becomes surface runoff. The runs below go through the packaged
adapter, on one 30-day case of the probe with 212.5 mm of rain over the scored
stretch:

| Run | Runoff | Evaporation | Runoff against PT1D, share of rain |
| --- | --- | --- | --- |
| PT1D | 158.3 mm | 51.0 mm | |
| PT1H | 185.8 mm | 34.4 mm | +13.0% |
| PT1H, `InfiltrationPot` not multiplied by `DtDay` | 156.8 mm | 48.8 mm | -0.7% |
| PT1H, each day's rain spread evenly over its hours | 158.1 mm | 51.0 mm | -0.08% |

Every run's budget closes to 6e-14 mm per step.

**`mass/area-invariance`: the cell, not the model.** The departure is the
harness's measure: the largest difference between the two runs over the
control's mean magnitude. The pairs are the probe's first gate seed
(66780378) at the stated 250 and 2500 km²:

| Cell | What the tenfold area changes | `mrro` departure | `channel` departure |
| --- | --- | --- | --- |
| packaged, 5 km representative cell | nothing | 0 | 0 |
| `.1`/`.2`, one cell of the catchment's area | cell and channel 15.8 km to 50 km | 4.6 | 7.8 |
| the same, with length held at 15.8 km | the cell area only | 4.4 | 7.5 |

- In both old pairs, evaporation, soil water and groundwater were identical.
  A component run of the first pair also found the hillslope's generated
  runoff identical.
- `.1` failed on all three gate seeds (`channel` 7.8, 14.0 and 7.2).
- The packaged cell passes on all three, with every variable identical.

**`mass/human-abstraction`.** See "A prescribed withdrawal" above: LISFLOOD
withdraws 17.0 to 17.4% of the prescription on the three gate seeds with the
sourced reserve, and 31.6 to 33.1% with none. It takes the groundwater share
in full, and from the channel only what one headwater cell holds above the
reserve.

**`mass/precipitation-counterfactual` outside the time limit.** The probe's
own cases, windows and criteria, run outside the 60 s limit: every criterion
passes on all three gate seeds, with container walls of 87 to 104 s. Each row
below is one seed's +20% variant, with its extra rain divided into what
evaporated, ran off and stayed stored:

| Gate seed | Extra rain | Evaporation | Runoff | Storage |
| --- | --- | --- | --- | --- |
| 872466880 | 1737 mm | 0.116 | 0.859 | 0.025 |
| 1379418994 | 1634 mm | 0.116 | 0.877 | 0.007 |
| 1886371108 | 1589 mm | 0.133 | 0.851 | 0.016 |

On the same seed, the partition matches to three decimals what earlier
adapters gave:
- seed 1886371108 on `.1`, the catchment-sized cell at `DtSecChannel`
  3600 s;
- seed 1379418994 on `.2`, at 21600 s.

Neither the cell nor the routing sub-step moves it. The spread between seeds
is the weather, and earlier comparisons across two seeds read it as a
sub-step effect.

**The two ERROR rows come from the host's speed.** Both probes score a ten-year
daily record (4015 rows with spinup) in a container with a 60 s budget, and the
harness stops a probe at its first timeout.
- Run outside the limit on this host, each `.4` variant took 96 to 110 s, two
  containers at a time.
- Left to finish, the first ten-year cases of a `.5` run took 87 and 92 s. The
  harness now kills a container at its 60 s budget, so a timed-out case no
  longer runs on into the next probe.
- Per step, a ten-year record runs at the speed of a 30-day case (21 to 23 ms).

Nothing in the model slows down; ten years of LISFLOOD's Python framework
under amd64 emulation simply need more than 60 s.

## Native re-run

**The two ERROR rows are provisional:** they come from an emulated host, and
this re-run replaces them.

Every row this package has archived was produced on an Apple-silicon host
running the amd64 image under emulation. Two probes score a ten-year daily
record in a container with a 60 s budget, `mass/precipitation-counterfactual`
and `mass/human-abstraction`. Under emulation a ten-year run takes about
90 s (87 and 92 s when left to finish), so their archived rows are ERROR
because of the host's speed alone. On an
x86-64 Linux host, from the repository root:

```bash
rm -f /tmp/lisflood-native.csv   # ht run --csv appends
docker build -t hydroturing/lisflood:5.0.0-onecell.5 -f models/lisflood/Dockerfile models/lisflood
ht verify-adapter --model lisflood
ht run --model lisflood --gate-seeds --csv /tmp/lisflood-native.csv --markdown
```

To replace the archived rows with the native ones, drop every `lisflood` row
from `models/result.csv` and append the new file's rows without its header,
leaving every other model's bytes as they are:

```bash
grep -v ',lisflood,' models/result.csv > /tmp/result.csv
tail -n +2 /tmp/lisflood-native.csv >> /tmp/result.csv
mv /tmp/result.csv models/result.csv
ht verify-adapter --model lisflood --csv models/result.csv
git diff origin/main -- models/result.csv   # only lisflood lines
```

Then rewrite everything that describes the emulated host:
- the verdict, its count and the ERROR narrative under "Result", and the
  time-budget paragraph under "Mechanisms checked with targeted runs";
- the ten-year timings in this README, in `model.yaml`'s closing comment and
  in the top-level README row;
- the share of `abstr` LISFLOOD withdraws and the residual, wherever they are
  quoted, from the native `mass/human-abstraction` row;
- the top-level README row and the three rows of `site/index.html`.

If the ten-year cases still exceed 60 s on the native host, the ERROR is the
model's own at its reference routing sub-step. The archived cases need about
1.5 times this host's speed to fit (87 and 92 s against 60 s).

## Running it

```bash
ht verify-adapter --model lisflood
ht run --model lisflood --gate-seeds --markdown
```

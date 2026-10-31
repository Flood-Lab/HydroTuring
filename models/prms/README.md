# PRMS 6.0.0 under HydroTuring

PRMS, the USGS Precipitation-Runoff Modeling System (Markstrom et al. 2015,
[tm6-B7](https://doi.org/10.3133/tm6B7)), from the official `prms_6.0.0_linux.zip`. Proposed by
Yan Chen as [Flood-Lab/HydroTuring#109](https://github.com/Flood-Lab/HydroTuring/issues/109).
CC0 1.0 per the release README, which the image carries at `/model/Readme_PRMS_6.0.0.pdf`.

## What runs

The release's own Linux executable, unmodified, as one HRU draining to one stream segment.
`ht_adapter.py` writes the control, parameter and climate-by-HRU files, runs `prms` once (twice on
a withdrawal case, see below), reads the statvar file and maps it onto `/io`. Nothing is
reimplemented.

USGS publishes the zip at a bare URL with no commit to name, so it is pinned by SHA-256
(`450c5b65...d00423`). The build refuses any other file, then runs the bundled Sagehen
`climate_hru` project and refuses to finish unless the `Total for run` line of the water budget
equals upstream's recorded `output-test`. Base image `debian:bookworm-slim` by digest. The image
is 45 MB and a decade of daily rows takes about 0.1 s.

Modules: `climate_hru` for temperature, precipitation, potential evaporation and solar radiation,
`transp_tindex`, `intcp`, `snowcomp`, `srunoff_smidx`, `soilzone`, `gwflow`, `muskingum`.

## What is reported

PRMS works in inches, so every value goes through `statvar`, which prints six decimals.

| Column | Native | Note |
| --- | --- | --- |
| `pr` | the forcing | echoed from the CSV as text |
| `evspsbl` | `basin_actet` | see the print quantum below |
| `mrro` | `basin_stflow_out` | what leaves the routing, so the water in the channel is in `channel` and not counted twice |
| `dis` | `mrro` over the case's `area_km2` | PRMS itself is given a fixed reference area, see below |
| `gwex` | `-(basin_gwsink + removed abstr)` | `gwsink_coef` is at its default of 0, so this is the withdrawal |
| `snm` | `snowmelt`, plus `net_rain` where `compute_infil` passes it around the pack | |
| `sbl` | `snow_evap`, plus `hru_intcpevap` when `intcp_form` is snow | |
| `mrso` | `basin_soil_moist + basin_ssstor` | capillary and gravity/preferential reservoirs |
| `snw` | `basin_pweqv` | |
| `canopy` | `basin_intcp_stor` | |
| `gw` | `basin_gwstor` | |
| `channel` | `basin_segment_storage` | PRMS's cumulative inflow minus outflow of the segment |

That answers the proposal review's second and third questions: `mrro` is separate from `dis`, and
both `gw` and `channel` are exposed. `mrro` is routed runoff. Upstream's pre-routing term
`basin_stflow_in` carries water that `channel` still holds.

**Print quantum.** `statvar` rounds to 1e-6 in (2.5e-5 mm), which is larger than the suite's 1e-6 mm
state-bound tolerance. Two adapter rules follow, and `run.json` counts the rows each one moves:

- `evspsbl` within one quantum of the row's `pet` is set to `pet`. PRMS holds actual at or below
  potential by construction (`avail_potet` in `soilzone`), and the print can read a quantum above.
- `canopy` above capacity by no more than one quantum is set to capacity.

Nothing else is moved, and a value further out is left to fail. Step residuals from rounding reach
1e-4 mm; the cumulative residual on the worst closure seed is -7.8e-4 mm out of 7708 mm of rain.

## What is not reported

- `basin_imperv_stor`: `hru_percent_imperv` is 0, so it is identically zero. Depression storage:
  `dprst_flag` is 0. Both are named under `not_reported` in `run.json`.
- `zero_stores` is `[canopy]` only on a case with `canopy_capacity_mm = 0`, where the HRU is given
  no cover. Otherwise it is empty, because PRMS has real snow and interception stores.
- Heat fluxes, `stage`: PRMS's water budget does not produce them, so the energy probes that need
  them and the stage probes are N/A (INCOMPLETE).

## Parameters

No calibration exists for the synthetic catchment. Anything not set is the value PRMS declares for
it in `declparam`. This is one documented point in parameter space.

| Parameter | Value | Source |
| --- | --- | --- |
| `soil_moist_max`, `sat_threshold` | `soil_capacity_mm` split 0.4457 / 0.5543 | case; the split is the ratio of Sagehen's HRU medians, 3.6759 and 4.5706 in |
| `srain_intcp`, `wrain_intcp`, `snow_intcp` | `canopy_capacity_mm` / `covden` | case, so the basin-mean store tops out at the capacity |
| `cov_type` | 3 (trees), or 0 (bare) with no canopy | PRMS default, 0 with no canopy |
| `covden_sum`, `covden_win` | 0.5 | PRMS default, 0 with no canopy |
| `rad_trncf` | 0.5 | PRMS default; 1 with no canopy (tm6-B7's fit, 0.9917 exp(-2.7557 covden) + 0.0083, at zero cover) |
| `tmax_allsnow` | `snow_threshold_degC` | case |
| `tmax_allrain_offset` | 1 C | PRMS's declared default |
| `hru_elev`, `hru_lat` | case | |
| `hru_area` | 1 km2 | reference, see below |
| `hru_percent_imperv` | 0 | no impervious area in the case |
| `gwflow_coef`, `ssr2gw_rate`, `slowcoef_lin`, `fastcoef_lin` | 0.015, 0.1, 0.015, 0.1 per day | PRMS defaults |
| `slowcoef_sq`, `fastcoef_sq` | 0.1, 0.8 | PRMS defaults |
| `snowinfil_max`, `settle_const`, `cecn_coef` | 2 in, 0.1, 5 | PRMS defaults |
| `carea_max`, `smidx_coef`, `smidx_exp` | 0.6, 0.005, 0.3 | PRMS defaults |
| `gwstor_init` | 2 in | PRMS default; soil starts empty at `soil_moist_init_frac` = 0 |
| `gwsink_coef`, `gwstor_min` | 0 | PRMS defaults |
| `transp_beg`, `transp_end`, `transp_tmax` | 1, 13, 1 | PRMS defaults, so transpiration runs all year |
| `K_coef`, `x_coef` | 1 h, 0.2 | PRMS defaults |
| `soilzone_aet_flag` | 1 | upstream's Sagehen control files |

**`hru_area` is a fixed 1 km2.** One HRU on one segment has no process that reads its area. PRMS
multiplies by it only to report volumes, in single precision, and divides again, so giving it the
case's area changes which way a value rounds at statvar's sixth decimal. Measured with the case's
area: 1e-6 in on one step of `basin_gwstor`, which puts `gw` 3.2e-6 relative off, over the probe's
1e-6. So `dis` alone takes the case's area, and `mass/area-invariance` passes by construction.

## Forcing

Precipitation and potential evaporation go in directly as depths per row, so `pet` is honoured and
none of PRMS's PET modules run. Two inputs are mocked from the row alone, with no calendar:

- `tmax` and `tmin` both equal `tas`, because the probes give one temperature. No diurnal range.
- `swrad` is the row's `pet` through Priestley-Taylor at a fixed 20 C to a net radiation, through
  an albedo of 0.23, in Langleys. Only `snowcomp` reads it. It is the same mock as `models/summa`.

The proposal suggested daily minimum and maximum temperature and internal PET and radiation
modules. The suite supplies `pet`, so using it is the more direct reading. `ddsolrad` is measured
below.

## Withdrawals

`abstr` goes through PRMS's own water-use module, as `gwr_transfer` from the groundwater reservoir
to an external location, in cfs per step. It is honoured, and it is capped.

`gwflow.f90:512` stops the run with `ERROR, not enough storage for transfer` rather than overdraw
the reservoir. So the adapter runs PRMS once, reads the recharge (a withdrawal from below cannot
change it), replays `gwflow`'s order of operations to find what the reservoir can give each day, and
runs PRMS again with that. What it could not supply is neither removed nor declared, and `run.json`
reports prescribed and removed volumes. The replay only looks backwards in time.

The water-use module has no sink on the soil zone or the canopy; those two only take gains. Its
other sources are stream segments, depression storage and lakes. The stream was not added because it
would change little, see the human-abstraction failure below.

## The step

**PRMS is daily and has no other step.** `prms_time.f90:115` stops with `ERROR, timestep < daily for
daily model`. The manifest declares PT1D, as the proposal does.

The resolution probe runs every model at its next finer step regardless, so hourly rows are fed one
per PRMS "day": depths per row in, fractions as 1-(1-k)^dt, amounts and first-order rates × dt,
`K_coef` / dt. The dates PRMS sees then advance a day per row, and `snowcomp`'s compiled half-day
energy terms stay per day. Transpiration does not read the date here (`transp_end` is 13). The effect
of the date on snow in that case was not measured, so the PT1H result is a weaker statement than
the daily ones, which is why the manifest does not declare it.

The proposal asked for a 30-day window; the manifest uses `window_days: full`, because a decade costs
0.1 s and the 30-day default would score less of the record.

## Evaluation

Gate seeds, all 37 probes, version `6.0.0-ht.3`: **20 PASS, 2 FAIL, 15 N/A** (12 INCOMPLETE,
3 INCOMPATIBLE). This is the full table:

| Probe | Verdict | Number |
| --- | --- | --- |
| `mass/catchment-closure` | PASS | residual 0.0000% of rain |
| `mass/snowpack-mass-closure` | PASS | residual 0.0001% of rain |
| `mass/extreme-event-closure` | PASS | worst event 5.6e-5 mm against 1e-3 mm |
| `mass/multi-decadal-drift` | PASS | -1.3e-4 mm over the last 1825 days |
| `mass/steady-state` | PASS | runoff 0.500 mm/day under 2.500 mm/day of rain |
| `mass/spinup-cycle-invariance` | PASS | worst departure 0.01% (`channel`) |
| `mass/dry-down` | PASS | 15 mm over 730 dry days |
| `mass/extreme-rain` | PASS | returns 1.00 of the 526 mm added |
| `mass/causality` | PASS | nothing moves before the storm; +0.998 of the rain after |
| `mass/antecedent-monotonicity` | PASS | 5.7 mm more from the 60 mm storm (0.10) |
| `mass/response-nonnegativity` | PASS | largest dip -1.8e-5 of the added 120 mm |
| `mass/runoff-bounds` | PASS | 0.48 of the rain, inside [0.11, 1.07] |
| `mass/precipitation-counterfactual` | PASS | +20%: evspsbl +0.201, mrro +0.772, storage +0.028 |
| `mass/warming-response` | PASS | warmer: mrro -0.306, evspsbl +0.312 per unit demand |
| `mass/phase-counterfactual` | PASS | 4.8% of the rain against a 5% limit |
| `mass/resolution-invariance` | PASS | 6.53% of the rain between PT1H and PT1D, worst of 3 seeds, limit 10% |
| `mass/area-invariance` | PASS | 0, by construction |
| `mass/time-origin-invariance` | PASS | 0 |
| `energy/pet-consistency` | PASS | 1.00 of potential wet, 0.26 dry |
| `momentum/routing-conservation` | PASS | channel non-negative |
| `mass/human-abstraction` | **FAIL** | 28.73% of the withdrawal unaccounted for, limit 5% |
| `mass/ungauged-basin-closure` | **FAIL** | `channel` reaches -0.000508 mm |

N/A (INCOMPLETE, a variable PRMS does not report or an input it does not consume): the six other
energy probes, `aquifer-recharge-ordering`, `groundwater-datum-invariance`,
`gw-sw-exchange-consistency`, `froude-regime`, `stage-discharge-monotonic`,
`uniform-flow-friction-consistency`. N/A (INCOMPATIBLE): `exchange-response` (`gwex` here is a
withdrawal, not head-driven, so `gwh` is not declared), `routing-lag-consistency`,
`wave-celerity-bounds`.

### `mass/human-abstraction`: the reservoir runs dry

Verified. Model property at the default parameters. On the worst seed the budget accounts for
271 mm of the 380 mm prescribed (`mrro` -270.1 mm, storage -0.8 mm); the three seeds leave 25.2%, 26.8%
and 28.7%. At `gwflow_coef` = 0.015 per day the withdrawal draws the reservoir to 0 in dry spells (`gw`
ranges 0.0 to 86 mm in the irrigated run), and PRMS has no withdrawal path from the soil, which still
holds water.

Replaying the natural run on that seed, over the whole record, gives this split of the 418 mm
prescribed: 287.8 mm from the reservoir, 8.2 mm more if the stream's remaining flow were taken too
(direct and interflow are small on those days), and 122 mm with no source. So the stream would move
the share by 2 points and not the verdict, and the adapter does not use it.

The failure is a statement about these parameters. Under Sagehen's, with `gwflow_coef` = 0.0037, the
reservoir stays above 0.97 mm and the probe passes at 0.44%.

### `mass/ungauged-basin-closure`: routing coefficients in single precision

Verified. Model property. Closure, the ET bound, variability and forcing fidelity all pass; only
`state_bounds` fails, on 5 of 12 seeds, with `channel` as low as -0.000508 mm.

`routing.f90:16` declares `C0`, `C1`, `C2` as `REAL`. At the defaults `K_coef` = 1 h, `x_coef` = 0.2
they are 0.23076925, 0.53846157, 0.23076925, which sum to 1 + 5.96e-8, so at steady flow Muskingum
releases 1 + 7.75e-8 of what enters. `basin_segment_storage` is the running sum of inflow minus
outflow (`muskingum.f90:358`), and it holds only about a twenty-fourth of a day's flow at K = 1 h, so
on a low-flow day the accumulated difference shows through. Over 6437 mm of throughflow that predicts
-0.000499 mm.

Replaying `muskingum.f90`'s hourly loop on the run's own inflow (seed 306340959):

| Coefficients | Channel minimum | Steps below zero | RMS against printed channel |
| --- | --- | --- | --- |
| PRMS's single precision | -0.000498 mm | 2413 of 4380 | 2.9e-7 in, the print quantum |
| Exact | +1e-13 mm | 0 | 1.2e-5 in |
| Printed by PRMS | -0.000508 mm | 2405 of 4380 | |

The amount is small and the bound is the probe's own, so this is reported as it stands. `channel` is
PRMS's own output and is not adjusted. At `K_coef` = 12 h the routing step is 12 h and the same ratios
round the other way (1 - 3.9e-8), so the running sum drifts up rather than down, the store reaches 17 mm,
and the probe passes. That is the sign of a rounding error, and it is not a reason to prefer that value.

## Passes that say little

- `mass/area-invariance`: trivial, because PRMS is given a fixed area.
- `momentum/routing-conservation`: at K = 1 h the channel holds about a twenty-fourth of a day's flow,
  so "within what a 15-day hydrograph can hold" is not a demanding bound here.
- `mass/time-origin-invariance`: the shift preserves month, day and day of year, so the date-reading
  parts of `snowcomp` see the same calendar.
- `mass/phase-counterfactual`: not trivial, but it passes at 4.8% against 5% and sits on the threshold,
  as the next section shows. A 10 C diurnal range moves it to 0.8%.

## How much each choice moves the result

A full evaluation of all 37 probes on the same gate seeds for each, at version `6.0.0-ht.3`, from
scratch copies of the tree that are not committed.

| Setting | PASS / FAIL | Verdicts that move | Other numbers |
| --- | --- | --- | --- |
| Evaluated | 20 / 2 | | abstraction 28.73%, resolution 6.53%, phase 4.8%, runoff 0.48 |
| **Sagehen parameters** | 20 / 2 | `human-abstraction` passes (0.44%); `antecedent-monotonicity` fails (+0.4 mm, 0.02 of the storm expected) | channel -0.000356 mm on 74 steps, resolution 4.72%, phase 5.0%, runoff 0.48 |
| `K_coef` = 12 h | 21 / 1 | `ungauged-basin-closure` passes | abstraction 28.73%, resolution 6.49% |
| Diurnal range 10 C | 20 / 2 | none | abstraction 43.89%, phase 0.8%, resolution 7.48%, runoff 0.49 |
| `ddsolrad` for `swrad` | 19 / 3 | `phase-counterfactual` fails at 5.1% | abstraction 28.64% |
| The case's area in PRMS | 19 / 3 | `area-invariance` fails, `gw` 3.2e-6 relative | |
| `soilzone_aet_flag` = 0 | 20 / 2 | none | abstraction 28.64%, runoff 0.49 |

The Sagehen set takes every value from upstream's `projects/sagehen/input/sagehen.params`, the median
over HRUs where it varies: `soil_rechr_max_frac` 0.568, `ssr2gw_rate` 0.0499, `ssr2gw_exp` 0.75,
`fastcoef_lin` 0.4, `pref_flow_den` 0.1, `smidx_coef` 0.0004, `gwflow_coef` 0.0037, `potet_sublim` 0.75,
`snowinfil_max` 2.75, `rad_trncf` 0.2331, `melt_force` 90.

So both failures, and two passes, depend on where in parameter space PRMS sits, and none of these
settings was chosen to change a verdict. The Sagehen `antecedent-monotonicity` failure is inferred to
come from `smidx_coef` being a twelfth of the default, which leaves the contributing area nearly
independent of antecedent soil moisture. That was not checked separately. Nor was the mechanism behind
the 43.89% with a diurnal range.

## Caveats

- Results describe PRMS 6.0.0 at its declared defaults, on one HRU, under mocked radiation and no
  diurnal range. They do not describe a calibrated NHM setup.
- PRMS's own `water_use` withdrawals are the only thing that removes `abstr`; the prescribed volume is
  not met where the reservoir is empty, and that is reported rather than forced.
- The proposal also asks how PRMS compares with SUMMA, SAC-SMA/Snow-17 and mHM. This package gives the
  verdicts and does not attempt that comparison.
- Hourly rows run on a calendar moving 24 times too fast, as above.

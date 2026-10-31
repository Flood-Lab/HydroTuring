# CFE under HydroTuring

CFE, the Conceptual Functional Equivalent of the National Water Model, written by
Fred Ogden and maintained by NOAA-OWP, packaged from
[NOAA-OWP/cfe](https://github.com/NOAA-OWP/cfe) at commit
`f9182dfb3d81c407c66c5b782114173893c30cf2`. Proposed by Xiaoyi Dong as
[Flood-Lab/HydroTuring#127](https://github.com/Flood-Lab/HydroTuring/issues/127).
Apache-2.0; the image carries the licence at `/model/LICENSE`.

## What runs

Upstream's shared library, built unmodified with its own CMake target
(`-DNGEN=ON`, `--target cfebmi`, the one option that needs no submodule).
`ht_cfe_driver.c` links against it and calls `initialize`, `set_value`, `update`,
`get_value` and upstream's `mass_balance_check`. It computes no hydrology. It reads from
upstream's `cfe_state_struct`, through upstream's header, what BMI does not export: the
two Nash cascade storages, the ET components, the infiltration and percolation terms and
the volume trackers. It writes one field, `time_step_size`, before `initialize`, because
`new_bmi_cfe` hard-codes 3600 and the config has no key for it; at the hourly step that
write is the same 3600. `ht_adapter.py` is the units layer, standard library only.

Base image `python:3.11-slim` pinned by digest (Debian 13); the source is pinned by commit
and the build refuses anything else. The evaluated image used cmake 3.31.6 and gcc 14.2 at
build time, is 45 MB, and needs no network at run time.

## What is reported

| Column | Native | Note |
| --- | --- | --- |
| `pr` | the forcing | echoed from the CSV as text |
| `evspsbl` | `ACTUAL_ET` | rain + retention depth + soil, over the row, in mm/day |
| `mrro` | `Q_OUT` | direct + lateral + baseflow released at the outlet, in mm/day |
| `dis` | `Q_OUT` over `area_km2` | CFE produces no other discharge; no probe that CFE can reach scores it |
| `mrso` | `SOIL_STORAGE` | absolute, mm |
| `gw` | `GW_STORAGE` | absolute, mm |
| `channel` | both Nash cascades | `nash_surface_params.nash_storage` + `nash_storage_subsurface`, from the struct, in mm |
| `snw`, `canopy` | identically zero | named under `zero_stores` in `run.json`; see below |

`mrro` is what leaves the routing, not the pre-routing `INFILTRATION_EXCESS`. That answers
the first question in the proposal review: the pre-routing term exists, but the water it
carries is in `channel` until it leaves, and counting it as runoff as well would count it
twice. The second question: the groundwater store and the routing stores are all real, and
are all reported, which goes beyond the proposal's `evspsbl`, `dis`, `mrso`.

## Zero stores and what is not reported

CFE has no snow module, all precipitation entering as rain and `tas` going unread, and no
interception store, `et_from_rainfall` taking PET straight from the row's rain with nothing held. So
`snw` and `canopy` are reported as identically zero and named under `zero_stores` in `run.json`. At
PT1H this reaches no probe: the closure family is daily, so it is N/A (INCOMPATIBLE) on the step
rather than on these columns. The proposal review expected `mrro` to open most of the mass suite. It
does not, for that reason.

There is no `gwex`: CFE's groundwater store exchanges nothing across the catchment boundary.
`abstr` is not read, because CFE has no withdrawal process.

## Parameters

No calibration exists for the synthetic catchment. Everything comes from upstream's
`configs/cfe_config_cat_87_pass.txt` at the pinned commit, the configuration the repository
ships for BMI forcing, except the soil depth, and two values that file leaves to its parser's
defaults.

| Key | Value | Unit | Meaning | Source |
| --- | --- | --- | --- | --- |
| `soil_params.depth` | capacity / 1000 / `smcmax` | m | soil depth | **from `soil_capacity_mm` in static.json**; upstream's 2.0 if absent |
| `soil_params.b` | 4.05 | - | Clapp-Hornberger exponent | upstream config |
| `soil_params.satdk` | 3.38e-6 | m/s | saturated conductivity | upstream config |
| `soil_params.satpsi` | 0.355 | m | saturated capillary head | upstream config |
| `soil_params.slop` | 0.01 | - | bottom-boundary gradient factor | upstream config |
| `soil_params.smcmax` | 0.439 | - | porosity | upstream config |
| `soil_params.wltsmc` | 0.066 | - | wilting point | upstream config |
| `soil_params.expon` | 1.0 | - | percolation outlet exponent | upstream config |
| `soil_params.expon_secondary` | 1.0 | - | lateral outlet exponent | upstream config |
| `max_gw_storage` | 0.25 | m | groundwater capacity | upstream config |
| `Cgw` | 1.8e-5 | m/h | groundwater outlet coefficient | upstream config |
| `expon` | 6.0 | - | groundwater exponential | upstream config |
| `gw_storage` | 0.125 | fraction | initial groundwater | upstream config |
| `alpha_fc` | 0.33 | - | field-capacity pressure fraction | upstream config |
| `soil_storage` | 0.585626 | fraction | initial soil water | upstream config |
| `K_nash_subsurface` | 0.03 | 1/h | lateral cascade constant | upstream config |
| `K_lf` | 0.01 | m/h | lateral outlet coefficient | upstream config |
| `nash_storage_subsurface` | 0, 0 | m | initial lateral cascade | upstream config |
| `surface_runoff_scheme` | NASH_CASCADE | | | upstream config |
| `N_nash_surface` | 2 | | surface cascade length | upstream config |
| `K_nash_surface` | 0.83089 | 1/h | surface cascade constant | upstream config |
| `nsubsteps_nash_surface` | 10 | | substeps per row | upstream config |
| `nash_storage_surface` | 0, 0 | m | initial surface cascade | upstream config |
| `surface_water_partitioning_scheme` | Schaake | | | upstream config, and fixed by the proposal |
| `Kinf_nash_surface` | 0.001 | 1/h | run-on infiltration | parser default (`bmi_cfe.c`) |
| `retention_depth_nash_surface` | 0.001 | m | retention in the first surface reservoir | parser default (`bmi_cfe.c`) |
| `refkdt` | 3.0 | - | Schaake constant | parser default (`bmi_cfe.c`) |

These are one catchment's example values, not a calibrated set for this one. The repository's
other example, `cfe_config_laramie_pass.txt`, differs in `b` (5.25) and `soil_params.expon` (6.0).

## The step

**The manifest declares PT1H only**, as the proposal does. Upstream runs at one step: 3600 s
is hard-coded with no config key, and `nash_cascade_surface` hard-codes `dt_h = 1.0`. Daily
probes are therefore N/A (INCOMPATIBLE), and that is most of them.

The resolution probe runs every model at the next finer step whatever the manifest says,
which for CFE is PT1M. So the adapter rescales upstream's hourly constants for any step other
than the hour, and at the hour leaves every value untouched:

| Constant | Rule |
| --- | --- |
| `Cgw`, `K_lf` (amounts per hour) | × dt in hours |
| `K_nash_subsurface` (fraction per hour) | 1 − (1 − k)^dt |
| `K_nash_surface`, `Kinf_nash_surface` (per hour against the hard-coded `dt_h`) | × dt in hours |
| `nsubsteps_nash_surface` | max(10, ⌈10 dt⌉), which keeps K/n below one |

`satdk`, the Schaake term and PET take the step from `time_step_size` themselves. Rain goes in
as the depth for the row, and fluxes come back as rates.

## Evaluation

Gate seeds, all 37 probes: **1 PASS, 0 FAIL, 36 N/A** (14 INCOMPLETE, 22 INCOMPATIBLE). The 22 are
the daily probes, together with the hourly probes that need inputs CFE does not consume; the 14 need
an output it does not have.

`mass/resolution-invariance` passes: `mrro` moves by 2.03% of the rain between PT1H and PT1M,
worst of three seeds, against a 10% limit.

A pass out of one scored probe says little about the model, and that is a limit of an hourly
model on a suite whose probes are nearly all daily, not a finding about CFE. The zero-store rule
changed no row's verdict here. The two sections below are the rest of what is known.

## CFE's mass accounting against what is reported

The adapter sums the reported columns and writes the result to `run.json` under
`water_budget` on every run, beside CFE's own `mass_balance_check` output.

On the harness's own cases:

| Case | P mm | ET mm | Q mm | ΔS mm | Residual mm | Steps off by > 1e-6 mm |
| --- | --- | --- | --- | --- | --- | --- |
| Resolution case, PT1H, 120 mm soil | 305.42 | 74.57 | 251.85 | −21.00 | −1.5e-13 | 0 of 960 |
| Same, PT1M | 305.42 | 74.46 | 252.28 | −21.32 | 2.8e-12 | 0 of 57600 |
| Same, PT1H, rain × 4 | 1221.67 | 74.84 | 1165.59 | −18.76 | −7.0e-13 | 0 of 960 |
| Closure case, PT1D (adapter's rescaling), 320 mm soil | 8937.51 | 4955.86 | 4013.49 | −31.84 | 3.3e-13 | 0 of 4015 |

CFE's global, precipitation, surface, soil and groundwater residuals are all below 1e-13 m on
these cases. Two differences between its accounting and these columns:

- **`volout` counts lateral flow when it enters the subsurface cascade**, where `mrro` counts
  it when it leaves. Measured: `volout − ET − Q` = 0.34327 mm, and the cascade held 0.34327 mm
  at the end of the run. The two agree once `channel` is included, which is why it is reported.
- **`mass_balance_check` leaves the subsurface cascade out** of its global term, for the same
  reason, and so closes without it.

### A water-making branch, not reached here

`cfe.c:191-201` reads `flux_perc_m`, which at that point still holds the *previous* row's
percolation. If it exceeds the current deficit, infiltration is set to the deficit and the
difference is added to the excess, whatever the rain was. On a rainless row that is water from
nowhere.

Verified with upstream's unmodified library. With the soil started full and no lateral outlet
(`soil_storage=1.0`, `K_lf=0`), five rainless rows each make exactly the previous row's
percolation, 0.12168 mm, 0.61 mm in all, and upstream's own check says
`WARNING: GLOBAL MASS BALANCE CHECK FAILED` with a residual of −6.084e-04 m.

**The benchmark configuration does not reach it.** It needs the soil within one row's
percolation of saturation, and `K_lf` drains the soil hard above field capacity: the most the
soil reached was 81.8 of 120 mm, and 185 of 320 mm. It did not fire at rain × 4, nor with
`satdk` raised to 1e-3 m/s, nor with both outlets off. It is worth telling upstream, not a
verdict about this configuration.

The branch above it, `cfe.c:163-170`, would overfill the soil, but it cannot fire with Schaake:
infiltration is `P·Ic/(P+Ic)`, below `Ic`, itself below the deficit. Checked by running with
both soil outlets off, where the soil peaked at 117.3 of 120 mm. Two more were read and not
exercised, so these are inferences: the groundwater overflow branch at `cfe.c:221-226` moves
the volume trackers without moving water (groundwater peaked at 31 of 250 mm), and the surface
cascade's `Q_out` is left stale when a last reservoir is exactly empty (`nash_cascade.c:79`),
which needs K/n ≥ 1.

## What daily probes would say (not archived)

This is a **diagnostic, not the manifest and not in the archive**: none of these rows are in
`result.csv`. Declaring `timestep: [PT1H, PT1D]` is one line, and the maintainer may want it, since
CFE's formulation has no structural objection to a daily row: Schaake is written for any step, and
`time_step_size` is a real field. But then every daily verdict rests on the adapter's rescaling, on a
step upstream does not offer, and the proposal claims only PT1H, so it stays undeclared. It was run
from an uncommitted copy of the tree, on the earlier adapter, before the zero-store rule, so the
closure family was N/A in it and its counts are not comparable with any model's current ones:

**11 PASS, 2 FAIL, 24 N/A** (22 INCOMPLETE, 2 INCOMPATIBLE).

| Probe | Result | Mechanism | |
| --- | --- | --- | --- |
| `momentum/routing-conservation` | FAIL, all 3 seeds, ratio up to 7.6 | The surface cascade's first reservoir keeps up to `retention_depth` (1 mm) out of routing. It stays in `channel` with no runoff to justify it, until ET takes it. With retention set to 0, no step exceeds the allowance. | verified on the first gate seed; the other two inferred |
| `mass/antecedent-monotonicity` | FAIL, 1 of 3 seeds, 0.0195 vs 0.02 | Schaake's infiltration barely moves with the deficit at a daily step: 33 mm more soil water before the storm adds 1.17 mm of direct runoff out of 60 mm. | verified |

Passing: `pet-consistency`, `area-invariance`, `causality`, `dry-down`, `extreme-rain`,
`phase-counterfactual`, `resolution-invariance`, `response-nonnegativity`, `runoff-bounds`,
`steady-state`, `warming-response`.

The routing result also shows a gap in the contract. Retention depth is ponded water, and there
is no store name for that. `channel` is the nearest one: CFE keeps it in the routing array, and
leaving it out would open the budget by exactly that water. So it stays there.

## How much each choice moves the result

On the resolution probe's hourly case, 40 days, first gate seed. "Step shift" is the runoff
moved between PT1H and PT1M over the whole case, in % of the rain, so it is smaller than the
probe's 2.03%, which is measured over the 7-day window.

| Choice | Q / P | ET / P | Step shift |
| --- | --- | --- | --- |
| **As packaged** | 0.825 | 0.244 | 0.14 |
| Upstream's 878 mm soil rather than 120 mm from the case | 0.890 | 0.254 | 0.05 |
| Initial soil 0.3 rather than 0.586 of capacity | 0.714 | 0.243 | 0.37 |
| Initial groundwater 0.5 rather than 0.125 of capacity | 1.056 | 0.244 | 0.10 |
| `K_lf` 0.001 rather than 0.01 | 0.787 | 0.251 | 0.83 |
| Xinanjiang with upstream's cat_87 values, not Schaake | 0.823 | 0.244 | 0.02 |

The initial groundwater matters most on a 40-day case because the store is slow; a Q/P above one
is that store releasing water it started with. On the 11-year closure case, run at PT1D with the
adapter's rescaling, the soil capacity matters most: Q/P is 0.45 at 320 mm and 0.33 at 878 mm. No
choice opens the budget.

GIUH routing was not measured: `ht_cfe_driver.c` reads the Nash arrays, which GIUH leaves null,
and would need the runoff queue added to `channel` first. The proposal does not ask for it.

## Caveats

- The soil capacity is met through depth, not porosity. Depth also sets the field-capacity
  threshold (70 mm at a 120 mm capacity, 185 mm at 320 mm) and the wilting point.
- The parameters are an example catchment's. No PASS or FAIL here is about CFE's skill.
- One PT1H case takes about 5 ms in the driver, and 0.15 s for 57,600 rows at PT1M. Repeat runs
  are bit-identical, and CFE draws no random numbers, so the request's seed is not used.

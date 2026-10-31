# VIC 5.1.0 under HydroTuring

VIC, the Variable Infiltration Capacity macroscale hydrologic model (Liang et al. 1994,
[doi:10.1029/94JD00483](https://doi.org/10.1029/94JD00483)), classic driver, from
[UW-Hydro/VIC](https://github.com/UW-Hydro/VIC) tag 5.1.0, commit
`14a371a8dbee7ea834152f8ac65e2d6d54a9c30a`. Proposed by Valmik Nahata as
[Flood-Lab/HydroTuring#60](https://github.com/Flood-Lab/HydroTuring/issues/60). VIC is MIT; the
sample parameters from [VIC_sample_data](https://github.com/UW-Hydro/VIC_sample_data) at
`3f2dcb4`, the commit the tag pins as its `samples/data` submodule, are GPL-2.0. The image carries
both licences at `/model/`.

**VIC 5.1.0 cannot be given a daily row, so every daily probe is N/A (INCOMPATIBLE) and only two
probes can score it.** The proposal lists PT1D. See [The step](#the-step).

## What runs

Upstream's `vic_classic.exe`, built from the tag with its own Makefile. The source is as tagged.
The build passes `LOG_LVL=30`, which the Makefile documents, and `CC="gcc -fcommon -O3"`: `-O3` is
the Makefile's own commented-out optimised line, and `-fcommon` is needed because VIC's headers
define globals (`LOG_DEST`, `funcd`) that GCC 10 and later refuse to link. The build checks the tag
still points at the commit and that it pins the sample data commit. Base image
`debian:bookworm-slim` by digest. The image is 44 MB.

`ht_adapter.py` writes the global parameter, soil, vegetation, constants and forcing files, runs VIC
once on one grid cell, and reads its ASCII output. It implements no hydrology. A case of 960 hourly
rows takes 0.07 s, and 57600 minute rows take 3 s.

Configuration: `FULL_ENERGY TRUE`, `QUICK_FLUX FALSE`, `FROZEN_SOIL FALSE`, ten thermal nodes, as in
upstream's `tests/science/global_param.classic.5.0_snotel.txt`. `BASEFLOW ARNO`, three soil layers,
no lakes, one snow band.

## The step

**Upstream does not support a daily row, and the reason is in three separate checks.** The docs
disagree with the code: `GlobalParam.md` says to set `MODEL_STEPS_PER_DAY` to 1 if `FULL_ENERGY` is
false, and `FAQ.md` lets the model and snow steps differ in water-balance mode. Daily
forcing was removed with MTCLIM in 5.0 ([GH#288](https://github.com/UW-Hydro/VIC/pull/288), in
`ReleaseNotes.md`), and `ForcingData.md` says the classic driver "requires subdaily forcings". The
proposal's daily step and minimum and maximum temperature describe VIC 4.

Each way of asking for a daily step was run on the unmodified executable, on 30 days of
`mass/catchment-closure`:

| `MODEL` | `SNOW` | `RUNOFF` | `FORCE` | `FULL_ENERGY` | Result |
| --- | --- | --- | --- | --- | --- |
| 1 | 1 | 4 | 1 | FALSE | exit 1, `get_global_param.c:738`: snow steps (1) < the minimum number of subdaily steps per day (4) |
| 1 | 4 | 4 | 1 | FALSE | exit 1, `get_global_param.c:979`: `FORCE_STEPS_PER_DAY must match SNOW_STEPS_PER_DAY` |
| 1 | 4 | 4 | 4 | FALSE | exit 1, `read_atmos_data.c:60`: with sub-daily forcing the model must run at the forcing's step (86400 s against 21600 s) |
| 1 | 24 | 24 | 24 | FALSE | exit 1, `read_atmos_data.c:60`, the same, against 3600 s |
| 1 | 4 | 4 | 4 | TRUE | exit 1, `read_atmos_data.c:60`, the same |
| 4 | 4 | 4 | 4 | FALSE | runs |
| 24 | 24 | 24 | 24 | FALSE | runs |

So the coarsest step the classic driver ran is six hours, in water-balance mode, and a daily row would have to be split
into at least four records before VIC could read it. `AGENTS.md` forbids that ("Do not aggregate
hourly rows to days inside the adapter, or split days into hours"), because the resolution probe
exists to measure the step, and because it would make the adapter invent a sub-daily weather.
Removing the check at `read_atmos_data.c:60` would be a change to upstream and would still need the
split. The image driver has the same `FORCE` = `SNOW` check (`vic_force.c:566`) and the same minimum
of 4 on the snow step; that is from reading it, and it was not run.

A daily model would need a different pin, such as a VIC 4 release, whose MTCLIM makes the split
internally. That is the maintainer's choice and was not attempted.

The manifest declares `[PT1H, PT15M, PT5M, PT1M]`, the steps in the suite that VIC accepts. The hour
is upstream's sample step. All four ran on the resolution case and closed, with the largest step
residual 1.1e-7 mm. VIC takes the step itself, so no rate constant is rescaled in the adapter:
`Dsmax` and the rest are per day in VIC's own input. Rows go in one to one.

What the daily probes say when a row is nevertheless held over sub-daily records is measured below, as
a diagnostic that is not for the archive.

## What is reported

| Column | Native | Note |
| --- | --- | --- |
| `pr` | the forcing | echoed from the CSV as text |
| `evspsbl` | `OUT_EVAP` | includes sublimation, and can be negative on net condensation |
| `mrro` | `OUT_RUNOFF + OUT_BASEFLOW` | |
| `hfls` | `OUT_LATENT + OUT_LATENT_SUB` | VIC reports the sublimation share separately |
| `hfss` | `OUT_SENSIBLE` | |
| `hfg` | `OUT_GRND_FLUX + OUT_DELTAH + OUT_FUSION` | VIC's flux leaves its thin surface layer at the first node; storage and phase change above it are added, as `AGENTS.md` asks |
| `rlus` | `OUT_LWDOWN - OUT_LWNET` | |
| `ts` | `OUT_SURF_TEMP` + 273.15 | |
| `mrso` | `OUT_SOIL_MOIST`, three layers summed | liquid and ice, absolute |
| `snw` | `OUT_SWE` | ground pack |
| `canopy` | `OUT_WDEW + OUT_SNOW_CANOPY` | |

## Is a store missing?

**No. VIC has no baseflow store and the classic driver has no routing store, so `mrso, snw, canopy` is
the whole of it.**

- **Baseflow.** `runoff.c` computes baseflow from the relative moisture of the bottom soil layer and
  takes it from that layer (`runoff.c:284-301`, ARNO form). That water is in `OUT_SOIL_MOIST`, so it
  is already in `mrso`. There is nothing below the third layer. `gw` is not reported, and is named
  under `not_reported`.
- **Routing.** The classic driver stops at grid-cell runoff, and routing is a separate offline tool
  (RVIC). `mrro` leaves when generated. `channel` and `dis` are not reported, so the probes that need
  them are N/A (INCOMPLETE). That is the correct result, and a zero `channel` would be false.
- **Pond.** `OUT_SURFSTOR` belongs to the lake module. It is identically 0 without lakes, and the
  adapter's own budget in `run.json` includes it.

Checked against the numbers: `P - E - Q - dS` over these three stores closes with the largest step
residual 1e-7 mm and a cumulative residual of 1e-7 mm over 305 mm of rain on the resolution case, and
3e-7 mm over 8937 mm in the 4015-day closure record held over records. Those are print rounding at ten
significant figures. A missing store would open the budget by what it holds. VIC's own
`OUT_WATER_ERROR` is at most 4e-14 mm.

`zero_stores` is `[canopy]` only on a case with `canopy_capacity_mm = 0`, which gets no vegetation
tiles (bare soil). Otherwise it is empty.

`sbl` is not reported, because VIC's vapour fluxes are net and go negative on deposition. `snm` is not
reported, because VIC has no output for liquid leaving the snow module with rain passed through.
`abstr` is not read: VIC has no withdrawal process. `pet` is not honoured as such, see below.

## Parameters

No calibration exists on the synthetic catchment. This is one documented point in parameter space:
upstream's Stehekin sample cell 86340 (`Stehekin_soil.txt`, `Stehekin_vegparam.txt`,
`Stehekin_veglib.txt`), with the case's capacities.

| Parameter | Value | Source |
| --- | --- | --- |
| `infilt`, `Ds`, `Dsmax`, `Ws`, `c` | 0.3, 0.1, 10.0 mm/day, 0.7, 2 | sample |
| `expt`, `Ksat`, `phi_s` | 12.676, 1257.1 mm/day, -99, all layers | sample |
| `bubble`, `quartz`, `bulk_density`, `soil_density` | 14.66, 0.693, 1569.1, 2619.1, all layers | sample |
| `Wcr_FRACT`, `Wpwp_FRACT` | 0.312, 0.184, all layers | sample |
| `rough`, `snow_rough` | 0.01, 0.03 m | sample |
| `dp`, `annual_prec`, `fs_active` | 4 m, 1318.5 mm, 0 | sample |
| layer depths | sample's 0.1, 0.3, 1.5 m × `soil_capacity_mm` / 761.7 | **from the case**; at 120 mm that is 0.016, 0.047, 0.236 m. 761.7 mm is the sample's porosity × depth |
| `init_moist` | sample's × the same scale | sample, scaled |
| `elev` | `elevation_m` | case |
| `lat` | `latitude_deg` | case |
| `off_gmt` | 0 | the record's clock is local |
| `avg_T` | mean `tas` over the first year of the record | case, no later rows |
| Vegetation | classes 1, 6, 7, 10, 11 at 0.643, 0.195, 0.122, 0.014, 0.027, monthly LAI from the sample | sample |
| `VEG_LAI_WATER_FACTOR`, `VEG_LAI_SNOW_MULTIPLIER` | scaled so the cell's peak canopy holds `canopy_capacity_mm` | **from the case** |
| `SNOW_MAX_SNOW_TEMP`, `SNOW_MIN_RAIN_TEMP` | `snow_threshold_degC` ± 0.5 | case, with VIC's default half-width |
| `EMISS_GRND`, `EMISS_VEG`, `EMISS_SNOW` | `eps`, where the case gives one | case; VIC's default is 1 |
| `WIND_H` | 10 m | sample |
| `RUNOFF_STEPS_PER_DAY`, `ATMOS_STEPS_PER_DAY` | max(24, model steps) | sample's hourly runoff; the driver's default atmos step |

## Forcing

VIC reads radiation, humidity, pressure and wind, and the probes give none of them except `rlds` on
some cases. These are mocked from the row alone, with no calendar, as in `models/summa`, and labelled
in `run.json`:

- `SWDOWN`: the row's `pet` through Priestley-Taylor at a fixed 20 C to a net radiation, through an
  albedo of 0.23.
- `LWDOWN`: the case's `rlds` where it gives one, otherwise the emission of a reference surface at air
  temperature, so no longwave deficit.
- `VP`: relative humidity 0.70 at the row's air temperature.
- `PRESSURE`: a standard atmosphere at the case's elevation.
- `WIND`: 2 m/s.

**VIC computes its own evaporative demand and never reads `pet`, except through `SWDOWN`.** So ET is
not bounded by the row's `pet`, and humidity moves it strongly (below). `rn` is not declared, because
VIC computes its own net radiation, so `energy/surface-energy-closure` is INCOMPATIBLE.

## Evaluation

Gate seeds, all 37 probes, version `5.1.0-14a371a.2`: **1 PASS, 1 FAIL, 35 N/A** (23 INCOMPATIBLE,
12 INCOMPLETE).

| Probe | Verdict | Number |
| --- | --- | --- |
| `mass/resolution-invariance` | PASS | `mrro` moves by 4.90%, 1.87%, 4.81% of the rain between PT1H and PT1M, limit 10%; `evspsbl` by 0.05% |
| `energy/radiation-consistency` | **FAIL** | 720 of 720 steps; worst residual -21.4 W m-2 against 2.1 allowed |

22 of the INCOMPATIBLE rows are daily probes. The other is `energy/surface-energy-closure`, for `rn`.

Two scored probes say little about VIC as a whole. That is the consequence of the step on a suite
that is nearly all daily, and it is no finding about the model.

### The ERROR, now resolved

The first evaluation had `mass/resolution-invariance` as FAIL (ERROR). The probe runs an hourly model at
PT1M. `ATMOS_STEPS_PER_DAY` was unset, the classic driver defaults it to 24, and 24 does not divide evenly
by 1440, so VIC stopped at `get_global_param.c:840`. That was an adapter fault. The adapter now sets it
to the model step below the hour, which leaves the hourly configuration unchanged.

### `energy/radiation-consistency`: VIC absorbs all of the sky

Verified. Model property. The criterion expects `rlus = eps sigma ts^4 + (1 - eps) rlds`. VIC's surface
has an emissivity but no longwave reflectivity: it emits at `EMISS_GRND` (`func_surf_energy_bal.c:632`)
and takes the whole downward flux in (`calc_surf_energy_bal.c:214`, `func_surf_energy_bal.c:636`).

Measured on one seed, `eps` 0.9814, 768 rows, bare soil:

| Candidate | Mean residual | Largest residual |
| --- | --- | --- |
| `rlus - eps sigma ts^4` | -0.027 W m-2 | 0.029 W m-2 |
| `rlus - eps sigma ts^4 - (1 - eps) rlds` | -7.03 W m-2 | 8.32 W m-2 |

The first is what VIC computes. The second is what the probe asks, and its mean equals the mean of
`(1 - eps) rlds`, 7.00 W m-2. The same holds for the report's first failing step: at `eps` 0.9519, `ts`
291.13 K and `rlds` 369, `eps sigma ts^4` is 387.75 against 387.72 reported, and `(1 - eps) rlds` is 17.75,
the gap to the 405.49 expected. That seed, 892092983, is also the worst, at 10.2 times the tolerance; the
other four reach 3.3 to 4.6 times.

The case has no canopy, so tiles are not involved, and `OUT_SURF_TEMP` and `OUT_RAD_TEMP` coincide there
to 1e-10 K. Computing `rlus` from `ts` and `eps` would pass and would report what VIC does not compute,
so the adapter reports `OUT_LWDOWN - OUT_LWNET`. Leaving emissivity at VIC's default of 1 would remove the
term and would make the entry not consume `eps`, which is INCOMPATIBLE.

## Diagnostic: daily rows held over 24 records

**Not for the archive.** These are full evaluations of all 37 probes on the same gate seeds, from scratch
copies of the tree with `timestep: [PT1D, PT1H]` and `VIC_HT_DAILY_RECORDS=24`: each daily row's values
are held for 24 records, precipitation split evenly, and the steps summed back. That is the
disaggregation this benchmark forbids, with a flat day, so it describes VIC under an invented
sub-daily weather. It also answers the proposal's own experiment.

| | PASS | FAIL | N/A |
| --- | --- | --- | --- |
| `FULL_ENERGY TRUE` | 15 | 8 | 14 |
| `FULL_ENERGY FALSE`, water outputs only | 15 | 5 | 17 |
| `TRUE`, alternative soil | 14 | 9 | 14 |

| Probe | Energy on | Water balance | Mechanism |
| --- | --- | --- | --- |
| `mass/catchment-closure` | FAIL | FAIL | Closure itself is 0.0000%. `et_plausible`: ET reaches -0.21 mm/day. Net condensation, inferred |
| `mass/multi-decadal-drift` | FAIL | FAIL | Closure and drift pass. `et_plausible`: ET is 1.47 and 1.45 of `pet`, which VIC does not read |
| `mass/ungauged-basin-closure` | FAIL | FAIL | Closure and bounds pass. `et_plausible`: 1.18 and 1.17 of `pet`, the same |
| `mass/human-abstraction` | FAIL | FAIL | 100% of the withdrawal: VIC has no withdrawal process and `abstr` is unread |
| `mass/steady-state` | FAIL | FAIL | `mrro` varies by 158% under constant weather. Inferred: the sample's monthly LAI, which keeps the run periodic. Not tested separately |
| `energy/latent-heat-et-consistency` | FAIL | N/A | `energy_closure` 15.1%, `flux_identity` on 399 of 3650 steps, where a pack is or could be present. VIC computes its own net radiation, so the applicability of `rn` here is open ([#160](https://github.com/Flood-Lab/HydroTuring/issues/160)). Not diagnosed further |
| `energy/evaporative-partition` | FAIL | N/A | `partition_shift` sum -331 against 5% of 2318. Same caveat |
| `energy/radiation-consistency` | FAIL | N/A | As above |

Both modes pass the other 15: `pet-consistency`, `antecedent-monotonicity`, `area-invariance` (trivial, VIC
never sees area), `causality`, `dry-down`, `extreme-event-closure`, `extreme-rain`, `phase-counterfactual`,
`precipitation-counterfactual`, `response-nonnegativity`, `runoff-bounds`, `spinup-cycle-invariance`,
`time-origin-invariance`, `warming-response`, `resolution-invariance`.

**The proposal's question.** Switching the energy balance on changes no verdict on any water probe here,
and very little of the water: on the closure case ET is 4763.0 and 4763.6 mm, runoff 4141.9 and 4141.2 mm,
and `steady-state` `mrro` varies 158.3% and 158.4%. The water budget closes in both. Two caveats on the
review's request for evidence that the flag isolates the energy solution: the configurations differ in
thermal nodes as well (10 against 3) and `QUICK_FLUX` is set only with energy on, and this is under flat,
mocked days, where a diurnal cycle cannot matter. It was not established from VIC's documentation, only
measured here.

The hourly water-balance evaluation scores one probe: `resolution-invariance` passes at 4.89%, against 4.90%
with energy on.

## How much each choice moves the result

Full evaluation under the alternative soil, which takes `infilt` 0.01, `Ds` 1e-4, `Dsmax` 3.05, `Ws` 0.93 from
the worked example in upstream's `docs/Documentation/Drivers/Classic/SoilParam.md`:

| Setting | Scored verdicts |
| --- | --- |
| Evaluated | resolution PASS 4.90%, radiation FAIL |
| **Alternative soil, hourly** | **resolution FAIL at 11.1%**, 6.8%, 8.0% on the three seeds, radiation FAIL |
| Water balance, hourly | resolution PASS 4.89%, radiation N/A |

So the resolution verdict depends on where in parameter space VIC sits. The cause of the larger movement
under the alternative soil was not diagnosed. Held over records, the alternative soil also fails
`antecedent-monotonicity` (+0.2 mm, 0.02 of the storm expected), and `ungauged-basin-closure` gains a runoff
ratio of 0.0064, outside [0.02, 0.98].

Targeted runs, not full evaluations. On the resolution case at PT1H, rain 305.4 mm:

| Setting | ET mm | Runoff mm | Storage change mm |
| --- | --- | --- | --- |
| Evaluated | 60.3 | 200.3 | 44.9 |
| Relative humidity 0.95 | 32.2 | 225.2 | 48.0 |
| Relative humidity 0.50 | 78.7 | 184.6 | 42.1 |
| Wind 4 m/s | 61.7 | 198.8 | 45.0 |
| Albedo 0.12 | 57.6 | 202.5 | 45.3 |
| Alternative soil | 60.1 | 179.0 | 66.3 |
| One grass tile | 35.2 | 223.0 | 47.3 |
| Water balance | 61.8 | 199.0 | 44.7 |

On the 4015-day closure case held over 24 records, rain 8937.5 mm:

| Setting | ET mm | Runoff mm |
| --- | --- | --- |
| Energy on | 4763.0 | 4141.9 |
| Water balance | 4763.6 | 4141.2 |
| Water balance, 4 records a day | 4547.8 | 4352.3 |
| Relative humidity 0.95 | 3331.2 | 5542.6 |
| Wind 4 m/s | 5065.9 | 3845.6 |
| Albedo 0.12 | 4694.6 | 4209.5 |
| Alternative soil | 6216.0 | 2643.0 |
| One grass tile | 4636.4 | 4268.8 |

Humidity moves ET by far the most of the mocks, because VIC builds its own demand from it. The record count
matters too, which is one more reason the split is not offered as a score: at 4 records a day runoff is 5%
higher than at 24.

## Caveats

- Two scored probes. Do not read the single PASS as a statement about VIC's conservation; the closure
  evidence is in the budget above and in the diagnostic, not in a verdict.
- Layer depths shrink with the case's capacity, to 1.6 cm for the top layer at 120 mm, while `dp` stays at 4 m.
  That is the capacity rule applied to a parameter set built for a much deeper column.
- Monthly LAI comes from the sample's parameters, so VIC does read the calendar through its own parameter file,
  though none of the mocks do.
- The mocks replace a real atmosphere. Radiation and humidity carry no diurnal cycle beyond what the row has.

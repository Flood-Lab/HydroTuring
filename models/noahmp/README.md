# Noah-MP under HydroTuring

Noah-MP (Niu et al. 2011, *JGR Atmospheres*,
[10.1029/2010JD015139](https://doi.org/10.1029/2010JD015139)), the land scheme of the
National Water Model, WRF and NASA LIS, in the refactored Fortran of He et al. Proposed by
Valmik Nahata as [Flood-Lab/HydroTuring#59](https://github.com/Flood-Lab/HydroTuring/issues/59).
UCAR open-source licence, for Noah-MP and for HRLDAS; the image carries both at
`/model/LICENSE.noahmp.txt` and `/model/LICENSE.hrldas.txt`, and is never pushed to a registry.

**Result, version `5.2.1-17751dc.3`: PASS on the 3 probes that could score it, N/A on the other
34.** Those 3 are all hourly, one of them is trivial for this model (see below), and none of them
is a water-budget closure. Noah-MP cannot be run at one step per daily row, so the closure family
and most of the suite cannot be put to it. Read the PASS with that in mind.

## What runs

| Piece | Pin |
| --- | --- |
| NCAR/noahmp | v5.2.1, `17751dcd7a2442a3c57138f530f523daf08dc1c5` |
| NCAR/hrldas, whose makefiles build it and which carries that commit as its submodule | `cd96df470220f7d7133cdbccd5f9c5355cf173e2` |
| Base image, both stages | `debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171` |
| Packages the build resolved from the bookworm archive | gfortran-12 12.2.0-14+deb12u1, libnetcdf-dev 1:4.9.0-3+b1, libnetcdff-dev 4.6.0+really4.5.4+ds-3, python3 3.11.2-1+b1, python3-numpy 1:1.24.2-1+deb12u1, python3-netcdf4 1.6.2-1+b1 |

The apt packages are not pinned individually; the table records what the evaluated image got. The
build checks both commits and refuses anything else. Upstream's source is compiled unmodified, in
single precision, with HRLDAS's serial gfortran options plus `-O2`. The image is 84 MB and uses no
network at run time.

`ht_noahmp_driver.F90` does what `hrldas.exe` does for one land point with the netCDF I/O removed:
`land_driver_ini`'s set-up calls, then for each row `land_driver_exe`'s forcing assignment and the
body of `NoahmpDriverMain`. It repeats those calls rather than calling `NoahmpDriverMain` so that
the column type is still in scope when `NoahmpMain` returns, which is the only way to read the
conductivity behind `hfg_bottom`. It computes no physics.

`verify_driver.py`, run inside the evaluated image, compares the driver with `hrldas.exe` over 3600
hourly steps of a winter into spring (peak SWE 222 and 237 mm), on croplands and on barren ground.
Across 16 series the worst difference is 3.8e-6 of the series' largest value. Thirteen series are
identical on every step on croplands and all but two on barren. The differences left are in `LH`,
the column's soil water and runoff, about the size of single-precision rounding in the sums.

## What is reported

| Column | Native | Note |
| --- | --- | --- |
| `pr` | the forcing | echoed as the text it arrived as |
| `evspsbl` | `EvapCanopyNet + Transpiration + EvapGroundNet` | net of dew and frost, so it can be negative |
| `mrro` | `RunoffSurface + RunoffSubsurface` | mm per step, converted to mm/day |
| `sbl` | `SublimSnowSfcIce + SublimCanopyIce` | gross sublimation |
| `snm` | `SnowBotOutflow`, the thin-snow pond terms, and rain with no snow layer | |
| `hfls` | `HeatLatentCanopy + HeatLatentGrd + HeatLatentTransp` | |
| `hfss` | `HeatSensibleSfc` | |
| `hfg` | `HeatGroundTot` | into the top of the snow-soil column, so under snow it is the snow surface |
| `rlus` | `RadLwEmitSfc` | downward longwave plus net longwave |
| `ts` | `TemperatureRadSfc` | Noah-MP inverts it from `RadLwEmitSfc` |
| `hfg_bottom` | `2 k1 (T1 - T2) / -z2` | Fourier flux through the base of the first soil layer, from the conductivity and end-of-step temperatures of the implicit solve; not a residual of the temperature change |
| `tsoil_layer` | `TemperatureSoilSnow(1)` | first soil layer, end of step |
| `mrso` | `SoilMoisture` over layer thickness | liquid and ice |
| `snw` | `SnowWaterEquiv` | |
| `canopy` | `CanopyLiqWater + CanopyIce` | |

Nothing is written as an invented zero, and `zero_stores` is empty: Noah-MP has both a snowpack and
a canopy. On a case whose `canopy_capacity_mm` is zero the column runs as bare ground (below), so
`canopy` comes out zero from the model's own arithmetic.

### Stores checked for the budget

Upstream's `BalanceWaterCheck` differences canopy, snow, soil, `WaterStorageAquifer` and
`WaterStorageWetland`. The question put to this package was whether a store was missing.

- **Aquifer: not reported, and correct under this configuration.** Upstream's example selects free
  drainage (`SUBSURFACE_RUNOFF_OPTION = 3`), under which `WaterStorageAquifer` is a constant 4900 mm
  (verified: minimum and maximum over the case are both 4900.0). A constant store reported as `gw`
  would be a store that does nothing, so it is named under `not_reported`. Under runoff options 1,
  upstream's groundwater scheme, the aquifer does move; see the alternative configuration below.
- **Ponding: not a store.** `PondSfcThinSnwMelt`, `PondSfcThinSnwComb` and `PondSfcThinSnwTrans` are
  reset every step and go into `SoilSfcInflow` in the same step. `WaterHeadSfc` only holds water
  under `WRF_HYDRO`, which this build does not define.
- **Wetland: off.** `WETLAND_OPTION = 0`.

With `mrso + snw + canopy` the budget closes without the aquifer: on a 40-day hourly case the worst
step residual is 1.1e-4 mm and the cumulative one -0.0016 mm of 305 mm of rain.

Also not reported: `dis` and `channel`, because Noah-MP has no routing. `gwex` is absent because
nothing crosses the boundary except drainage, which is in `mrro`. `abstr` is not read. `lwsnl` and
`csnow` could be built from Noah-MP's layers but are not emitted: the only probe asking for them
also needs a daily step and a prescribed `rn`, so it could not score the model anyway.

## Why PT1D is not declared

`NOAH_TIMESTEP` is a real setting, and upstream's example uses 1800 s under 3600 s forcing. The
adapter runs one Noah-MP step per row, so a daily row means `NOAH_TIMESTEP = 86400`. That does not
work, for two independent reasons found in the model, neither of them in the adapter.

**1. Upstream's water-budget check stops the run, within weeks, on all 8 records tried (verified).**
`SNOW_RET_FAC` is 5.0e-5 per second, and `SnowpackHydrology` releases bottom-layer liquid explicitly
as `SnowLiqReleaseFac * MainTimeStep * excess`. At 86400 s that factor is 4.32, so the release
overdraws the layer. In an instrumented prototype the summed snow liquid starts the failing step
negative (-5.16 mm on the first one), and `SoilSnowWaterPhaseChange`'s `max(0, total - ice)` raises
the pack's mass by 4.82 mm, against the 4.819 mm Noah-MP reports. Every error is a gain, and all 43
over the four records traced occur with snow on the ground. Lowering `SNOW_RET_FAC` to 1.0e-5, purely as a test, removes
them at every step tried: the worst error at 86400 s falls from 12.4 mm to 0.001 mm.

**2. With that check out of the way, the top soil layer's temperature blows up (existence verified,
explanation inferred).** It alternates between implausible values until `EnergyMain` stops on
negative emitted longwave, on days 81, 123, 280 and 293 of the four records traced. On the two also
run with the lower `SNOW_RET_FAC` it stops on the same days, 280 and 293. The likely cause is the
explicit coupling between the surface balance and the soil. A rough estimate of the gain,
h dt / (C dz), with h about 9 W m-2 K-1 and a 0.1 m layer, is about 3 at 86400 s, past the limit of 2
for stability, and about 1.6 at 43200 s, where the ladder stays stable; the estimate itself was not
tested separately.

The ladder below holds each daily row for k equal Noah-MP steps, with the budget STOPs removed, over
the first 400 days of `mass/catchment-closure` seed index 2. It is a prototype, not the adapter.

| Step | Days reached | Budget errors over 0.1 mm | Largest error | Top soil temperature |
| --- | --- | --- | --- | --- |
| 3600 s to 21600 s | 400 | 0 | 0.002 mm | 266 to 301 K |
| 28800 s | 400 | 7 | 1.25 mm | 266 to 301 K |
| 43200 s | 400 | 19 | 10.1 mm | 265 to 302 K |
| 86400 s | 280 | 10 | 12.4 mm | 169 to 450 K, then stops |
| 86400 s, `SNOW_RET_FAC` 1.0e-5 | 280 | 0 | 0.001 mm | 169 to 450 K, then stops |

`mass/warming-response` seed index 0 behaves the same way, reaching about -15 to 614 K before it
stops. Eleven single changes of physics option (compaction, time scheme twice, frozen soil,
supercooled water, snow conductivity, runoff, precipitation partition, snow cover, drag, surface
resistance) only moved the date of the budget stop, on six records: none of the 66 runs finished, and
the latest stop was row 317. Changing a table parameter to get past the stop would be
tuning, and would still not get past the second failure.

Holding a daily row for several hourly Noah-MP steps does run. It is not done here, because AGENTS.md
rules out splitting days into hours to reach the step a model prefers. HRLDAS's own route to a Noah
step shorter than the forcing step also reads the next forcing record to interpolate
(`READFORC_HRLDAS`), which would break causality, so a held row would not be upstream's behaviour
either. Whether a model's own step setting should count as an exception is a question for the
maintainers. The evidence for it is under "If daily rows were held", clearly not as a verdict.

## Parameters

No calibration exists on the synthetic catchment. These results describe one documented point in
parameter space.

| Quantity | Value | Source |
| --- | --- | --- |
| Physics options | every option in `namelist.hrldas_example` unedited; among them dynamic vegetation 4, runoff 3/3, snow compaction 2, time scheme 3 | upstream's example |
| `NOAH_TIMESTEP`, `FORCING_TIMESTEP` | the row's step | one step per row |
| `SOIL_TIMESTEP` | 0, the main step | upstream's example |
| All other parameters | `NoahmpTable.TBL`, unedited | upstream |
| Land class | MODIS 12, croplands | upstream's `bondville.dat` |
| Soil class | STAS 8, silty clay loam, porosity 0.464 | upstream's `bondville.dat` |
| Maximum and minimum vegetation | 96 and 1 percent | upstream's `bondville.dat` |
| Initial soil moisture | 0.298, 0.294, 0.271, 0.307 | upstream's `bondville.dat`, unless the case sets one |
| Initial LAI | 2.0 | upstream's `bondville.dat` |
| Layers | 0.10, 0.30, 0.60, 1.00 m | upstream's example, reshaped as below |
| Column depth | `soil_capacity_mm` / porosity | **static.json**; top layer stays 0.1 m where it can, the rest shared 0.3 : 0.6 : 1.0 |
| `CH2OP` | `canopy_capacity_mm` / class's largest monthly LAI + SAI | **static.json** |
| Bare ground, MODIS 16, no vegetation | when `canopy_capacity_mm` is 0 | **static.json** |
| `EG`, soil emissivity | `eps` | **static.json**, table's 0.97 otherwise |
| First layer, `CSOIL`, initial column | `soil_layer_depth_m`, `soil_heat_capacity_areal` through Noah-MP's mixture rule, `soil_temperature_initial`, `soil_water_content_initial` | **static.json** |
| Latitude | `latitude_deg`, for the phenology's hemisphere; Bondville's 40.01 otherwise | **static.json** |
| Deep and initial soil temperature | `soil_temperature_initial`, or the first year's mean air temperature | **static.json**, or the record's first year and nothing later |

The soil class is fixed, so a case whose `soil_porosity` is not 0.464 stops the adapter rather than
running on a porosity it did not ask for; no case the model can be scored on does this. A small
capacity gives an odd column: `soil_capacity_mm` of 120 gives layers 0.1, 0.025, 0.050 and
0.083 m, a top layer thicker than those below it.

Noah-MP's mixture rule makes the first layer's capacity move with water content, so on the soil heat
probe it spreads 0.2 to 0.8 percent around the case's value.

### Mock forcing

Built from each row alone, with no calendar, and labelled in `run.json` under `mock_inputs`.
`rsds`, `rlds`, `sfcWind`, `huss` and `ps` go in directly wherever a case supplies them.

| Input | Mock | Source |
| --- | --- | --- |
| Net radiation target | the probe's `rn`; otherwise `pet` times the Priestley-Taylor (alpha 1.26) conversion at a fixed 20 C | `models/summa` |
| Shortwave, longwave | split so a reference surface at air temperature, albedo 0.23 and Noah-MP's soil emissivity has that net radiation | `models/summa` |
| Humidity, wind, pressure | 70 percent, 2 m/s at 10 m, standard atmosphere at `elevation_m` | `models/summa` |
| Cosine of the zenith angle | 0.5 where the row has shortwave, 0 where it has none | this package: HRLDAS computes it from the clock, and Noah-MP discards shortwave at zero |

A row labelled t forces the step ending at t + dt, as HRLDAS's forcing stamps do. The clock reaches
Noah-MP only through the phenology's monthly LAI and SAI.

`pet` is declared under `uses_forcing`, but it is only read when a case supplies neither `rn` nor
both `rsds` and `rlds`. No probe that can score this model depends on it. A daily step, if ever
allowed, would make this worth revisiting.

### Prescribed `rn`, and the maintainer's question on #59

`rn` is **not declared**, so `energy/surface-energy-closure` is N/A (INCOMPATIBLE). Declaring it would
mean accepting the supplied net radiation as the energy input to the reported budget, and Noah-MP
does not: it computes net radiation from its own albedo, emissivity and surface temperature, and has
no prescribed-radiation mode. The adapter uses `rn` to build a sky and nothing holds the budget to it.

As asked on the issue, both are kept: `run.json` carries `mean_prescribed_rn_w_m2`,
`mean_noahmp_net_radiation_w_m2` and `mean_abs_prescribed_rn_minus_noahmp_net_radiation_w_m2` beside
the budget terms, and the developer switch `NOAHMP_HT_DIAG=1` writes `rn_prescribed`, `rn_target` and
`rn_noahmp` for every row. `rsds` and `rlds` are used directly where a probe supplies them.

Applying that probe's criterion by hand, on all 5 of its seeds, separates the two possibilities the
maintainer asked about. This is a check and not a verdict.

| Compared with `hfls + hfss + hfg` | Blocks failing, of 28 | Worst block |
| --- | --- | --- |
| the probe's `rn` | 20 to 24 | 6.7 to 8.5 W m-2 against 2 allowed, at night |
| Noah-MP's own net radiation | 0 | 0.40 to 1.47 W m-2 |
| Noah-MP's own, with precipitation-advected heat credited | 0 | 3e-6 W m-2 |

So the gap is an input mismatch and the energy budget itself closes. Mean `rn` is 76 to 82 W m-2 and
Noah-MP's own is 70 to 76, with an albedo near 0.14 where the reference surface assumes 0.23. Noah-MP's
own budget also carries `HeatPrecipAdvSfc` and `HeatCanStorageChg`, which have no contract column;
`run.json` reports their means. On a vegetated surface they account for the residual left against
Noah-MP's own net radiation, by upstream's `EnergyBalanceError`, which stays below 1e-4 W m-2.

The same runs, on a bare snow-free surface, show `hfls / evspsbl` fixed at 2.5104e6 J/kg, 2.75 percent
above a temperature-dependent latent heat at these temperatures.

## Verdicts

Run through the harness on the gate seeds, with `window_days: full`, on image digest `e554f1152eaa`.

| Probe | Verdict | Number | Mechanism |
| --- | --- | --- | --- |
| `energy/radiation-consistency` | PASS | worst step 0.01 of tolerance (-0.035 W m-2 against 2.66) | **Trivial.** `EnergyMain` computes `TemperatureRadSfc` by inverting this identity from `RadLwEmitSfc` (verified in source). Only two things are tested: that `eps` reaches `EmissivitySfc`, and the time labelling |
| `energy/soil-heat-storage-consistency` | PASS | worst phase 0.035 W m-2 against 1 | Implicit solver's own consistency: the flux uses its conductivity and end-of-step temperatures. The capacity was set from the case. Inferred |
| `mass/resolution-invariance` | PASS | worst `mrso_end` 1.857 percent of rain, limit 10 | Not trivial: hourly against minute, volumes move by about 1 percent of rain. Says nothing about a daily step, which is where Noah-MP fails |

**N/A (INCOMPATIBLE), 24.** Hourly step does not cover the probe's daily step: `energy/evaporative-partition`,
`energy/latent-heat-et-consistency`, `energy/pet-consistency`, `mass/antecedent-monotonicity`,
`mass/area-invariance`, `mass/catchment-closure`, `mass/causality`, `mass/dry-down`,
`mass/extreme-event-closure`, `mass/extreme-rain`, `mass/human-abstraction`,
`mass/multi-decadal-drift`, `mass/phase-counterfactual`, `mass/precipitation-counterfactual`,
`mass/response-nonnegativity`, `mass/runoff-bounds`, `mass/snowpack-mass-closure`,
`mass/spinup-cycle-invariance`, `mass/steady-state`, `mass/time-origin-invariance`,
`mass/ungauged-basin-closure`, `mass/warming-response`, `momentum/routing-lag-consistency`. `rn` not
consumed: `energy/surface-energy-closure`.

**N/A (INCOMPLETE), 10.** Variables Noah-MP does not report here: `energy/snowmelt-energy-water`,
`mass/aquifer-recharge-ordering`, `mass/exchange-response`, `mass/groundwater-datum-invariance`,
`mass/gw-sw-exchange-consistency`, `momentum/froude-regime`, `momentum/routing-conservation`,
`momentum/stage-discharge-monotonic`, `momentum/uniform-flow-friction-consistency`,
`momentum/wave-celerity-bounds`.

The proposal hoped to roughly double what the energy law can say. As packaged, it scores none of the
four energy probes it named, because three are daily and the fourth requires a prescribed `rn`.

An earlier version came back ERROR on `mass/resolution-invariance`: 40 days of minute rows with 16
columns at 12 significant digits wrote 9.85 MB against the probe's 8 MB limit. That was the adapter's
doing. Output is now 8 significant digits, about what single precision resolves, and 7.54 MB. Nine
digits still gave 8.12 MB. Keeping `window_days: full` was preferred to shrinking the month that probe
integrates over.

## What the choices move

Direct runs of the evaluated image on one seed of `mass/resolution-invariance` (212.5 mm of rain in
the scored rows), with the probe's own criteria applied by hand to the other two probes. Each row
changes one choice.

| Choice | `mrro` hourly, mm | `evspsbl` hourly, mm | `mrro`, hourly against minute, % of rain | `evspsbl`, same |
| --- | --- | --- | --- | --- |
| as evaluated | 155.5 | 60.7 | 0.94 | 1.14 |
| grassland, MODIS 10 | 166.1 | 49.9 | 0.95 | 0.12 |
| reference albedo 0.15 | 158.8 | 57.2 | 0.42 | 0.87 |
| relative humidity 0.5 | 146.6 | 70.9 | 2.29 | 2.03 |
| wind 4 m/s | 141.6 | 76.5 | 3.74 | 3.37 |
| cosine of zenith angle 0.3 | 157.5 | 58.0 | 0.65 | 1.03 |
| runoff options 1, with the aquifer | 158.1 | 62.5 | 0.10 | 0.05 |
| snow compaction option 1 | 155.5 | 60.7 | 0.94 | 1.14 |

Wind and humidity move the partition the most, up to 14 mm of runoff here. No verdict moves under any
row: all stay inside the 10 percent limit, and the radiation and soil heat probes pass under every
one. Snow compaction changes nothing because no snow forms on that record.

### Alternative configuration: upstream's groundwater scheme

The variant `noahmp_gw` ran runoff options 1 with `WaterStorageAquifer` reported as `gw`, through the
harness on all probes. It is not part of this package.

- `energy/radiation-consistency`: PASS. `energy/soil-heat-storage-consistency`: PASS.
- `mass/resolution-invariance`: PASS, worst `mrro` 1.061 percent of rain. This variant has one more
  column, so it wrote 7 significant digits to fit the 8 MB limit.
- The other 34 probes stay N/A, as they do for the evaluated package: 3 PASS and no verdict moved.
- The aquifer moves between 4898.3 and 5000.0 mm on the 40-day case, and the budget with it included
  closes to -0.089 mm of 305 mm. Nearly all of that is probably the first step, because the variant
  took the aquifer's initial value from the end of step 1 (inferred).

Free drainage stays the evaluated configuration because it is upstream's example. Reporting the
aquifer under runoff options 1 would add a `gw` store but no additional scored probe: the groundwater
exchange probes also need river-exchange terms and forcing this model does not take.

### If daily rows were held

This is **not an evaluation of the package and these are not verdicts.** The variant `noahmp_held`
broke the contract on purpose, holding each daily row for 24 hourly Noah-MP steps, flat within the
day, to measure what the rule costs and what a maintainer's exception would open. The 37 rows came out
14 PASS, 9 FAIL, 2 ERROR, 2 INCOMPATIBLE and 10 INCOMPLETE.

- PASS: `energy/pet-consistency`, `energy/radiation-consistency`, `energy/soil-heat-storage-consistency`,
  `mass/area-invariance`, `mass/causality`, `mass/extreme-rain`, `mass/phase-counterfactual`,
  `mass/resolution-invariance`, `mass/response-nonnegativity`, `mass/runoff-bounds`,
  `mass/snowpack-mass-closure`, `mass/spinup-cycle-invariance`, `mass/time-origin-invariance`,
  `mass/warming-response`.
- FAIL, none of them diagnosed, so none is a finding about Noah-MP: `energy/evaporative-partition`,
  `energy/latent-heat-et-consistency`, `mass/antecedent-monotonicity`, `mass/catchment-closure`,
  `mass/dry-down`, `mass/human-abstraction`, `mass/precipitation-counterfactual`, `mass/steady-state`,
  `mass/ungauged-basin-closure`.
- ERROR: `mass/extreme-event-closure` and `mass/multi-decadal-drift`. The binary files filled the
  container's 64 MiB `/tmp`: `forcing.bin` and the truncated `output.bin` add up to 67.0 MB in both
  cases (verified by the arithmetic). That is this variant's own limit. The evaluated adapter's largest
  case, 57600 rows, needs 26.7 MB.

If the maintainers allow this route, every FAIL above has to be diagnosed before it means anything.

## Limits

- Three scored probes, one trivial. This package says little about whether Noah-MP conserves water at
  catchment scale, since the closure probes cannot be put to it. Noah-MP's own budget check, at 0.1 mm
  per step, is the only water-closure evidence here, along with the residuals quoted above.
- One documented point in parameter space, a single column, on mock radiation, humidity and wind.
- Single-precision physics, with rounding about 1e-4 mm per step in the budget.
- A row's forcing is flat over the step, and the zenith angle is mocked, so there is no solar geometry.
- Cost: 0.11 s of wall time for 960 hourly rows and 1.1 s for 57600 minute rows.
- `NOAHMP_HT_OPTIONS`, `NOAHMP_HT_VEG`, `NOAHMP_HT_SOIL`, `NOAHMP_HT_ALBEDO`, `NOAHMP_HT_RH`,
  `NOAHMP_HT_WIND`, `NOAHMP_HT_COSZ` and `NOAHMP_HT_DIAG` are developer switches. None is set in the
  evaluated image.

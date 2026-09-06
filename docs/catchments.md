# Real catchments, generated weather

*Proposal, 2026-09-05. Data and a prototype generator are in the repository;
the probes have not been moved onto them yet.*

## The problem this fixes

Every probe today runs on one invented catchment: seven ad-hoc numbers in
`static.json` (a 250 km² area, a 320 mm soil, a 2 mm canopy) and a weather
generator whose climate is nobody's. Two consequences showed up in the first
two evaluations:

- **The adapter decides the verdict.** A submitted model reads dozens of
  catchment attributes the probe does not provide, so each adapter invents
  them. On δHBV 2.0 one reasonable guess for a seasonality index moved the
  runoff's coefficient of variation from 1.1 to 0.07 and two probe verdicts
  with it. Whatever the model is doing, that number was ours.
- **Physics violation and distribution shift are indistinguishable.** The
  invented climate's PET seasonality sits outside anything in CONUS. A model
  trained on real basins that misbehaves there may be wrong about physics or
  may simply be far from its training data, and the probe cannot tell.

The fix is to make the catchment real and keep only the weather synthetic.

## Four catchments

Each is a gauged CAMELS-US basin with its full HydroATLAS row (196
attributes: monthly precipitation, temperature, potential and actual
evaporation, snow cover and soil water; soils, lithology, land cover,
terrain, human footprint) plus Caravan's climate indices, area and
coordinates, copied verbatim from Caravan (Kratzert et al. 2023). Nothing in
`catchments/*.json` is invented.

| Id | Basin | Area | P | PET | PET/P | Snow share | Runoff | What it exercises |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `alpine-snow` | South Fork Payette River at Lowman, ID (`camels_13235000`) | 1156 km², 2090 m | 548 mm | 878 mm | 1.6 | 0.69 | 265 mm | a snowpack that holds half the year's water and releases it in a spring pulse; a summer with no supply and full demand |
| `maritime-rain` | Tilton River near Cinebar, WA (`camels_14236200`) | 363 km², 669 m | 2102 mm | 825 mm | 0.4 | 0.00 | 1322 mm | energy-limited: winter-peaked rain on full forest, runoff two thirds of P, storms of 300 mm |
| `humid-karst` | North Fork River near Tecumseh, MO (`camels_07057500`) | 1456 km², 333 m | 1090 mm | 1253 mm | 1.15 | 0.00 | 414 mm | 100 % karst: a large, slow groundwater store; summer demand above supply; the dry-down and steady-state probes' hardest case |
| `lowland-agricultural` | Fish Creek near Crystal, MI (`camels_04115265`) | 96 km², 269 m | 833 mm | 962 mm | 1.15 | 0.19 | 290 mm | flat, 45 % crops, 24 % lakes, water table at 75 cm: small runoff with a seasonal snowpack and a strong human imprint |

A fifth type is wanted and could not be filled from the data available on
this machine: a **semi-arid** basin with PET/P above 2 and ephemeral flow
(a CAMELS-US Arizona or Texas basin, or a CAMELS-BR or CAMELS-AUS one for a
tropical wet-dry regime). The Caravan archive is 12.5 GB; from a local copy,
`scripts/catchment_from_caravan.py <attributes dir> <gauge_id> <name>
"<description>"` writes the file in one command. The eight basins shipped
in the Caravan tutorial data, which is where these four come from, are all
humid or snow-fed.

## Weather that matches the catchment

`src/hydroturing/weather.py` generates daily precipitation, temperature and
potential evaporation from the catchment's own HydroATLAS monthly
climatology. The seed decides which days are wet, how much falls, and how
the temperature anomaly runs; the climate is not the seed's to change:

- **Precipitation.** A two-state wet/dry chain (storms persist) with a gamma
  depth whose expectation reproduces each month's HydroATLAS total.
- **Temperature.** Monthly means interpolated to the day, plus an AR(1)
  anomaly.
- **Potential evaporation.** Monthly totals interpolated to the day and
  scaled with the anomaly, so a warm spell raises demand as any
  temperature-based PET does.

Over 40 generated years the annual totals land within 1–3 % of HydroATLAS
and every monthly mean within 0.8 °C; monthly precipitation and PET within
16 % and 17 % of their targets, which is sampling noise on a heavy-tailed
process plus the mid-month interpolation, and can be tightened. Wet-day
frequency and storm persistence are not in HydroATLAS; they default to
values inside the range Caravan's `high_prec_freq`/`low_prec_freq` span for
these basins and can be set per catchment.

The consequence for a submitted model is that the attributes it reads and
the forcing it receives describe the same place. δHBV's `aridity`, `meanP`,
`ETPOT_Hargr`, `snowfall_fraction` and the seasonality indices are now
read, not derived; Google's 84 Caravan attributes are the row itself. The
adapter becomes a renaming, and the sensitivity table in
`models/dhbv2/README.md` becomes unnecessary.

## What changes in the contract

1. `static.json` becomes the flattened real row: `area_km2`,
   `latitude_deg`, `longitude_deg`, then every Caravan index as
   `caravan.<name>` and every HydroATLAS field as `hydroatlas.<name>`
   (`weather.static_attributes`). Published names, so an adapter maps by
   name and never guesses.
2. The bucket parameters the reference models need (soil capacity, canopy
   capacity, degree-day factor, baseflow coefficient, snow threshold) leave
   `static.json`. The reference bucket derives them from the HydroATLAS row
   in its own adapter, exactly as a submitted model must, and documents the
   mapping. The gate then proves the probe on the same inputs every model
   sees.
3. Probe bounds that today read `soil_capacity_mm` read a HydroATLAS-derived
   quantity instead: root-zone storage from `swc_pc_syr` and soil depth, or
   simply `hydroatlas.swc_pc_syr`-scaled capacity, stated in the probe.
   `state_bounds` then bounds a store against something the model was
   actually told.
4. `case.catchments: [alpine-snow, maritime-rain, humid-karst,
   lowland-agricultural]` in `probe.yaml`; every seed runs on every listed
   catchment and all must pass. A probe whose expectation needs a
   particular regime lists only those: warming-response on the two
   snow-free basins, so timing does not confound the sign; dry-down on all
   four, with karst as the slow case.
5. Windows and reports carry the catchment id; `result.csv` gains a column.

## What it costs

Four catchments multiply the model runs by four. δHBV runs the full record
in 4 s per case; Google takes 20 s on a 30-day window, so a seven-probe
evaluation goes from about 10 to about 40 minutes locally and the CI
timeout has to rise or the probe set has to be split per job. The gate
runs fourteen reference models on four catchments and stays under a minute.

## What is deliberately not proposed

Real observed weather. Feeding a real basin's real forcing would let a
model that has seen Caravan recognise the record, and the unseen-sample
guarantee is the benchmark's strongest property. The climatology is
public; the sequence is not.

# δHBV 2.0 under HydroTuring

The differentiable HBV behind the MHPI group's CONUS-scale water model
(Song et al. 2025, *WRR*, [10.1029/2024WR038928](https://doi.org/10.1029/2024WR038928)),
packaged from the NextGen module [mhpi/dhbv2](https://github.com/mhpi/dhbv2)
at v0.5.4 (`6d4bb68`) with the published daily weights `hbv_2_ep100.pt`
from the MHPI release bucket. Submitted as
[Flood-Lab/HydroTuring#2](https://github.com/Flood-Lab/HydroTuring/issues/2).

An LSTM writes three HBV parameters (`parBETA`, `parK0`, `parBETAET`) day
by day from the forcing and the catchment attributes; an MLP writes the
other thirteen, and the two unit-hydrograph parameters, from the
attributes alone; the HBV bucket model runs with them, three components in
parallel, averaged. It is the first submission that carries explicit
stores and an explicit evaporation, so it is the first that can be scored
on closure rather than left N/A (INCOMPLETE).

## Licence

Code and weights are under the PSU Non-Commercial Software License. This
evaluation is the non-commercial research use it permits. The image
carries the licence at `/model/LICENSE`; do not push it to a public
registry.

## What the adapter reports

| Column | What it is |
| --- | --- |
| `pr` | the forcing, echoed |
| `evspsbl` | HBV actual evaporation (`AET_hydro`) |
| `mrro` | routed streamflow (`streamflow`) |
| `dis` | the same over the catchment area, m3/s |
| `snw` | `SNOWPACK + MELTWATER`: the snowpack and the liquid water it holds |
| `mrso` | `SM`: the soil moisture box |
| `gw` | `SUZ + SLZ`: the upper and lower groundwater boxes |
| `channel` | cumulative unrouted minus cumulative routed flow: runoff generated but still inside the unit hydrograph |
| `canopy` | zero: HBV has no interception store |

States are averaged over the three HBV components, as the model averages
its fluxes. `canopy` is reported as identically zero because the model has
no such store; that is a statement about its structure, not a value
invented to satisfy a column, and it is recorded as such in `run.json`.

## Inputs

The model takes precipitation, temperature and potential evaporation,
which is what the probe generates, so no forcing has to be mocked. The
probe's `pet` goes in directly; the module's Hargreaves routine is not
used. The adapter never resamples: rows go in at the step of the case.

### Units at a step other than a day

The networks were trained to read forcing in mm per day and to write HBV
parameters in daily units. At an hourly step the adapter therefore does
three things, none of which changes anything at the daily step:

- the networks are fed rates in mm per day, as they were trained;
- the HBV core is fed depths per row, rate × step length, and its fluxes
  are divided back to rates;
- every parameter with a time in its units is put in the units of the
  step before the core uses it: `parK0`, `parK1`, `parK2`, `parC` (per-day
  fractions) become `1 − (1 − k)^dt`; `parPERC`, `parCFMAX`, `parRT`
  (per-day amounts) become `amount × dt`; the unit hydrograph keeps its
  shape in days by scaling its time constant by `1/dt` and its length from
  15 days to `15/dt` steps. `parFC`, `parUZL`, `parBETA`, `parLP`, `parTT`,
  `parCFR`, `parCWH`, `parBETAET`, `parAC` carry no time and are untouched.

Without the scaling the stores drained twenty-four times too fast and the
integrated runoff moved by 216–482 % of the rain between the hourly and
daily step; that number measured the adapter's units, not the model. With
it, under `.5`, the worst volume differs by 2.9, 11.3 and 12.1 % of the rain
on the resolution probe's three gate seeds. On the worst seed that is `mrro`,
84 mm hourly against 54 mm daily of 248 mm of rain, with evaporation at 1.7 %
and `mrso_end` at 0.5 %. What is left is the LSTM's own
recurrence, whose memory was learned with one step meaning one day. Measured
on the first gate seed before `.5` and not repeated since: holding the hourly
physics but feeding it the daily run's parameter series repeated 24× brings both volumes within 1 % of the daily run, and the
hourly LSTM writes a lower `parBETAET` (0.29 against 0.32 normalised) and
a higher `parBETA` (0.38 against 0.30), which is where the extra
evaporation comes from. That is the model's step dependence, and it is
what the probe now reports.

Of the 28 attributes the parameterisation network reads:

- **Derived from the first year of the record** (a climatology the model
  has seen, through which nothing later can reach back): `meanP`,
  `ETPOT_Hargr`, `aridity`, `meanTa`, and `snowfall_fraction` as the share
  of precipitation falling below the snow threshold, a sub-daily record being
  grouped into whole days for that one attribute, since `.5`.
- **From `static.json`**: `uparea` is the catchment area, which also
  drives the regional groundwater term (below). `meanelevation` comes from
  `elevation_m` when present, otherwise the training mean.
- **Zero**: `glaciers`, `permafrost`.
- **Training mean** (zero after standardisation): soil textures, porosity,
  permeability, slope, NDVI, free water, because a synthetic lumped
  catchment has no counterpart for them; and `seasonality_P`,
  `seasonality_PET`, `snow_fraction`, because the definitions the training
  set used could not be established. That last choice matters and is
  recorded here rather than hidden. A first version of the adapter derived
  the seasonalities as the amplitude of a fitted annual sine relative to
  the mean, which put the synthetic climate's PET seasonality at 1.07
  against a training range of 0.40 to 0.67, and the parameter network
  answered by pushing `parRT` to 19 mm/day, `parFC` to its 1000 mm cap and
  the runoff to a near-constant. On the closure probe's first gate seed:

  | Attributes | Q/P | AET/PET | runoff CV | residual / P | `parRT` | `parAC` | `parFC` |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | sine-fit seasonalities and snow-day share included | 0.63 | 0.98 | 0.07 | −58% | 18.8 | 1658 | 1000 |
  | as shipped (those three at the training mean) | 1.03 | 0.76 | 1.14 | −74% | 6.6 | 58 | 673 |
  | every attribute at the training mean | 0.81 | 0.71 | 1.06 | −48% | 4.3 | 99 | 915 |

  The runoff's variability depends on that guess; the size of the budget
  residual does not. Whatever the attributes, the model's outputs exceed
  its inputs by half to three quarters of the precipitation.

## Structure read from the checkpoint

The `config.yaml` first shipped with the weights declared four HBV components
and no routing; the one republished on 2026-10-31 declares three, and still no
routing. The weights disagree: the LSTM's output layer has 9 units
(three dynamic parameters × three components) and the MLP's has 41
(thirteen static parameters × three, plus two unit-hydrograph
parameters). Loaded as configured, the checkpoint does not fit. The
adapter reads the number of components and the presence of routing from
the checkpoint's own shapes, uses those, and records both what it found
and what the configuration declared in `run.json`. Nothing else is
overridden. The module's own daily example is marked "coming soon", so it
is likely nobody has loaded the daily weights through this path before.

## What to expect from the budget

HBV 2.0 in this package differs from textbook HBV in one term that
matters here. Every day the lower groundwater box receives

    LF = parRT * clamp((Ac - parAC) / 1000, -1, 1)          (Ac < 2500 km2)

where `Ac` is the upstream area and `parRT` (0–20 mm/day) and `parAC`
(0–2500 km2) are learned. It is the model's representation of regional
groundwater exchange between unit basins, and it is a source or a sink
that no reported flux or store accounts for. With `parRT` and `parAC`
static, it is a constant. On the closure probe's 250 km2 catchment it is a
gain of one to two mm every day, half to three quarters of the
precipitation over the record, and its size changes with the catchment
area handed in and with the attributes above. The closure criterion
measures exactly this residual. The unit hydrograph itself is normalised
to sum to one and conserves what enters it.

## The soil capacity

The probes' `soil_capacity_mm` (320 mm on the closure probe) is in `static.json`
but this model does not read it: HBV takes its capacity from the learned `parFC`.
`mrso` runs above that number on most probes, 494-668 mm on the closure probe. The
upper bound is N/A for this model, and no store runs away, so `state_bounds` passes
on what it can still check. It is not a claim that `mrso` stays below `parFC`, which
the contract has no way to report.

## Running it

```bash
ht verify-adapter --model dhbv2
ht run --model dhbv2 --csv models/result.csv
```

The image is built from the pinned commit and the weights URL at build
time; the container runs with no network. One lumped catchment over the
full ten-year record takes about five seconds, most of it loading the
400 MB of MLP weights.

## The weights changed on 2026-10-31

The URL carries no version. MHPI republished the archive that day with a new
`config.yaml`. The checkpoint's file date (2026-06-30) predates this model's
first evaluation, but no hash of the earlier files was kept, so what else
changed cannot be said. The adapter needed
one line, a default for the `multimodel_type` key the new config no longer
sets. SUMMA reproduced its archived closure numbers the same day, which points
at the files rather than the harness, and three of this model's verdicts moved
(`mass/extreme-rain` to PASS; `mass/resolution-invariance` and
`mass/response-nonnegativity` to FAIL). Version `0.5.4-hbv2ep100.4` is that
state, and the Dockerfile now pins its three files by SHA-256.

MHPI then removed the archive from the bucket, while its setup docs still point
to the URL, so a build with no cached copy of the download fails until the
archive returns. The `.4` rows were produced from a copy taken out of the local
image, whose hashes match the pins, with the download step swapped for a `COPY`.
If the archive returns with other files, the hash check fails the build, and
that is a new version to evaluate.

## Version `.5`

`0.5.4-hbv2ep100.5` changes the adapter only. `snowfall_fraction` is defined over
days and was computed per row, so an hourly and a daily run of the same weather
handed the parameter network different attributes. Daily rows are unchanged. On
`mass/resolution-invariance` the worst reading moved from `mrso_end` at 41.9 % of
the rain to `mrro` at 12.1 %, still over the 10 % limit, so the verdict stands.

## Version `.6`

`0.5.4-hbv2ep100.6` also changes the adapter only. The zero-rain guard on
`snowfall_fraction` checked the whole record while the fraction divided by the
whole days within it, so rain falling only in a last partial day gave the network
a NaN. The guard now checks those days. No row changed.

Archived on the gate seeds under Docker on Apple Silicon at `66118d8`: **FAIL
(VIOLATION), 15 PASS, 6 FAIL, 16 N/A**.

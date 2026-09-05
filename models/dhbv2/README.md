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
on closure rather than declared INCOMPLETE.

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
used. The forcing arrives as rates in mm per day at the case's step; the
adapter multiplies by the step length to give the model the depth per row
it wants and divides its fluxes back. It never resamples.

Of the 28 attributes the parameterisation network reads:

- **Derived from the first year of the record** (a climatology the model
  has seen, through which nothing later can reach back): `meanP`,
  `ETPOT_Hargr`, `aridity`, `meanTa`, and `snowfall_fraction` as the share
  of precipitation falling below the snow threshold.
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

The `config.yaml` shipped with the weights declares four HBV components
and no routing. The weights disagree: the LSTM's output layer has 9 units
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

## Running it

```bash
ht verify-adapter --model dhbv2
ht run --model dhbv2 --csv models/result.csv
```

The image is built from the pinned commit and the weights URL at build
time; the container runs with no network. One lumped catchment over the
full ten-year record takes about five seconds, most of it loading the
400 MB of MLP weights.

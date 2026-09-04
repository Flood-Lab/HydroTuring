# google_flood_forecast

Google's operational flood forecasting model, as released in
[google-research/flood-forecasting](https://github.com/google-research/flood-forecasting)
(OpenHydroNet) at commit `828dfc5`, with the published
`google-floodhub-settings-55-epochs` weights. Architecture: the mean-embedding
forecast LSTM ([Gauch et al. 2025](https://hess.copernicus.org/articles/29/6221/2025/));
the lineage is [Nearing et al. 2024](https://doi.org/10.1038/s41586-024-07145-1).
Requested in [issue #1](https://github.com/Flood-Lab/HydroTuring/issues/1).

## What it emits

Streamflow, and nothing else. `mrro` is the model's own prediction in mm/day;
`dis` is that depth times the catchment area, in m3/s. There is no
evapotranspiration, no storage, and the adapter invents none, so on every
budget probe the verdict is **FAIL (INCOMPLETE)**: the model has not violated
conservation, it has declined to be falsifiable on it. That is the honest
outcome for every streamflow-only model, and the reason the report carries a
reason next to the verdict.

## How the probe's forcing reaches the model

The model was trained on four weather products (HRES, GraphCast, IMERG, CPC)
and 84 Caravan/HydroATLAS attributes. A probe hands it daily precipitation,
temperature and potential ET for one lumped catchment. The adapter maps one
onto the other; where the probe has nothing to give, the input is mocked
from the forcing and labelled as such, which is acceptable for the
synthetic mass-balance test:

| Model input | Given |
| --- | --- |
| GraphCast, IMERG, CPC precipitation; GraphCast temperature | the forcing's `pr` and `tas`, standardised with each product's own training statistics |
| HRES precipitation and temperature | the forcing's `pr` and `tas`, as above |
| HRES net solar radiation, net thermal radiation, surface pressure | **mock inputs**, allowed for the synthetic test, derived from the forcing and nothing else: net radiation inverted from the probe's potential evaporation with Priestley-Taylor at the given temperature, net longwave from air temperature at 70 percent humidity and cloudiness on wet days, net solar as the remainder, pressure from the elevation the static attributes assume. Nothing follows the calendar, so constant weather gives constant inputs and a storm cannot reach earlier rows. `run.json` lists them under `mock_inputs`. `GFF_MOCK_HRES=0` marks the product missing instead, through the model's own NaN-aware mean over products |
| 14 climate attributes (`p_mean`, `pet_mean`, `aridity`, `frac_snow`, `moisture_index`, `seasonality`, high/low precipitation frequency and duration, annual P, PET, aridity index, mean temperature) | derived from the first year of the forcing the model is given, using Caravan's definitions. The first year only, so that nothing later in the record reaches back through the attributes: the causality probe adds a storm in year two and requires year one to be untouched |
| the other 70 attributes (land cover, terrain, soils, human footprint, ...) | the training mean, i.e. zero after standardisation: the mock catchment is an average Caravan basin at the probe's latitude |

Each day is its own forecast issue with a 365-day hindcast window, run from a
fresh state exactly as the operational model does. The reported value is the
day-0 member of that day's forecast: the prediction for the issue day given
everything up to and including it. Days early in the record use the history
that exists. The point prediction is the exact median of the model's own
CMAL mixture, found with the package's deterministic quantile search
rather than by drawing samples: the same quantity its tester estimates
from 7500 draws, without the estimator's noise, so that identical days give
identical answers and a stress probe sees the model rather than its
sampler. `run.json` records all of this per run.

## What each probe variable reaches inside the model

A probe that perturbs one forcing variable and holds the others needs to
know where each one goes. This is the complete map for this adapter:

| Probe variable | Model inputs it reaches | Held fixed by |
| --- | --- | --- |
| `pr` | the four precipitation inputs (HRES, GraphCast, IMERG, CPC); the cloudiness in the mocked longwave and solar radiation; the precipitation attributes `p_mean`, `aridity`, `moisture_index`, `seasonality`, `frac_snow`, `high_prec_*`, `low_prec_*`, `pre_mm_syr`, `ari_ix_sav` | a perturbation of `tas` or `pet` |
| `tas` | the two temperature inputs (HRES, GraphCast); the mocked net thermal and net solar radiation, through air temperature and the slope of the vapour pressure curve; `frac_snow` and `tmp_dc_syr` | a perturbation of `pr` |
| `pet` | the mocked net solar radiation, as the energy Priestley-Taylor needs to produce that demand; the demand attributes `pet_mean`, `aridity`, `moisture_index`, `seasonality`, `pet_mm_syr`, `ari_ix_sav` | a perturbation of `pr` |
| `static` | only `area_km2`, to convert depth to discharge; every other attribute the model wants is derived from the first year of the forcing above or held at its training mean | everything |

So the warming probe, which raises `tas` and `pet` together over identical
rain, reaches the model through its two temperature products and its
mocked radiation; the attributes, taken from the unperturbed first year,
stay fixed, as does every precipitation input. Attributes fixed is the
usual way an LSTM is run under a changed climate, and the reason its
response to warming is weaker than when the attributes were allowed to
follow the warmed record.

## Evaluation window

The submission asked for no particular window, so the manifest takes the
default for a daily model: the largest 30-day flood event of the record,
located by the probe's reference model, with the full spinup in front of it
(395 rows). One case takes about 20 s on 8 CPU cores; the full ten-year
record would take several minutes per seed and does not fit the probe's time
budget.

## Results (2026-09-03, suite 0.1.0)

| Check | Outcome |
| --- | --- |
| `ht verify-adapter` on `mass/catchment-closure`, gate seed 598896396 | contract OK: 395 rows in 19 s, columns `mrro`, `dis`, no missing products |
| `ht run` on `mass/catchment-closure` | FAIL (INCOMPLETE): does not report `pr`, `evspsbl`, `mrso`, `snw`, `canopy` |
| `ht run` on `mass/resolution-invariance` | FAIL (VIOLATION): runoff volume differs by 66% of precipitation between hourly and daily runs of the same month (17%, 40%, 66% on the three gate seeds) |
| `ht run` on `mass/warming-response` | FAIL (VIOLATION), marginal: the sign is right both ways on every seed, but under warming runoff falls by only 0.07 to 0.15 of the added demand (0.61 on one seed), under the 0.10 the probe requires on three of five seeds |
| `ht run` on `mass/causality` | PASS: nothing changes before the added storm, to floating point, and runoff answers it by 0.32 of the added rain |
| `ht run` on `mass/dry-down` | FAIL (VIOLATION): 224 mm drain in two rainless years, inside the 322 mm bound and never rising, but as a near-constant 0.3 mm/day that does not decay (runoff CV 0.02); a third dry year would break the bound |
| `ht run` on `mass/steady-state` | PASS: settles to 0.63 mm/day of runoff under 2.5 mm/day of rain, with zero variation in the last year |
| `ht run` on `mass/extreme-rain` | FAIL (VIOLATION): the response flattens past twice the largest storm; the ladder returns 0.52, 0.04 and -0.01 of the rain added on its three rungs, 0.07 overall |

Archived in [`models/result.csv`](../result.csv).

What the model did on the event it was handed (seed 598896396, 2006-03-18 to
2006-04-16, the reference bucket's largest 30-day flood), from
[`event_window_seed598896396.csv`](event_window_seed598896396.csv). The
second model column is the same adapter with the HRES product marked missing
instead of mocked, kept as a sensitivity check:

| | precipitation | reference_bucket | google_flood_forecast, all products | google_flood_forecast, HRES missing |
| --- | --- | --- | --- | --- |
| window total (mm) | 189 | 212 | 100 | 184 |
| peak (mm/day) | 55 | 49 | 12 | 29 |
| daily correlation with reference_bucket | | | 0.74 | 0.65 |

With every product present the model reproduces the timing of each
rain-driven peak and recession, but returns about half the rain that fell as
runoff and a quarter of the exact model's peaks, and it does not see the
late-March melt pulse that the bucket's snowpack produces. Dropping the HRES
product roughly doubles its runoff and makes it flash on winter
precipitation the bucket stores as snow. The model is that sensitive to a
product the probe cannot supply, which is the main caveat on any number
here: what is measured is the released weights under mock radiation and
pressure and average-catchment attributes, not the operational system.

![forcing, reference snowpack and runoff around the scored event](event_window_seed598896396.png)

### Response to warming and cooling

The warming probe holds the rain and shifts the air temperature by three
degrees over a full year, warmer in one variant and cooler in another,
with potential evaporation following it, and reaches this model through
the two temperature products and the mocked radiation (see the map
above). The attributes stay at the unperturbed first year's values, which
is how an LSTM is normally run under a changed climate.

The sign is the physical one on every seed and in both directions: runoff
falls when demand rises and rises when it falls. The size is small. Per
unit of demand, runoff falls by 0.07 to 0.15 under warming (0.61 on one
seed) and rises by 0.10 to 0.21 under cooling, against 0.20 and 0.24 for
the exact reference bucket on the same weather. The probe asks for at
least 0.10, and three of five seeds fall short under warming, so the
verdict is FAIL by a margin the report shows. An earlier version of this
adapter let the attributes follow the warmed record, and the model then
answered with 0.14 to 0.48: most of its apparent climate sensitivity was
carried by the attributes, not by the weather it was shown day to day.

### Stress probes

Four probes push the model where no training record goes; only runoff is
needed, so all four score it.

- **Causality** passes exactly. A 40 mm storm added in the second year
  changes nothing before its day, to floating point, and adds 0.32 of its
  water to the runoff afterwards. The adapter earns part of this: its
  attributes come from the first year only and its mocked inputs from the
  forcing alone, so nothing in the record reaches back.
- **Steady state** passes exactly. Under constant weather the model settles
  to 0.63 mm/day of runoff on 2.5 mm/day of rain and does not move in the
  last year. An earlier adapter that estimated the median by sampling
  showed 4 to 8 percent of day-to-day noise here; the deterministic
  quantile search removed it, and what remained was a model that settles.
- **Dry-down** fails, and the way it fails is the finding. Over two years
  without rain the model runs off 224 mm, inside the 322 mm the catchment
  could have held and never rising, but as a near-constant 0.3 mm/day that
  does not decay: a baseflow floor learned from its training basins. The
  exact bucket drains 55 mm and fades. A third dry year would put the
  model over the bound.
- **Extreme rain** fails in the classic way. Scaling the largest storm to
  twice its size returns 0.52 of the added rain as runoff; scaling it to
  five and ten times returns 0.04 and then nothing (-0.01). The response
  flattens at the edge of what the model has seen, on 526 mm of added
  rain a catchment with a 320 mm soil cannot absorb.

### Dependence on the step

The resolution probe serves one month of the same weather at the day and
at the hour and runs the model at both. The adapter feeds the rows as
given, so at the hourly step the model's 365-row hindcast covers fifteen
days and each hourly rate is read as a daily depth. Integrated over the
month, its runoff comes to 290 mm from the daily record and 126 mm from
the hourly one, on 248 mm of rain (seed 1760592044). That is the size of
the model's dependence on the step it was trained at, as a number rather
than a declaration: a model whose arithmetic assumed nothing about the
step would return the same volume from either record, as the reference
bucket does to within 1 percent.

## Reproduce

```bash
ht verify-adapter --model google_flood_forecast --csv models/result.csv   # builds the image, ~20 s per case
ht run --model google_flood_forecast --gate-seeds --csv models/result.csv # INCOMPLETE, without starting the container
ht verify-adapter --model google_flood_forecast --window full             # full record; exceeds the 60 s budget
```

To reproduce the HRES-missing column, rebuild with `ENV GFF_MOCK_HRES=0` in
the Dockerfile (or edit `MOCK_HRES` in the adapter) and run the contract
check again.

The image pins the source revision and the weight files; nothing is fetched
at run time.

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
| HRES net solar radiation, net thermal radiation, surface pressure | **mock inputs**, allowed for the synthetic mass-balance test: FAO-56 extraterrestrial radiation at the probe's latitude, attenuated on wet days, with albedo 0.23; net longwave from air temperature at 70 percent humidity; pressure from the elevation the static attributes assume. Their means land inside the training distribution (156 vs 147 W/m2, -59 vs -63 W/m2, 92.5 vs 93.0 kPa). `run.json` lists them under `mock_inputs`. `GFF_MOCK_HRES=0` marks the product missing instead, through the model's own NaN-aware mean over products |
| 14 climate attributes (`p_mean`, `pet_mean`, `aridity`, `frac_snow`, `moisture_index`, `seasonality`, high/low precipitation frequency and duration, annual P, PET, aridity index, mean temperature) | derived from the forcing the model is given, using Caravan's definitions |
| the other 70 attributes (land cover, terrain, soils, human footprint, ...) | the training mean, i.e. zero after standardisation: the mock catchment is an average Caravan basin at the probe's latitude |

Each day is its own forecast issue with a 365-day hindcast window, run from a
fresh state exactly as the operational model does. The reported value is the
day-0 member of that day's forecast: the prediction for the issue day given
everything up to and including it. Days early in the record use the history
that exists. The point prediction is the median of the model's own CMAL
mixture samples, seeded from the request so a case reproduces exactly.
`run.json` records all of this per run.

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
| `ht verify-adapter` on `mass/catchment-closure`, gate seed 598896396 | contract OK: 395 rows in 20 s, columns `mrro`, `dis`, no missing products |
| `ht run` on `mass/catchment-closure` | FAIL (INCOMPLETE): does not report `pr`, `evspsbl`, `mrso`, `snw`, `canopy` |
| `ht run` on `mass/resolution-invariance` | FAIL (VIOLATION): runoff volume differs by 58% of precipitation between hourly and daily runs of the same month (11%, 38%, 58% on the three gate seeds) |

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

### Dependence on the step

The resolution probe serves one month of the same weather at the day and
at the hour and runs the model at both. The adapter feeds the rows as
given, so at the hourly step the model's 365-row hindcast covers fifteen
days and each hourly rate is read as a daily depth. Integrated over the
month, its runoff comes to 259 mm from the daily record and 115 mm from
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

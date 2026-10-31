# MTS-LSTM under HydroTuring

The Multi-Timescale LSTM of [Gauch et al. (2021)](https://doi.org/10.5194/hess-25-2045-2021),
proposed in [issue #126](https://github.com/Flood-Lab/HydroTuring/issues/126). One network with a
daily and an hourly branch, run from the paper's pinned `neuralhydrology` with its published
weights. Nothing was trained here, by the maintainer's decision.

| Pin | Value |
| --- | --- |
| Code | [neuralhydrology](https://github.com/neuralhydrology/neuralhydrology) at `9626578d6013aa08864caf29d9520321cdcafc24`, BSD-3-Clause, the revision the proposal names |
| Weights | [Zenodo 18953458](https://doi.org/10.5281/zenodo.18953458) (concept 4071885), `models-mtslstm.tar.gz`, run `mtslstm_seed333_1208_060947`, `model_epoch030.pt` and `config.yml`. CC-BY-4.0 |
| Scaler | **Not this run's own**: `smtslstm_noregularization_seed555_1408_051556/train_data/train_data_scaler.p` from `models-smtslstm.tar.gz`. See below |
| Checks | sha256 of each member, recorded from our own first download. The archives are 22 and 6 GB and only the needed members are streamed, so Zenodo's archive md5 is not checked |
| Base images | `python:3.12.15-slim` by digest (torch 2.5.1 CPU, numpy 1.26.4, pandas 2.2.3, xarray 2026.9.0); `python:3.8.20-slim` by digest for one build stage |
| Image | 451 MB, no network at run time, never pushed to a registry |

## The scaler is borrowed, and how we know that is safe

This run's own scaler was never published. The complete listing of `models-mtslstm.tar.gz` has
40 runs, the ten `mtslstm_seed*` among them, and every `train_data` directory is empty. The
statistics therefore come from a sister sMTS-LSTM run on the same data and period. Its configuration
is not identical: hourly sequence length 72 rather than 336, and the scaler of a third run, the
naive hourly LSTM, differs from it by up to 3.1% on the input and target statistics. Nothing about this
scaler can be assumed from the configs, so it was measured.

With this scaler, the driver reproduces **this checkpoint's own** published simulation
(`test_results.p` of `mtslstm_seed333`) on all 6 of the 8 downloaded basins that the simulation
contains, over 457 days and 10968 hours of real NLDAS forcing and CAMELS attributes:

| Basin | Daily branch, RMSE (mm/h) | Largest | Hourly branch, RMSE (mm/h) | Largest | Published volume, daily and hourly branch (mm/day) |
| --- | --- | --- | --- | --- | --- |
| 02102908 | 7.7e-9 | 7.0e-8 | 5.4e-8 | 2.6e-7 | 1.021, 1.069 |
| 03026500 | 1.3e-8 | 1.0e-7 | 6.9e-8 | 1.2e-6 | 1.518, 1.498 |
| 05525500 | 1.4e-8 | 1.2e-7 | 7.5e-8 | 1.3e-6 | 1.160, 1.078 |
| 06477500 | 6.0e-9 | 2.5e-8 | 8.3e-9 | 5.6e-8 | 0.037, 0.048 |
| 07142300 | 5.8e-9 | 2.2e-8 | 1.0e-8 | 8.5e-8 | 0.067, 0.062 |
| 08202700 | 4.9e-9 | 1.8e-8 | 9.4e-9 | 1.1e-7 | 0.063, 0.070 |

That is single-precision agreement, which a scaler off even slightly would not give, so this scaler
is numerically indistinguishable from the run's own in everything the network reads. It is still a
measurement on 6 basins, not the original file. The image was opened on a compute node and hashed:
`model_epoch030.pt`, `config.yml` and `scaler.json` are byte-identical to the validated copies.

The pickle holds 2020 xarray objects that no longer open on a current stack, as the Zenodo record's
own 2026 note says, so `export_scaler.py` reads it during the build under the paper's pins
(numpy 1.18.5, pandas 1.0.5, xarray 0.16.0) and writes JSON. That is the Python 3.8 stage.

## What runs

`MTSLSTM` from upstream, unmodified, built from the run's own `config.yml`. Its branches are called
directly, because `MTSLSTM.forward` needs both resolutions in one sample and the `/io` contract
serves one step per case. The table above is what shows this matches upstream's own output.

The step in the request selects the branch, which answers the proposal review's questions:

- **PT1D:** the daily LSTM and head on the daily rows. That branch never sees hourly inputs.
- **PT1H:** the hourly LSTM and head on the hourly rows, started from the state the daily branch
  hands over. Each day of rows is one sample and its last 24 hours are kept, as upstream's tester does.
- **Any other step:** refused, since there is no branch for it.
- **`mrro` is reported** alongside `dis`, and what the verdict adds over a single-timescale LSTM is
  that the resolution probe compares two trained heads, each on the rows it was built for.

The proposal lists PT1H; the manifest says `[PT1D, PT1H]`. The checkpoint has both heads, daily
first makes every daily probe answerable without resampling, and hourly first would pair PT1H with
PT1M on the resolution probe, a step with no branch.

### One place a model input is aggregated

AGENTS.md says not to aggregate rows to reach a step the model prefers. At PT1H the daily branch
still has to read something, and upstream's own pipeline gives it the daily mean of the hourly record
(`BaseDataset`, `resample(freq).mean()`). The adapter does the same, for that branch alone and only
for whole days before the hourly window; no row is reported from it. This is flagged for the
maintainer to judge. `MTSLSTM_DAILY_BRANCH=none` skips it and transfers a zero state instead: the
worst seed goes from 49.7% to 62.2% on `resolution-invariance`, and the verdict does not move.

## What is reported

`mrro` is the head's output, de-standardised, with negatives set to zero as `clip_targets_to_zero`
does, converted from mm/hour to mm/day. `dis` is that depth times `area_km2`.

The network has no stores and no evaporation, so `states: []`. The zero-store rule does not apply,
because a discharge-only network's storage is unknown rather than absent. There is no `zero_stores`,
and `run.json` says so under `not_reported`. `pr` is not echoed and no flux is filled in. `abstr`,
`gwh` and `q_in` are not read and not declared.

## How the forcing reaches the network

Only the row's own `pr`, `tas` and `pet` are used, with no calendar. Nothing here is calibrated.
Statistics are the scaler's. Rows with less history use the history that exists.

| Network input | Given | Source |
| --- | --- | --- |
| `total_precipitation`, `potential_evaporation` | `pr / 24`, `pet / 24`, in mm/hour | direct, in the training units |
| `temperature` | `tas` | direct |
| `shortwave_radiation` | Priestley-Taylor (1.26) inverted from `pet`, plus FAO-56 net longwave with clearness 0.93 dry and 0.60 wet (1 mm/day), over 1 - albedo 0.23; mean over the row, W/m2 | mock |
| `longwave_radiation` | Brutsaert clear-sky emissivity at 70% relative humidity, W/m2 | mock |
| `specific_humidity` | from `tas` at 70% relative humidity, kg/kg | mock |
| `pressure`, `convective_fraction`, `potential_energy`, `wind_u`, `wind_v` | training mean, zero after standardisation | nothing in the row speaks to them |
| `p_mean`, `pet_mean`, `aridity`, `frac_snow`, `high_prec_freq`, `high_prec_dur`, `low_prec_freq`, `low_prec_dur` | derived from the first year of the record, or all of it when shorter | Addor et al. (2017) definitions |
| the other 18 attributes, `area_gages2` among them | training mean | no input to derive them from |

`area_km2` reaches `dis` only.

### An adapter fault found and fixed (0.1.0 to 0.1.1)

The first adapter applied the day-based indices to rows, so at PT1H the same weather described a
different catchment: `low_prec_freq` 339 against 201 days per year, `high_prec_dur` 0.09 against 1.0
days (seed 746687816). The CAMELS indices are defined over days, so they are now computed over whole
days at any step, and no input the network reads per row comes from that. It was the largest single
number in the first report: `mass/resolution-invariance` went from 197.6% to 49.7%, the hourly month on
the worst seed from 54 mm to 456 mm. The verdict did not change, and daily probes are unaffected, to
the digit.

## Results

From `ht run` on the gate seeds, under Apptainer, image `mtslstm_0.1.1-9626578-seed333_73d5ebca9a05`,
every row at version `0.1.1-9626578-seed333`. **8 PASS, 3 FAIL, 26 N/A** (24 INCOMPLETE, 2 INCOMPATIBLE),
no ERROR. The full evaluation took 17 minutes, of which the image build was 12.

| Probe | Verdict | Number | Mechanism | Status |
| --- | --- | --- | --- | --- |
| `mass/dry-down` | FAIL | Weekly runoff rises by 5.9% of its first block | Seasonal response to temperature and demand, below | verified |
| `mass/resolution-invariance` | FAIL | 49.7% of rain (27.6, 21.0); limit 10% | The two branches diverge mainly on attributes from a 40-day record, below | verified |
| `mass/warming-response` | FAIL | Runoff falls by 0.084 per unit of demand on 1 of 5 seeds; 0.10 required. The sign is right on every seed | Attributes follow the warmed spinup and push the other way, below | verified |
| `mass/antecedent-monotonicity` | PASS | +1.8, +2.1 and +6.1 mm of the 60 mm storm; 1.2 mm required | The smallest margins are thin | |
| `mass/area-invariance` | PASS | 0 departure | Trivial, below | |
| `mass/causality` | PASS | Nothing moves before the storm; +0.197 of the added rain after | Partly by design, below | |
| `mass/extreme-rain` | PASS | Rises on every rung; 0.24 of the 526 mm returned at the top | **Fragile**: fails on seed 111 and in three sensitivity runs | |
| `mass/phase-counterfactual` | PASS | At most 4.7% of the rain; limit 5% | `tas` does reach the network, so this is not trivial | |
| `mass/response-nonnegativity` | PASS | No dip beyond numerical noise | | |
| `mass/runoff-bounds` | PASS | Runoff is 0.12 to 0.18 of the rain | | |
| `mass/steady-state` | PASS | 0.397 mm/day on 2.5 mm/day of rain, no variation | Trivial, below | |

N/A (INCOMPLETE), 24 probes: all energy probes, every closure-family probe, the groundwater
exchange probes, `snowpack-mass-closure`, `routing-conservation`, `froude-regime`,
`stage-discharge-monotonic` and `uniform-flow-friction-consistency`. N/A (INCOMPATIBLE):
`routing-lag-consistency` and `wave-celerity-bounds`, whose inputs it does not consume.

### Resolution

This is the question the proposal asks.

| Seed (rain in the month) | PT1D, mm | PT1H, mm | Difference, % of rain | PT1H, each day's mean rate for 24 rows | Attributes at the mean: difference, % of rain |
| --- | --- | --- | --- | --- | --- |
| 746687816 (212 mm) | 206 | 147 | 27.6 | 178 | 2.1 |
| 1253639930 (228 mm) | 504 | 456 | 21.0 | 426 | 15.2 |
| 1760592044 (248 mm) | 272 | 148 | 49.7 | 191 | 8.4 |

On real basins with real inputs, the published simulations of the two branches agree well (last
column of the first table: within 0.09 mm/day on the wetter basins). Here they do not, and most of
the difference is the configuration. This probe's record is 40 days, shorter than a year, so the
attributes describe one wet month (`p_mean` 7 to 9 mm/day, `aridity` 0.22 to 0.25). On seed
1253639930 no day reaches the high-precipitation threshold, which gives `high_prec_freq` and
`high_prec_dur` of zero, outside the training range, and both branches then return about twice the
rain that fell or more.

At the training mean the difference falls to 2 to 15%, still over the limit on one seed. Within-day
structure moves the hourly volume by 30 to 43 mm. The 10-day spinup also fills little of the
daily branch's 365-row window: repeating the record in front moves the daily month from 206, 504
and 272 mm to 163, 466 and 211 mm.

So the probe measures the branches outside anything they trained on, and this FAIL is largely a
property of the configuration with a smaller model part underneath it. It should not be read as
evidence against the paper's cross-timescale consistency on real forcing. Deriving from whatever
record exists follows `google_flood_forecast` and `dhbv2`. Whether a record under a year should yield
attributes at all is a question for the maintainer. The verdict is FAIL either way.

### Warming

| Seed | Control runoff, mm | Warmer, as run | Warmer, attributes frozen at the control | Attributes alone |
| --- | --- | --- | --- | --- |
| 732983943 | 123 | **-0.084** | -0.136 | +0.057 |
| 1239936057 | 182 | -0.151 | -0.181 | +0.018 |
| 613308752 | 190 | -0.164 | -0.152 | -0.007 |
| 1746888171 | 282 | -0.257 | -0.262 | +0.015 |
| 106356638 | 287 | -0.289 | -0.309 | +0.033 |

The evaluation window moves the 365 spinup rows inside the warmed years, so on all five seeds the
first year the adapter sees is itself warmed (up to 3 degC). `pet_mean` and `aridity` then rise. On
the failing seed the network answers a more arid description with **more** runoff (+0.057), and the
weather's own response is cut to below the limit. With attributes frozen at the control's values
every seed passes in both directions, so this verdict turns on the attribute route. Seed 111 passes
as run. Holding the mock forcing at its means leaves that route alone, and runoff then rises under
warming, by +0.087 per unit.

This follows AGENTS.md, which says to take a climatology from the first year of the record, so it
is not an adapter fault. The probe's text does say the variants enter the window from the same state,
which the windowed spinup contradicts.

### Dry-down

Weekly runoff rises and falls against temperature through rows with no rain at all (correlation
-0.49, -0.76, -0.73 by seed). With `tas` and `pet` held constant over the dry stretch, the worst rise
on the three seeds falls from 8.6, 9.8 and 5.9% to 0.00, 0.04 and 0.00%, and runoff recedes smoothly
(0.160 to 0.111, 0.159 to 0.106, 0.237 to 0.195 mm/day). So there is a real recession with a learned
seasonal signal on top. The report names 5.9% because it picks the seed with the largest drained
volume, which is not the seed with the largest rise. Totals are 86 to 141 mm in two years against a
322 mm bound. A property of the model.

### Trivial passes

- `mass/area-invariance`: `area_gages2` is held at the mean, so area enters only `dis`. Feeding the
  case's area to that attribute instead fails the probe by a relative departure of 21.
- `mass/steady-state`: with a fresh state per row, constant weather gives identical windows after 365
  rows, so the output cannot vary.
- The causality criterion's "nothing moves": a forward LSTM with first-year attributes and row-only
  mocks cannot see a later storm. That it answers the storm afterwards is a real result.

## Attribute sensitivity

Prototype runs of the same adapter under the real harness on the gate seeds. The all-derived column
matches the archived rows on every verdict and number. Each cell gives the verdict and the report's
number, which for a PASS may be one seed's rather than the worst.

| Probe | All derived | Precipitation at mean | Demand at mean | Snow at mean | Everything at mean |
| --- | --- | --- | --- | --- | --- |
| antecedent-monotonicity, mm | PASS +6.1 | PASS +7.1 | PASS +6.4 | PASS +7.3 | PASS +8.8 |
| area-invariance | PASS | PASS | PASS | PASS | PASS |
| causality, share answered | PASS 0.197 | PASS 0.282 | PASS 0.227 | PASS 0.220 | PASS 0.425 |
| dry-down, rise | FAIL 5.9% | FAIL 4.4%, also near-constant | FAIL 3.4%, also near-constant | FAIL 12.6% | FAIL 2.6%, also near-constant |
| extreme-rain, returned at the top | PASS 0.24 | **FAIL** 0.10 | PASS 0.24 | **FAIL** 0.09 | PASS 0.26 |
| phase-counterfactual | PASS 4.7% | PASS 5.0% | PASS 4.4% | PASS 1.2% | PASS 0.5% |
| resolution-invariance | FAIL 49.7% | FAIL 19.2% | FAIL 35.9% | FAIL 37.6% | FAIL 15.2% |
| response-nonnegativity | PASS | PASS | PASS | PASS | PASS |
| runoff-bounds, runoff ratio | PASS 0.18 | PASS 0.40 | PASS 0.22 | PASS 0.19 | PASS 0.38 |
| steady-state, mm/day | PASS 0.397 | PASS 0.447 | PASS 0.434 | PASS 0.230 | PASS 0.528 |
| warming-response | FAIL -0.084 | PASS | PASS | FAIL -0.091 | PASS |
| **PASS / FAIL** | **8 / 3** | **8 / 3** | **9 / 2** | **7 / 4** | **9 / 2** |

Dry-down and resolution fail in every column. Warming follows the demand and precipitation groups,
as above. Extreme rain is the verdict least settled: holding precipitation or snow at the mean turns
it into a FAIL.

`MTSLSTM_MOCK=mean`, every input but precipitation at the training mean, gives 7 / 4: dry-down
(near-constant, cv 0.046), extreme-rain (0.07), resolution (91.3%) and warming (+0.087) fail.
`MTSLSTM_AREA=static` was checked on three probes: area-invariance fails, causality and runoff-bounds
pass.

## Another checkpoint

Results describe one documented point: one repetition, not the ten-member ensemble the paper
reports.

| Checkpoint | PASS / FAIL | What moves |
| --- | --- | --- |
| `mtslstm_seed333` (evaluated) | 8 / 3 | |
| `mtslstm_seed111` | 8 / 3 | Warming passes; extreme-rain fails (0.08 returned); dry-down 10.6%; resolution 82.9% |

Dry-down and resolution fail on both. The count is the same but one verdict swaps for another, so
the warming and extreme-rain verdicts belong to this checkpoint at this setting, not to the
architecture.

## Against `google_flood_forecast` and `ealstm`

All three are discharge-only LSTMs run from published weights with mocked radiation and attributes
derived from the first year. The Google rows are the latest archived for each probe, from 2026-09,
and were not re-run.

| Probe | google_flood_forecast | ealstm | mtslstm |
| --- | --- | --- | --- |
| antecedent-monotonicity | PASS | FAIL | PASS |
| area-invariance | PASS | PASS | PASS |
| causality | PASS | FAIL (`non_degenerate` only) | PASS |
| dry-down | FAIL (near-constant floor) | FAIL (seasonal rises) | FAIL (seasonal rises on a recession) |
| extreme-rain | FAIL (0.07 returned) | PASS (0.31) | PASS (0.24) |
| phase-counterfactual | FAIL (8.4%) | PASS (3.8%) | PASS (4.7%) |
| resolution-invariance | FAIL (65.9%) | FAIL (60.5%) | FAIL (49.7%) |
| response-nonnegativity | FAIL (0.176 mm/day dip) | PASS | PASS |
| runoff-bounds | PASS | FAIL (1 seed of 5) | PASS |
| steady-state | PASS | PASS | PASS |
| warming-response | FAIL (-0.088) | FAIL (-0.069) | FAIL (-0.084) |

These differ in inputs as well as architecture (84, 27 and 26 attributes; four, one and one weather
products), so a row-by-row difference is not a statement about architecture. The Google adapter still
computes the day-based indices on rows at PT1H (`climate_attributes`, line 210), the fault fixed here,
which moved this model by 148 percentage points. Its 65.9% may carry the same effect. That was not
measured.

## Limits

- One repetition, seed 333 at epoch 30, with a sister run's scaler. The paper reports an ensemble.
- A synthetic catchment with mocked radiation and humidity, five inputs and 18 of 26 attributes at the
  training mean. This measures the published network on inputs it never saw, not its skill on CAMELS
  basins.
- Budget questions cannot be put to it: 26 of 37 probes are N/A.
- The sensitivity table and the targeted checks ran in a prototype environment, with the same adapter
  and upstream code as the image, and its numbers match the image's rows.
- `mtslstm_seed111` uses the same scaler. Its numbers here were not re-checked against its own
  published simulation; only seed 333's were.

## Reproduce

```bash
sbatch --nodes=1 --chdir="$PWD" --job-name=ht-mtslstm --wait scripts/hpc/run_model.sbatch mtslstm
```

The sensitivity columns use `MTSLSTM_HOLD=precipitation|demand|snow|all`, `MTSLSTM_MOCK=mean`,
`MTSLSTM_AREA=static` and `MTSLSTM_DAILY_BRANCH=none`, and `MTSLSTM_RUN_DIR` for another checkpoint.
None of these is set in the image.

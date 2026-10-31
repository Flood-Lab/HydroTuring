# EA-LSTM under HydroTuring

**This proposal is not accepted.** [Issue #57](https://github.com/Flood-Lab/HydroTuring/issues/57)
carries the label `revision`, not `accepted`. The maintainer asked for this package anyway,
and CONTRIBUTING.md's normal route was not followed.

The Entity-Aware LSTM of [Kratzert et al. (2019)](https://doi.org/10.5194/hess-23-5089-2019),
run from the paper's own code with the paper's published weights. Nothing was trained here,
by the maintainer's decision.

| Pin | Value |
| --- | --- |
| Code | [kratzert/ealstm_regional_modeling](https://github.com/kratzert/ealstm_regional_modeling) at `d11815824793737f28f125d55477e42f8900a78e`, Apache-2.0 |
| Weights | [HydroShare 83ea5312635e44dc824eeb99eda12f06](https://www.hydroshare.org/resource/83ea5312635e44dc824eeb99eda12f06/), run `run_1906_1004_seed111`, `model_epoch30.pt`: EA-LSTM, NSE loss, repetition seed 111. CC BY (the resource's own rights statement) |
| Archive check | md5 `224985a9c4d0932e607f4c9fe8d61ee8`, the checksum HydroShare states; the build stops on a mismatch |
| Base image | `python:3.12.15-slim` by digest; torch 2.5.1 (CPU), numpy 1.26.4, pandas 2.2.3 |
| Image | 378 MB, no network at run time, never pushed to a registry |

## How this differs from the proposal

The proposal names `neuralhydrology` v1.13.0 and says that retraining under it is acceptable.
Published weights only rules that out, because these weights belong to the paper's code, so:

- The code is the paper's repository, not `neuralhydrology`, and its licence is Apache-2.0,
  not the BSD-3-Clause the proposal lists.
- The window is the checkpoint's 270 rows, from its `cfg.json`, not the proposal's 365.
- The matched comparison with MC-LSTM (#56) that the proposal argues for, same library and same
  training run, is **not** what this package gives. An MC-LSTM trained under v1.13.0 would
  differ from this model in code, training period and forcing handling as well as in
  architecture.

## What runs, and the scaler

`main.Model` and `papercode.ealstm.EALSTM` from upstream, unmodified, built from the run's own
`cfg.json`, as `main.evaluate` builds them. The image was opened on a compute node and hashed:
`model_epoch30.pt`, `cfg.json` and `attributes.db` are byte-identical to the members of the
published archive.

The statistics are upstream's own:

- Inputs and output use `SCALER`, hard-coded in `papercode/datautils.py`.
- Static attributes use the mean and standard deviation over the 531 basins of the run's
  `attributes.db`, which is how `main.evaluate` computes them.

**Reproduction.** With real Maurer forcing and real CAMELS attributes, the driver reproduces
the published `ealstm_seed111.p` on all 7 basins checked (730 days each): RMSE at most
2.3e-6 mm/day, largest single difference 3.5e-5 mm/day, against a standard deviation of up to
8.3 mm/day. So what follows comes from the published network under the inputs described here,
not from the driver.

## What is reported

`mrro` is the network's prediction de-standardised by upstream's `rescale_features`, with
negatives set to zero as `evaluate_basin` does, in mm/day. `dis` is that depth times `area_km2`.

The network has no stores and no evaporation, so `states: []`. The zero-store rule does not
apply, because a discharge-only network's storage is unknown rather than absent. There is no
`zero_stores`, and `run.json` says so under `not_reported`. `pr` is not echoed and no flux is
filled in. `abstr`, `gwh` and `q_in` are not read and not declared.

## How the forcing reaches the network

Each row is predicted many-to-one from a fresh state with a window of 270 rows ending on that
row, as `main.evaluate` does. Rows with less history behind them use the history that exists.
Only the row's own `pr`, `tas` and `pet` are used, with no calendar. Nothing here is calibrated.

| Network input | Given | Source |
| --- | --- | --- |
| PRCP | `pr` | direct |
| Tmax, Tmin | `tas` plus and minus half the difference of the training means (6.67 degC) | mock |
| Vp | saturation vapour pressure at Tmin, in Pa | mock, FAO-56 eq. 48 |
| SRAD | Priestley-Taylor (1.26) inverted from `pet`, plus FAO-56 net longwave with clearness 0.93 dry and 0.60 wet (1 mm/day), over 1 - albedo 0.23, doubled to a daylight mean under a fixed 12-hour day | mock. The training mean of 372 W/m2 is a daylight mean; the true day length follows the calendar, so it cannot be used |
| `p_mean`, `pet_mean`, `aridity`, `frac_snow`, `high_prec_freq`, `high_prec_dur`, `low_prec_freq`, `low_prec_dur` | derived from the first year of the record, or all of it when shorter | Addor et al. (2017) definitions |
| the other 19 attributes, `p_seasonality` and `area_gages2` among them | training mean, zero after standardisation | no definition, or no input, to derive them from |

Pressure in the SRAD inversion uses the training mean elevation. `area_km2` reaches `dis` only.

### An adapter fault found and fixed (0.1.0 to 0.1.1)

The first adapter applied the day-based indices to rows, so at PT1H the same weather described a
different catchment: `low_prec_freq` 339 against 201 days per year, `high_prec_dur` 0.09 against
1.0 days (seed 746687816). The CAMELS indices are defined over days, so they are now computed
over whole days at any step. This is the only place rows are grouped, and no input the network
reads per row comes from it. It moved `mass/resolution-invariance` from 70.1% to 60.5% and
changed no verdict. Daily probes are unaffected, to the digit.

## Results

From `ht run` on the gate seeds, under Apptainer, image `ealstm_0.1.2-d118158-seed111_d48f2acbf59b`,
every row at version `0.1.2-d118158-seed111`. **5 PASS, 6 FAIL, 26 N/A** (24 INCOMPLETE, 2 INCOMPATIBLE),
no ERROR. The full evaluation took 7 minutes including the 3-minute image build.

| Probe | Verdict | Number | Mechanism | Status |
| --- | --- | --- | --- | --- |
| `mass/antecedent-monotonicity` | FAIL | +0.8 mm of a 60 mm storm on the wetter catchment (+1.3, +6.5 on the others); 1.2 mm required | Low-runoff regime, below | verified |
| `mass/causality` | FAIL | Causality itself passes, nothing moves before the storm. `non_degenerate`: runoff ratio 0.0169 (0.0387, 0.0281); 0.02 required | Low-runoff regime, below | verified |
| `mass/dry-down` | FAIL | Weekly runoff rises by 19.0% of its first block | Seasonal response to temperature and demand, below | verified |
| `mass/resolution-invariance` | FAIL | 60.5% of rain (32.3, 37.4); limit 10% | No step in the network's arithmetic, plus attributes from a 40-day record, below | verified |
| `mass/runoff-bounds` | FAIL | 972 mm ran off where at least 1134 mm had to; 1 seed of 5 | Low-runoff regime, below | verified |
| `mass/warming-response` | FAIL | Runoff falls by 0.069 and 0.051 per unit of demand on 2 of 5 seeds; 0.10 required. The sign is right on every seed | Low-runoff regime, below | verified |
| `mass/area-invariance` | PASS | 0 departure | Trivial, below | |
| `mass/extreme-rain` | PASS | Rises on every rung; 0.31 of the 526 mm returned at the top | | |
| `mass/phase-counterfactual` | PASS | At most 3.8% of the rain | `tas` does reach the network, so this is not trivial | |
| `mass/response-nonnegativity` | PASS | Largest dip 4.4e-4 of the added 120 mm | | |
| `mass/steady-state` | PASS | 0.259 mm/day on 2.5 mm/day of rain, no variation | Trivial, below | |

N/A (INCOMPLETE), 24 probes: all energy probes, every closure-family probe, the groundwater
exchange probes, `snowpack-mass-closure`, `routing-conservation`, `froude-regime`,
`stage-discharge-monotonic` and `uniform-flow-friction-consistency`. N/A (INCOMPATIBLE):
`routing-lag-consistency` and `wave-celerity-bounds`, whose inputs it does not consume.

### One mechanism behind four failures: the low-runoff regime

The failures on causality, runoff bounds, antecedent wetness and warming are one finding. Derived
from the first year, the attributes describe a fairly dry basin, inside the training range apart
from `high_prec_dur`, which sits just under its minimum of 1.075 on most seeds (1.00 to 1.06).
For example, on causality seed 1016628858: `p_mean` 1.59 (z -1.3), `aridity` 1.70 (z +1.4),
`low_prec_freq` 306 (z +1.6). Given those, the network returns very little water.

Clipping is not the cause: at most 5% of rows are clipped on these probes, and the unclipped
volume is the same to 1 mm. The mock forcing is not the cause either: with it held at the
training means, causality still fails (0.0160). The precipitation attributes are:

| Check, on the failing seed | All derived | With attributes moved |
| --- | --- | --- |
| Causality, runoff ratio | 0.0169 | `p_mean` at the mean 0.035; `aridity` at the mean 0.037 |
| Runoff bounds, runoff ratio | 0.105 | `p_mean` at the mean 0.200; `low_prec_freq` at the mean 0.148 |
| Antecedent, wet and dry after the storm | 3.4 and 2.6 mm | Precipitation group at the mean: 10.3 and 5.9 mm |
| Warming, seeds by control runoff | 30, 67, 101, 213, 248 mm | Response -0.051, -0.069, -0.184, -0.237, -0.260: it shrinks with the runoff there is |

On runoff bounds the first year is also drier than the decade it stands for: `p_mean` 1.91 mm/day
against 2.54 over the record. The network has no budget through which later rain could correct
that. This is a property of the published network at this attribute setting, not an adapter
fault: the definitions are the published ones and were not chosen to move a result.

On warming, the attributes do not stay fixed. The evaluation window moves the 365 spinup rows
inside the warmed years, so the first year the adapter sees is itself warmed on all five seeds
(up to 3 degC; `pet_mean` rises by 0.08 to 0.39 mm/day). Frozen at the control's values, the two
failing seeds still fail, at -0.042 and -0.060, so this verdict does not depend on it.

### Dry-down

Weekly runoff rises and falls with the season, against temperature (correlation -0.97, -0.97,
-0.95 by seed), through rows that have no rain at all. With `tas` and `pet` held constant over
the dry stretch, the worst rise on the three seeds falls from 19.0, 36.6 and 32.8% to 0.04, 0.16
and 2.53%, the last only in the first block. The report names 19.0% because it picks the seed
with the largest drained volume, which is not the seed with the largest rise. The volumes are
small, 7 to 13 mm in two years against a 322 mm bound, so a 19% rise is about 0.006 mm/day. A
property of the model.

### Resolution

| Seed (rain in the month) | PT1D, mm | PT1H, mm | Difference, % of rain | PT1H, each day's mean rate for 24 rows | Attributes at the mean: difference, % of rain |
| --- | --- | --- | --- | --- | --- |
| 746687816 (212 mm) | 295 | 167 | 60.5 | 215 | 5.1 |
| 1253639930 (228 mm) | 253 | 180 | 32.3 | 322 | 6.1 |
| 1760592044 (248 mm) | 310 | 217 | 37.4 | 258 | 10.7 |

The network has no step in its arithmetic. At PT1H the window covers 11 days instead of 270 and
each hourly rate is read as a daily one, and removing the structure inside the day moves the
hourly volume by 41 to 142 mm, which is the model's own dependence on the step.

Most of the size comes from the configuration. This probe's record is 40 days, shorter than a
year, so the attributes describe one wet month (`p_mean` 7 to 9 mm/day, `aridity` 0.22 to 0.25),
and on seed 1253639930 no day reaches the high-precipitation threshold at all, which gives
`high_prec_freq` and `high_prec_dur` of zero, outside the training range. At the training mean
the difference is 5 to 11%, still over the limit on one seed.

The 10-day spinup also leaves most of the 270-row window empty: repeating the record in front
moves the daily month from 295, 253 and 310 mm to 196, 299 and 275 mm.

Deriving from whatever record exists follows `google_flood_forecast` and `dhbv2`. Whether a
record under a year should yield attributes at all is a question for the maintainer. The verdict
is FAIL either way.

### Trivial passes

- `mass/area-invariance`: `area_gages2` is held at the mean, so area enters only `dis`. Feeding
  the case's area to that attribute instead fails the probe by a relative departure of 21.
- `mass/steady-state`: with a fresh state per row, constant weather gives identical windows after
  270 rows, so the output cannot vary.
- The causality criterion itself: a forward LSTM with first-year attributes and row-only mocks
  cannot see a later storm. Its 0 departure is the design, not a measurement of restraint.

## Attribute sensitivity

Prototype runs of the same adapter under the real harness on the gate seeds. The all-derived
column matches the archived rows on every verdict and number. Each cell gives the verdict and
the report's number, which for a PASS may be one seed's rather than the worst.

| Probe | All derived | Precipitation at mean | Demand at mean | Snow at mean | Everything at mean |
| --- | --- | --- | --- | --- | --- |
| antecedent-monotonicity | FAIL +0.8 mm | PASS +10.5 | FAIL +0.8 | FAIL +1.0 | PASS +9.8 |
| area-invariance | PASS | PASS | PASS | PASS | PASS |
| causality | FAIL ratio 0.017 | PASS | PASS | FAIL ratio 0.014 | PASS |
| dry-down, rise | FAIL 19.0% | FAIL 5.5% | FAIL 16.4% | FAIL 14.1% | FAIL 10.3% |
| extreme-rain, returned at the top | PASS 0.31 | PASS 0.31 | PASS 0.30 | PASS 0.31 | PASS 0.32 |
| phase-counterfactual | PASS 3.8% | PASS 2.8% | PASS 1.9% | PASS 3.8% | PASS 2.7% |
| resolution-invariance | FAIL 60.5% | FAIL 16.0% | FAIL 61.7% | FAIL 49.2% | FAIL 10.7% |
| response-nonnegativity | PASS | PASS | PASS | PASS | **FAIL** 0.457 mm/day dip |
| runoff-bounds | FAIL 972 mm | PASS | PASS | FAIL 1079 mm | PASS |
| steady-state, mm/day | PASS 0.259 | PASS 0.609 | PASS 0.300 | PASS 0.277 | PASS 0.764 |
| warming-response | FAIL -0.069 | PASS | PASS | FAIL -0.030 | PASS |
| **PASS / FAIL** | **5 / 6** | **9 / 2** | **8 / 3** | **5 / 6** | **8 / 3** |

The precipitation group carries most of it, and removing all attributes is not simply better:
an average basin passes the low-runoff probes and then fails `response-nonnegativity`, which
the derived configuration passes. Dry-down and resolution fail in every column. Steady runoff
under the same rain ranges from 0.26 to 0.76 mm/day, a factor of three, on attributes alone.

`EALSTM_MOCK=mean`, holding SRAD, Tmax, Tmin and Vp at their training means, gives 6 / 5:
antecedent passes, and the other five failures stay. `EALSTM_AREA=static` was checked on three
probes: area-invariance then fails, and causality and runoff-bounds keep failing.

## Other checkpoints

Results describe one documented point: one repetition of the paper's ensemble, not the ensemble.

| Checkpoint | PASS / FAIL | What moves |
| --- | --- | --- |
| `run_1906_1004_seed111`, NSE (evaluated) | 5 / 6 | |
| `run_1906_1005_seed222`, NSE | 6 / 5 | Causality and runoff-bounds pass; extreme-rain now fails (runoff falls by 6 mm between the x5 and x10 rungs); resolution 100.8% |
| `run_1906_2336_seed111`, MSE | 8 / 3 | Antecedent, causality, runoff-bounds and warming pass; phase-counterfactual fails (9.7%); resolution 119.9%, dry-down 50.0% |

Dry-down and resolution fail on all three. The low-runoff verdicts are not stable across
repetitions, so treat them as properties of this checkpoint at this setting, not of the
architecture.

## Against `google_flood_forecast`

Both are discharge-only LSTMs run from published weights with mocked radiation and attributes
derived from the first year, so the comparison is of the same kind. It is not of the same
vintage: the Google rows are the latest archived for each probe, from 2026-09, and were not
re-run.

| Probe | google_flood_forecast | ealstm |
| --- | --- | --- |
| antecedent-monotonicity | PASS | FAIL |
| area-invariance | PASS | PASS |
| causality | PASS | FAIL (`non_degenerate` only) |
| dry-down | FAIL (near-constant floor, cv 0.02) | FAIL (seasonal rises) |
| extreme-rain | FAIL (flattens, 0.07 returned) | PASS (0.31) |
| phase-counterfactual | FAIL (8.4%) | PASS (3.8%) |
| resolution-invariance | FAIL (65.9%) | FAIL (60.5%) |
| response-nonnegativity | FAIL (0.176 mm/day dip) | PASS |
| runoff-bounds | PASS | FAIL (1 seed of 5) |
| steady-state | PASS | PASS |
| warming-response | FAIL (-0.088) | FAIL (-0.069) |

Two caveats on reading this. The Google model takes 84 attributes and four weather products, this
one 27 and one, so a row-by-row difference is not a statement about architecture. And the Google
adapter still computes the day-based indices on rows at PT1H (`climate_attributes`, line 210), the
fault fixed here, so its 65.9% may carry the same effect. That was not measured.

## Limits

- One repetition, seed 111 at epoch 30. The paper reports an ensemble.
- A synthetic catchment under mock radiation and humidity, with 19 of 27 attributes at the
  training mean. This measures the published network on inputs it never saw, not its skill on
  CAMELS basins.
- Budget questions cannot be put to it: 26 of 37 probes are N/A.
- PT1D is the only step declared, the only one upstream supports. PT1H is run because the
  resolution probe runs it regardless.

## Reproduce

```bash
sbatch --nodes=1 --chdir="$PWD" --job-name=ht-ealstm --wait scripts/hpc/run_model.sbatch ealstm
```

The sensitivity columns use `EALSTM_HOLD=precipitation|demand|snow|all`, `EALSTM_MOCK=mean` and
`EALSTM_AREA=static`, and `EALSTM_RUN_DIR` for another checkpoint. None of these is set in the image.

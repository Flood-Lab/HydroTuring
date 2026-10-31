# xaj

The Xin'anjiang model (Zhao 1992), a conceptual rainfall-runoff model used for flood
forecasting in China since the 1980s: three-layer evaporation, saturation-excess
runoff over a distributed tension-water capacity, a free-water store separating
surface, interflow and groundwater runoff, and lag-recession channel routing.
Proposed in issue #97.

## Provenance

The code is `hydromodel/models/xaj.py` from
[OuyangWenyu/hydromodel](https://github.com/OuyangWenyu/hydromodel) v0.4.0, commit
`89d7a8ed1d72ce4fffbbd9897490b089382ecbac`, GPL-3.0. It is run unmodified, not
reimplemented, with its defaults: `name="xaj"`, `source_type="sources"`,
`source_book="HF"`.

The image copies that file and the three modules it imports (`model_config.py`,
`param_utils.py`, `unit_hydrograph.py`) byte for byte, under empty package `__init__`
files, and installs numpy, scipy, numba, llvmlite and pyyaml at the versions
upstream's `uv.lock` records. It does not `pip install` the package, for two reasons.
The package's `__init__` imports `hydroutils` only to get a cache path, and that
brings in async-retriever, cartopy and a data stack forward simulation never uses.
And `uv.lock` cannot be installed as shipped, because it resolves `hydrodataset` from
a sibling checkout outside the repository. The base image is pinned by digest.

`xaj()` returns streamflow and evaporation for every row but only the final stores.
So `mrro` and `evspsbl` are what `xaj()` itself returned, divided by the row length,
and the stores come from driving the same module's `generation`, `sources` and
`linear_reservoir` row by row. The two are compared on every run, the adapter refuses
to report if they differ by more than 1e-9 mm per row, and `run.json` records the
difference.

## Stores

Catchment-average, in mm, weighted as upstream weights the runoff they feed.

| Column | What it is |
| --- | --- |
| `mrso` | `(1 - IM) * (WU + WL + WD + S * FR)`: the three tension layers, and the free water, which is `S` deep over the runoff-producing fraction `FR` |
| `gw` | the interflow and groundwater linear reservoirs. Upstream carries each as its outflow `Q` alone; the water that stands for is `Q * C / (1 - C)`, which balances upstream's recursion exactly. Both are subsurface, so the interflow reservoir is not in `channel` |
| `channel` | the channel routing: runoff inside the lag plus the recession store, `sum(last L inflows) + Q * CS / (1 - CS)`, which balances row by row to 1e-14 mm |

## Parameters

No calibration exists for a synthetic catchment. Nothing below was chosen by looking
at a verdict. The values are `configs/example_xaj_params.yaml` at the pinned commit,
the file the proposal runs, which upstream describes as "typical calibrated values
for demonstration purposes".

| Parameter | Meaning | Value | At a finer step |
| --- | --- | --- | --- |
| `K` | ratio of potential evaporation to the forcing's | 0.75 | - |
| `B` | exponent of the tension-water capacity curve | 0.25 | - |
| `IM` | impervious fraction | 0.06 | - |
| `UM`, `LM`, `DM` | tension-water capacity of each layer | 18, 80, 95 mm, then scaled so they sum to `soil_capacity_mm`: 29.8, 132.6, 157.5 mm at 320 mm | - |
| `C` | deep evaporation coefficient | 0.12 | - |
| `SM` | free-water capacity | 45 mm, kept, as `models/sacsma_snow17` keeps its lower free-water stores | - |
| `EX` | exponent of the free-water capacity curve | 1.3 | - |
| `KI`, `KG` | free-water outflow to interflow, to groundwater | 0.35, 0.45 per day | upstream's own conversion in `sources5mm` |
| `CS` | channel recession | 0.5 per day | `CS^dt` |
| `L` | lag | 5.5, which upstream's `int()` makes 5 days | `int(L) / dt` rows, so the same time |
| `CI`, `CG` | interflow, groundwater recession | 0.65, 0.985 per day | `C^dt` |

Initial stores are upstream's for `warmup_length=0`: each layer half full, `S` at half
`SM`, `FR` 0.1, both reservoir outflows 0.1. Those two outflows are in mm per row
upstream, so the adapter passes them through `xaj()`'s `initial_states` as 0.1 mm per
day, which is the identity at a daily step. See "A failure that was the adapter's".

## Choices, and how much each moves

Measured on the gate records of `mass/catchment-closure` (runoff ratio, E/PET, budget
residual) and `mass/resolution-invariance`. These are measurements outside the
harness; only the first row went through it. The default's ratios here are the mean
of five seeds.

| Configuration | Runoff ratio | E/PET | Residual, % of rain | mrro moved PT1D to PT1H, worst % of rain |
| --- | --- | --- | --- | --- |
| **default** | 0.513 | 0.562 | -2.435 | 4.48 |
| `UM`, `LM`, `DM` as published, capacity ignored | 0.560 | 0.507 | -2.133 | 4.28 |
| `K = 1` | 0.437 | 0.646 | -2.720 | 4.03 |
| `IM = 0.01`, upstream's lower bound | 0.492 | 0.562 | -0.406 | 5.11 |
| rate constants and lag left per row | 0.513 | 0.562 | -2.435 | 35.61 |

`K` sets `energy/pet-consistency` directly: wet-soil evaporation is `K` times the
forcing's potential, and the probe's floor is 0.7, so 0.75 passes with little margin.
That is the value upstream ships, not one chosen for the probe.

**Routing.** Upstream offers two routes, and the choice decides one verdict.

| Option | Runoff ratio | Mean record peak, mm/day | `routing-conservation`, by the criterion's formula |
| --- | --- | --- | --- |
| `name="xaj"`, `sources` (**used**) | 0.513 | 24.5 | FAIL |
| `name="xaj"`, `sources5mm` | 0.513 | 18.0 | not measured |
| `name="xaj_mz"`, `sources` | 0.513 | 19.5 | pass on all three seeds |

`xaj` is used because it is the function's default, and because the parameter file the
proposal runs is keyed `xaj` and names `CS` and `L`. The proposal's command pairs that
file with `example_simulate_config.yaml`, which says `name: xaj_mz`; from reading
`run_xaj_simulate.py`, that lookup would raise, though this was not run. Under
`xaj_mz`, `CS` and `L` would be reinterpreted as the kernel's `A` and `THETA`.

Other choices:

- **Daily only in the manifest.** PT1D is the step in the proposal. The adapter also
  accepts PT1H, because the resolution probe runs every model at its next finer step.
  Rows are fed as they arrive and never resampled.
- **`dis`.** `mrro` converted through `area_km2`. XAJ's routing has no hydraulics, so
  this says nothing independent of `mrro`.
- **`window_days: full`.** A ten-year record takes about two seconds.

## Zero stores and what is not reported

`snw` and `canopy` are reported as identically zero and named under `zero_stores` in
`run.json`: the model has no snow module, so all precipitation is treated as liquid whatever the
temperature, and it has no interception store. That is what lets the budget probes score it.

`gwex` and `stage` are listed under `not_reported` in `run.json`. The model computes none of these:

| Variable | Why |
| --- | --- |
| `gwex` | no declared exchange with the outside |
| `stage` | no channel geometry or hydraulics; a stage here would be written by the adapter, not the model |
| `hfls`, `hfss`, `hfg` and the rest of the energy terms | no energy balance |

## Result

**18 PASS, 3 FAIL, 16 N/A**, on the gate seeds, under Apptainer (`--runner apptainer`), on the
whole record of every probe (`window_days: full`).

| Probe | Result |
| --- | --- |
| `energy/pet-consistency` | PASS: 0.75 of potential when wet, 0.41 when dry, 0.58 over the record |
| `mass/antecedent-monotonicity` | PASS: 0.21 of the 60 mm storm |
| `mass/area-invariance` | PASS: exactly 0; see below |
| `mass/catchment-closure` | PASS: 2.51% of rain against 5%; see "The budget" |
| `mass/causality` | PASS: 0 before the storm, +0.995 of the added rain after |
| `mass/dry-down` | PASS: 13 mm over 730 dry days |
| `mass/extreme-event-closure` | **FAIL (VIOLATION)**: 133 of 1316 events fail, the worst 0.249 mm on 0.008 mm of rain; see "The budget" |
| `mass/extreme-rain` | PASS: returns 1.00 of the 997 mm added |
| `mass/human-abstraction` | **FAIL (VIOLATION)**: all 380 mm unaccounted for; see below |
| `mass/multi-decadal-drift` | PASS: storage changes 0 mm over the final 1825 days; closure 1.67% |
| `mass/phase-counterfactual` | PASS: exactly 0; see below |
| `mass/precipitation-counterfactual` | PASS: of +1589 mm, 0.122 evaporates, 0.868 runs off, 0.016 stays |
| `mass/resolution-invariance` | PASS: 4.48% of rain against a 10% limit |
| `mass/response-nonnegativity` | PASS: no dip |
| `mass/runoff-bounds` | PASS: 0.55 of the rain, inside [0.12, 1.06] |
| `mass/spinup-cycle-invariance` | PASS: 0.00% between spin-up lengths |
| `mass/steady-state` | PASS: 1.000 mm/day under 2.5 mm/day of rain |
| `mass/time-origin-invariance` | PASS: exactly 0; closure 3.22% |
| `mass/ungauged-basin-closure` | PASS: 3.89% of rain against 5%, the closest this model comes; ET at most 0.750 of potential |
| `mass/warming-response` | PASS: mrro -0.308 and -0.334 per unit of demand |
| `momentum/routing-conservation` | **FAIL (VIOLATION)**: 56.6 mm in the channel against an allowance of 0.2 mm |

Some of these say nothing about the physics. `mass/phase-counterfactual` and
`mass/time-origin-invariance` are exactly 0 because the model never reads temperature
or the date. `mass/area-invariance` is exactly 0 because `mrro` never reads the area.
`state_bounds` on `snw` and `canopy` holds on a column of zeros.

### The failure: `mass/human-abstraction`

Structural, verified. XAJ has no withdrawal term and the adapter never reads `abstr`, so the
natural and irrigated runs are identical: `evspsbl`, `mrro` and storage differ by exactly 0 on all
three seeds, and the whole 380 mm is the residual. AGENTS.md says a model that never reads the
column fails this way, and it does.

### The failure: `momentum/routing-conservation`

The criterion allows `channel <= 15 days * the largest mrro of the preceding 30 days`.

**Mechanism, verified.** `L` is a pure delay of 5 days: release does not begin to
respond until 5 rows after the runoff was generated. When a storm follows a dry spell,
the channel fills while `mrro` is still the flow from before the storm. On seed
1431599080:

| Date | Rain | Into the channel | `mrro` | `channel` |
| --- | --- | --- | --- | --- |
| 2003-01-29 | 0.00 | 0.008 | 0.008 | 0.048 |
| 2003-01-30 | 48.00 | 5.102 | 0.008 | 5.141 |
| 2003-02-01 | 64.65 | 12.779 | 0.008 | 24.730 |
| 2003-02-03 | 48.00 | 18.237 | 0.008 | 56.579 |
| 2003-02-04 | 48.00 | 25.803 | 2.555 | 79.828 |
| 2003-02-10 | 0.00 | 2.141 | 26.040 | 48.193 |

**This is the model's routing, not the adapter's report.** The store is real water:
`channel` balances `inflow - mrro` row by row to 1e-14 mm, stays positive, and drains
(0.29 mm at the end of the record). It never exceeds 0.65 of fifteen days of the peak
of what is actually flowing into it. The verdict follows the lag and nothing else:

| `L` | 0 | 1 | 2 | 3 | 5 |
| --- | --- | --- | --- | --- | --- |
| seed 1431599080 | pass | FAIL, 3 rows | FAIL, 7 rows | FAIL, 11 rows | FAIL, 25 rows |

`L` was not changed, because that would be tuning the model to the probe. The FAIL
stands as the harness reported it.

**For the maintainers.** The criterion reads the feed from released `mrro`, which by
definition has not yet seen water held in a lag. On these records XAJ fails at every
whole-row delay tried, from one row up, although its store balances in each case. By
the same arithmetic any router with such a delay should fail when a storm follows a
dry spell, but only XAJ was tested. Whether that is the intended meaning of the probe
is a question for its authors; no probe was edited.

### A failure that was the adapter's

The first evaluation failed `mass/resolution-invariance` at 23.8% of rain. That was not
the model. Upstream starts both reservoirs at an outflow of 0.1 in mm per row, and at
an hourly step the groundwater store that implies is 158 mm rather than 6.6 mm, which
then drained into a 40-day record. With the initial outflows in the units of the step
it is 4.48%, and the same run with the old value reproduces 23.80% exactly on the seed
the harness named. What is left looks like timing rather than loss: with the same daily
weather spread evenly over the hours, runoff falls by 4.0 to 7.7% of rain and `gw` and
`channel` end higher by 4.7 to 7.8%.

### The budget, and the failure on `mass/extreme-event-closure`

Over a whole record the budget misses by 1.67 to 3.89% of rain on the probes that score it, always
in the same direction, which passes the 5% limit on every one of them. Event by event it fails: 107
to 133 events on each of the five seeds of `mass/extreme-event-closure` miss their allowance of 5%
of the event's rain, or 0.001 mm.

The impervious fraction is booked inconsistently. Runoff from it is `IM * PE`, where
`PE` is rain net of evaporation, and the tension layers take the whole of `PE - R`
with no `1 - IM` on it. No reading of the stores closes:

| Stores read as | Residual, seed 598896396 | What it equals | Rows |
| --- | --- | --- | --- |
| weighted by `1 - IM`, as reported | -202.925 mm | `-IM * (E - P)`: evaporation reported on the impervious share with no water to draw on | dry |
| unweighted | -201.137 mm | `-IM * (PE - R)`: impervious runoff on top of the same rain infiltrating | wet |

The first line predicts each residual to under 0.001 mm on all five seeds, and the
residual scales with `IM`: at `IM = 0.01` it is -0.406%.

**Mechanism on the event probe, verified on its own five cases.** Each of the 1245 to 1316 scored
events' residual equals `-IM * (E - P)` summed over the event's rows where `E > P`, to 6e-5 mm. So a
failing event is a rainy day too small to pay for the day's evaporation: every failing event holds
under 2.7 mm of rain, with medians of 0.22 to 0.31 mm by seed. The worst, on seed 999454640, is
0.008 mm of rain on 2019-08-12 under 4.16 mm of evaporation, giving 0.249 mm. A storm's own rain sets
a far larger allowance, so no big event fails, and the largest residual of any event is 0.33 mm. This
is the same term as the whole-record miss, but only a small part of it: about 96% of that miss, 378
to 399 mm of 393 to 414 mm, falls on rainless rows, which belong to no event. It is XAJ's bookkeeping,
not the adapter's: `evspsbl` is upstream's and is compared with `xaj()` on every run.

 The largest single-row
imbalance is `IM` times the largest evaporation, 0.264 mm at 4.4 mm/day. The reported
`evspsbl` is upstream's and has not been corrected.

Two things that could also lose water were checked. The free-water cap,
`ss = min(ss, SM)`, can destroy water when `FR` falls sharply, but it bound on no row
of the five closure records or of the resolution records at either step. The clip on
the tension layers removes 0.001 mm in ten years.

### The first rows of the routing

Upstream passes the first `L` rows straight through and then lags those same rows
again, releasing them twice, which creates water: 3.7 mm on row 5 of a 365-day spinup
(routing-conservation seed 1431599080), 2.7 to 15.4 mm on the resolution probe's
10-day spinup. Every other scored probe has a 365-day spinup. The `channel` column is
zero while it passes them through and steps up on row `L`, so this shows in the
reported series rather than being hidden. Row `L` is inside every spinup, but on the
resolution probe's 10-day spinup the recession tail of that second release does reach
the scored window, at both steps. Its size there was not measured separately.

### N/A

- INCOMPLETE, outputs it lacks: `exchange-response` (`gwex`), `snowpack-mass-closure`
  (`snm`), `aquifer-recharge-ordering`, `groundwater-datum-invariance`,
  `gw-sw-exchange-consistency` (`gw_sw_exchange`), `froude-regime`,
  `stage-discharge-monotonic`, `uniform-flow-friction-consistency` (`stage`), and
  the six energy probes that need radiative or heat outputs.
- INCOMPATIBLE: `routing-lag-consistency` (does not consume channel lengths; its routing
  constants are not derived from geometry) and `wave-celerity-bounds` (PT1H, `q_in`).

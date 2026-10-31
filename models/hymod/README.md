# hymod

HYMOD is a five-parameter conceptual rainfall-runoff model: a
probability-distributed soil store (Moore 1985) whose excess is split between
three quick-flow linear reservoirs in series and one slow reservoir. Proposed in
issue #133.

## Provenance

The code is `hydromodel/models/hymod.py` from
[OuyangWenyu/hydromodel](https://github.com/OuyangWenyu/hydromodel) v0.4.0, commit
`89d7a8ed1d72ce4fffbbd9897490b089382ecbac`, GPL-3.0. It is run unmodified, not
reimplemented.

The image copies that file and the two modules it imports (`model_config.py`,
`param_utils.py`) byte for byte, under empty package `__init__` files, and installs
numpy, scipy, numba, llvmlite and pyyaml at the versions upstream's `uv.lock`
records. It does not `pip install` the package, for two reasons. The package's
`__init__` imports `hydroutils` only to get a cache path, and that brings in
async-retriever, cartopy and a data stack forward simulation never uses. And
`uv.lock` cannot be installed as shipped, because it resolves `hydrodataset` from a
sibling checkout outside the repository. The base image is pinned by digest.

`hymod()` returns the streamflow series but only the final stores, and it never
returns evaporation: the second item it returns, named `et`, is the effective
rainfall. So:

- `mrro` is what `hymod()` itself returned, divided by the row length.
- The stores come from driving the same module's `excess` and `linres` row by row.
  The two series are compared on every run, the adapter refuses to report if they
  differ by more than 1e-9 mm per row, and `run.json` records the difference. It has
  been exactly 0.
- `evspsbl` is not a residual of the budget, which would close the budget by
  construction. `excess` is called twice on each row, with the row's demand and
  with none, and the difference between the two end states is what its evaporation
  line removed. An error in the runoff accounting would not be absorbed by it.

## Parameters

No calibration exists for a synthetic catchment. Nothing below was chosen by looking
at a verdict.

| Parameter | Value | Source |
| --- | --- | --- |
| `cmax` | `soil_capacity_mm * (bexp + 1)`: 448.6 mm at the usual 320 mm | The catchment. `cmax` is the largest point capacity and the store's mean capacity is `cmax / (bexp + 1)`, so the stated capacity sets it through that. Without `soil_capacity_mm` the published 273.302 is used. |
| `bexp` | 0.402 | published set, below |
| `alpha` | 0.706 | published set, below |
| `ks` | 0.031 per day | published set, below; `1 - (1 - k)^dt` at a finer step |
| `kq` | 0.675 per day | published set, below; `1 - (1 - k)^dt` at a finer step |

The published set is the SCE-UA calibration in the
[`sceua` package documentation](https://sceua.readthedocs.io/en/latest/examples/calibration/):
USGS 01666500, GridMET forcing, NWIS streamflow 2000-2019, objective NSE. All five
values lie inside upstream's parameter ranges.

This is a documentation example, not a peer-reviewed table, and the limits of that
should be plain. The paper upstream's docstring cites (Quan et al. 2015,
doi:10.5194/piahs-368-180-2015) reports only parameter ranges, with `Kq` and `Ks`
as residence times in days, and no best-fit set. Vrugt et al. (2003) for the Leaf
River could not be retrieved, and its values are not quoted from memory. The `sceua`
example cites the same paper and uses the same parameter names; its model code was
not read line by line against upstream's.

Initial stores are upstream's: all empty. The harness's spinup does the rest.

## Choices, and how much each moves

Measured on the gate records of `mass/catchment-closure` (runoff ratio, E/PET, budget
residual) and `mass/resolution-invariance`. These are measurements outside the
harness; only the first row went through it. The default's ratios here are the mean
of five seeds.

| Configuration | Runoff ratio | E/PET | Residual, % of rain | mrro moved PT1D to PT1H, worst % of rain |
| --- | --- | --- | --- | --- |
| **default** | 0.439 | 0.615 | 0.000 | 2.05 |
| `cmax` as published, capacity ignored | 0.501 | 0.548 | 0.000 | 1.28 |
| midpoints of upstream's ranges, capacity used | 0.491 | 0.558 | 0.000 | 2.00 |
| midpoints of upstream's ranges, capacity ignored | 0.580 | 0.460 | 0.000 | 1.99 |
| `ks`, `kq` left per row | 0.439 | 0.615 | 0.000 | 13.31 |

The parameter set moves the partition by up to 0.14 of the rain and does not touch
the budget. Rescaling the two rate constants is what keeps the step dependence small:
left per row, the movement would be over the probe's 10% limit, and that would be the
adapter's units rather than the model.

Other choices:

- **Daily only in the manifest.** PT1D is the step in the proposal. The adapter also
  accepts PT1H, because the resolution probe runs every model at its next finer step.
  Rows are fed as they arrive and never resampled.
- **`mrso`, `gw`, `channel`.** The soil store, the slow reservoir, and the three
  quick reservoirs summed. Every store is reported, so this mapping does not affect
  a budget, only which store a criterion reads.
- **`dis`.** `mrro` converted through `area_km2`. HYMOD has no hydraulics, so this
  says nothing independent of `mrro`.
- **`window_days: full`.** A ten-year record takes about two seconds.

## Zero stores and what is not reported

`snw` and `canopy` are reported as identically zero and named under `zero_stores` in
`run.json`: the model has no snow module, so all precipitation enters the soil store as liquid
whatever the temperature, and it has no interception store. That is what lets the budget probes
score it.

`gwex` and `stage` are listed under `not_reported` in `run.json`. The model computes none of these:

| Variable | Why |
| --- | --- |
| `gwex` | no exchange with the outside |
| `stage` | no channel geometry or hydraulics; a stage here would be written by the adapter, not the model |
| `hfls`, `hfss`, `hfg` and the rest of the energy terms | no energy balance |

## Result

**20 PASS, 1 FAIL, 16 N/A**, on the gate seeds, under Apptainer (`--runner apptainer`), on the
whole record of every probe (`window_days: full`).

| Probe | Result |
| --- | --- |
| `energy/pet-consistency` | PASS: 0.97 of potential when wet, 0.37 when dry, 0.61 over the record |
| `mass/antecedent-monotonicity` | PASS: 0.23 of the 60 mm storm |
| `mass/area-invariance` | PASS: exactly 0; see below |
| `mass/catchment-closure` | PASS: 0.0000% of rain, largest row 1.6e-13 mm; ET at most 0.65 of potential |
| `mass/causality` | PASS: 0 before the storm, +0.747 of the added rain after |
| `mass/dry-down` | PASS: 6 mm over 730 dry days |
| `mass/extreme-event-closure` | PASS: worst event residual 9.2e-14 mm against 0.001 mm |
| `mass/extreme-rain` | PASS: returns 1.00 of the 997 mm added |
| `mass/human-abstraction` | **FAIL (VIOLATION)**: all 380 mm unaccounted for; see below |
| `mass/multi-decadal-drift` | PASS: storage changes 0 mm over the final 1825 days |
| `mass/phase-counterfactual` | PASS: exactly 0; see below |
| `mass/precipitation-counterfactual` | PASS: of +1589 mm, 0.149 evaporates, 0.830 runs off, 0.021 stays |
| `mass/resolution-invariance` | PASS: 2.05% of rain against a 10% limit |
| `mass/response-nonnegativity` | PASS: largest dip -9.9e-06 of the storm |
| `mass/runoff-bounds` | PASS: 0.49 of the rain, inside [0.12, 1.06] |
| `mass/spinup-cycle-invariance` | PASS: 0.00% between spin-up lengths |
| `mass/steady-state` | PASS: 0.909 mm/day under 2.5 mm/day of rain |
| `mass/time-origin-invariance` | PASS: exactly 0; see below |
| `mass/ungauged-basin-closure` | PASS: 0.0000% of rain over 12 catchments; ET at most 0.891 of potential |
| `mass/warming-response` | PASS: mrro -0.203 and -0.258 per unit of demand |
| `momentum/routing-conservation` | PASS: within 0.14 of what a 15-day hydrograph can hold |

### The failure: `mass/human-abstraction`

Structural, verified. HYMOD has no withdrawal term and the adapter never reads `abstr`, so the
natural and irrigated runs are identical: `evspsbl`, `mrro` and storage differ by exactly 0 on all
three seeds, and the whole 380 mm is the residual. AGENTS.md says a model that never reads the
column fails this way, and it does.

### Passes that say little

- `mass/phase-counterfactual` and `mass/time-origin-invariance` are exactly 0 because the model
  never reads temperature or the date. `mass/area-invariance` is exactly 0 because `mrro` never
  reads the area.
- `state_bounds` on `snw` and `canopy` holds on a column of zeros, and on `mrso` by construction:
  `excess` cannot take the store past `cmax / (bexp + 1)`, which is the stated capacity.

### The budget

The closure family now scores the model and it closes on every probe, to about 2e-13 mm on any row.
The harness flags this `suspicious_exact`, as it must for any budget closed to rounding. It is not a
residual: `evspsbl` is what the evaporation line of `excess` removed, measured by a second call
without demand, and the runoff and the stores come from the row-by-row drive that the adapter checks
against `hymod()`. An error in the runoff accounting would not be absorbed by it.

### N/A

- INCOMPLETE, outputs it lacks: `exchange-response` (`gwex`), `snowpack-mass-closure`
  (`snm`), `aquifer-recharge-ordering`, `groundwater-datum-invariance`,
  `gw-sw-exchange-consistency` (`gw_sw_exchange`), `froude-regime`,
  `stage-discharge-monotonic`, `uniform-flow-friction-consistency` (`stage`), and
  the six energy probes that need radiative or heat outputs.
- INCOMPATIBLE: `routing-lag-consistency` (does not consume channel lengths; its routing
  constants are not derived from geometry) and `wave-celerity-bounds` (PT1H, `q_in`).

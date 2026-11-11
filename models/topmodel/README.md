# TOPMODEL under HydroTuring

TOPMODEL (Beven and Kirkby 1979, *Hydrol. Sci. Bull.* 24,
[10.1080/02626667909491834](https://doi.org/10.1080/02626667909491834)) in Wouter Buytaert's R
package, built from [`ICHydro/topmodel@f60752c`](https://github.com/ICHydro/topmodel/tree/f60752ced41ffac47778ffb6616ae01cae56456e),
the commit "update to 0.7.5" of 2022-11-15. The repository has no tags, so this commit is what
"V0.7.5" in [#117](https://github.com/Flood-Lab/HydroTuring/issues/117) resolves to. The package
says it is an implementation of Beven's 1995 Fortran.

Licence: GPL (>= 2). The source the image was built from stays in the image, under
`/opt/topmodel-src`.

## How the model is run

The package's C core runs with no equation changed, but not through the package's R function.
`topmodel()` returns `Q`, `qo`, `qs`, `S`, `fex` and `Ea`, then frees the root-zone deficit, the
unsaturated store and the water in the delay function, and a budget cannot be closed without them.
The only other way to get them would be a residual, which is what the contract forbids. So
`ht_driver.c` is the same loop as `c_topmodel()` over the same unmodified functions, compiled
against the pinned sources, and it copies those stores out before they are freed. `check_driver.R`
runs at image build and fails it if the driver's six outputs differ from `topmodel()`'s by more than
1e-12 on the package's example at 15 minutes, one hour and one day.

## What the adapter reports

| Column | What it is |
| --- | --- |
| `pr` | the forcing, echoed as text |
| `evspsbl` | `Ea`, actual evapotranspiration, taken from the root zone |
| `mrro` | `Q`, the flow at the outlet after the delay function |
| `dis` | `mrro × area_km2 / 86.4`, m3/s |
| `mrso` | the root zone: `Srmax - Srz`, over the core's class weights |
| `gw` | the unsaturated store `Suz` plus the saturated zone, `datum - S_mean` |
| `channel` | the delay function's contents: what was on its way at the start, plus `qo + qs` so far, less `Q` so far |
| `snw`, `canopy` | identically zero; see below |

All fluxes are mm/day at every step and all stores are absolute, in mm. Every row is reported,
spinup included.

**Zero stores.** `snw` and `canopy` are reported as identically zero and named under `zero_stores`
in `run.json`: TOPMODEL has no snow module, so all precipitation is liquid, and no interception
store, rain going straight to the root zone. That is what lets the budget probes score it.

**Not reported, because the model has none of them**, and named under `not_reported` in `run.json`:
`gwex` (nothing crosses the boundary but `Q` and `Ea`), `snm`, every energy flux, `stage`.

**Why `Suz` and the saturated zone are `gw`.** `mrso` is held to `soil_capacity_mm` on the probes
that bound it, and the root zone is the store that capacity describes. `Suz` is capped by the local
deficit rather than by any soil capacity, and sits below the root zone.

**The datum.** TOPMODEL's saturated zone is a deficit below saturation with no floor, and a store
has to be reported as an amount. It is reported on a datum of ten times the parameter `m`, 213 mm, the depth at which
transmissivity has fallen to e⁻¹⁰ of its value at saturation and the store has effectively stopped
draining. Only differences enter a budget. The absolute value does reach two bounds, because
`mass/runoff-bounds` and `mass/dry-down` add the first `gw` to the capacity. It loosens them by
about 100 mm without being able to change their verdicts: dry-down drains 29 mm against a 320 mm
capacity alone, and the runoff ratio of 0.43 is far inside the bound with or without it.

At PT1D `channel` is identically zero on most probes, and that is a computed result: the delay
function's 5 km at 1199 m/h is shorter than a day, so all of it arrives inside the row. So
`momentum/routing-conservation` passes without being tested on the default geometry.

## Parameters and the topographic index

There is no calibration on a synthetic catchment, no DEM behind the probes, and nothing in
static.json describes a topographic index.

**The topographic-index distribution is the package's own example catchment:** `huagrahuma$topidx`,
16 classes. It is not the probe catchment's, because the probe catchment has none, and it comes from
a small Andean catchment far from the probes' 250 km2. **So the index distribution here is an
assumption and not a property of the synthetic catchment.** The verdicts below are for this
catchment's distribution. Its sensitivity is measured further down.

| Parameter | Value | Where it came from |
| --- | --- | --- |
| `qs0`, initial subsurface flow | 3.168e-5 m/h | `huagrahuma$parameters` |
| `lnTe`, log of mean surface transmissivity | −0.5991 (m²/h) | `huagrahuma$parameters` |
| `m`, transmissivity decay | 0.02130 m | `huagrahuma$parameters` |
| `Sr0`, initial root-zone deficit | 0.00263 m | `huagrahuma$parameters` |
| `Srmax`, root-zone capacity | `soil_capacity_mm` (320 mm on most probes); 868 mm when absent | static.json, as `models/flex_lumped` and `models/sacsma_snow17` do; the fallback is `huagrahuma$parameters` |
| `td`, unsaturated delay | 2.85 h/m | `huagrahuma$parameters` |
| `vch`, channel velocity outside the catchment | 1000 m/h | `huagrahuma$parameters`, and it acts on nothing here because the delay table starts at zero distance |
| `vr`, channel velocity inside the catchment | 1199 m/h | `huagrahuma$parameters` |
| `k0`, surface conductivity | 9.361e-3 m/h | `huagrahuma$parameters` |
| `CD`, capillary drive | 0.7236 m | `huagrahuma$parameters` |
| `dt` | 24 or 1 h | the case's step |
| delay table | 0, 1000, 3000, 5000 m at cumulative area 0, 0.1, 0.5, 1 | `huagrahuma$delay` |

The data set's documentation cites Buytaert et al. (2005, *Hydrol. Process.* 19) and does not say
how these values were obtained. Its `delay` is documented with the columns the other way round; the
C code and the package's worked example read column 1 as distance and column 2 as area, and that is
what is used.

**Geometry.** When a case supplies `main_channel_length_km` and `centroid_channel_length_km`, they
replace the delay table: half the area within the centroid's distance and all of it within the main
channel's. TOPMODEL's delay function is natively distance against area, so this is its own input
rather than a mapping onto something else. `area_km2` sets `dis` and nothing else.

**Steps.** The core takes `dt` in hours and scales `lnTe`, `vr` and the unsaturated delay by it
itself, so no parameter is rescaled by the adapter. Rates are turned into metres per row and back,
and no row is resampled. Only PT1D and PT1H are declared, because those are what have been run. The
parameters came from a 15-minute model, and the arithmetic does not restrict it to these two.

## Conservation

The closure family scores the model and passes it: 0.06 to 0.40% of rain over a record, against 5%,
and on `mass/extreme-event-closure` no event of about 1300 per seed misses, the worst at 0.32 of its
allowance. The terms behind that were measured off the harness on the first gate seed of
`mass/catchment-closure`, over 11 years and 8938 mm of rain:

**The model does not close, by a small amount: P − ET − Q − ΔS = −13.64 mm, −0.15% of the rain**, so
it reports slightly more water than it was given. Per step the residual reaches 0.28 mm on storm
days. Two terms account for it exactly, and `run.json` reports them:

- **−16.40 mm from saturation runoff.** The core takes 431.13 mm out of the stores as excess and
  reports 447.54 mm as `qo`. At the class where saturation ends, `core_topmodel.c` writes
  `Aatb_r * ex[j-1] / (ex[j-1] - ex[j]) * ex[j-1] / 2`, using the *node* weight where the fully
  saturated branch uses the *interval* area. The water removed and the water reported are therefore
  two different integrals, and they differ by half the excess times half the difference of
  neighbouring class areas.
- **+2.77 mm from the class weights.** The core's weights sum to 1 − A₀/2 and make the area whole
  only if the first class has zero area, as `make.classes()` builds it. The shipped `huagrahuma$topidx`
  has A₀ = 0.00062, so 0.031% of each row's rain enters no class.

The residual is far inside the closure limit. The discontinuity behind it matters more than its
size, as the next section shows. The closure rows are not flagged `suspicious_exact`, because this
budget does not close to rounding.

## Result on the gate seeds

**19 PASS, 3 FAIL, 15 N/A**, under Apptainer, on the whole record of every probe.

| Probe | Result |
| --- | --- |
| `energy/pet-consistency` | PASS: 1.00 of potential when wet, 0.47 when dry, 0.70 over the record |
| `mass/antecedent-monotonicity` | **FAIL (VIOLATION)**: see below |
| `mass/area-invariance` | PASS: exactly 0; `mrro` never reads the area |
| `mass/catchment-closure` | PASS: 0.19% of rain; ET at most 0.739 of potential |
| `mass/causality` | PASS: 0 before the storm, +1.004 of the added rain after |
| `mass/dry-down` | PASS: 29 mm over 730 dry days |
| `mass/extreme-event-closure` | PASS: worst event 0.070 mm against 0.219 mm allowed; 0 of 1299 fail |
| `mass/extreme-rain` | PASS: returns 1.00 of the 997 mm added |
| `mass/human-abstraction` | **FAIL (VIOLATION)**: all 380 mm unaccounted for |
| `mass/multi-decadal-drift` | PASS: storage changes 0 mm over the final 1825 days; ET 0.995 of potential |
| `mass/phase-counterfactual` | PASS: exactly 0; the model never reads temperature |
| `mass/precipitation-counterfactual` | PASS: of +1737 mm, 0.151 evaporates, 0.848 runs off, 0.004 stays |
| `mass/resolution-invariance` | PASS: 0.83% of rain, on `mrso_end` |
| `mass/response-nonnegativity` | PASS: largest dip -1.9e-05 of the storm; see below for why this depends on the index table |
| `mass/runoff-bounds` | PASS: 0.43 of the rain, inside [0.11, 1.07] |
| `mass/spinup-cycle-invariance` | PASS: 0.00% between spin-up lengths |
| `mass/steady-state` | PASS: 0.993 mm/day under 2.5 mm/day of rain, `gw` still approaching its level |
| `mass/time-origin-invariance` | PASS: exactly 0; the model never reads the date |
| `mass/ungauged-basin-closure` | PASS: 0.40% of rain over 12 catchments; ET at most 0.980 of potential |
| `mass/warming-response` | PASS: mrro -0.315 and -0.379 per unit of demand |
| `momentum/routing-conservation` | PASS, trivially: `channel` is zero at PT1D on this geometry |
| `momentum/routing-lag-consistency` | **FAIL (VIOLATION)**: see below |

Passes that say little: `area-invariance`, `phase-counterfactual`, `time-origin-invariance` and
`routing-conservation`, for the reasons in the table, and `state_bounds` on `snw` and `canopy`, which
holds on a column of zeros.

N/A: 14 INCOMPLETE, the outputs it lacks: `exchange-response` (`gwex`), `aquifer-recharge-ordering`, `groundwater-datum-invariance`, `gw-sw-exchange-consistency`
(`gw_sw_exchange`), `froude-regime`, `stage-discharge-monotonic`,
`uniform-flow-friction-consistency` (`stage`), `snowpack-mass-closure` (`snm`), and the six energy
probes that need radiative or heat outputs. 1
INCOMPATIBLE: `momentum/wave-celerity-bounds`, which needs `q_in` and channel geometry the model does
not take.

### `mass/human-abstraction`

Structural, verified. TOPMODEL has no withdrawal term and the adapter never reads `abstr`, so the
natural and irrigated runs are identical: `evspsbl`, `mrro` and storage differ by exactly 0 on all
three seeds, and the whole 380 mm is the residual. AGENTS.md says a model that never reads the
column fails this way, and it does.

### `mass/antecedent-monotonicity`: FAIL on 2 of 3 seeds

> the storm of 60 mm ran off +0.0 mm more on the wetter catchment; at least 0.02 of the storm is expected

Runoff on the two catchments is identical to every printed digit (3.3944 mm on both). **This is the
root zone's threshold.** Rain fills the root-zone deficit first, and nothing reaches `Suz`, the
saturated zone or `qo` until the deficit is gone, even in a class already saturated at the surface.
The wet run's extra 120 mm therefore stays entirely in the root zone: `mrro` does not differ by
anything while it falls and `gw` differs by exactly 0 at the storm. On the two failing seeds the
storm is smaller than the deficit in both runs (dry 191 and wet 105 mm; dry 181 and wet 92 mm), so
neither produces any runoff from it. On the third the wet deficit is 56.7 mm, under the 60 mm storm,
and it passes with 0.09. The wet run does remember the water, in a root zone about 85 mm fuller,
but runoff cannot see that until the deficit closes.

### `momentum/routing-lag-consistency`: FAIL on all three seeds

> medium lag 0.00 d is outside [0.06, 2.75] d; xlarge lag 10.00 d is outside [1.17, 7.17] d; variant
> 'small' has no measurable storm response

Three separate things, all checked on the probe's cases:

- **No response at all on two seeds.** The same threshold: the 50 mm storm goes into a root-zone deficit
  of about 100 mm, and every variant's `mrro` keeps receding without noticing it.
- **Lag 0 at medium.** TOPMODEL has no hillslope delay, only `d / vr` along the channel, and 74% of a
  39 km channel's kernel arrives within the row at a daily step.
- **Lag 10 d at xlarge.** `vr` of 0.33 m/s comes from a small Andean catchment and is applied to a
  320 km channel, 11.1 days end to end, against an expected 3.3. With 1 m/s that violation goes
  (as flex_lumped's celerity anchors it), and the other two remain, so the verdict does not change.

## How far each choice moves the results

Each alternative was re-run on the gate seeds, on the 14 probes that scored the model before the
zero-store rule. The closure family was not part of this comparison.

| Alternative | What changes |
| --- | --- |
| **`Srmax` = 868 mm**, the package's own, ignoring static.json | **Considerably worse.** Five further FAILs, on `mass/area-invariance`, `causality`, `extreme-rain`, `response-nonnegativity` and `warming-response`: the deeper root zone takes almost every storm, the runoff ratio on the warming case falls to 0.0196 and runoff stops correlating with rain |
| **The index from the package's DEM**, by `sinkfill`, `topidx` and `make.classes` as in its example, 16 classes with A₀ = 0 | **One verdict flips:** `mass/response-nonnegativity` now FAILS, 0.209 mm/day below the control. See below |
| `vr` = 1 m/s | The xlarge lag violation goes. The verdict on that probe does not change, and nothing else moves |
| `k0` very high, no infiltration excess | No verdict moves. Infiltration excess is 0 mm on the closure case at this `k0` |

**So the PASS on `mass/response-nonnegativity` depends on the table, and should be read that way.**
A sweep of one step's rain over a fixed state shows why: at each class boundary `qo` switches branch
and *drops*. On the shipped table, 0.05 mm more rain lowers `qo` by up to 0.32 mm (11.656 to 11.341 mm at
63.45 mm of rain, while the water actually removed as excess keeps rising). On the DEM table it is
0.54 mm. This is the same discontinuity as the −16.40 mm above. Both tables have it, and whether a
probe's storm happens to land across a jump is what differs. It is the model's arithmetic and not the
adapter's: the sweep ran through the driver, whose `qo` the build check holds to `topmodel()`'s.

`Srmax` is the parameter the most results depend on, and it comes from the catchment description
rather than from anything calibrated for TOPMODEL. The rule of taking capacities from static.json makes
this model look better than its own example's value would.

## Packaging

`rocker/r-ver:4.4.2` is pinned by digest, the package by commit, and jsonlite 1.8.9 by SHA-256 from
CRAN's archive. Apt's `git`, `curl` and `ca-certificates` in the throwaway source stage are not
pinned, and nothing from that stage reaches the image except the checked sources. The adapter is
written in R. It uses no random numbers, so the request's seed plays no part, and no statistic of the
record: each row depends only on the rows before it.

One ten-year case of 4015 rows takes about 0.7 s through the container, most of it starting R.

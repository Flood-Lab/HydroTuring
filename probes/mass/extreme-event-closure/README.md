# mass/extreme-event-closure

Implementation for [proposal #31](https://github.com/Flood-Lab/HydroTuring/issues/31).
In the final scored year, rainfall events are overlapped until at least one
1-, 3- or 7-day depth exceeds its fitted 100-year target. The resulting
storms are placed from May 1 with one dry day between them. An independent
100-year synthetic record supplies the targets. Water closure is checked
on every complete wet event, including groups that remain below target.

## Water-budget assertion

For each complete precipitation event B:

```text
R_B = sum_B[(pr + gwex - evspsbl - mrro) * dt_days]
      - (S_after_B - S_before_B)                         [mm]
```

`pr` is supplied precipitation. `gwex` is positive into the catchment and
zero when absent. S includes every reported water store, including `gw`
and `channel`; its endpoints are the states after the preceding step and
after the final wet step. `sbl` is already included in ET, and `dis`
represents runoff in different units, so neither is an additional loss.

Every event must satisfy `abs(R_B) / P_B <= 0.05`, using its own supplied
rainfall P_B, the project's 5% engineering tolerance and no absolute floor.
The largest ratio is reported. Separate events cannot cancel or dilute
one another's errors; signed errors within an event are integrated together.

Events are maximal `pr > 0` runs separated by a dry step (one day here).
They are identified from full forcing and scored after spinup when dry
neighbors are observed on both sides. Events crossing spinup or touching
record edges are skipped and reported. A case with no complete event fails.
All complete events in ordinary and transformed periods are scored,
independent of construction labels or target exceedance.

Only wet steps enter the event budget. Delayed runoff remains in routing
storage at event end; the recession tail need not be appended.
Whole-window `closure`, `state_bounds`, `non_degenerate` and
`forcing_fidelity` are separate checks. Event diagnostics report dates,
budget terms, endpoint stores, residuals and verdict. Non-finite budget data
reaching the criterion fail. Full event diagnostics are retained for the
worst seed, with aggregate results for others.

## Weather, adaptive overlap and storm placement

The climate equations and static attributes match `mass/catchment-closure`:
seasonal rainfall occurrence, gamma wet-day depths, seasonal temperature
with persistent weather variability, and temperature/day-length PET.
Snow remains possible. Each seed spawns two PCG64 streams through NumPy
`SeedSequence`: child 0 for model weather, child 1 for calibration.
Calibration does not consume the model-weather stream.
See [NumPy's documentation](https://numpy.org/doc/stable/reference/random/parallel.html).

The 4,015 daily rows contain 365 spinup days, nine ordinary 365-day scored
blocks and the final 365-day transformation block. These are row-count
blocks, not calendar years. Only rows `[3650:4015)` are transformed;
the first 3,650 rows remain unchanged. The model runs continuously with
all 3,650 scored days, even when it requests a shorter window.

Within the final block:

1. Take complete original wet events chronologically. Align their starts
   and add their daily depths, padding shorter events with zeros.
2. After adding each whole event, calculate the group's maximum 1-, 3- and
   7-day depths. Close the group at the first strict exceedance of **any**
   corresponding Q100 depth, then restart with the next unused source.
   A first source already above target closes as a singleton.
3. If sources run out before a target is exceeded, combine all remaining
   events into one final unmet group, including a singleton.
4. Pack all groups, including the unmet group, in chronological order.
   Start on the actual calendar date May 1 within the final block, then place each next group after
   exactly one dry day: `next_start = previous_stop + 1`, with exclusive stops.

The anchor is 2010-05-01, row 3773 (zero-based), in the current record.
It is found from the dates, not a fixed day offset within a 365-day block.

For these construction depths, rain is zero outside the isolated group's
hyetograph. A duration longer than the event therefore receives its total
rainfall. Neighboring events cannot help it exceed a target. Comparison
uses six-decimal rainfall depths and strict `>`; equality does not stop
accumulation. The OR rule exceeds at least one marginal duration-specific
100-year level; it does not establish a joint 100-year event frequency.

Overlap first forms each group at its first source's original start.
Packing copies these hyetographs, clears their old positions, and writes
them at the new dates within the final block. Every eligible source is used once and whole;
rainfall is not rescaled. Boundary-crossing and unfinished record-edge
events remain unchanged. Total rain, record length, dates, temperature and
PET are preserved. Packing changes the gaps between storms without changing
their isolated depths or target status; this clustering is an additional
perturbation. If rainfall is insufficient, target exceedance is not guaranteed.
Temperature and PET retain their values on each calendar date; temperatures
are not manually warmed and negative-reference parameters are not tuned.

Host-only `_event_id`, `_event_start`, `_event_end` and `_regime`
describe construction. They do not determine closure eligibility.

## Calibration and rainfall reports

Each seed's child 1 generates 100 synthetic 365-day years under the same
rainfall law, plus six prefix days. Annual maxima at 1-, 3- and 7-day
durations include dry gaps and year-crossing windows, assigned by ending
day. A separate GEV fit uses unbiased sample L-moments at each duration,
Hosking's shape `k = -xi`, and 100 deterministic bisection steps.
No extra dependency or data download is needed. See the
[Hosking lmom reference](https://cran.r-universe.dev/lmom/doc/manual.html)
and [NOAA Atlas 14 methods](https://www.weather.gov/media/owp/oh/hdsc/docs/Atlas14_Volume12.pdf)
for context; this probe does not implement NOAA's regional frequency analysis.

Q100 supplies the construction targets from the 100-year calibration record.
Reports provide return levels at 2, 5, 10, 20, 50, 100, 200 and 500 years, and estimated return
periods for constructed groups and the modified year's maximum depths.
The latter use windows ending anywhere in the final block, including up to
six preceding days, dry gaps or separate groups. With one-day gaps, a
3- or 7-day window can include adjacent storms and greatly exceed the
isolated-group statistic that stops construction.

Estimates describe the original synthetic climate, not external training
limits or the occurrence frequency of deliberately constructed storms.
Values above 100 years are flagged as extrapolations; no confidence
intervals are computed. Clustered 7-day maxima produce extremely large,
poorly constrained tail estimates or exceed the fitted upper endpoint;
these are not precise recurrence claims.
Upper-endpoint or unrepresentable-tail exceedances receive null and a status.
Independent DDF crossings are flagged without smoothing; seed 2013358868
has crossings at 200 and 500 years.

## Current gate results

All four physical references (`reference_bucket`, `flex_lumped`,
`flex_topo`, `sacsma_snow17`) pass every criterion on all five seeds.
Construction counts and the negative reference's event-closure outcomes
are:

| Seed | Groups above target / total | Packed first–last wet date | Negative reference |
| --- | ---: | --- | --- |
| 492502526 | 5 / 6 | 2010-05-01–2010-05-24 | FAIL |
| 999454640 | 6 / 7 | 2010-05-01–2010-05-20 | FAIL |
| 1506406754 | 5 / 6 | 2010-05-01–2010-05-18 | FAIL |
| 2013358868 | 6 / 7 | 2010-05-01–2010-05-20 | FAIL |
| 372827335 | 6 / 7 | 2010-05-01–2010-05-24 | FAIL |

Of 33 constructed groups, 28 exceed the target and consume 6–18 source
events each. Each seed has one final unmet group, which remains in the
closure test. Packed sequences span 18–24 days including their dry gaps.
These counts describe the five gate seeds, not a guarantee for arbitrary rainfall.
Their constructed wet days have temperatures from 5.586 to 23.196 °C;
none is below freezing. This observation does not guarantee snow-free
conditions for arbitrary seeds.

The sole negative, `reference_in_sample`, fails on all five seeds;
its aggregate failure is solely `event_water_closure`. It silently removes
surface runoff above 55 mm/day rainfall, capped by available surface runoff.
A large rainfall return period alone does not guarantee a budget failure.
The defect requires high rainfall and surface overflow on the same day.
Failures in either ordinary or transformed periods are valid; the criterion
does not require the transformation to cause the defect.

As a supplementary check of the existing `state_bounds` safeguard,
`reference_cheater` was also evaluated on all five gate seeds. It invents
storage to balance independently chosen fluxes: both whole-window and
event water closure pass on all five seeds, while `state_bounds` fails on
all five. This reuses the project's existing counterexample to show why
closure alone does not establish physical realism. `reference_in_sample`
remains the sole declared negative in this probe's acceptance gate.

## Reproduction

From the repository root:

```bash
python probes/mass/extreme-event-closure/generate.py --gate-seeds --output-dir runs/extreme-event-closure/may-ddf-100
pytest -q tests/test_extreme_event_generator.py tests/test_extreme_event_ddf.py
pytest -q tests/test_event_water_closure.py tests/test_event_water_references.py
ht validate
ht gate --probe mass/extreme-event-closure
ht run --model reference_in_sample --probe mass/extreme-event-closure --gate-seeds --json runs/extreme-event-closure/may-ddf-100/negative-reference.json
ht run --model reference_cheater --probe mass/extreme-event-closure --gate-seeds --json runs/extreme-event-closure/may-ddf-100/cheater-delivery-audit.json
```

The generator CLI writes:

- `report.json`: full calibration, thresholds, `anchor_row`/`anchor_time`, actual placement bounds,
  pre-packing `original_start`/`original_stop` and unchanged source-event indices.
- `events.csv`: each group's actual dates, duration-specific depth, threshold, source
  count, exceedance and stopping status.
- `modified_year.csv`: final-year maximum depths and estimated return periods.
- `ddf.csv`: return levels and crossing/extrapolation flags.

For one seed, use `--seed 20260912` instead of `--gate-seeds`.
Diagnostics reside in `forcing.attrs["rainfall_diagnostics"]` on the host
and are absent from model forcing, static inputs and requests. Standard
`ht run --json` reports budget scores but does not export these rainfall
attributes; use the generator CLI for the rainfall report.

## Runtime

The probe PR template requires a run under one minute on a two-core runner.
`max_runtime_s: 60` is a separate timeout for each adapter invocation.
The complete `ht gate --probe mass/extreme-event-closure` workload was timed
three times on Windows 11, an i9-13900K and Python 3.12.14, restricting the
benchmark process and its descendants to logical CPUs 0 and 1 and capping
OMP, OpenBLAS, MKL and NumExpr threads at two.

The gate includes imports, 25 case generations with DDF calibration,
25 adapter calls (four physical references and one negative, each on five
seeds), I/O and scoring. All three runs passed in **15.739, 15.479 and
15.636 seconds**, with the slowest below 60 seconds. Generating just the five
seeds took 0.052–0.054 seconds, or 0.598–0.599 seconds including fresh-process
startup and imports. These measurements establish the local two-logical-CPU
result; they do not establish the runtime on a GitHub CI runner.

## Contribution boundary

The planned PR contains this probe's three core files, the
`event_water_closure` implementation and registration, and focused tests.
Calibration stays in `generate.py`; no new reference model or dependency
is introduced.

Authorship is also recorded in the repository-wide CONTRIBUTORS and CITATION
files. By the contributor's requested division of work, repository-wide README,
ROADMAP, site and model-archive updates will be
requested from the maintainer at integration time. Documentation, reporting
and archive checks may remain failing until those updates are completed.
This probe's README and author information remain part of the contribution.

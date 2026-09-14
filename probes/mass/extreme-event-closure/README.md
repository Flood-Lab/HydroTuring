# mass/extreme-event-closure

Implementation for [proposal #31](https://github.com/Flood-Lab/HydroTuring/issues/31).
The intuition is to test whether water-balance closure generalizes to extreme
rainfall that a model may not have encountered during training. Closure
learned as a statistical regularity may fail under unfamiliar forcing.
Extremeness is defined against the synthetic baseline climate, without
assuming knowledge of a submitted model's actual training distribution.

Of twenty scored years, select the year whose total rainfall is closest to
the median annual total. Overlap its rainfall events until at least one
1-, 3- or 7-day depth exceeds its fitted 100-year target. Place the resulting
storms from May 1 with one dry day between them. An independent 100-year
synthetic record supplies the targets. Water closure is checked on every
complete wet event across all twenty years, including groups below target.

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

Every event must satisfy `abs(R_B) <= max(0.05 * P_B, 0.001 mm)`, using its
own supplied rainfall P_B. The 0.001 mm absolute allowance protects tiny
events against finite output precision; the project's 5% engineering
tolerance dominates for P_B >= 0.02 mm. No rainfall denominator is replaced
and no small event is skipped. The allowance is a numerical precision budget,
not a claim that every adapter precision is supported. Five-seed regressions
round the conservative bucket's reported fluxes and stores to four decimals
or float32, while retaining the exact precipitation echo. All pass with this
allowance; four-decimal output fails on all five seeds with the allowance off.

The scalar score is the largest `abs(R_B) / allowed_residual_mm`, with a
threshold of 1. Event diagnostics retain the actual precipitation, signed
residual, rain-normalized residual and allowed millimetres. The worst event
is selected by the fraction of its allowance used, so a rescued tiny event
cannot hide a genuine failure elsewhere. Separate events cannot cancel or dilute
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
reaching the criterion fail. Every event is evaluated; the report retains
the worst event, the twenty worst failing events (or fewer), total and failed
counts, the omitted-failure count, and p50/p90/p95/p99/max summaries of event
rainfall, relative residual and residual/allowance. Ties retain chronological
order. These diagnostics are retained for the worst seed by the harness.

Like whole-window closure, `event_water_closure` checks consistency between
supplied precipitation and reported fluxes and stores; it cannot independently
establish that the reported storage is physically authentic. A model can
fabricate storage changes to absorb a budget residual. The inherited
`state_bounds` configuration limits soil and canopy storage and requires
nonnegative snow storage, but sets no snow upper bound or bounds on `gw`
and `channel`, even though these stores enter the budget when reported.
It therefore does not rule out storage fabrication. This probe targets
event-scale budget errors illustrated by `reference_in_sample`; broader
storage validation would require separately justified physical constraints.

## Weather, adaptive overlap and storm placement

The climate equations and static attributes match `mass/catchment-closure`:
seasonal rainfall occurrence, gamma wet-day depths, seasonal temperature
with persistent weather variability, and temperature/day-length PET.
Snow remains possible. Each seed spawns two PCG64 streams through NumPy
`SeedSequence`: child 0 for model weather, child 1 for calibration.
Calibration does not consume the model-weather stream.
See [NumPy's documentation](https://numpy.org/doc/stable/reference/random/parallel.html).

The 7,665 daily rows contain 365 spinup days and twenty 365-day scored
blocks. These are row-count blocks, not calendar years. Compute each block's
original precipitation total and the median of the twenty totals, then
select the block closest to that median. For an even number of blocks,
the median is the mean of the two middle totals. Equal distances are broken
by choosing the earliest chronological block; this can select either middle
total, not necessarily the lower one. Totals and distances are compared in
integer micro-mm at the forcing's six-decimal precision so numerical rounding
does not resolve a mathematical tie.

Only the selected block's rainfall is transformed. Spinup and the other
nineteen blocks remain unchanged. The model runs continuously through all
7,300 scored days, even when it requests a shorter window. The selected
year's share of total scored rainfall is reported; it is not assumed to be
5%. Extending the record and selecting a median-wet year do not by themselves
guarantee that a defective model passes whole-window closure.

Within the selected block:

1. Take complete original wet events chronologically. Align their starts
   and add their daily depths, padding shorter events with zeros.
2. After adding each whole event, calculate the group's maximum 1-, 3- and
   7-day depths. Close the group at the first strict exceedance of **any**
   corresponding Q100 depth, then restart with the next unused source.
   A first source already above target closes as a singleton.
3. If sources run out before a target is exceeded, combine all remaining
   events into one final unmet group, including a singleton.
4. Pack all groups, including the unmet group, in chronological order.
   Start on the actual calendar date May 1 within the selected block, then
   place each next group after exactly one dry day:
   `next_start = previous_stop + 1`, with exclusive stops.

The selected year and its May 1 anchor can differ between seeds. The anchor
is found from dates, not a fixed day offset within a 365-day block.

For these construction depths, rain is zero outside the isolated group's
hyetograph. A duration longer than the event therefore receives its total
rainfall. Neighboring events cannot help it exceed a target. Comparison
uses six-decimal rainfall depths and strict `>`; equality does not stop
accumulation. The OR rule exceeds at least one marginal duration-specific
100-year level; it does not establish a joint 100-year event frequency.

Overlap first forms each group at its first source's original start.
Packing copies these hyetographs, clears their old positions, and writes
them at the new dates within the selected block. Every eligible source is used once and whole;
rainfall is not rescaled. Boundary-crossing and unfinished record-edge
events remain unchanged. Total rain, record length, dates, temperature and
PET are preserved. Packing changes the gaps between storms without changing
their isolated depths or target status; this clustering is an additional
perturbation. If rainfall is insufficient, target exceedance is not guaranteed.
Temperature and PET retain their values on each calendar date. May placement
does not impose a temperature threshold or guarantee snow-free conditions.

Construction spans and selected-year information live only in
`forcing.attrs["rainfall_diagnostics"]`. The forcing has exactly `time`,
`pr`, `tas` and `pet`; there are no underscore annotation columns to impose
an additional window restriction. The explicit 7,300-day minimum remains
because both the twenty-year aggregate and all events are part of this case.

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
The latter use windows ending anywhere in the selected block, including up to
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
Independent DDF crossings are flagged without smoothing.

## Negative reference

The probe reuses the existing
[`reference_in_sample`](../../../models/reference_in_sample/ht_adapter.py).
It runs the bucket model's snow, interception, soil and runoff processes,
then removes an unreported quantity from surface runoff:

```text
loss = min(surface_runoff, 0.55 * max(P - 55, 0))
```

Daily precipitation P is in mm/day and echoed unchanged. The lost water is
neither stored nor declared as an exchange. The fixed 55 mm/day cutoff and
55% excess-loss fraction are the existing reference's illustrative defect
parameters; this probe does not change them. They stand in for unfamiliar
forcing rather than identify a measured training boundary. The adapter
reads no event labels, selected-year information, calibration statistics or
future forcing.

The defect needs high precipitation and available surface runoff on the
same day; a large rainfall return period alone does not guarantee a loss.
Such losses can pass the whole-window 5% tolerance when normalized by twenty
years of rainfall, while exceeding 5% of an individual event's rainfall.
This is passage within an engineering tolerance, not exact whole-window
conservation. Losses can occur in any year, so both outcomes must be
measured on the actual cases. Errors in ordinary events are valid failures
too; the criterion does not require the rearrangement to cause the defect.

This is an existing process model with an injected defect, not a newly
constructed model without hydrological physics. It demonstrates a
limitation of aggregate closure and the additional discrimination from an
event budget, without claiming to reproduce a trained model's extrapolation.

## Current gate results

`reference_bucket`, `flex_lumped`, `flex_topo` and `sacsma_snow17` pass
every criterion on all five fixed gate seeds. The sole declared negative,
`reference_in_sample`, passes whole-window closure, forcing fidelity,
state bounds and non-degeneracy, but fails `event_water_closure` on all five.
Both residual columns below are percentages of their respective rainfall
denominators. The worst negative-reference events are large enough that the
5% relative tolerance, rather than the absolute allowance, applies.

| Seed | Selected scored year | Selected share of 20-year rain (%) | Groups above Q100 / total | Whole-window residual (%) | Worst event residual (%) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 492502526 | 11 | 5.120 | 6 / 7 | 0.490 | 26.276 |
| 999454640 | 16 | 4.920 | 6 / 7 | 0.626 | 26.476 |
| 1506406754 | 3 | 5.116 | 7 / 7 | 0.609 | 24.379 |
| 2013358868 | 4 | 5.041 | 6 / 7 | 0.772 | 27.686 |
| 372827335 | 1 | 5.100 | 6 / 7 | 0.611 | 32.108 |

The 35 constructed groups include 31 above-target groups and four unmet
tails. All storms start on May 1 in their selected year; the sequences end
between May 23 and May 27. Constructed wet-day temperatures are positive
on these five seeds, with a minimum of 3.783 degrees C. These are observed
properties of the gate cases, not guarantees for every seed.

Selected-year losses are approximately 81--123 mm, while the other
nineteen years contribute approximately 0--13 mm. Both contributions enter
the measured whole-window residual. The result demonstrates dilution under
the aggregate 5% tolerance, not exact conservation or a proof that selecting
a median-wet year always hides a leak. A five-seed regression test requires
aggregate closure and the other safeguards to pass on the same runs whose
event closure fails; the gate alone only requires the declared failure.

### Ordinary-weather comparison

The same seed's unmodified `generate_baseline` weather produces the following
worst event residuals for `reference_in_sample`. Aggregate closure passes
on all five ordinary runs as well as all five transformed runs.

| Seed | Ordinary event residual (%) | Ordinary event verdict | Transformed event verdict |
| --- | ---: | --- | --- |
| 492502526 | 3.220 | PASS | FAIL |
| 999454640 | < 0.000001 | PASS | FAIL |
| 1506406754 | 8.927 | FAIL | FAIL |
| 2013358868 | 13.802 | FAIL | FAIL |
| 372827335 | < 0.000001 | PASS | FAIL |

Thus the transformation exposes the defect on three additional seeds and
raises the worst event residual on every seed. It does not establish that
the reference fails exclusively under constructed weather or reproduce a
trained model's extrapolation. Both ordinary failures are one-day events:

| Seed | Date | Rain (mm) | Residual (mm) | 1-day return period (years) | 3-day return period | 7-day return period |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1506406754 | 2002-03-22 | 65.657204 | 5.861462 | 6.20 | 2.23 | 1.13 |
| 2013358868 | 2010-05-14 | 73.425792 | 10.134186 | 5.43 | 2.46 | 1.24 |

These estimates use the same independent 100-year calibration and isolated
event convention as construction: zero rain outside each one-day event,
so its 1/3/7-day depths are equal. They are marginal fitted return periods,
not a joint event frequency or estimates of surrounding multi-storm windows.
Neither event meets the construction's Q100 target. The fixed 55 mm cutoff
corresponds to 1-day return periods of 2.45 and 1.81 years for these two seeds;
their fitted 1-day Q50 values are 84.63 and 100.55 mm, respectively. No change
to the existing reference's cutoff is made in this contribution.

The comparison is pinned in `tests/test_event_water_references.py`.
The frequency estimates can be reproduced with `generate_baseline(seed)`,
`calibrate_ddf(seed)`, `event_duration_maxima(event_rain)` and
`gev_return_period(calibration["fits"][str(duration)], depth)` from `generate.py`.

## Reproduction

From the repository root:

```bash
python probes/mass/extreme-event-closure/generate.py --gate-seeds --output-dir runs/extreme-event-closure/median20-ddf-100
pytest -q tests/test_extreme_event_generator.py tests/test_extreme_event_ddf.py
pytest -q tests/test_event_water_closure.py tests/test_event_water_references.py
ht validate
ht gate --probe mass/extreme-event-closure
ht run --model reference_in_sample --probe mass/extreme-event-closure --gate-seeds --json runs/extreme-event-closure/median20-ddf-100/negative-reference.json
```

The generator CLI writes:

- `report.json`: full calibration, annual rainfall totals, median and selected
  year, its share of total scored rain, thresholds, `anchor_row`/`anchor_time`,
  actual placement bounds, pre-packing `original_start`/`original_stop` and
  unchanged source-event indices.
- `events.csv`: each group's actual dates, duration-specific depth, threshold, source
  count, exceedance and stopping status.
- `modified_year.csv`: selected-year maximum depths and estimated return periods.
- `ddf.csv`: return levels and crossing/extrapolation flags.

For one seed, use `--seed 20260912` instead of `--gate-seeds`.
Diagnostics reside in `forcing.attrs["rainfall_diagnostics"]` on the host
and are absent from model forcing, static inputs and requests. Standard
`ht run --json` reports budget scores but does not export these rainfall
attributes; use the generator CLI for the rainfall report.

## Runtime

The probe PR template requires a run under one minute on a two-core runner.
`max_runtime_s: 300` is a separate timeout for each adapter invocation;
the longer limit gives containerized models room to process 7,665 rows.
`max_output_mb: 4.0` preserves output-file headroom at this record length.
The complete twenty-year `ht gate --probe mass/extreme-event-closure`
workload was measured three times on Windows 11, an i9-13900K and Python
3.12.14. The benchmark process and its descendants were restricted to
logical CPUs 0 and 1, with OMP, OpenBLAS, MKL and NumExpr thread caps of two.
The initial revision's three gates passed in **29.578, 29.177 and 29.272 seconds**.
After adding the numerical floor and compact diagnostics, a repeat under the
same restrictions passed in **29.879 seconds**; generating all five seeds
took 0.068 seconds excluding imports. Source hashes and child affinity were
verified again for this revision.

The timing includes imports, 25 case generations with DDF calibration,
25 adapter calls (four physical references and one negative, each on five
seeds), I/O and scoring. Generating just the five seeds took
0.067--0.078 seconds, or 0.601--0.648 seconds including process startup and
imports. Child affinity and unchanged source hashes were verified. These
are local two-logical-CPU measurements; GitHub runner timing remains to be
verified rather than inferred from the adapter timeout.

## Integration and remaining container evaluations

The PR contains this probe's three core files, the
`event_water_closure` implementation and registration, and focused tests.
Calibration stays in `generate.py`; no new reference model or dependency
is introduced. Generated rainfall and detailed JSON reports are not committed; measured
model summaries are appended to the shared CSV archive.

Authorship is recorded in CONTRIBUTORS and CITATION. README, ROADMAP,
the three site languages, the flowchart, the criterion guide and pinned
report headlines are updated with the contribution. The archive includes
new five-seed results for the three physical models and an honest N/A for
the discharge-only Google model. This Windows development machine has
neither Docker nor WSL. Five submitted-model evaluations therefore remain
pending; their previous standings are explicitly marked as incomplete for
the expanded suite rather than filled with guessed outcomes.

On a Docker-capable machine, from the repository root after checking out this
PR, install the project (`python -m pip install -e '.[dev]'`) and run:

```bash
ht run --model dhbv2 --probe mass/extreme-event-closure --gate-seeds --csv models/result.csv
ht run --model wflow_sbm --probe mass/extreme-event-closure --gate-seeds --csv models/result.csv
ht run --model summa --probe mass/extreme-event-closure --gate-seeds --csv models/result.csv
ht run --model cwatm --probe mass/extreme-event-closure --gate-seeds --csv models/result.csv
ht run --model lisflood --probe mass/extreme-event-closure --gate-seeds --csv models/result.csv
```

The runner builds each image from its committed Dockerfile when necessary.
Inspect any ERROR before treating it as an archived model outcome. Once the
five results are available, update their README/site standings from the
archive and run `pytest -q tests/test_docs_in_sync.py`. Until then, the two
standing-completeness checks and the archive-completeness check remain open.
This probe's README and author information remain part of the contribution.

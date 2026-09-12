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
denominators; the tolerance is 5% in each case.

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

The 217 focused tests pass. The full local suite reports 623 passed,
5 skipped and 7 known documentation/archive synchronization failures.

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
`max_runtime_s: 60` is a separate timeout for each adapter invocation.
The complete twenty-year `ht gate --probe mass/extreme-event-closure`
workload was measured three times on Windows 11, an i9-13900K and Python
3.12.14. The benchmark process and its descendants were restricted to
logical CPUs 0 and 1, with OMP, OpenBLAS, MKL and NumExpr thread caps of two.
All three gates passed in **29.578, 29.177 and 29.272 seconds**.

The timing includes imports, 25 case generations with DDF calibration,
25 adapter calls (four physical references and one negative, each on five
seeds), I/O and scoring. Generating just the five seeds took
0.067--0.078 seconds, or 0.601--0.648 seconds including process startup and
imports. Child affinity and unchanged source hashes were verified. These
are local two-logical-CPU measurements; GitHub runner timing remains to be
verified rather than inferred from the adapter timeout.

## Contribution boundary

The planned PR contains this probe's three core files, the
`event_water_closure` implementation and registration, and focused tests.
Calibration stays in `generate.py`; no new reference model or dependency
is introduced. Generated rainfall and result files are not committed.

Authorship is also recorded in the repository-wide CONTRIBUTORS and CITATION
files. By the contributor's requested division of work, repository-wide README,
ROADMAP, site and model-archive updates will be
requested from the maintainer at integration time. Documentation, reporting
and archive checks may remain failing until those updates are completed.
This probe's README and author information remain part of the contribution.

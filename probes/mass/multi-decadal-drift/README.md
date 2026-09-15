# Multi-decadal storage drift

A budget can close exactly while a model hides a persistent runoff reporting
error in invented storage. This probe asks whether that storage stays physical
over fifty years and whether total reported storage keeps drifting under a
repeated climate.

## Case

A seed generates a five-year block of warm daily weather. One wet day is
sampled within each consecutive three-day group (including the final partial
group). The seasonal phase, rain depth, and daily temperature/PET anomalies
are seeded. Temperature stays between 9.5 and 20.5 degrees C and PET between
0.65 and 1.35 mm/day. Repeat that
block for five years of spinup and fifty scored years (20,075 input rows).
There is no imposed secular trend in the block totals. The case is periodic,
not realistic stochastic weather or a statistically stationary process.

Model years contain 365 days. Gregorian timestamps advance continuously;
seasonality follows the forcing index, not Gregorian day-of-year. The warm
case avoids snow accumulation and isolates soil/canopy capacity. Each seed
generates different rain dates, rain depths, temperatures and PET, reproducibly,
without a committed forcing dataset.

The minimum evaluation window is 18,250 days, preventing the default submitted
model flood-event window from discarding the long horizon. Model states run
continuously through spinup; the reporter's bias is not reset at scoring time.

## Criteria and units

At the daily step, the budget residual is

`r[t] = P[t] + gwex[t] - ET[t] - Q[t] - (S[t] - S[t-1])` in mm.

`S` includes every reported water store, including groundwater and routing
storage where present. `gwex` includes declared external exchange; the supplied
precipitation remains the denominator. Closure requires
`abs(sum(r)) / sum(P) <= 0.05`. Rain is positive in aggregate, so no denominator
floor is needed. The initial storage is the last spinup state, not the first
scored state.

The first discriminating criterion is `state_bounds`: soil stays in [0, 320]
mm, canopy in [0, 2] mm, and the remaining reported stores stay nonnegative.
The existing numerical slack is 1e-6 mm. Extra groundwater is not folded into
the soil capacity. FLEX uses case capacities (area-weighted for FLEX-Topo);
SAC-SMA scales its upper/tension-zone capacities and reports lower free water
as groundwater. No universal finite groundwater capacity is imposed here.

The second discriminating criterion is `total_storage_drift`. It sums every
reported water store and compares the final repeated five-year block endpoint
with the storage immediately before that block. The pass condition is

`abs(delta_S_final_block) <= 2e-4 * P_final_block + 1e-6 mm`.

The allowance is relative to supplied precipitation integrated over that final
block. For roughly 3,300 mm of rain it is about 0.66 mm, allowing small late
equilibration while catching a persistent reporting bias above about
0.0004 mm/day in any reported store.

The existing ET plausibility, non-degeneracy, and forcing-fidelity criteria
guard against trivial compensating partitions and altered drivers.

## Discrimination

`reference_slow_drift` runs the exact bucket internally, but reports
`Q' = Q - epsilon` and `soil' = soil + epsilon * elapsed_days`, with
`epsilon = 0.012 mm/day`. The deficit and storage increment cancel in closure.
The bias is independent of seed, probe identity, and requested run length.
It does not feed back into the bucket's internal state.

Over fifty scored years, the added storage grows by 219 mm, in addition to
21.9 mm accumulated during spinup. Over a ten-year scored prefix it grows by
only 43.8 mm. Tests require the short prefix to stay within bounds but fail
the total-storage trend check, and the long run to fail both `state_bounds`
and `total_storage_drift`, while both close their budgets and the negative
control's runoff remains nonnegative in this case.

The gate requires reference_bucket, flex_lumped, flex_topo, and sacsma_snow17
to pass. `reference_gw_slow_drift` parks the same kind of compensating bias in
groundwater, where no finite capacity is imposed, and must fail
`total_storage_drift`. Merely extending the precipitation-normalized closure
check does not amplify a constant fractional leak; endpoint storage drift under
repeated forcing is the added scientific target.

The final forcing strengthens the seasonal rain and PET cycles relative to
the issue's preliminary experiment, so correct models satisfy the existing
anti-degeneracy checks without lowering their thresholds. Five- and ten-year
warmups are both tested. FLEX-Topo still adjusts slightly after five years;
the probe permits this bounded adjustment, and the tests compare late-cycle
means rather than requiring the first scored cycle to be fully equilibrated.

## Limits

This probe cannot detect every slow leak. A model can fabricate bounded states
and compensating fluxes below the final-block tolerance and still pass; passing
is not proof of physical understanding. Legitimate equilibration, seasonal
storage, and declared groundwater exchange are not themselves failures. The
trend criterion deliberately scores only the final repeated block rather than
requiring the whole sequence of block changes to be monotone or zero, so very
slow transients below the stated allowance are permitted. The detection floor
is case-dependent, not a universal minimum detectable trend.

Moving the reporting bias into `gw` or `channel` no longer escapes the probe,
because `total_storage_drift` sums every reported water store. Adding an
arbitrary groundwater capacity is still avoided.
The negative reporter can produce negative runoff in drier, unrelated cases;
the tests explicitly exclude that alternative failure on this forcing.

## Reproduce

Re-run the gate after changes to forcing or upstream model implementations.
The focused tests cover independent seeds, reproducibility, prefix consistency,
longer spinup and the documented detection limits. Results from the original
PR do not establish results for the revised generator.

The revised gate passes all five fixed seeds for all four physical baselines;
the soil slow reporter fails `state_bounds` and `total_storage_drift`, and the
groundwater slow reporter fails `total_storage_drift`. The three evaluated
physical models have regenerated PASS archive rows, and Google has N/A
(INCOMPLETE).

The previous local dhbv2 daemon failure has been removed from the published
archive and replaced with a real five-seed Docker evaluation: FAIL (VIOLATION)
on `state_bounds` only. The worst seed reports soil storage of 488-501 mm
against the 320 mm capacity. This is an existing capacity mismatch, not evidence
that dhbv2 develops a secular drift; its other four criteria pass.

Each container emitted all 20,075 rows within the 300-second run limit, under
the manifest's two-CPU/four-GB limits on local Docker Desktop. Adapter-reported
wall times ranged from 66.65 to 121.16 seconds; result CSVs were
3,090,736-3,110,230 bytes, below 10 MB. These are local measurements, not a
GitHub-hosted runner benchmark, and exclude image construction.

Evaluation of the new upstream submitted models (wflow_sbm, summa, cwatm and
lisflood) remains required before merge. Missing evaluations are not filled
with fabricated PASS, FAIL or N/A rows.

```sh
ht validate
ht gate --probe mass/multi-decadal-drift
pytest -q tests/test_multi_decadal_drift.py
```

Author: Bing Li (@hiter-joe). Proposal: https://github.com/Flood-Lab/HydroTuring/issues/13.

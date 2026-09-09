# Multi-decadal storage drift

A budget can close exactly while a model hides a persistent runoff reporting
error in invented storage. This probe asks whether that storage stays physical
over fifty years, not whether a fitted storage trend is exactly zero.

## Case

A seed generates a five-year block of warm daily weather. Rain arrives every
third day with a seasonal envelope and seeded amplitude; temperature remains
between 10 and 20 degrees C, and PET between 0.7 and 1.3 mm/day. Repeat that
block for five years of spinup and fifty scored years (20,075 input rows).
There is no imposed secular trend in the block totals. The case is periodic,
not realistic stochastic weather or a statistically stationary process.

Model years contain 365 days. Gregorian timestamps advance continuously;
seasonality follows the forcing index, not Gregorian day-of-year. The warm
case avoids snow accumulation and isolates soil/canopy capacity. Each seed
generates different rain, reproducibly, without a committed forcing dataset.

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

The discriminating criterion is `state_bounds`: soil stays in [0, 320] mm,
canopy in [0, 2] mm, and the remaining reported stores stay nonnegative.
The existing numerical slack is 1e-6 mm. Extra groundwater is not folded into
the soil capacity. FLEX uses case capacities (area-weighted for FLEX-Topo);
SAC-SMA scales its upper/tension-zone capacities and reports lower free water
as groundwater. No universal finite groundwater capacity is imposed here.

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
only 43.8 mm. Tests require the short prefix to stay within bounds and the
long run to fail specifically on `state_bounds`, while both close their
budgets and the negative control's runoff remains nonnegative in this case.

The gate requires reference_bucket, flex_lumped, flex_topo, and sacsma_snow17
to pass. Merely extending the precipitation-normalized closure check does
not amplify a constant fractional leak; eventual storage-bound violations
are the added scientific target.

The final forcing strengthens the seasonal rain and PET cycles relative to
the issue's preliminary experiment, so correct models satisfy the existing
anti-degeneracy checks without lowering their thresholds. Five- and ten-year
warmups are both tested. FLEX-Topo still adjusts slightly after five years;
the probe permits this bounded adjustment, and the tests compare late-cycle
means rather than requiring the first scored cycle to be fully equilibrated.

## Limits

This probe cannot detect every slow leak or a trend that remains within valid
bounds. A model can fabricate bounded states and compensating fluxes and still
pass; passing is not proof of physical understanding. Legitimate equilibration,
seasonal storage, and declared groundwater exchange are not themselves failures.
The negative reporter can produce negative runoff in drier, unrelated cases;
the tests explicitly exclude that alternative failure on this forcing.

## Reproduce

The local gate passes all five fixed seeds for all four physical baselines,
and the slow reporter fails only `state_bounds`. The focused tests also cover
three independent seeds, reproducibility, prefix consistency and longer spinup.
The archive records the three evaluated physical models as PASS and Google as
INCOMPLETE (required budget variables absent). The local dhbv2 attempt is
recorded as ERROR because Docker was unavailable, not as a physical violation;
its container evaluation remains to be rerun before final acceptance.

```sh
ht validate
ht gate --probe mass/multi-decadal-drift
pytest -q tests/test_multi_decadal_drift.py
```

Author: Bing Li (@hiter-joe). Proposal: https://github.com/Flood-Lab/HydroTuring/issues/13.

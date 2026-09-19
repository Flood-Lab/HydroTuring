# mass/routing-network-closure
Reach-by-reach mass closure in a branching network.

## What it asserts

This is a component-level probe for physical, hybrid, and learned hydrological
models that expose an independently driven, reach-resolved routing state. It
does not test rainfall-runoff generation. Outlet-only rainfall-runoff models
and models without an explicit routing store are N/A, not failures.

Two prescribed headwater inflows enter reaches A and B. Their outflows meet at
zero-storage junction J and become the inflow to reach C.

```text
q_in,A -> reach A -> q_out,A --\
                                  J -> q_in,C -> reach C -> q_out,C
q_in,B -> reach B -> q_out,B --/
```

For each reach `i` and daily interval `t`, and for the junction,

```text
r_i,t = (q_in,i,t - q_out,i,t) * dt - (V_i,t - V_i,t-1)      [m3]
r_J,t = (q_out,A,t + q_out,B,t - q_in,C,t) * dt              [m3]
```

`q_in` and `q_out` are interval-mean discharge in `m3/s`; `V` is
interval-end river-channel storage in `m3`; and `dt = 86,400 s`. A row's
`time` is the start of its interval. For the first scored interval,
`V_i,t-1` is read from the final spin-up row, so no scored interval is lost.

Accept when the sum of absolute interval residuals is under **5 percent of
incoming volume**, separately for every reach and for the junction:

```text
sum_t(abs(r_i,t)) / sum_t(q_in,i,t * dt)                   <= 0.05
sum_t(abs(r_J,t)) / sum_t((q_out,A,t + q_out,B,t) * dt)    <= 0.05
```

Absolute interval residuals prevent loss in one interval from cancelling
creation in another. Scoring reaches separately prevents equal-and-opposite
errors in A and B from cancelling basin-wide.

The denominator can approach zero if a model routes essentially nothing. A
case therefore fails as degenerate when incoming volume is below 1 mm over
the relevant contributing area: 60,000 `m3` at A, 90,000 `m3` at B, and
150,000 `m3` at C and J. The floor is an eligibility test; it does not replace
the actual denominator.

### Proposed routing contract extension

The current scalar `/io` contract has no reach dimension. This proposal adds
two prescribed forcing columns and one long-format output while leaving
`result.csv` unchanged:

```text
/io/input/forcing.csv
time,pr,tas,pet,q_in_A,q_in_B

/io/output/routing.csv
time,reach_id,q_in,q_out,channel_storage
```

There is exactly one routing row per forcing time and declared reach.
`reach_id` must match `static.json.routing_network.reaches`; `(time, reach_id)`
must be unique; and all values must be finite. `q_in` and `q_out` are in
`m3/s`; `channel_storage` is the reach's river storage alone in `m3`, excluding
land, soil, groundwater, floodplain, and overland stores.

For A and B, `q_in` is the external inflow actually accepted by the routing
component and must reproduce the prescribed `q_in_A` and `q_in_B` within

```text
abs(q_reported - q_forcing)
<= 1e-9 m3/s + 1e-6 * abs(q_forcing).
```

For C, `q_in` is the routing component's independently reported upstream
inflow. It is not constructed by the harness or adapter as
`q_out,A + q_out,B`: doing that would make the junction criterion true by
construction and prevent it from distinguishing a reach-internal budget error
from a junction-transfer error.

The proposed declarations are:

```yaml
# probe.yaml
requires:
  routing: [q_in, q_out, channel_storage]
  forcing: [q_in_A, q_in_B]
  static: [routing_network]

# model.yaml
emits:
  fluxes: [pr, evspsbl, mrro, dis, gwex]
  states: [mrso, snw, canopy, channel]
  routing: [q_in, q_out, channel_storage]
needs_forcing: [pr, tas, pet, q_in_A, q_in_B]
needs_static: [routing_network]
```

Both probe and model schemas, request/output validation, compatibility checks,
and the variable table must be extended for these declarations. A model without
the table is `N/A (INCOMPLETE)`; a model that does not declare consumption of
both headwater series and the topology is `N/A (INCOMPATIBLE)`.

## The case

`generate.py` produces, from a seed:

- a fixed Y network with contributing areas of 60, 90, and 150 `km2`;
- 365 spin-up days with intermittent prescribed inflows to both headwaters;
- a 90-day scored window with one independent A event and one independent B
  event in the first 14 days, followed by at least 70 inflow-free days;
- `q_in_A` and `q_in_B` in `m3/s`, converted exactly from generated depth and
  the two headwater areas;
- standard `pr/tas/pet` columns retained for contract compatibility, with both
  precipitation and PET set to zero; and
- no channel evaporation, lateral channel inflow, groundwater-channel
  exchange, withdrawal, or inter-basin transfer.

Across seeds 0--299, each scored headwater event supplies 30.0--65.0 mm of
inflow-equivalent depth and daily headwater peaks span 5.23--33.45 `m3/s`.
Zero precipitation isolates routing from runoff generation. Zero PET is equally
important: Wflow otherwise removes open-water evaporation from the river, which
would be a legitimate sink missing from the stated three-term reach budget.

The case is intentionally recognisable as a routing-component test. The
contract and topology already declare its applicability; seeded event timing
and volume prevent memorisation of the numerical answer. Storage is not reset
at the scoring boundary, and end-of-window storage remains in every budget.

Both `period_days` and `min_window_days` are supported by the current schema.
Full `ht validate` will still fail until `requires.routing`, `emits.routing`,
the routing output path, and the `routing_*` criteria are implemented.

## Why these criteria

Closure alone is trivially satisfiable. Each additional criterion closes a
specific escape.

| Criterion | The cheat it closes off |
| --- | --- |
| `routing_input_fidelity` | rescaling or replacing the prescribed headwater inflows |
| `routing_reach_closure` | a local sink, reach-to-reach compensation, or temporal cancellation |
| `routing_junction_closure` | water lost or created during transfer through J |
| `routing_recession_bound` | retaining routed water indefinitely after its hydrograph has receded |
| `routing_state_bounds` | negative discharge/storage or storage above declared capacity |

An instantaneous zero-storage router is allowed. It may have no resolvable lag
at PT1D, but it does not violate mass conservation. Timing and physically
realistic lag belong to routing-timing probes rather than this mass probe.

The recession criterion is the per-reach analogue of the merged
`momentum/routing-conservation` criterion. For each reach,

```text
peak_i,t = max q_in,i over the preceding lookback_days

allowed_i,t = max_lag_days * peak_i,t * 86,400
              * (1 + tolerance)
              + min_allowance_mm * area_i * 1,000

V_i,t <= allowed_i,t
```

The 70-day dry tail is part of the test mechanism: after 30 days the event peak
leaves the lookback window and the proportional allowance falls to zero. The
remaining absolute allowance prevents harmless dead storage or numerical traces
from being labelled a leak. `0.05 mm` is a provisional starting value borrowed
from the merged scalar-routing probe; it must be recalibrated against native
per-reach Wflow storage before the acceptance gate is reported.

`reference_network_stuck_router` retains 10 percent of routed water, proving
the criterion catches an obvious fault. `reference_network_leaky_router`
retains 1 percent, constraining how far the absolute allowance may be relaxed.
Both report the retained water as storage, so reach closure still passes and
the recession criterion must do the catching.

The exact positive is `reference_network_exact`. The independent physical
positive is the existing `wflow_sbm`, version-bumped and extended from its
current one-cell schematisation to three river cells. Its adapter must:

1. map A and B to Wflow's external river-inflow forcing;
2. build A and B draining to C in the generated static maps;
3. report native interval-mean river inflow and discharge and interval-end
   river storage for all three cells; and
4. write `routing.csv` without reconstructing C inflow from the other reported
   columns.

The four standard lumped baselines lack reach-indexed routing output and are
N/A; they neither pass nor fail this probe. The Wflow gate result and calibrated
absolute allowance must be archived before merge.

| Model | Required failure |
| --- | --- |
| `reference_network_input_rescale` | `routing_input_fidelity` |
| `reference_network_compensation` | `routing_reach_closure` |
| `reference_network_junction_loss` | `routing_junction_closure` |
| `reference_network_temporal_shift` | `routing_reach_closure` |
| `reference_network_stuck_router` | `routing_recession_bound` |
| `reference_network_leaky_router` | `routing_recession_bound` |
| `reference_network_negative` | `routing_state_bounds` |

## Reproducing a failure

Every run records its seed. To re-run one:

```bash
ht run --model <name> --probe mass/routing-network-closure --seed <n>
```

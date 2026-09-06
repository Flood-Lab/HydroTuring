# momentum/routing-conservation

`channel` was added to the contract so that a model with a unit hydrograph
or a channel store could report the water it holds in transit and close
its budget. This probe asks whether that store behaves like one. Two
assertions: it is never negative, and at any step it holds at most
`max_lag_days` (15) times the largest runoff rate of the preceding month,
which is the most a hydrograph of that length can have accumulated from
the flow feeding it. A kernel that does not sum to one, or a store that
accumulates, breaks the second within a few storms.

It sits under the momentum law because the channel is where momentum is
conserved, and conservation of the water in it is the shadow that every
routing scheme casts on the mass budget. A model that does not route reports
its channel store as identically zero, as the exact bucket does, and passes
trivially; the two FLEX models, which route through a triangular lag, are
the real test.

| Criterion | Asserts |
| --- | --- |
| `routing_conservation` | the channel store is non-negative and never exceeds 15 days of the largest recent runoff |
| `non_degenerate` | runoff varies with the weather |

Must-fail: `reference_stuck_router`, the bucket with a three-day kernel that
sums to 0.9.

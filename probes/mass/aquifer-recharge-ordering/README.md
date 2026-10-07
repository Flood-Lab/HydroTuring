# mass/aquifer-recharge-ordering

Two runs of the same passive aquifer under the same seeded river stage differ
only by a direct `gw_recharge` pulse of 5 mm/day on the first two scored days.
The question is the groundwater comparison principle (#135): water added to an
aquifer whose outflow does not fall as its storage rises can never leave it
with less water, can never make the exchange less favourable to the river, and
has to be accounted for.

## The case

One aquifer under 1 km² with specific yield 0.2, its head at the river's
stage, draining to the river through a conductance of 1000 m²/day, the value
`mass/groundwater-datum-invariance` uses. The record is 18 scored days after
2 days of spinup. At this conductance the pulse raises the exchange toward the
river by about 0.8 mm over the record, eight times the tolerance, so a reversal
of it is visible; at 100 m²/day it was 0.08 mm and no reversal could be seen.

## `aquifer_recharge_ordering`

Let `I_n` be the cumulative added recharge, `D_n` the cumulative incremental
exchange toward the river (control minus perturbed `gw_sw_exchange`, which is
positive into the aquifer), `B_n` the cumulative change in a declared
`gw_boundary` (positive into the aquifer) and `ΔS_n` the paired storage
difference. Three one-sided errors are scored:

* `E_S = max_n [-ΔS_n]_+`, the largest fall of storage below the control;
* `E_Q = Σ_i [-d_i]_+ Δt`, the exchange increment sent the wrong way; and
* `E_B = max_n |ΔS_n + D_n - B_n - I_n|`, the worst running partition residual.

`E_B` largely restates `groundwater_balance`, which this probe also checks on
both runs: if each run closes its own budget, the paired one closes to within
the two runs' tolerances, which can add to slightly more than `ε`. The probe's
independent content is `E_S` and `E_Q`.

All three share one tolerance, `ε = max(0.01 · I_*, 1e-3 mm, min(1e-5 · S_max,
0.05 mm))`, where `I_*` is the total added recharge (10 mm, so ε = 0.1 mm) and
`S_max` the largest absolute storage either run reports. The scale comes from
the intervention rather than the background storage; the capped third term
absorbs rounding in a large absolute `gw` without letting a datum buy unlimited
slack, as in `groundwater_balance`.

## A declared boundary

A model with a boundary on the aquifer, a GHB or a well, reports it as
`gw_boundary`. Both `groundwater_balance` and the paired budget credit it into
the aquifer. A steady boundary cancels between the runs; a head-dependent one
drains part of the pulse, and the order still holds because that flux too does
not fall as the head rises.

## What this probe does not catch

The order is non-strict. A zero exchange increment passes, as it must for an
exchange branch that is genuinely insensitive, so a model whose exchange
ignores its own aquifer head stores the whole pulse, never releases it, and
passes. Asking that the water drain back out needs the aquifer's response
time. For a single store that is `S_y · A / C`, from the case's own
attributes, but for a distributed aquifer it depends on a transmissivity and a
geometry the case does not supply: MODFLOW 6's 10 × 10 grid at K = 1 m/day
drains on a time scale of thousands of days, and a floor built on
`S_y · A / C` failed it. `tests/test_recharge_ordering.py` pins the gap, so a
recovery check added later shows up there as a change of scope.

## Baselines

The four standard catchment references have no direct aquifer recharge, so
`must_pass` is the exact reference alone. MODFLOW 6, a physical model with this
process, is archived as passing in `models/result.csv`.

| Reference | What it does | What catches it |
| --- | --- | --- |
| `reference_recharge_order_exact` | the exact daily solution of one aquifer cell draining to the river | nothing, it must pass |
| `reference_recharge_blind` | ignores `gw_recharge`, so its own budget no longer closes | `aquifer_recharge_ordering` (`E_B`), and `groundwater_balance` |
| `reference_recharge_overshoot` | sends 120% of a recharge day to the river | `aquifer_recharge_ordering`: `E_S` at about 20 times the tolerance, and `E_Q` at about 1.5 as the river refills the deficit |

A reversed exchange increment that keeps each run's budget closed trips `E_Q`
alone. No reference fakes it; `tests/test_recharge_ordering.py` builds it from
the exact reference and checks that `E_Q` is the arm that fails.

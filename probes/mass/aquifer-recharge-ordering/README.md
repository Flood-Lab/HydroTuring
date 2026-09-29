# mass/aquifer-recharge-ordering

The paired case adds only direct `gw_recharge` to the same passive aquifer. The
`aquifer_recharge_ordering` criterion evaluates storage ordering, the sign of
the incremental exchange, and the cumulative partition identity using the
added recharge volume as its scale. A plateau with zero incremental exchange
is valid; a recharge-blind or nonmonotone update is not.

The exact reference supplies the deterministic gate. The independent MODFLOW 6 baseline is enabled for paired perturbations and is archived separately; the four standard catchment references do not expose a direct aquifer recharge interface.

## Criterion and tolerance

For each scored step, let `I_n` be cumulative added recharge and `D_n` be
cumulative incremental exchange (positive when the added recharge sends more
water to the river). The criterion records three one-sided errors:

* `E_S = max([-ΔS_n]_+)`, the largest storage decrease;
* `E_Q = sum([-d_i]_+ · Δt_i)`, the total exchange increment in the wrong
  direction; and
* `E_B = max(|ΔS_n + D_n + B_n - I_n|)`, the worst cumulative partition
  residual, where `B_n` is the cumulative declared boundary-source change.

The common tolerance is `ε = max(0.01 · I_* , 1e-3 mm, f · S_max)` where
`I_*` is total added recharge, `S_max` is the largest absolute storage in
either run, and `f · S_max` is a capped storage precision floor (`f = 1e-5`,
cap `0.05 mm`). Every error must be no larger than `ε`. This makes the scale
come from the perturbation itself while accounting for storage rounding without
allowing a large absolute datum to buy unlimited tolerance.

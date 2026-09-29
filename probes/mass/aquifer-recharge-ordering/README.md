# mass/aquifer-recharge-ordering

The paired case adds only direct `gw_recharge` to the same passive aquifer. The
`aquifer_recharge_ordering` criterion evaluates storage ordering, the sign of
the incremental exchange, and the cumulative partition identity using the
added recharge volume as its scale. A plateau with zero incremental exchange
is valid; a recharge-blind or nonmonotone update is not.

The exact reference supplies the deterministic gate. The independent MODFLOW 6 baseline is enabled for paired perturbations and is archived separately; the four standard catchment references do not expose a direct aquifer recharge interface.

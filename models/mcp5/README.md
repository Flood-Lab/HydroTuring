# mcp5: MCPModel5, the three-node distributed-input MCP network

MN_SALO^DI(3) of Wang & Gupta (2025, *WRR* 61(12), e2025WR040178), proposed in issue #93. Three
stores that each take a learned share of the rain, with gates that read all three stores through
learned matrices, and a discharge that is the sum of what they release. Trained on the Leaf River.

Source: [YuanHWang/5MCPmodelnet_HydroTuring](https://github.com/YuanHWang/5MCPmodelnet_HydroTuring)
at `5944b5f23adb7265a84fd10d62a8bd5501e25208`, Apache-2.0 for code and weights. The checkpoint ships
in the repository; nothing was trained for this evaluation. The repository's directory `MCPmodel5`
holds a script named `evaluate_model4_mn3_distributed_input.py` (and `MCPmodel4` one named
`evaluate_model5_...`); its README's table agrees with the directories, and so does this package.

```
Q_i,t = oo_i,t * C_i,t        E_i,t = olc_i,t * C_i,t,  olc_i,t = min(ol_i,t, PET_t / C_i,t)
C_i,t+1 = C_i,t - Q_i,t - E_i,t + a_i P_t
oo_i,t = s_o,i * sigmoid(b_o,i + sum_j n_j,t W_o[j,i])         n_j,t = (C_j,t - c_mean) / c_std
ol_i,t = s_l,i * sigmoid(sum_j n_j,t W_l[j,i] + b_l,i + v_i (PET_t - 2.9086) / 1.8980)
Q_t = sum_i Q_i,t
```

## Every parameter and where it came from

From `MCPmodel5/best_model_epoch2999.pt` unless stated, nodes 0-2.

| Parameter | Value |
| --- | --- |
| `weight_r_yom` | 0.7315, 4.0560, 0.2645 |
| `weight_r_ylm` | 1.4161, 1.2485, 1.6540 |
| `weight_r_yfm` | 4.0528, 3.6728, 4.1752 |
| derived s_o | 0.0326, 0.5740, 0.0182 |
| derived s_l | 0.0646, 0.0346, 0.0730 |
| `bias_b0_yom` | -9.4491, -7.8631, -7.2375 |
| `bias_b0_ylm` | -3.8526, -4.1780, -3.8965 |
| `weight_b2_ylm` (v) | 0.6256, 2.0152, -0.2742 |
| `wi`, softmax to the input shares a | 0.7232, 1.1101, 1.2942; a = 0.2357, 0.3471, 0.4172 |
| `c_mean`, `c_std` | 412.9139085, 77.49943867: defaults of the evaluation script |
| PET normalisation | 2.9086, 1.8980: hard-coded in the upstream forward |
| initial stores | 0: upstream forward |
| `spinLen`, `traintimeLen` | 0, row count (upstream: 1095, 8400): only slice the observations for a diagnostic the adapter does not read; `y_obs` is zeros |

`weight_b1_yom` (W_o), rows j, columns i:

```
  1.9651  -7.7861   1.2787
  3.4783   2.3773   3.4863
  5.0414   5.6674   8.1685
```

`weight_b1_ylm` (W_l), rows j, columns i:

```
  0.5234   0.1750  -0.0392
  0.2066   1.9736  -1.9869
  0.0100   0.0020   2.1946
```

The checkpoint loads with `strict=True`. Run as the evaluation script runs it on the Leaf River, the
adapter's model gives NSE 0.863 / 0.841 / 0.828 on the train / selection / test flags. The repository
publishes no number to compare that with.

## What is reported, and the one choice in it

`pr` (echoed), `mrro` (the upstream output Q), `evspsbl`, in mm/day; `mrso`, in mm, at the end of the row.

**Upstream sums only the discharge.** It returns evaporation and storage per node and never combines
them. The adapter sums both. The input shares sum to one, so the nodes' own updates add up to one
catchment budget under a plain sum and no other. Note that the shares divide the rain but not the
demand: each node's cap is the whole PET, as upstream writes it, so the sum of the nodes' evaporation is
not capped at PET. That is the model's and is reported as it comes; it is the `pet-consistency` failure
below.

Upstream's `c_n[t]` is each store *before* row t. The harness differences a state against the row's
fluxes, so row t reports the stores after it, read from `c_n[t+1]`: the model is given one extra row of
no weather at the end of the record and that row's outputs are thrown away.

## Not reported, and why

Named under `not_reported` in `run.json`. The model has none of these, and none is filled in.

| Output | Why not |
| --- | --- |
| `snw` | no snow store; temperature is never read |
| `canopy` | no interception store |
| `gw`, `channel` | three undifferentiated stores, reported together as `mrso`; no routing |
| `dis` | the model knows no catchment area |
| `gwex` | each node conserves mass; nothing crosses the boundary but Q and E |
| energy fluxes | no energy balance |

**What that does to the standing.** `mass/catchment-closure` and seven other budget probes lack only `snw`
and `canopy` from this model, so this model is N/A (INCOMPLETE) on the probe that asks. Off the harness, the reported
columns close: at most 9.0e-5 mm in any row and 1.2e-4 mm cumulative over 5805 mm of rain (2000 rows,
PT1D). That number is not an archived verdict.

## Choices, and how much each moves the results

**Spin-up.** Upstream repeats the first water year three times in front of the record, from empty stores.
The adapter does the same with the first 365 days of the rows it is given, or all of them when the record
is shorter. Those rows run before row 0 and are never reported. On daily probes that is three times the
case's own 365-day spin-up, nothing later than it, so it is causal. The table is the prototype under the
same harness and gate seeds. This model moves more than the others.

| Repeats | `dry-down` | `resolution-invariance` | `response-nonnegativity` |
| --- | --- | --- | --- |
| 0 | PASS | PASS (`evspsbl` 0.9% of rain) | FAIL (0.331 mm/day dip) |
| 1 | FAIL (11.8%) | PASS (`mrro` 9.4%) | PASS (dip 7e-4 of the storm) |
| 3 (used) | FAIL (4.8%) | **FAIL** (`mrso_end` 10.7%) | FAIL (0.143 mm/day dip) |

Those three rows depend on the spin-up. Their mechanisms below are the model's, but whether they cross a
limit on these seeds is sensitive, so treat the verdicts as not robust. A first version cycled a short
record until it filled three years, inventing about 2.6 m/yr of rain; it was dropped.

**Step scaling.** The upstream class knows one step. The resolution probe runs it at PT1H regardless, so
below a day the adapter runs **its own copy of the upstream loop** (`forward_step_scaled`), not the
repository's code. It differs in three ways: rain is a depth per row, each per-day gate becomes
1 - (1 - k)^dt, and the PET cap is a depth per row. The gates still read PET in mm/day. At dt = 1 the copy
matches the upstream forward to single-precision round-off (at most 1.2e-4 mm on a store over 3000 rows).
Every daily row goes through the upstream class. The manifest declares PT1D only, the step in the proposal.
Scaling the total outflow fraction and splitting it, instead of each gate, leaves the resolution numbers
unchanged to 0.01% of rain.

**Precision.** The model runs in float32, as upstream does; see the closure figure above.

**Window.** `window_days: full`, since a record takes seconds. On the one probe where this differs,
`warming-response` passes either way.

**Seed.** Recorded and unused; the model is deterministic.

## Verdict

Gate seeds, Apptainer, 2026-10-30: **FAIL (VIOLATION), 7 of 12 scored**, 25 N/A (24 INCOMPLETE,
1 INCOMPATIBLE: `routing-lag-consistency` needs catchment attributes the model does not read).

| Probe | Result | Detail |
| --- | --- | --- |
| energy/pet-consistency | **FAIL** | 1.14 of demand when wet |
| mass/antecedent-monotonicity | PASS | +9.8 mm on a 60 mm storm (0.16 of it) |
| mass/area-invariance | PASS | trivially: the model reads no area |
| mass/causality | PASS | nothing moves before the storm |
| mass/dry-down | **FAIL** | `mrro` rises 4.8% in a week with no rain; depends on the spin-up |
| mass/extreme-rain | **FAIL** | 1.11 of the added rain comes back |
| mass/phase-counterfactual | PASS | trivially: the model reads no temperature |
| mass/resolution-invariance | **FAIL** | `mrso_end` 10.7% of rain; depends on the spin-up |
| mass/response-nonnegativity | **FAIL** | `mrro` 0.143 mm/day below control; depends on the spin-up |
| mass/runoff-bounds | PASS | runoff 0.41 of rain |
| mass/steady-state | PASS | settles at 1.238 mm/day |
| mass/warming-response | PASS | -0.24 and -0.29 per unit of demand |

## The failures

All verified with targeted runs. Three come from the negative entry of W_o: node 1 carries the flow in dry
periods, and column 1 of W_o reads node 0's store at -7.79.

**pet-consistency: three caps, one atmosphere.** Each node is capped at the whole PET, and the catchment's
evaporation is their sum. Each node stays at or under PET (largest ratio 1.00002), while the sum exceeds it on
35% of rows and reaches 3.0 times PET at most, and 1.14-1.17 of demand over the wettest fifth. The shares that
divide the rain were never applied to the demand, so this is structural.

**dry-down: a gate opens as another store drains.** At the worst block the stores fall by -3.5, -4.4 and -1.6
mm. Through column 1 of W_o that is +0.108 on the logit, a gate 11% wider, against a store 1.2% lower: predicted
+10%, measured +10% on node 1's flow.

**response-nonnegativity: the same term.** At the dip node 0 is 3.9 mm higher, node 1 14 mm higher and node 2 6
mm lower than in the control. Through column 1 that is -0.39 + 0.43 - 0.44 on the logit, and node 1's `oo` goes
from 0.0013 to 0.0009. Node 1's flow drops by 0.173 mm/day with the other two nearly unchanged.

**extreme-rain: stored water is released with the storm.** On the failing rung `mrro` rises by 1.11 of the 342 mm
added and evaporation by 0.30, and the store at the end of the record is 0.41 of it lower. The budget holds, so the
excess was already in the catchment; that much is measured. That it leaves through gates the storm opened is inferred
from W_o rather than isolated.

**resolution-invariance: the two steps enter the scored month in different states.** Volumes over the month agree
(`mrro` within 2.2% of rain, `evspsbl` within 1.1%). What exceeds the limit, by 0.7 of a percent, is `mrso_end`, an
absolute state: 737 against 720 mm entering the month and 879 against 856 mm leaving it, on the first seed. Most of
that gap was there before the month began.

## Packaging

The image is built like `models/mcp1`: `python:3.10-slim` at a fixed digest, the repository's torch 2.9.0 and numpy
2.2.2 with torch's dependencies pinned, and only `MCPBRNN_lib_tools`, the checkpoint and `LICENSE` from the pinned
commit. No network at run time. A case takes about 2 s plus container start. `MCP_SPINUP_REPEATS` exists only for
the table above; nothing sets it in the container.

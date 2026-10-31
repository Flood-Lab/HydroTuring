# mcp4: MCPModel4, the five-node distributed-state MCP network

MN_SALO^DS(5) of Wang & Gupta (2025, *WRR* 61(12), e2025WR040178), proposed in issue #92. Five
stores that each take the whole rain, with gates that read all five stores through learned
matrices, and a discharge that is their outflows averaged with learned weights. Trained on the
Leaf River.

Source: [YuanHWang/5MCPmodelnet_HydroTuring](https://github.com/YuanHWang/5MCPmodelnet_HydroTuring)
at `5944b5f23adb7265a84fd10d62a8bd5501e25208`, Apache-2.0 for code and weights. The checkpoint
ships in the repository; nothing was trained for this evaluation. The repository's directory
`MCPmodel4` holds a script named `evaluate_model5_mn5_distributed_state.py` (and `MCPmodel5` one
named `evaluate_model4_...`); its README's table agrees with the directories, and so does this
package.

```
Q_i,t = oo_i,t * C_i,t        E_i,t = olc_i,t * C_i,t,  olc_i,t = min(ol_i,t, PET_t / C_i,t)
C_i,t+1 = C_i,t - Q_i,t - E_i,t + P_t
oo_i,t = s_o,i * sigmoid(b_o,i + sum_j n_j,t W_o[j,i])         n_j,t = (C_j,t - c_mean) / c_std
ol_i,t = s_l,i * sigmoid(sum_j n_j,t W_l[j,i] + b_l,i + v_i (PET_t - 2.9086) / 1.8980)
Q_t = sum_i w_i Q_i,t
```

The off-diagonal entries of W are where this model differs from a set of independent buckets:
each node's gate depends on the other nodes' stores, some of them negatively.

## Every parameter and where it came from

From `MCPmodel4/best_model_epoch2707.pt` unless stated, nodes 0-4.

| Parameter | Value |
| --- | --- |
| `weight_r_yom` | 1.0148, 2.7074, 1.7768, 1.0264, 0.5027 |
| `weight_r_ylm` | 1.9386, 1.1083, 1.4025, 2.8937, 2.0237 |
| `weight_r_yfm` | 3.6950, 2.0844, 2.9386, 3.3782, 4.2822 |
| derived s_o | 0.0552, 0.5752, 0.2048, 0.0556, 0.0203 |
| derived s_l | 0.1391, 0.1162, 0.1408, 0.3600, 0.0927 |
| `bias_b0_yom` | -6.7707, -8.4439, -10.4903, -7.2332, -8.9809 |
| `bias_b0_ylm` | -4.1241, -3.8294, -3.7684, -2.5371, -3.9311 |
| `weight_b2_ylm` (v) | 0.2640, -1.5267, 1.5553, 0.9854, 1.0158 |
| `wo`, softmax to the output weights w | -1.1211, 1.4589, 3.1454, -0.4491, 4.4026; w = 0.0030, 0.0390, 0.2109, 0.0058, 0.7413 |
| `c_mean`, `c_std` | 412.9139085, 77.49943867: defaults of the evaluation script |
| PET normalisation | 2.9086, 1.8980: hard-coded in the upstream forward |
| initial stores | 0: upstream forward |
| `spinLen`, `traintimeLen` | 0, row count (upstream: 1095, 8400): only slice the observations for a diagnostic the adapter does not read; `y_obs` is zeros |

`weight_b1_yom` (W_o), rows j, columns i:

```
  6.2119   0.5939   2.5891   3.2116   4.0589
 -0.0453  -0.4540  -1.1554   1.6635  -0.4156
  1.2079   1.3747   2.0074   0.9486   1.5755
 -4.8944   1.2210  -3.7073   6.9855  -2.8261
  0.7541  -0.0086   0.6903  -0.4124   0.8666
```

`weight_b1_ylm` (W_l), rows j, columns i:

```
  2.2411   0.3170   2.0809   0.0535   1.0496
 -0.6950   0.4236  -0.4171   0.2381   0.1260
 -0.2252   0.5481   1.2105  -0.5419   0.0475
  0.0933   1.7003   0.7142   0.8956  -0.0457
 -0.7078   0.1734  -0.2774  -0.3251   1.1879
```

The checkpoint loads with `strict=True`. Run as the evaluation script runs it on the Leaf River,
the adapter's model gives NSE 0.927 / 0.891 / 0.906 on the train / selection / test flags, the
best of the five. The repository publishes no number to compare that with.

## What is reported, and the one choice in it

`pr` (echoed), `mrro` (the upstream output Q), `evspsbl`, in mm/day; `mrso`, in mm, at the end of
the row.

**Upstream combines only the discharge.** It returns evaporation and storage per node and never
combines them. The adapter uses the output weights for both. Every node takes the whole rain and
the weights sum to one, so P = sum_i w_i P, and this is the only combination under which the
nodes' own updates add up to one catchment budget. The weights then act as area fractions. Equal
weights, the obvious alternative, would not close, and would move evaporation on
`pet-consistency` from 0.97/0.33 (wet/dry, first seed, 0.74 overall) to 0.63/0.22 (0.47 overall).

Upstream's `c_n[t]` is each store *before* row t. The harness differences a state against the
row's fluxes, so row t reports the store after it, read from `c_n[t+1]`: the model is given one
extra row of no weather at the end of the record and that row's outputs are thrown away.

## Zero stores, and what is not reported

`zero_stores` in `run.json`: `snw` (no snow store; temperature is never read) and `canopy` (no
interception store). Both are reported as identically zero, so the budget probes score the model.

Named under `not_reported`. The model has none of these, and none is filled in.

| Output | Why not |
| --- | --- |
| `gw`, `channel` | five undifferentiated stores, reported together as `mrso`; no routing |
| `dis` | the model knows no catchment area |
| `gwex` | each node conserves mass; nothing crosses the boundary but Q and E |
| energy fluxes | no energy balance |

## Choices, and how much each moves the results

**Spin-up.** Upstream repeats the first water year three times in front of the record, from empty
stores. The adapter does the same with the first 365 days of the rows it is given, or all of them
when the record is shorter. Those rows run before row 0 and are never reported. On daily probes
that is three times the case's own 365-day spin-up, nothing later than it, so it is causal. The table is the adapter run directly under the same harness and gate seeds, on every probe that scores this model, after the zero-store rule; at three repeats that run reproduces every archived verdict.

| Repeats | Result |
| --- | --- |
| 0 | `resolution-invariance` PASS (`mrso_end` 0.5% of rain); no other verdict changes |
| 1 | `resolution-invariance` PASS (`mrro` 2.9%); no other verdict changes |
| 3 (used) | `resolution-invariance` **FAIL** (`mrso_end` 15.5%) |

So this model's `resolution-invariance` FAIL depends on the spin-up; see the failure below.

A first version cycled a short record until it filled three years. On the resolution probe that
fed the model a 40-day stretch of 7-8 mm/day for three years, about 2.6 m/yr, and `mrso_end` then
differed by 47.6% of rain, because the stores split into two regimes: some runs reached node 1 near
2200 mm and node 0 near 120 mm, others stayed at 370-460 mm, and which one depended on the step.
The Leaf River never takes node 1 above about 1000 mm. That forcing was the adapter's invention, not
the record's, so it was dropped; the bistability it exposed is still the model's.

**Step scaling.** The upstream class knows one step. The resolution probe runs it at PT1H regardless,
so below a day the adapter runs **its own copy of the upstream loop** (`forward_step_scaled`), not the
repository's code. It differs in three ways: rain is a depth per row, each per-day gate becomes
1 - (1 - k)^dt, and the PET cap is a depth per row. The gates still read PET in mm/day. At dt = 1 the
copy matches the upstream forward to single-precision round-off (at most 5.8e-4 mm on a store over
3000 rows). Every daily row goes through the upstream class. The manifest declares PT1D only, the step
in the proposal. Scaling the total outflow fraction and splitting it, instead of each gate, moves the
resolution numbers by at most 0.03% of rain, so the form of the rule is not what is being measured.

**Precision.** The model runs in float32, as upstream does; see the closure figure under the verdict.

**Window.** `window_days: full`, since a record takes seconds. `warming-response` was also checked on the default window and passes either way.

**Seed.** Recorded and unused; the model is deterministic.

## Verdict

Gate seeds, Apptainer, 2026-10-30, version `5944b5f-epoch2707.2`: **FAIL (VIOLATION), 7 of 20 scored**, 17 N/A
(16 INCOMPLETE, 1 INCOMPATIBLE: `routing-lag-consistency` needs catchment attributes the model does not
read). The INCOMPLETE rows ask for output this model does not have: heat fluxes and radiation (six energy probes), river-aquifer exchange, `snm`, `gwex`, `channel`, `stage` and `dis`.

**`resolution-invariance` is not robust.** Its FAIL holds at three spin-up repeats and is a PASS at zero and at one
(table above). Read that row as undecided.

`closure` passes on all eight probes that score it, worst residual 4.5e-7 of the rain, which is single-precision
round-off in updates that conserve mass by construction. Trivial passes for this model: `area-invariance` (reads no
area), `phase-counterfactual` (reads no temperature) and `time-origin-invariance` (reads no calendar), as are the
bounds on the zero `snw` and `canopy`.

Five of the 13 FAILs are `state_bounds` and nothing else.

| Probe | Result | Detail |
| --- | --- | --- |
| energy/pet-consistency | PASS | 0.99 of demand when wet, 0.39 when dry |
| mass/antecedent-monotonicity | **FAIL** | -2.3 mm: the wetter catchment runs off less |
| mass/area-invariance | PASS | trivially: the model reads no area |
| mass/catchment-closure | **FAIL** | `state_bounds`: `mrso` 436-840 mm against 320 |
| mass/causality | PASS | nothing moves before the storm |
| mass/dry-down | **FAIL** | `mrro` rises 11.4% in a week with no rain |
| mass/extreme-event-closure | **FAIL** | `state_bounds`: `mrso` 413-936 mm against 320; event closure passes |
| mass/extreme-rain | **FAIL** | 1.17 of the added rain comes back |
| mass/human-abstraction | **FAIL** | 380 mm of 380 unaccounted for; `state_bounds`: `mrso` 424-833 mm against 320 |
| mass/multi-decadal-drift | **FAIL** | `state_bounds`: `mrso` 645-695 mm against 320; `et_plausible`: 1 + 4.2e-8 of demand, round-off; no drift |
| mass/phase-counterfactual | PASS | trivially: the model reads no temperature |
| mass/precipitation-counterfactual | **FAIL** | `state_bounds`: `mrso` 433-821 mm against 320; the partition passes |
| mass/resolution-invariance | **FAIL** | `mrso_end` 15.5% of rain; **depends on the spin-up** |
| mass/response-nonnegativity | **FAIL** | `mrro` 0.614 mm/day below control after added rain |
| mass/runoff-bounds | PASS | runoff 0.37 of rain |
| mass/spinup-cycle-invariance | **FAIL** | `state_bounds`: `mrso` 429-743 mm against 320; the invariance passes |
| mass/steady-state | **FAIL** | `mrro` still moving 11.7% in year 3 of constant weather |
| mass/time-origin-invariance | PASS | trivially: the model reads no calendar |
| mass/ungauged-basin-closure | **FAIL** | `state_bounds` on 12 of 12 seeds: `mrso` 456-762 mm against 124 |
| mass/warming-response | PASS | -0.17 and -0.21 per unit of demand |

## The failures

Verified with targeted runs unless a sentence says inferred. All are the model; the adapter's arithmetic contributes
only to `et_plausible`, below, which changes no verdict.

**state_bounds, on seven probes: the stores have no capacity, and their gates open far above the bound.** The probes take
the bound from `soil_capacity_mm` in `static.json`: 320 mm, or a draw of 124-682 mm on `ungauged-basin-closure`. Nothing
in the model can use it. The checkpoint's tensors are gate scales, biases, weights and output weights, no store has an
upper limit, and the adapter never opens `static.json`. The stores run wherever the learned gates let out the runoff the
climate requires. Scaling the five stores together from their mean proportions until the reported `mrso` is 320 mm, a
state the model does not visit, gives a release of 0.002 mm/day, and 4.5 mm/day at 700 mm, while the probe's records need
0.68-0.88 mm/day. So on `catchment-closure` the reported store spans 426-840 mm across the seeds, with node means of
362-607 mm on the first, and on the Leaf River, where it was trained, 438-954 mm (mean 648), above 320 mm on every day. On
`ungauged-basin-closure` its peak is 639-782 mm whatever the draw, and no seed passes, not even the draws of 679 and 682
mm. The weighted sum is the adapter's aggregate, the only one that closes (above).

The adapter could not honestly do otherwise. There is no capacity to hand `soil_capacity_mm` to, and writing it into
the normalisation or the gates would make a different model from the one trained. Subtracting a datum would not fit
either: the reported store's range on `catchment-closure` is 378-404 mm wide on each seed, wider
than the bound. Moving the water above the bound into `gw` would invent a partition the model does not have. So this is
a finding about the model, as with `parFC` in `models/dhbv2`, except that here not even a learned capacity exists. It
does not depend on the spin-up: the same probe fails the same way at zero, one and three repeats.

**human-abstraction: the model has no withdrawal term.** Structural, and verified from the row and the
code: the adapter reads `pr` and `pet` and nothing else, the paired runs are identical (`evspsbl`, `mrro` and
storage all +0.0 mm), and the whole 380 mm goes unaccounted for. AGENTS.md scores this as a failure, not as N/A.
Subtracting the column in the adapter would add a process the published model does not have, so nothing was added.

**multi-decadal-drift, `et_plausible`: round-off at the cap, not excess evaporation.** In this probe's mild climate (PET
about 1 mm/day, rain 1.8) all five nodes sit at their PET cap on every scored row, so evaporation equals demand, and the
cumulative ratio is between 1 - 9.8e-9 and 1 + 4.2e-8 on the five seeds. The criterion allows 1 + 1e-9, so four fail. Two
single-precision sources add up to that. One is the upstream `C * (ol - relu(ol - PET / C))`. The other is the adapter's,
which combines the nodes with upstream's float32 softmax weights, and those sum to 1 + 1.65e-8. Renormalising the weights
in double still leaves two seeds over, at +2.6e-8 and +3.1e-9. The same equations in float64 reach 1 + 4.4e-16 with no row
above PET. The model stays in float32, as upstream runs it, and `state_bounds` fails all five seeds of this probe
regardless. Do not read this row as evaporation exceeding demand.

Four of the other six come from the negative entries of W_o. Node 2 carries
most of the outflow in dry periods (w = 0.21), and column 2 of W_o reads node 3 at -3.71 and node 1 at
-1.16: when those stores rise, node 2's gate closes, and when they fall, it opens.

**dry-down: a gate opens as another store drains.** At the worst block the stores all fall, by -4.3, -0.1,
-9.3, -11.4 and -5.7 mm. Through column 2 of W_o that adds up to +0.115 on the logit, a gate 12% wider,
which beats the store's 1.8% fall: predicted +10.2% on node 2's flow, measured +10.2%.

**antecedent-monotonicity: more water closes the main gate.** With 120 mm of extra rain, node 3 is 70 mm
higher and node 1 110 mm higher at the storm. Node 2's logit falls by about 2.3, its gate by about a tenth,
and it releases 12.3 mm less in the 30 days, which at w = 0.21 accounts for the -2.3 mm. Extra runoff on the
three seeds is -2.29, -1.01 and +3.71 mm against +1.2 required.

**response-nonnegativity: the same term, later.** At the dip, 88 days after the storm on the first seed, the
added rain has left node 1 200 mm lower and nodes 0, 3 and 4 higher, which shuts the gates: `oo` on node 3 is
0.0070 against 0.0253 in the control. Flow falls on every node, by -0.61 mm/day in total. Added water ends
up lowering the flow.

**extreme-rain: stored water is released with the storm.** On the failing rung `mrro` rises by 1.17 of the 58
mm added, evaporation by 0.00, and the store at the end of the record is 0.17 of it lower. The budget holds,
so the excess is water that was already in the catchment; that much is measured. That it leaves through gates the
storm opened, nodes 2 and 4 carrying the rise, is inferred from W_o rather than isolated.

**steady-state: it settles, slowly and not monotonically.** Under the probe's constant weather the daily map
has a stable fixed point (stores 482, 96, 413, 513, 763 mm; Q 0.504 mm/day), with eigenvalue magnitudes 0.9465,
0.9865 and 0.9976, the first two oscillating pairs. It is not a limit cycle. In the probe's three years Q falls
to a minimum of about 0.39 mm/day late in year 2 and then climbs back toward 0.504, so the last year travels as
far as the one before (111% on the worst seed) and the criterion's slow-settling allowance, which asks for 70%,
does not apply. Continuing the run past the probe, each year's range in Q is 5.3%, 2.4%, 1.0%, 0.4%, 0.2% of its
level, and zero to four digits within about twelve more years. The verdict is correct on its terms; the model takes
about a decade to stop moving.

**resolution-invariance: the two steps enter the scored month in different states.** Volumes over the month
agree (`mrro` within 5.2% of rain, `evspsbl` within 0.13%). What exceeds the limit is `mrso_end`, an absolute
state: 738 against 760 mm entering the month and 752 against 785 mm leaving it, on the first seed. Starting
both steps from one shared state, the month itself differs by 1.1-3.3% of rain. So the difference accumulates
over the spin-up, and with no spin-up in front the verdict is PASS (table above).

## Packaging

The image is built like `models/mcp1`: `python:3.10-slim` at a fixed digest, the repository's torch 2.9.0 and
numpy 2.2.2 with torch's dependencies pinned, and only `MCPBRNN_lib_tools`, the checkpoint and `LICENSE` from
the pinned commit. No network at run time. A case takes about 2-3 s plus container start. `MCP_SPINUP_REPEATS`
exists only for the table above; nothing sets it in the container.

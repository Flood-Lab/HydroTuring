# mcp3: MCPModel3, the MA_5 three-node architecture

MA_5 of Wang & Gupta (2024, *WRR* 60(10), e2024WR037224), proposed in issue #91. A primary store
takes the rain and has three outflows; a surface node and a groundwater node route two of them;
what those two release is the discharge. Trained on the Leaf River.

Source: [YuanHWang/5MCPmodelnet_HydroTuring](https://github.com/YuanHWang/5MCPmodelnet_HydroTuring)
at `5944b5f23adb7265a84fd10d62a8bd5501e25208`, Apache-2.0 for code and weights. The
checkpoint ships in the repository; nothing was trained for this evaluation.

```
Qm_t = oo_t * C_t     R_t = oog_t * C_t      E_t = olc_t * C_t,  olc_t = min(ol_t, PET_t / C_t)
C_{t+1} = C_t - Qm_t - R_t - E_t + P_t
Qs_t = os_t * S_t     S_{t+1} = S_t - Qs_t + Qm_t
Qg_t = og_t * G_t     G_{t+1} = G_t - Qg_t + R_t
Q_t  = Qs_t + Qg_t
```

Each gate is a learned scale times a sigmoid in its own store; `ol` reads PET instead.

## Every parameter and where it came from

From `MCPmodel3/model_epoch619.pt` unless stated.

| Parameter | Value | Note |
| --- | --- | --- |
| primary `weight_r_yom`, `_yom_gw`, `_ylm`, `_yfm` | 2.3979, 0.7041, 1.9713, 2.5880 | softmax to 0.3283 (surface), 0.0604 (recharge), 0.2143 (loss), 0.3971 (retained) |
| primary `bias_b0_yom`, `weight_b1_yom` | -5.1488, 2.7403 | surface gate |
| primary `bias_b0_yom_gw`, `weight_b1_yom_gw` | -1.5712, 3.6248 | recharge gate |
| primary `bias_b0_ylm`, `weight_b2_ylm` | -3.3837, 0.9848 | loss gate |
| surface `weight_r_yom`, `weight_r_yfm` | 0.3160, 1.5886 | scale 0.2188 |
| surface `bias_b0_yom`, `weight_b1_yom` | -11.0208, 0.3028 | |
| groundwater `weight_r_yom`, `weight_r_yfm` | 0.7780, -0.2251 | scale 0.7317 |
| groundwater `bias_b0_yom`, `weight_b1_yom` | -2.9328, 4.3223 | |
| `c_mean`, `c_std` | 412.9139085, 77.49943867 | defaults of `evaluate_model3_mca5.py` |
| `SW_mean`, `SW_std` | 3.051248249, 3.365261828 | same script |
| `GW_mean`, `GW_std` | 1.87176748, 2.085571043 | same script |
| `Ini_C` | 1.6871063 | same script: initial groundwater node |
| PET normalisation | 2.9086, 1.8980 | hard-coded in the upstream forward |
| initial primary and surface stores | 0 | upstream forward |
| `spinLen`, `traintimeLen` | 0, row count (upstream: 1095, 8400) | only slice the observations for a diagnostic the adapter does not read; `y_obs` is zeros |

The checkpoint loads with `strict=True`. Run as the evaluation script runs it on the Leaf River,
the adapter's model gives NSE 0.769 / 0.730 / 0.758 on the train / selection / test flags. The
repository publishes no number to compare that with.

## What is reported

`pr` (echoed), `evspsbl` (the constrained loss E), `mrro` (the total Q), in mm/day; and every
store, in mm, at the end of the row:

| Column | Store | Why there |
| --- | --- | --- |
| `mrso` | primary store C | the proposal's soil store |
| `gw` | groundwater node G | upstream's `MCPBRNN_GWVariant_Routing`, fed by recharge |
| `channel` | surface node S | upstream's `MCPBRNN_SW_Variant_Routing`: generated runoff not yet released |

`gw` and `channel` are additions to the proposal's list. A store left out would open the budget
by exactly the water in it. `channel` is this probe's definition of S even though S turns out
to be very slow (below); filing it under `gw` would only dodge `routing-conservation`.

Upstream's `c_n[t]` is each store *before* row t. The harness differences states against the
row's fluxes, so row t reports the store after it, read from `c_n[t+1]`: the model is given one
extra row of no weather at the end of the record and that row's outputs are thrown away.

## Zero stores, and what is not reported

`zero_stores` in `run.json`: `snw` (no snow store; temperature is never read) and `canopy` (no
interception store). Both are reported as identically zero, so the budget probes score the model.

Named under `not_reported`. The model has none of these, and none is filled in.

| Output | Why not |
| --- | --- |
| `dis` | the model knows no catchment area |
| `gwex` | every node conserves mass; nothing crosses the boundary but Q and E |
| energy fluxes | no energy balance |

## Choices, and how much each moves the results

**Spin-up.** Upstream repeats the first water year three times in front of the record. The
adapter does the same with the first 365 days of the rows it is given, and none when the record is shorter than a year: it has no year to repeat, and all of it may be scored, so replaying it would set row 0's store from weather that comes later. Those rows run before row 0 and are never reported. On daily probes that is
three times the case's own 365-day spin-up, nothing later than it, so it is causal. The table is the adapter run directly under the same harness and gate seeds, on every probe that scores this model, after the zero-store rule; at three repeats that run reproduces every archived verdict.

| Repeats | Result |
| --- | --- |
| 0 (this probe) | `resolution-invariance` PASS (`mrso_end` 9.2% of rain); `dry-down` still FAIL, at 17.0% rather than 3.7%; no verdict changes |
| 1 | `resolution-invariance` **FAIL** (`mrro` 10.2%); no other verdict changes |
| 3 (other probes) | `resolution-invariance` PASS (`mrro` 6.7%) |

The `resolution-invariance` row above is from before a short record stopped being replayed. That probe's record is
40 days, all of it weather the model is later scored on, so replaying it set row 0's store from what comes later. It
now runs as the 0 row. Every other probe has a full first year and still runs as the 3 row.

**Step scaling.** The upstream classes know one step. The resolution probe runs them at PT1H
regardless, so below a day the adapter runs **its own copy of the upstream loops**
(`forward_step_scaled`, `routing_step_scaled`), not the repository's code. They differ in three
ways: rain is a depth per row, each per-day gate becomes 1 - (1 - k)^dt, and the PET cap is a
depth per row. The gates still read PET in mm/day. At dt = 1 the copy matches the upstream
forward to single-precision round-off (at most 1.1e-4 mm on the stores over 3000 rows). Every
daily row goes through the upstream classes. The manifest declares PT1D only, the step in the
proposal.

**Precision.** The model runs in float32, as upstream does; see the closure figure under the verdict.

**Window.** `window_days: full`, since a record takes seconds. `warming-response` was also checked on the default window and passes either way.

**Seed.** Recorded and unused; the model is deterministic.

## Verdict

Gate seeds, Apptainer, 2026-10-31, version `5944b5f-epoch619.4`: **FAIL (VIOLATION), 11 of 21 scored**, 16 N/A
(15 INCOMPLETE, 1 INCOMPATIBLE: `routing-lag-consistency` needs catchment attributes the model does not
read). The INCOMPLETE rows ask for output this model does not have: heat fluxes and radiation (six energy probes), river-aquifer exchange, `snm`, `gwex`, `stage` and `dis`.

`closure` passes on all eight probes that score it, worst residual 4.2e-7 of the rain, which is single-precision
round-off in updates that conserve mass by construction. Trivial passes for this model: `area-invariance` (reads no
area), `phase-counterfactual` (reads no temperature) and `time-origin-invariance` (reads no calendar), as are the
bounds on the zero `snw` and `canopy`.

Six of the ten FAILs are `state_bounds` and nothing else.

| Probe | Result | Detail |
| --- | --- | --- |
| energy/pet-consistency | PASS | 0.95 of demand when wet, 0.49 when dry |
| mass/antecedent-monotonicity | **FAIL** | +0.1 mm on a 60 mm storm; 1.2 mm expected |
| mass/area-invariance | PASS | trivially: the model reads no area |
| mass/catchment-closure | **FAIL** | `state_bounds`: `mrso` 135-480 mm against 320 |
| mass/causality | PASS | nothing moves before the storm |
| mass/dry-down | **FAIL** | `mrro` rises 3.7% in a week with no rain |
| mass/extreme-event-closure | **FAIL** | `state_bounds`: `mrso` 68-546 mm against 320; event closure passes |
| mass/extreme-rain | PASS | returns 1.00 of 997 mm |
| mass/human-abstraction | **FAIL** | 380 mm of 380 unaccounted for; `state_bounds`: `mrso` 80-500 mm against 320 |
| mass/multi-decadal-drift | **FAIL** | `state_bounds`: `mrso` 349-393 mm against 320; no drift |
| mass/phase-counterfactual | PASS | trivially: the model reads no temperature |
| mass/precipitation-counterfactual | **FAIL** | `state_bounds`: `mrso` 123-474 mm against 320; the partition passes |
| mass/resolution-invariance | PASS | `mrso_end` 9.2% of rain |
| mass/response-nonnegativity | PASS | no dip |
| mass/runoff-bounds | PASS | runoff 0.39 of rain |
| mass/spinup-cycle-invariance | **FAIL** | `state_bounds`: `mrso` 165-426 mm against 320; the invariance passes |
| mass/steady-state | PASS | settles at 0.846 mm/day |
| mass/time-origin-invariance | PASS | trivially: the model reads no calendar |
| mass/ungauged-basin-closure | **FAIL** | `state_bounds` on 8 of 12 seeds: `mrso` 211-446 mm against 124 |
| mass/warming-response | PASS | -0.32 and -0.43 per unit of demand |
| momentum/routing-conservation | **FAIL** | `channel` 62.4 mm against a 1.9 mm bound |

## The failures

All verified with targeted runs. All are the model, not the adapter.

**state_bounds, on seven probes: the primary store has no capacity, and has to climb past the bound to pass water.** The
probes take the bound from `soil_capacity_mm` in `static.json`: 320 mm, or a draw of 124-682 mm on
`ungauged-basin-closure`. Nothing in the model can use it. The checkpoint's tensors are gate scales, biases and weights,
the store has no upper limit, and the adapter never opens `static.json`. At 320 mm the surface and recharge gates release
0.02 and 0.05 mm/day, at 400 mm 0.5 and 2.5, and they half-open at 559 and 447 mm, while the probe's records need
0.65-0.90 mm/day. The store's mean (274-304 mm on `catchment-closure`) is under the bound, unlike the other MCP models,
but it spans 78-480 mm across the seeds and exceeds 320 mm on up to about half the rows. On the Leaf River, where it was
trained, it spans 65-577 mm (mean 322), above 320 mm on 60% of days. On `ungauged-basin-closure` its peak is 408-446 mm
whatever the draw, and the four seeds that pass are the ones that drew 545, 656, 679 and 682 mm. The probes put no upper
bound on `gw` and `channel`.

The adapter could not honestly do otherwise. There is no capacity to hand `soil_capacity_mm` to, and writing it into
the normalisation or the gates would make a different model from the one trained. Subtracting a datum would not fit
either: the store's range on `catchment-closure` is 345-389 mm wide on each seed, wider than the
bound. Moving the water above the bound into `gw` would invent a partition the model does not have. So this is
a finding about the model, as with `parFC` in `models/dhbv2`, except that here not even a learned capacity exists. It
does not depend on the spin-up: the same probe fails the same way at zero, one and three repeats.

**human-abstraction: the model has no withdrawal term.** Structural, and verified from the row and the
code: the adapter reads `pr` and `pet` and nothing else, the paired runs are identical (`evspsbl`, `mrro` and
storage all +0.0 mm), and the whole 380 mm goes unaccounted for. AGENTS.md scores this as a failure, not as N/A.
Subtracting the column in the adapter would add a process the published model does not have, so nothing was added.

**routing-conservation: the surface node is a near-dead pool.** Its gate is
0.2188 * sigmoid(-11.02 + 0.3028 (S - 3.05) / 3.37), which at S = 60 mm releases 6e-4 of the store
a day, a residence of about 1700 days, and at S = 30 mm releases 4e-5. On the three seeds S sits
at means of 41-55 mm while median runoff is 0.02-0.04 mm/day, so it holds far more than 15 days
of recent peak flow. The same holds on the model's own catchment: on the Leaf River S ranges over
58-205 mm (mean 74) and exceeds 15 times the preceding 30 days' peak on 68% of rows. So this is
what training produced, not a result of the synthetic climate. The node passes a storm through
only once S is high, so it behaves as a threshold store rather than a 15-day hydrograph.

**dry-down: the groundwater node keeps filling after the rain stops.** On the failing seed G
rises from 1.27 to 1.38 mm over the first three rainless weeks while the primary store falls from
311 to 299 mm, because recharge still exceeds release, and the node's steep gate (4.32 per 2.09 mm)
turns that into runoff rising from 0.025 to 0.030 mm/day. It is a delay longer than the probe's
7-day allowance, and small: 3.7% against a 2% limit, the worst of the three seeds.

**antecedent-monotonicity: the gates are shut and the surface node absorbs what passes.** The
storm falls on a primary store of 122-190 mm (dry) or 208-285 mm (wet), below the surface gate's
opening, and S barely moves (63.58 against 63.61 mm on the first seed). Extra runoff in 30 days is
+0.02, +0.07 and +2.79 mm against 1.2 mm required. The third seed passes because its store reaches
the gate.

## Packaging

The image is built like `models/mcp1`: `python:3.10-slim` at a fixed digest, the repository's
torch 2.9.0 and numpy 2.2.2 with torch's dependencies pinned, and only `MCPBRNN_lib_tools`, the
checkpoint and `LICENSE` from the pinned commit. No network at run time. A case takes about 2.5 s
plus container start. `MCP_SPINUP_REPEATS` exists only for the table above; nothing sets it in
the container.

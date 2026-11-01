# mcp2: MCPModel2, the single-node MCP with state-dependent mass relaxation

MC{O_sigma L_sigma^con M_sigma_r^R} of Wang & Gupta (2024, *WRR* 60(4), e2023WR036461),
proposed in issue #90. MCPModel1's store and gates, plus a learned term that adds water below
a learned store level and removes it above. Trained on the Leaf River.

Source: [YuanHWang/5MCPmodelnet_HydroTuring](https://github.com/YuanHWang/5MCPmodelnet_HydroTuring)
at `5944b5f23adb7265a84fd10d62a8bd5501e25208`, Apache-2.0 for code and weights. The
checkpoint ships in the repository; nothing was trained for this evaluation.

```
Q_t = oo_t * C_t             oo_t = s_o * sigmoid(b_o + w_o * (C_t - c_mean) / c_std)
E_t = olc_t * C_t            olc_t = min(ol_t, PET_t / C_t)
                             ol_t = s_l * sigmoid(b_l + w_l * (PET_t - 2.9086) / 1.8980)
M_t = ov_t * |C_t - 500 b|   ov_t = min(r * tanh(s * (C_t / 500 - b)), 1 - oo_t - olc_t)
C_{t+1} = C_t - Q_t - E_t + P_t - M_t
```

## Every parameter and where it came from

| Parameter | Value | Source |
| --- | --- | --- |
| `weight_r_yom`, `weight_r_ylm`, `weight_r_yfm` | 1.4279, 2.3210, 4.0176 | `MCPmodel2/model_epoch3.pt`; softmax to s_o = 0.0596, s_l = 0.1457, retained 0.7947 |
| `bias_b0_yom` (b_o), `weight_b1_yom` (w_o) | -7.2404, 2.6522 | checkpoint |
| `bias_b0_ylm` (b_l), `weight_b2_ylm` (w_l) | -3.3472, 1.7004 | checkpoint |
| `weight_r_yvm` | -3.0175 | checkpoint; sigmoid to r = 0.0466 |
| `weight_s_yvm` | -2.0513 | checkpoint; exp to s = 0.1286 |
| `bias_b0_yrm` (b) | 1.0974 | checkpoint; the level 500 b = 548.7 mm |
| `c_mean`, `c_std` | 475.0253683, 89.86101728 | defaults of `evaluate_model2_mass_relax.py` |
| `scale_mr` | 500 | hard-coded in the upstream forward |
| PET normalisation | 2.9086, 1.8980 | hard-coded in the upstream forward |
| initial store | 0 | upstream forward |
| `spinLen`, `traintimeLen` | 0, row count (upstream: 1095, 8400) | only slice the observations for a diagnostic the adapter does not read; `y_obs` is zeros |

The checkpoint loads with `strict=True`. Run as the evaluation script runs it on the Leaf
River, the adapter's model gives NSE 0.760 / 0.772 / 0.764 on the train / selection / test
flags. The repository publishes no number to compare that with.

## What is reported

`pr` (echoed), `evspsbl` (the constrained loss E), `mrro` (Q), `gwex` (-M), in mm/day; `mrso`
(the store C at the end of the row), in mm.

**`gwex` is an addition to the proposal's list**, which names only `evspsbl`, `mrro` and `mrso`.
The relaxation moves water across the model's boundary and nothing else accounts for it, so
it is the exchange AGENTS.md says to declare, positive into the catchment. It is computed from
the model's own returned gate and learned level, `-gate_ov * |c_n - 500 b|`, not from a budget
residual. Leaving it out would hide the term this model exists to test.

Upstream's `c_n[t]` is the store *before* row t. The harness differences a state against the
row's fluxes, so row t reports C_{t+1}, read from `c_n[t+1]`: the model is given one extra row
of no weather at the end of the record and that row's outputs are thrown away.

## Zero stores, and what is not reported

`zero_stores` in `run.json`: `snw` (no snow store; temperature is never read) and `canopy` (no
interception store). Both are reported as identically zero, so the budget probes score the model.

Named under `not_reported`. The model has none of these, and none is filled in.

| Output | Why not |
| --- | --- |
| `gw`, `channel` | one undifferentiated store, no routing |
| `dis` | the model knows no catchment area |
| energy fluxes | no energy balance |

## Choices, and how much each moves the results

**Spin-up.** Upstream repeats the first water year three times in front of the record, from an
empty store. The adapter does the same with the first 365 days of the rows it is given, and none when the record is shorter than a year: it has no year to repeat, and all of it may be scored, so replaying it would set row 0's store from weather that comes later. Those rows run before row 0 and are never reported. On
daily probes that is three times the case's own 365-day spin-up, nothing later than it, so it
is causal. The table is the adapter run directly under the same harness and gate seeds, on every probe that scores this model, after the zero-store rule; at three repeats that run reproduces every archived verdict.

| Repeats | Result |
| --- | --- |
| 0 (this probe) | `resolution-invariance` **FAIL** (`mrso_end` 15.4% of rain); no other verdict changes |
| 1 | `resolution-invariance` **FAIL** (`mrso_end` 21.1%); no other verdict changes |
| 3 (other probes) | `resolution-invariance` PASS (`evspsbl` 5.9%) |

The `resolution-invariance` row above is from before a short record stopped being replayed. That probe's record is
40 days, all of it weather the model is later scored on, so replaying it set row 0's store from what comes later. It
now runs as the 0 row. Every other probe has a full first year and still runs as the 3 row.

**Step scaling.** The upstream class knows one step. The resolution probe runs it at PT1H
regardless, so below a day the adapter runs **its own copy of the upstream loop**
(`forward_step_scaled`), not the repository's code. It differs in these ways: rain is a depth
per row, each per-day gate becomes 1 - (1 - k)^dt (the signed relaxation gate on its magnitude,
keeping its sign), and the PET cap is a depth per row. The gates still read PET in mm/day. At
dt = 1 the copy matches the upstream forward to single-precision round-off (at most 2.4e-4 mm
on the store over 3000 rows). Every daily row goes through the upstream class. The manifest
declares PT1D only, the step in the proposal.

**Precision.** The model runs in float32, as upstream does; see the closure figure under the verdict.

**Window.** `window_days: full`, since a record takes seconds. `warming-response` was also checked on the default window and passes either way.

**Seed.** Recorded and unused; the model is deterministic.

## Verdict

Gate seeds, Apptainer, 2026-10-31, version `5944b5f-epoch3.4`: **FAIL (VIOLATION), 9 of 20 scored**, 17 N/A
(15 INCOMPLETE; 2 INCOMPATIBLE: `routing-lag-consistency` needs catchment attributes the model does not read, and
`exchange-response` needs `gwh`, not declared because the relaxation reads the store rather than a head).
The INCOMPLETE rows ask for output this model does not have: heat fluxes and radiation (six energy probes), river-aquifer exchange, `snm`, `channel`, `stage` and `dis`.

`closure` passes on all eight probes that score it, `gwex` counted as a source, worst residual 8.7e-7 of the rain,
which is single-precision round-off. Trivial passes for this model: `area-invariance` (reads no area),
`phase-counterfactual` (reads no temperature) and `time-origin-invariance` (reads no calendar), as are the bounds on
the zero `snw` and `canopy`.

Five of the ten FAILs are `state_bounds` and nothing else.

| Probe | Result | Detail |
| --- | --- | --- |
| energy/pet-consistency | **FAIL** | 0.75 of demand when dry, 0.75 when wet |
| mass/antecedent-monotonicity | **FAIL** | +0.3 mm on a 60 mm storm; 1.2 mm expected |
| mass/area-invariance | PASS | trivially: the model reads no area |
| mass/catchment-closure | **FAIL** | `state_bounds`: `mrso` 297-730 mm against 320 |
| mass/causality | PASS | nothing moves before the storm |
| mass/dry-down | **FAIL** | `mrso` rises 2.0% in a week with no rain |
| mass/extreme-event-closure | **FAIL** | `state_bounds`: `mrso` 233-840 mm against 320; event closure passes |
| mass/extreme-rain | PASS | returns 0.94 of 526 mm |
| mass/human-abstraction | **FAIL** | 380 mm of 380 unaccounted for; `state_bounds`: `mrso` 257-726 mm against 320 |
| mass/multi-decadal-drift | **FAIL** | `state_bounds`: `mrso` 585-627 mm against 320; no drift |
| mass/phase-counterfactual | PASS | trivially: the model reads no temperature |
| mass/precipitation-counterfactual | **FAIL** | responses account for 0.83-0.89 of the change in rain; `state_bounds`: `mrso` 268-709 mm against 320 |
| mass/resolution-invariance | **FAIL** | `mrso_end` 15.4% of rain: the stores end the month apart, from the 10-day spin-up alone |
| mass/response-nonnegativity | PASS | no dip |
| mass/runoff-bounds | PASS | runoff 0.37 of rain, declared exchange +0.06 of rain |
| mass/spinup-cycle-invariance | **FAIL** | `state_bounds`: `mrso` 322-665 mm against 320; the invariance passes |
| mass/steady-state | PASS | settles at 1.114 mm/day |
| mass/time-origin-invariance | PASS | trivially: the model reads no calendar |
| mass/ungauged-basin-closure | **FAIL** | `state_bounds` on 10 of 12 seeds: `mrso` 379-692 mm against 124 |
| mass/warming-response | PASS | -0.38 and -0.49 per unit of demand |

## The failures

All verified with targeted runs. All are the model, not the adapter.

**state_bounds, on seven probes: the store has no capacity, and its gate opens far above the bound.** The probes take
the bound from `soil_capacity_mm` in `static.json`: 320 mm, or a draw of 124-682 mm on `ungauged-basin-closure`. Nothing
in the model can use it. The checkpoint's ten tensors are gate scales, biases and weights and the relaxation's rate and
level, the store has no upper limit, and the adapter never opens `static.json`. The store runs wherever the learned
output gate lets out the runoff the climate requires. At 320 mm that gate releases 0.0001 mm/day, half-opens at 720 mm,
and the probe's records need 0.63-0.86 mm/day. The relaxation pulls the same way, adding water below its learned level
of 548.7 mm. So on `catchment-closure` the store spans 232-730 mm across the seeds, and on the Leaf River, where it was
trained, 214-808 mm (mean 543), above 320 mm on 97% of days. On `ungauged-basin-closure` its peak is 634-692 mm whatever
the draw, and the two seeds that pass are the ones that drew 679 and 682 mm.

The adapter could not honestly do otherwise. There is no capacity to hand `soil_capacity_mm` to, and writing it into
the normalisation or the gates would make a different model from the one trained. Subtracting a datum would not fit
either: the store's range on `catchment-closure` is 433-475 mm wide on each seed, wider than the
bound. Moving the water above the bound into `gw` would invent a partition the model does not have. So this is
a finding about the model, as with `parFC` in `models/dhbv2`, except that here not even a learned capacity exists. It
does not depend on the spin-up: the same probe fails the same way at zero, one and three repeats.

**human-abstraction: the model has no withdrawal term.** Structural, and verified from the row and the
code: the adapter reads `pr` and `pet` and nothing else, the paired runs are identical (`evspsbl`, `mrro` and
storage all +0.0 mm), and the whole 380 mm goes unaccounted for. AGENTS.md scores this as a failure, not as N/A.
Subtracting the column in the adapter would add a process the published model does not have, so nothing was added. The relaxation reads only the store, so it does not stand in for one.

**precipitation-counterfactual: the share the criterion cannot see is the relaxation.** Beyond `state_bounds`,
`counterfactual_response` fails every variant on all three seeds: evaporation, runoff and storage account for 0.83-0.89
of the rain added or removed, against 1 within 0.10. More rain raises the store, and the store-dependent relaxation then
adds less. The change in the declared `gwex` is -0.11 to -0.17 of the change in rain, and with it the sum is 1.000 to
within 2e-6 on all nine pairs. `models/mcp1`, which has no relaxation, accounts for 1.000 on the same record. So no water
is lost: it is in a column this criterion does not count, although `closure` and the runoff ratio on the same suite do.
What stands about the model is that a learned exchange takes 11-17% of any change in rain. Relabelling it would hide the
term, so it stays. `state_bounds` fails this probe either way.

**dry-down: the relaxation makes water.** Without rain the store falls below the learned level
of 548.7 mm, where the term is a source that grows as the store falls. At the worst block on
the three seeds the store is 147-161 mm, `gwex` is +1.8 to +1.9 mm/day, and evaporation is only
0.4-1.0 mm/day. So the store climbs with no rain falling. The rise is the term the model was
built with, measured in its own declared flux.

**pet-consistency: the relaxation keeps the dry regime from being dry.** Evaporation is
min(ol C, PET), so it depends on the store only through C. In the driest fifth of the store's
range the relaxation is adding +0.56 to +0.67 mm/day (-0.08 to -0.10 in the wettest, +0.14 to
+0.19 overall). On the failing seed the ratio to demand is 0.74 wet and 0.75 dry, against a
required contrast of 0.05. Switching the term off in memory, a diagnosis and not a
configuration anyone should run, gives 0.74 and 0.66, which clears it, and the other two seeds
widen the same way.

**antecedent-monotonicity: the output gate is shut, as in MCPModel1.** The storm falls on a
store of 292-381 mm (dry) or 404-495 mm (wet). The gate opens higher than that, so extra runoff
in 30 days is +0.10, +0.29 and +4.17 mm against 1.2 mm required. The third seed passes because
its store reaches the gate.

## Packaging

The image is built like `models/mcp1`: `python:3.10-slim` at a fixed digest, the repository's
torch 2.9.0 and numpy 2.2.2 with torch's dependencies pinned, and only `MCPBRNN_lib_tools`, the
checkpoint and `LICENSE` from the pinned commit. No network at run time. A case takes about 2 s
plus container start. `MCP_SPINUP_REPEATS` exists only for the table above; nothing sets it in
the container.

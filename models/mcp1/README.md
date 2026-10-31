# mcp1: MCPModel1, the single-node Mass-Conserving Perceptron

MC{O_sigma L_sigma^con} of Wang & Gupta (2024, *WRR* 60(4), e2023WR036461), proposed
in issue #89. One store and two gated outflows, trained on the Leaf River.

Source: [YuanHWang/5MCPmodelnet_HydroTuring](https://github.com/YuanHWang/5MCPmodelnet_HydroTuring)
at `5944b5f23adb7265a84fd10d62a8bd5501e25208`, Apache-2.0 for code and weights. The
checkpoint ships in the repository; nothing was trained for this evaluation.

```
Q_t = oo_t * C_t                  oo_t = s_o * sigmoid(b_o + w_o * (C_t - c_mean) / c_std)
E_t = olc_t * C_t                 olc_t = min(ol_t, PET_t / C_t)
C_{t+1} = C_t - Q_t - E_t + P_t   ol_t = s_l * sigmoid(b_l + w_l * (PET_t - 2.9086) / 1.8980)
```

Row t's outflows read the store as it stood before row t's rain, so rain reaches the
outflows one row later.

## Every parameter and where it came from

| Parameter | Value | Source |
| --- | --- | --- |
| `weight_r_yom`, `weight_r_ylm`, `weight_r_yfm` | 0.8, 1.9518, 3.5941 | `MCPmodel1/model_epoch27.pt`; softmax to s_o = 0.0488, s_l = 0.1542, retained 0.7970 |
| `bias_b0_yom` (b_o), `weight_b1_yom` (w_o) | -10.0975, 2.7342 | checkpoint |
| `bias_b0_ylm` (b_l), `weight_b2_ylm` (w_l) | -3.6252, 0.8528 | checkpoint |
| `c_mean`, `c_std` | 412.9139085, 77.49943867 | defaults of `evaluate_model1_single_node.py` |
| PET normalisation | 2.9086, 1.8980 | hard-coded in the upstream forward (`MCNZoo.py`) |
| initial store | 0 | upstream forward |
| `spinLen`, `traintimeLen` | 0, row count (upstream: 1095, 8400) | only slice the observations for a diagnostic the adapter does not read; `y_obs` is zeros |

The checkpoint loads with `strict=True`. Run the way the evaluation script runs it on the
Leaf River, the adapter's model gives NSE 0.766 / 0.764 / 0.768 on the train / selection /
test flags. The repository publishes no number to compare that with.

## What is reported

`pr` (echoed), `evspsbl` (the constrained loss E), `mrro` (Q), in mm/day; `mrso` (the store
C at the end of the row), in mm. The proposal's own list.

Upstream's `c_n[t]` is the store *before* row t. The harness differences a state against the
row's fluxes, so row t reports C_{t+1}, read from `c_n[t+1]`: the model is given one extra row
of no weather at the end of the record and that row's outputs are thrown away. No term is
computed as a residual.

## Zero stores, and what is not reported

`zero_stores` in `run.json`: `snw` (no snow store; temperature is never read) and `canopy` (no
interception store). Both are reported as identically zero, so the budget probes score the model.

Named under `not_reported`. The model has none of these, and none is filled in.

| Output | Why not |
| --- | --- |
| `gw`, `channel` | one undifferentiated store, no routing |
| `dis` | the model knows no catchment area |
| `gwex` | the store update conserves mass; nothing crosses the boundary but Q and E |
| energy fluxes | no energy balance |

## Choices, and how much each moves the results

**Spin-up.** Upstream repeats the first water year three times in front of the record, from an
empty store. The adapter does the same with the first 365 days of the rows it is given, or all
of them when the record is shorter. Those rows run before row 0 and are never reported. On
daily probes that is three times the case's own 365-day spin-up, nothing later than it, so it
is causal. The table is the adapter run directly under the same harness and gate seeds, on every probe that scores this model, after the zero-store rule; at three repeats that run reproduces every archived verdict.

| Repeats | Result |
| --- | --- |
| 0 | `resolution-invariance` PASS (`mrso_end` 6.0% of rain); no other verdict changes |
| 1 | `resolution-invariance` **FAIL** (`mrso_end` 13.2%); no other verdict changes |
| 3 (used) | `resolution-invariance` PASS (`mrro` 3.2%) |

So this model's `resolution-invariance` verdict depends on the spin-up. Treat that row as
not robust. A first version cycled a short record until it filled three years, feeding the
model about 2.6 m/yr of rain that no record carried; it was dropped for all five MCP models
because it invented forcing (3.2% here, and see `models/mcp4`).

**Step scaling.** The upstream class knows one step. The resolution probe runs it at PT1H
regardless, so below a day the adapter runs **its own copy of the upstream loop**
(`forward_step_scaled`), not the repository's code. It differs in three ways: rain is a depth
per row, each per-day gate becomes 1 - (1 - k)^dt, and the PET cap is a depth per row. The
gates still read PET in mm/day. At dt = 1 the copy matches the upstream forward to
single-precision round-off (at most 1.2e-4 mm on the store over 3000 rows). Every daily row
goes through the upstream class. The manifest declares PT1D only, the step in the proposal.

**Precision.** The model runs in float32, as upstream does; see the closure figure under the verdict.

**Window.** `window_days: full`, since a record takes seconds. `warming-response` was also checked on the default window and passes either way (-0.35/-0.51 per unit of demand on the full
record, -0.68/-0.72 on the default window).

**Seed.** Recorded and unused; the model is deterministic.

## Verdict

Gate seeds, Apptainer, 2026-10-31, version `5944b5f-epoch27.3`: **FAIL (VIOLATION), 12 of 20 scored**, 17 N/A
(16 INCOMPLETE, 1 INCOMPATIBLE: `routing-lag-consistency` needs catchment attributes the model does not
read). The INCOMPLETE rows ask for output this model does not have: heat fluxes and radiation (six energy probes), river-aquifer exchange, `snm`, `gwex`, `channel`, `stage` and `dis`.

**`resolution-invariance` is not robust.** Its PASS turns to FAIL at one spin-up repeat, and back at zero (table
above). Read that row as undecided.

`closure` passes on all eight probes that score it, worst residual 3.3e-7 of the rain, which is single-precision
round-off in an update that conserves mass by construction. Trivial passes for this model: `area-invariance` (reads
no area), `phase-counterfactual` (reads no temperature) and `time-origin-invariance` (reads no calendar), as are the
bounds on the zero `snw` and `canopy`.

Six of the eight FAILs are `state_bounds` and nothing else.

| Probe | Result | Detail |
| --- | --- | --- |
| energy/pet-consistency | PASS | 0.93 of demand when wet, 0.59 when dry |
| mass/antecedent-monotonicity | **FAIL** | +0.4 mm on a 60 mm storm; 1.2 mm expected |
| mass/area-invariance | PASS | trivially: the model reads no area |
| mass/catchment-closure | **FAIL** | `state_bounds`: `mrso` 337-725 mm against 320 |
| mass/causality | PASS | nothing moves before the storm |
| mass/dry-down | PASS | nothing rises |
| mass/extreme-event-closure | **FAIL** | `state_bounds`: `mrso` 243-850 mm against 320; event closure passes |
| mass/extreme-rain | PASS | returns 1.00 of 997 mm |
| mass/human-abstraction | **FAIL** | 380 mm of 380 unaccounted for; `state_bounds`: `mrso` 303-722 mm against 320 |
| mass/multi-decadal-drift | **FAIL** | `state_bounds`: `mrso` 571-624 mm against 320; no drift |
| mass/phase-counterfactual | PASS | trivially: the model reads no temperature |
| mass/precipitation-counterfactual | **FAIL** | `state_bounds`: `mrso` 329-707 mm against 320; the partition passes |
| mass/resolution-invariance | PASS | `mrro` 3.2% of rain; **depends on the spin-up** |
| mass/response-nonnegativity | PASS | no dip |
| mass/runoff-bounds | PASS | runoff 0.36 of rain |
| mass/spinup-cycle-invariance | **FAIL** | `state_bounds`: `mrso` 382-662 mm against 320; the invariance passes |
| mass/steady-state | PASS | settles at 0.887 mm/day |
| mass/time-origin-invariance | PASS | trivially: the model reads no calendar |
| mass/ungauged-basin-closure | **FAIL** | `state_bounds` on 10 of 12 seeds: `mrso` 438-686 mm against 124 |
| mass/warming-response | PASS | -0.35 and -0.51 per unit of demand |

## The failures

All verified with targeted runs. All are the model, not the adapter.

**state_bounds, on seven probes: the store has no capacity, and its gate opens far above the bound.** The probes take
the bound from `soil_capacity_mm` in `static.json`: 320 mm, or a draw of 124-682 mm on `ungauged-basin-closure`. Nothing
in the model can use it. The checkpoint's seven tensors are gate scales, biases and weights, the store has no upper
limit, and the adapter never opens `static.json`. The store runs wherever the learned output gate lets out the runoff
the climate requires. At 320 mm that gate releases under 0.001 mm/day, half-opens at 699 mm, and the probe's records
need 0.56-0.81 mm/day. So on `catchment-closure` the store spans 265-725 mm across the seeds, and on the Leaf River, where
it was trained, 239-814 mm (mean 549), above 320 mm on 98% of days. On `ungauged-basin-closure` its peak is 618-686 mm
whatever the draw, and the two seeds that pass are the ones that drew 679 and 682 mm.

The adapter could not honestly do otherwise. There is no capacity to hand `soil_capacity_mm` to, and writing it into
the normalisation or the gates would make a different model from the one trained. Subtracting a datum would not fit
either: the store's range on `catchment-closure` is 388-437 mm wide on each seed, wider than the
bound. Moving the water above the bound into `gw` would invent a partition the model does not have. So this is
a finding about the model, as with `parFC` in `models/dhbv2`, except that here not even a learned capacity exists. It
does not depend on the spin-up: the same probe fails the same way at zero, one and three repeats.

**human-abstraction: the model has no withdrawal term.** Structural, and verified from the row and the
code: the adapter reads `pr` and `pet` and nothing else, the paired runs are identical (`evspsbl`, `mrro` and
storage all +0.0 mm), and the whole 380 mm goes unaccounted for. AGENTS.md scores this as a failure, not as N/A.
Subtracting the column in the adapter would add a process the published model does not have, so nothing was added.

**antecedent-monotonicity: the output gate is shut at the stores this climate produces.**
Verified with a per-seed run. The gate is a steep sigmoid in the store:

| C (mm) | 400 | 500 | 550 | 600 | 650 | 700 | 800 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Q (mm/day) | 0.001 | 0.022 | 0.14 | 0.86 | 4.8 | 17 | 38 |

The storm falls in midsummer on a store of 340-407 mm (dry) or 442-510 mm (wet), and leaves
it at 396-463 and 497-565 mm. Both runs stay below the gate's opening, so little runs off in
30 days, and the wet run's extra 100 mm is still in the store when the window closes. Extra
runoff is +0.15, +0.44 and +4.70 mm on the three seeds, against 1.2 mm required; the third
passes because its store reaches the gate. The model's memory is real but expresses itself
after the window. It was trained on about 1470 mm/yr of Leaf River rain; the probe gives
about 880.

## Packaging

The image is 275 MB: `python:3.10-slim` at a fixed digest, the repository's torch 2.9.0 and
numpy 2.2.2 with torch's dependencies pinned, and only `MCPBRNN_lib_tools`, the checkpoint and
`LICENSE` from the pinned commit. The base image's pip cannot read the PyTorch index's package
names, so pip is pinned first. No network at run time. A case takes about 1.5-2 s plus
container start. `MCP_SPINUP_REPEATS` exists only for the table above; nothing sets it in the
container.

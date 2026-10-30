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

## Not reported, and why

Named under `not_reported` in `run.json`. The model has none of these, and none is filled in.

| Output | Why not |
| --- | --- |
| `snw` | no snow store; temperature is never read |
| `canopy` | no interception store |
| `gw`, `channel` | one undifferentiated store, no routing |
| `dis` | the model knows no catchment area |
| `gwex` | the store update conserves mass; nothing crosses the boundary but Q and E |
| energy fluxes | no energy balance |

**What that does to the standing.** `mass/catchment-closure` and seven other budget probes
lack only `snw` and `canopy` from this model, so a model built to conserve mass is N/A (INCOMPLETE) on the
probe that asks. Off the harness, the reported columns close: at most 8e-5 mm in any row and
1e-3 mm cumulative over 5805 mm of rain (2000 rows, PT1D), which is single-precision
round-off. That number is not an archived verdict.

## Choices, and how much each moves the results

**Spin-up.** Upstream repeats the first water year three times in front of the record, from an
empty store. The adapter does the same with the first 365 days of the rows it is given, or all
of them when the record is shorter. Those rows run before row 0 and are never reported. On
daily probes that is three times the case's own 365-day spin-up, nothing later than it, so it
is causal. The table is the prototype under the same harness and gate seeds.

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

**Precision.** The model runs in float32, as upstream does; see the closure figure above.

**Window.** `window_days: full`, since a record takes seconds. On the one probe where this
differs, `warming-response` passes either way (-0.35/-0.51 per unit of demand on the full
record, -0.68/-0.72 on the default window).

**Seed.** Recorded and unused; the model is deterministic.

## Verdict

Gate seeds, Apptainer, 2026-10-30: **FAIL (VIOLATION), 11 of 12 scored**, 25 N/A (24 INCOMPLETE,
1 INCOMPATIBLE: `routing-lag-consistency` needs catchment attributes the model does not read).

| Probe | Result | Detail |
| --- | --- | --- |
| energy/pet-consistency | PASS | 0.93 of demand when wet, 0.59 when dry |
| mass/antecedent-monotonicity | **FAIL** | +0.4 mm on a 60 mm storm; 1.2 mm expected |
| mass/area-invariance | PASS | trivially: the model reads no area |
| mass/causality | PASS | nothing moves before the storm |
| mass/dry-down | PASS | nothing rises |
| mass/extreme-rain | PASS | returns 1.00 of 997 mm |
| mass/phase-counterfactual | PASS | trivially: the model reads no temperature |
| mass/resolution-invariance | PASS | `mrro` 3.2% of rain; depends on the spin-up, above |
| mass/response-nonnegativity | PASS | no dip |
| mass/runoff-bounds | PASS | runoff 0.36 of rain |
| mass/steady-state | PASS | settles at 0.887 mm/day |
| mass/warming-response | PASS | -0.35 and -0.51 per unit of demand |

## The failure

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

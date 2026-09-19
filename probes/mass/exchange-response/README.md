# mass/exchange-response

## What it asserts

`gwex` is the contract's confession channel. A model whose catchment gains or
loses water across its boundary — regional groundwater, an inter-basin transfer,
a prescribed withdrawal — declares that flux, and its budget closes over a term
it admits to instead of hiding it in the residual. `AGENTS.md` puts it plainly:
*declared, it is a source in the budget and the budget can close; hidden, it is
the residual.*

Nothing yet says what a declared exchange may look like. This probe does.

Write the budget **without** the exchange:

```
R0 = pr - evspsbl - mrro - d(mrso + snw + canopy + gw + channel)/dt   [mm/day]
```

so that the budget the suite checks is `R = R0 + gwex`, and `closure` requires
`|Σ R Δt| <= 0.05 Σ pr Δt`.

**A model that closes has `R ≈ 0` and therefore `gwex ≈ −R0` identically.**
Measured on a CREST implementation over a ten-year record, `max |R|` per step is
**2.3e-13 mm/day**. So comparing the declared flux with the residual it closes is
the same number twice, and neither the *magnitude* of `gwex` nor its timing on
its own can separate an honest exchange from an invented one. **A discriminator
therefore has to use information that is not algebraically fixed by closure.**
This probe uses a counterfactual: the same weather, the same model, the same
seed, and one change to the one driver the case prescribes.

### The assertion

The case supplies an external hydraulic head, `gwh` (metres), as a visible
forcing column, and runs the record three times — with the head as given,
**raised by 0.5 m**, and **lowered by 0.5 m** — with every other column
byte-identical. The shift begins with the first scored step; **the spinup is
identical across the three runs** (see *Timing* below for why that matters).

The probe **requires** `gwh`: a model that does not list it in `needs_forcing`
or `uses_forcing` is INCOMPATIBLE and not scored. Declaring it is a semantic
opt-in, stated in `AGENTS.md`:

> `gwh` is the prescribed external hydraulic head associated with the model's
> declared `gwex`; positive `gwex` is into the catchment. A model declaring
> `gwh` asserts that its external exchange responds monotonically to this
> hydraulic potential.

For such a model, with the integrated exchange over the scored record and the
gross movement of the control run,

```
G     = Σ_t gwex_t · Δt                  [mm]
TV(g) = Σ_t |gwex_t − gwex_{t−1}| · Δt   [mm, the control run's total variation]

G(raised)  − G(control)  >=  +max(s · TV(g), ε)
G(lowered) − G(control)  <=  −max(s · TV(g), ε)      s = 3e-4,  ε = 1e-8 · max(G_gross, 1 mm)
```

Raising the external head must bring more water in, and lowering it less, **by
at least a small share of how much the exchange the model itself declared moves
from day to day.** `ε` is for floating point and nothing else.

### Why a share, and why of the variation

An earlier draft asked for no minimum at all, on the argument that any absolute
floor would fail a genuine boundary of small conductance. That is true of an
absolute floor and false of a share: for a general-head boundary `Q = C (H − h)`
the response and everything about the exchange scale with `C` together, so a
ratio does not. Without a share, a model that keeps the accounting sink and adds
a head term at a millionth of the controls' conductance answers the head with
the right sign and passes; `reference_token_exchange` is that model.

The first share was taken against the **gross** exchange, `Σ|gwex|`, with the
derivation that in the fast limit the response is `S·Δh`, the gross is
`S·TV(H)`, storativity cancels too, and the ratio is bounded below by
`Δh / TV(H) ≈ 2e-3`. A native MODFLOW aquifer confirmed that limit exactly. It
also showed what the derivation had assumed: that the boundary is the aquifer's
*only* source. On an ordinary losing catchment — half the runoff recharging the
aquifer and draining out through the same boundary — the gross becomes the
throughput, of the order of the recharge and independent of `S`, while the
response stays `S·Δh`. The ratio then falls with storativity, and at
`S = 1 mm/m`, the top of the usual confined range, an exact, monotone,
budget-closing boundary measured **2.4e-4** and failed. That was a false FAIL,
not a disclosed cheat: every millimetre of that gross was the same boundary.

The gross is the wrong denominator because a steady throughput inflates it
without moving. The **total variation** `TV(g) = Σ|gwex_t − gwex_{t−1}|` does
not have that defect: throughput at a steady mean adds little to how much the
exchange moves day to day, and what it does add is the throughput's own
variation, which is bounded. On the same losing catchment, sweeping storativity
at a fixed ten-day time constant on `reference_recharge_exchange`:

| S (mm/m) | Δ raised (mm) | gross (mm) | TV (mm) | resp / gross | **resp / TV** |
| --- | --- | --- | --- | --- | --- |
| 10 | +3.5 | 1878 | 283 | 1.5e-3 | **1.0e-2** |
| 3 | +1.19 | 1772 | 111 | 6.3e-4 ✗ | **1.0e-2** |
| 1 | +0.50 | 1767 | 66 | 2.4e-4 ✗ | **6.3e-3** |
| 0.3 | +0.15 | 1767 | 56 | 8.2e-5 ✗ | **2.6e-3** |
| 0.1 | +0.05 | 1767 | 54 | 2.8e-5 ✗ | **9.2e-4** |

The gross floor at `1e-3` failed everything from `S = 3 mm/m` down — the whole
confined range. The variation floor at `3e-4` passes to `S = 0.1 mm/m`, a tenth
of the confined top, with a factor of three in hand there. The token cheat sits
at **1.1e-6**, three hundred times under the floor, because its variation is
the accounting error's, which is large, and its response is the token's, which
is not.

### What the floor excludes

No physical bound is claimed for the ratio. The variation has a floor of its
own, set by how the throughput varies (about 54 mm here), while the response
falls linearly with storativity, so a boundary of **very small storativity that
also carries a strongly varying throughput** falls under `3e-4` eventually — on
this catchment at `S = 0.03 mm/m`, a thirtieth of the confined top, where three
of five gate seeds pass (5/5 at 0.05, 0/5 at 0.02, checked independently). That
is a named class of false FAIL, and it belongs here rather than in the cheat
paragraph. It is far narrower than the class the gross floor excluded, and a
model in it can be told from the cheat by its diagnostics — a response that is
exactly `S·Δh` and a variation that is the recharge's — but the criterion does
not make that call.

The other limit is the mixture: a `gwex` that combines a small head-driven
part with a large, strongly varying unrelated one has a variation the unrelated
part inflates and a share that is honest but small. From outside it is
indistinguishable from a cheat that reads the head, keeps the sink, and sizes
its token term to `3e-4` of the sink's variation. Disclosed, and the price of
not asking the contract to split the head-driven part of `gwex` into its own
variable.

### Timing

The shift begins with the first scored step, not with the spinup, and this is
load-bearing. A responsive aquifer answers a shift by admitting the water that
raises its own head to meet the new external one — and then by nothing. An
earlier draft applied the shift from the start of spinup; a native MODFLOW 6.7.0
confined aquifer behind a general-head boundary, run on this probe's head series
with a five-day time constant, admitted its extra 5 mm during spinup, raised its
internal head by 0.5 m, and showed a paired difference of nearly zero inside the
scored window. It failed on all five seeds for having responded correctly. With
the shift at the scored start and the spinup identical, the same model gives
+5/−5 mm and passes. The criterion checks the spinup identity and refuses a
case that breaks it.

### What is reported and deliberately not gated

Everything two earlier drafts gated on: the magnitude-weighted reversal
fraction of the exchange, its count of significant reversals inside the
rainless windows, the prescribed head's turning points there, and the exchange's
gross and net size against precipitation. All of it is in the report; none of it
decides anything. Between two turning points of the external head the flux
`h_ext − h_catchment` crosses zero as often as the internal head moves, and in a
real aquifer the internal head does move, so no bound on reversal follows from
the external head; and a physical aquifer behind an oscillating boundary
reverses its flux with no rain at all, so no bound follows from rainlessness
either. Both were tried and both were wrong.

## The case

`generate.py` produces, from a seed and a variant, ten scored years of daily
weather behind 365 days of spinup, on the closure probe's climate: roughly
**843 mm/yr precipitation** and **759 mm/yr potential ET** (aridity ~0.9, the
energy-limited side of the Budyko curve, with a real seasonal snowpack).
Gamma-distributed depths on a seasonal occurrence probability; temperature an
annual cycle plus an AR(1) anomaly; potential ET temperature-driven with a
daylength factor.

Three features are specific to this probe:

*   **The external head `gwh`**: a 10 m datum, a 0.15 m annual cycle and a
    0.35 m cycle of 20 days, so the head turns about every ten days. The
    `raised` and `lowered` variants add and subtract 0.5 m **from the first
    scored step on**; the spinup is byte-identical across the three, and so is
    every other column, and the criterion checks both.
*   **Two rainless stretches, 90 and 75 days, labelled `_regime = dry`**, in
    the last two scored years and in summer. They decide nothing; the reversal
    diagnostics are computed inside them so a reviewer can see how an exchange
    behaves where direct precipitation forcing is absent. Columns beginning
    with an underscore are probe annotation, stripped before the model sees
    the forcing.
*   **The whole record scored.** `min_window_days` is 3650: an integrated
    response is a statement about the record, and storage carries across it.

## Why these criteria

| Criterion | The cheat it closes off |
| --- | --- |
| `closure` | A precondition, not the point. The declared exchange has to be carrying the closure, and it has to be **the water that actually moved** — see the finding below. Without it a model could declare a number no store ever paid for. |
| `exchange_response` | The probe's own question: a model that declares the prescribed head and whose exchange does not answer it, by at least a share of itself that any head-driven boundary clears, is absorbing its budget under a provenance it cannot back. Also bounds a declared loss by the water available to lose, in every variant. |
| `forcing_fidelity` | Stops the exchange closing a budget against a rescaled driver. |
| `state_bounds` | The cheater's original escape, still closed here. |
| `non_degenerate` | A declared exchange can manufacture runoff from nothing; the partition still has to be non-trivial. `non_degenerate` already counts `gwex` as supply in the runoff ratio. |

### The models that pass closure while doing no physics

`reference_noise_sink` is the exact bucket, conserving water internally, which
misreports its evaporation by a seeded ±30% and declares the difference as a
groundwater exchange. It declares that it consumes the prescribed head, and
never reads it. Its reported budget closes to **0.0000% of precipitation**; it
passes `closure`, `forcing_fidelity`, `state_bounds` and `non_degenerate`, and
across the rest of the suite it passes every budget probe — closure, the
counterfactuals, the invariances, the bounds — and fails only where its ±30%
evaporation noise shows as noise: `mass/steady-state`, where evaporation still
varies by 190% of its level under constant weather, and depending on the seed
`energy/pet-consistency`. The point stands with the exceptions named: no probe
that asks about the *budget* catches it. Because the harness
holds the model seed fixed across variants and nothing it computes reads `gwh`,
its exchange is identical with the head raised and lowered: **a response of
exactly 0 both ways.**

`reference_token_exchange` is the same model with one line added: a head-driven
trickle a million times weaker than the positive controls, folded into the soil
and into the declared exchange so the budget still closes exactly. It answers
the head with the right sign — +0.0018 mm raised, −0.0018 mm lowered — and is
caught because that is **1.5e-6 of its own gross exchange**, where a
head-driven boundary of any conductance answers with at least about 2e-3.

Both are `reference_cheater` moved one column over. The cheater solves for
*storage* as whatever balances the budget and is caught by `state_bounds`,
because invented storage cannot stay physical. Until now nothing asked whether
an invented *flux* answers the driver it claims.

## Calibration

Across **all five gate seeds** of `mass/exchange-response`, full ten-year
record. Only models that declare `gwh` are scored; the four physical models and
every submitted model in the repository decline it and are N/A.

| Model | Verdict | Δ raised (mm) | Δ lowered (mm) | gross (mm) | TV (mm) | resp / gross | **resp / TV** | required (mm) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `reference_driven_exchange` | **PASS** | **+3650** | **−3333 to −3447** | 1702 | 511 | 2.0 | **6.5–6.8** | 0.15 |
| `reference_evolving_exchange` | **PASS** | **+5** | **−5** | 1474 | 465 | 3.4e-3 | **1.1e-2** | 0.14 |
| `reference_recharge_exchange` | **PASS** | **+0.50** | **−0.34 to −0.44** | 1584–1767 | 52–66 | 2.0e-4 to 2.8e-4 | **5.4e-3 to 8.6e-3** | 0.016–0.020 |
| `reference_noise_sink` | **FAIL** | 0 | 0 | 1124–1218 | 1589–1754 | 0 | **0** | 0.48–0.53 |
| `reference_token_exchange` | **FAIL** | +0.0018 | −0.0018 | 1124–1218 | 1589–1754 | 1.5e-6 | **1.0e-6 to 1.1e-6** | 0.48–0.53 |

Three positive controls, because two were not enough. `reference_driven_exchange`
holds its internal head constant and answers with a sustained flux, exactly
`2 mm/m/day × 0.5 m × 3650 d = +3650 mm`; it would pass under any timing and
any floor, and on its own it hid both the spinup problem and the token escape.
`reference_evolving_exchange` lets its internal head move, `S dh/dt = C (H − h)`
with `S = 10 mm/m` and a five-day time constant; it is in equilibrium when
scoring begins and answers with `S·Δh = 5 mm` — the number a native MODFLOW run
gave — and it is the control that fails if the shift begins during spinup.
`reference_recharge_exchange` is the evolving boundary on a losing catchment,
`S = 1 mm/m`, half the runoff recharging the aquifer and draining back out
through the same boundary; its response is `S·Δh = 0.5 mm` against about
1700 mm of throughput, exact and monotone in the head, and it is the control
that pins the normalisation: a share of the gross fails it on every seed and a
share of the variation does not.

The floor of `3e-4` is set from the storativity sweep above: three times under
the value at `S = 0.1 mm/m`. Over a wider draw — the five gate seeds plus seeds
0–99, 105 paired runs per control, run independently at review — the margins
are a little narrower than the gate seeds alone suggest and are the figures to
quote: `reference_recharge_exchange` spans **3.95e-3 to 1.36e-2** (minimum 13×
above the floor), `reference_token_exchange` **1.02e-6 to 1.24e-6** (floor 242×
above it), `reference_evolving_exchange` is 1.076e-2 to four figures on every
seed since the weather does not enter it, and nothing lands on the wrong side
on any seed. The floor sits about four times above the geometric centre of that
gap, spending margin against the cheat to buy room for honest low-storativity
boundaries; for a conservation benchmark a false FAIL of a model that closes its
budget exactly is the more damaging error. The noise sink's variation is larger
than its gross because the accounting error changes sign almost every day.

### A finding from building it: the nominal term is not the applied term

An adapter that reports a model's exchange from its parameters rather than its
states can over-declare it. HBV 2.0's `LF` is a deterministic function of
drainage area and two parameters, but hydrodl2 applies it as
`SLZ = clamp(SLZ + LF, min=0)`: when the lower groundwater box is empty the
withdrawal is silently truncated, so the nominal value is not the water that
moved. Declaring the nominal `LF` over-declared the exchange by nearly a factor
of two — 34.1% of precipitation against an applied 18.1% — and `closure` failed
at **16.06%** until the applied value was recovered from the reported states and
fluxes. That is why `closure` is kept as a precondition rather than assumed:
`closure` requires the declared number to be the water that actually moved, and
`exchange_response` requires the water that actually moved to answer its
driver. Either alone is escapable.

## What this probe does not claim

**It scores no model in the current field.** No submitted model declares
`gwh` — the column did not exist — and none of the physical baselines has a
head-driven exchange to declare. All are N/A, INCOMPATIBLE, and their standings
are unchanged by this probe. That is the same footing as
`energy/radiation-consistency` on a model that does not consume `rlds`: a
criterion that did not judge a model does not add to its count. `dhbv2`, whose
`gwex` is a learned regional term with no external head behind it, is N/A here
and continues to fail `state_bounds` on the other budget probes for its learned
field capacity.

**It does not close the mixed-exchange escape, and it excludes one narrow
physical class.** See *What the floor excludes* above.

**It does not bound the size of an exchange.** `Σ|gwex|/Σpr` and `|Σgwex|/Σpr`
are reported and never gated. Gross bidirectional movement is not bounded by
precipitation, and how much external water a basin may legitimately depend on
is not a question a conservation probe can settle.

## What this adds to the contract, and how it relates to `mass/gw-sw-exchange-consistency`

**`gwh` is a new optional forcing column**, added the way `mass/human-abstraction`
added `abstr`, with its semantic stated in `AGENTS.md`; the `gwex` output
contract is unchanged. **`RunResult.model`** is a new optional field, set by
`Runner.run`, so a criterion can read what the model declared; every other
criterion ignores it. The criterion treats a `RunResult` with no manifest — the
repository's way of building one by hand in a test — as one to judge, not one to
excuse, and `tests/test_exchange_response.py` asserts that the gate fires in
that case.

**This probe and `mass/gw-sw-exchange-consistency` (#79) are complementary, not
overlapping.** They differ in the variable, the mechanism and the cheat:

| | `mass/gw-sw-exchange-consistency` | `mass/exchange-response` |
| --- | --- | --- |
| Variable | `gw_sw_exchange`, river–aquifer, between two stores *inside* the catchment | `gwex`, across the catchment *boundary* |
| Criteria | single-run: components sum to the net, recharge + exchange = Δ`gw`, signs | paired: the same weather with the prescribed head raised and lowered |
| Question | does the model's exchange bookkeeping agree with itself? | does the model's exchange answer the driver it declares? |
| Broken model | `reference_exchange_sign_error` — balance exact, one component's sign wrong | `reference_noise_sink` — bookkeeping exact, closure 0.0000%, exchange does not answer the head |

Its README says of its own checks that they are *"bookkeeping, not physics,
and hold exactly for any model that means it"*; this probe exists because
bookkeeping can be exact and mean nothing. `reference_noise_sink` would pass
every criterion of #79 if it reported those variables, since its accounts are
perfect; `reference_exchange_sign_error` is N/A here, since it reports no
`gwex`. Neither probe's negative control is caught by the other. The earlier
#41 proposal from which #79 grew included a stage counterfactual; the merged
probe does not perturb anything, so the mechanisms no longer meet.

The two also fit together on one model. #79 introduces `gw_boundary`, a
boundary term acting on the aquifer alone — a GHB, a well, a regional
exchange — and says that where it also crosses the catchment boundary it is
reported in `gwex` as well. A general-head boundary on a catchment aquifer is
exactly that: reported as `gw_boundary`, it is held to the aquifer balance by
#79; reported as `gwex` with `gwh` declared, it is held to its driver by this
probe. The same flux is checked from two sides, and no budget adds the two.

## Reproducing a failure

Every run records its seed. To re-run one:

```bash
ht run --model reference_token_exchange --probe mass/exchange-response --seed 1541944022
ht run --model reference_evolving_exchange --probe mass/exchange-response --seed 1541944022
ht gate --probe mass/exchange-response
```

## Proposal

[#76](https://github.com/Flood-Lab/HydroTuring/issues/76). The prescribed
driver, the conditional gate, the reversing positive control, the move from a
reversal bound to a paired counterfactual, the shift's timing, the second and
third positive controls, and the share floor's numerator and denominator all
follow from the reviews there and on the pull request, which supplied the
MODFLOW counter-examples, the token-term cheat and the losing-catchment case.

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
G       = Σ_t gwex_t · Δt          [mm]
G_gross = Σ_t |gwex_t| · Δt        [mm, control run]

G(raised)  − G(control)  >=  +max(s · G_gross, ε)
G(lowered) − G(control)  <=  −max(s · G_gross, ε)      s = 1e-3,  ε = 1e-8 · max(G_gross, 1 mm)
```

Raising the external head must bring more water in, and lowering it less, **by
at least a thousandth of the exchange the model itself declared.** `ε` is for
floating point and nothing else.

### Why a share, and not an absolute minimum

An earlier draft asked for no minimum at all, on the argument that any absolute
floor would fail a genuine boundary of small conductance. That is true of an
absolute floor and false of a share, and the distinction is the whole design.
For a general-head boundary `Q = C (H − h)`:

```
response = C · Δh · T                 gross = C · Σ|H − h| dt
response / gross = Δh · T / Σ|H − h| dt         — independent of C
```

Both halves scale with conductance, so their ratio does not. A boundary a
million times weaker than the positive control has a response a million times
smaller *and a gross exchange a million times smaller*, and its ratio is
unchanged: about **2** for this head series against a constant internal head.

Against an internal head that **moves**, the ratio is smaller but still bounded
below. With `S dh/dt = C (H − h)`, a shift is answered by the water that lifts
the internal head, `S · Δh`, and in the fast limit the gross exchange is
`S · TV(H)`, the storativity times the head's total variation. Both scale with
`S`, neither with `C`, and the ratio tends to `Δh / TV(H)` — about **2e-3** for
this head series. That is the physical minimum a head-driven exchange of *any*
conductance and *any* storativity can show, and the floor of `1e-3` sits under
it.

What falls under the floor is the model this probe exists to catch: a declared
exchange that is the day's accounting error with a token head term added. Its
response scales with the token's conductance; its gross does not, because the
gross is the error. `reference_token_exchange` measures **1.5e-6**, six orders
of magnitude under the floor. Without the share, that model passes: its answer
to the head has the right sign, and a signed test is all an earlier draft
asked for.

### What the share does not close

A model whose declared `gwex` mixes a genuine head-driven part with a large
unrelated one — a deep loss, a withdrawal, a `SIDE`-type baseflow share — has a
gross that the unrelated part inflates and a share that is honest but small.
From outside it is indistinguishable from the token cheat. The floor is placed
at a thousandth, a factor of two under the physical minimum, precisely to
tolerate a mixture of that order; a cheat that sizes its token term to a
thousandth of its sink passes, and is disclosed here as the price of not asking
the contract to split the head-driven part of `gwex` into its own variable.

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
it passes every one of the 28 probes currently in the suite. Because the harness
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

| Model | Verdict | Δ raised (mm) | Δ lowered (mm) | gross (mm) | response / gross | required (mm) | reversals, dry | sign changes, record |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `reference_driven_exchange` | **PASS** | **+3650** | **−3333 to −3447** | 1702 | **1.96–2.03** | 1.7 | 15 | 360 |
| `reference_evolving_exchange` | **PASS** | **+5** | **−5** | 1474 | **3.4e-3** | 1.5 | 3 | 365 |
| `reference_noise_sink` | **FAIL** | 0 | 0 | 1124–1218 | **0** | 1.1–1.2 | 69–75 | 1644–1741 |
| `reference_token_exchange` | **FAIL** | +0.0018 | −0.0018 | 1124–1218 | **1.5e-6 to 1.6e-6** | 1.1–1.2 | 69–75 | 1749–1831 |

Two positive controls, because one is not enough. `reference_driven_exchange`
holds its internal head constant and answers with a sustained flux, exactly
`2 mm/m/day × 0.5 m × 3650 d = +3650 mm`; it would pass under any timing and
with any floor up to a ratio of about 2, and on its own it hides both of the
things this revision fixed. `reference_evolving_exchange` lets its internal head
move, `S dh/dt = C (H − h)` with `S = 10 mm/m`, `C = 2 mm/m/day` and a time
constant of five days; it is in equilibrium when scoring begins and answers with
`S · Δh = 5 mm` — the number the MODFLOW run gave — over the first weeks of the
scored record, and with nearly nothing thereafter. It is the control that fails
if the shift is applied during spinup, and the one the floor is calibrated
against: its ratio of 3.4e-3 is the small response a real aquifer gives, 3.4×
above the floor and close to the fast-limit minimum of `Δh / TV(H) ≈ 2e-3`.

The floor of `1e-3` is therefore set from the physics of the head series, not
from either control: half the fast-limit minimum, so that no head-driven
boundary can fall under it, and three orders of magnitude above the token cheat.

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

**It does not close the mixed-exchange escape.** See *What the share does not
close* above.

**It does not bound the size of an exchange.** `Σ|gwex|/Σpr` and `|Σgwex|/Σpr`
are reported and never gated. Gross bidirectional movement is not bounded by
precipitation, and how much external water a basin may legitimately depend on
is not a question a conservation probe can settle.

## What this adds to the contract, and how it relates to #41

**`gwh` is a new optional forcing column**, added the way `mass/human-abstraction`
added `abstr`, with its semantic stated in `AGENTS.md`; the `gwex` output
contract is unchanged. **`RunResult.model`** is a new optional field, set by
`Runner.run`, so a criterion can read what the model declared; every other
criterion ignores it. The criterion treats a `RunResult` with no manifest — the
repository's way of building one by hand in a test — as one to judge, not one to
excuse, and `tests/test_exchange_response.py` asserts that the gate fires in
that case.

**This probe is close to #41 in mechanism.** Both perturb a hydraulic driver
and ask for a directionally correct exchange response. #41 perturbs
surface-water stage and is framed around gross GW→SW and SW→GW fluxes not yet in
`/io`; this probe perturbs the external head of the existing catchment-boundary
`gwex` channel and requires only that existing net flux. If the maintainers
would rather have one generic hydraulic-driver response criterion with two probes
on it, the implementation here factors that way.

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
reversal bound to a paired counterfactual, the shift's timing, the second
positive control and the share floor all follow from the reviews there and on
the pull request, which supplied the MODFLOW counter-examples and the token-term
cheat.

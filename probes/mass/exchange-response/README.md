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
`|sum R| <= 0.05 * sum pr`.

**A model that closes has `R ≈ 0` and therefore `gwex ≈ -R0` identically.**
Measured on a CREST implementation over a ten-year record, `max |R|` per step is
**2.3e-13 mm/day**. So comparing the declared flux with the residual it closes is
the same number twice, and neither the *magnitude* of `gwex` nor its timing on
its own can separate an honest exchange from an invented one. **A discriminator
therefore has to use information that is not algebraically fixed by closure.**
This probe uses a counterfactual: the same weather, the same model, the same
seed, and one change to the one driver the case prescribes.

### The assertion

The case supplies an external head, `gwh` (metres), as a visible forcing column,
and runs the record three times — with the head as given, **raised by 0.5 m**,
and **lowered by 0.5 m** — with every other column byte-identical. For a model
that **declares it consumes that column** in `needs_forcing` or `uses_forcing`:

> raising the external head must bring more water into the catchment over the
> record, and lowering it must bring less.

`gwex` is positive into the catchment. With the integrated exchange over the
scored record and the gross movement of the control,

```
G       = Σ_t gwex_t · Δt          [mm]
G_gross = Σ_t |gwex_t| · Δt        [mm, control run]
ε       = 1e-8 · max(G_gross, 1 mm)

G(raised)  − G(control)  >  +ε
G(lowered) − G(control)  <  −ε
```

That is the whole gate. It asserts nothing about how the exchange moves in time,
only that it answers the one quantity that was varied, in the direction physics
requires. **The tolerance is for floating point and nothing else.** No physical
minimum is asked for: a boundary with a very small conductance responds to a
0.5 m shift by a very small amount and is entirely physical, so the only thing
excluded is a response that is absent or wrongly signed. An earlier draft
required one percent of the gross exchange; that was not a law, and a model with
`C → 0` would have failed it while doing exactly what a general-head boundary
does.

**For a model that does not declare the head**, nothing is decided by the
response; its response is reported with the reason named. A model that declares
no exchange, or one under `1e-4` of the record's precipitation, returns
**`PASS`, reason "negligible exchange"** — never an unscored state. A declared
loss that exceeds the water available to lose fails, for every model and in
every variant, because that is a law and not a judgement.

### What is reported and deliberately not gated

Everything two earlier drafts of this probe gated on: the magnitude-weighted
reversal fraction of the exchange, its count of significant reversals inside
the rainless windows, the prescribed head's turning points there, and the
exchange's gross and net size against precipitation. All of it is in the
report; none of it decides anything. The reasons are in the next section,
because they are the reasons this probe has the shape it has.

### Why a counterfactual, and not a statement about reversal

The first draft gated on reversal alone inside rainless windows, on the argument
that with precipitation suppressed nothing drives a boundary exchange back and
forth. A native MODFLOW 6.7.0 run showed that to be false: a confined aquifer
behind a general-head boundary, no recharge, no ET, no pumping, water balance
closed to 7e-10 mm — and a boundary head of `10.5 + 0.5 sin(2πt/4)` m. Its
flux reversed every two days and would have failed. A regional boundary is
driven by conditions outside the basin, and local rainlessness constrains
nothing about them.

The second draft prescribed the head and gated on reversal *relative to it*: no
more reversals of the exchange than turning points of the head, on the claim
that between two turning points of the external head the difference
`h_ext − h_catchment` can cross zero at most once, *whatever the catchment's
head is*. That claim is false. It holds only if the catchment's head is
constant. If the internal head moves — `h_ext = t`, `h_catchment = t + 0.1
sin(10t)` is enough — the difference reverses as often as the internal head
does, and in any real aquifer the internal head does move: it is what MODFLOW
solves for. Reversal frequency is not a law of head-driven exchange unless the
internal head is constrained too, and the case does not prescribe the internal
head.

So the hard condition is now a controlled response. Paired runs, identical in
everything but the external head, and the question is whether the exchange
answered the change in the direction physics requires. That is a statement the
case can actually make, because `gwh` is the only quantity varied between the
paired runs — the model's internal states and heads respond and differ, and that
response is what is measured. It is the footing `mass/warming-response`,
`mass/human-abstraction` and `mass/precipitation-counterfactual` already stand
on.

### What declaring `gwh` means

`needs_forcing` and `uses_forcing` establish, mechanically, that a model *reads*
`gwh`. That alone would not license the gate: a model might read the column for
some other purpose, and its exchange might be driven by `gwh` together with its
own groundwater, river stage, soil water and routing, or not at all. So the
contract makes the declaration a **semantic opt-in**, stated in `AGENTS.md`:

> `gwh` is the prescribed external hydraulic head associated with the model's
> declared `gwex`; positive `gwex` is into the catchment. A model declaring
> `gwh` asserts that its external exchange responds monotonically to this
> hydraulic potential.

That is the semantic of a general-head boundary, `Q = C (H_ext − h)`, and it is
the only thing the gate uses. It asks for no functional form and no magnitude:
whatever else drives a head-driven exchange, raising the external head cannot
make less water come in. A model whose exchange does not respond to the head it
declares has claimed a provenance it cannot back, and that is the cheat. A model
whose exchange is not head-driven should not declare `gwh`, and is then not held
to it; its response is reported.

## The case

`generate.py` produces, from a seed and a variant, ten scored years of daily
weather behind 365 days of spinup, on the closure probe's climate so that a
model passing that probe is not asked to work in a new climate as well as to
account for its exchange: roughly **843 mm/yr precipitation** and **759 mm/yr
potential ET** (aridity ~0.9, the energy-limited side of the Budyko curve, with
a real seasonal snowpack). Gamma-distributed depths on a seasonal occurrence
probability; temperature an annual cycle plus an AR(1) anomaly; potential ET
temperature-driven with a daylength factor.

Three features are specific to this probe:

*   **The external head `gwh`**: a 10 m datum, a 0.15 m annual cycle and a
    0.35 m cycle of 20 days. The `raised` and `lowered` variants add and
    subtract 0.5 m, of the order of the head's own swing, so a head-driven
    exchange answers unmistakably; every other column is byte-identical across
    the three variants, and the criterion checks that it is.
*   **Two rainless stretches, 90 and 75 days, labelled `_regime = dry`**, in
    the last two scored years and in summer. They no longer decide anything;
    the reversal diagnostics are computed inside them so that a reviewer can
    see how an exchange behaves where direct precipitation forcing is absent.
    Columns beginning with an underscore are probe annotation: the harness
    strips them before the model sees the forcing.
*   **The whole record scored.** `min_window_days` is 3650: an integrated
    response is a statement about the record, and storage carries across it.

## Why these criteria

| Criterion | The cheat it closes off |
| --- | --- |
| `closure` | A precondition, not the point. The declared exchange has to be carrying the closure, and it has to be **the water that actually moved** — see the finding below. Without it a model could declare a number no store ever paid for. |
| `exchange_response` | The probe's own question: a model that declares the prescribed head and whose exchange does not answer it has claimed a provenance it cannot back. Also bounds a declared loss by the water available to lose, in every variant. |
| `forcing_fidelity` | Stops the exchange closing a budget against a rescaled driver. |
| `state_bounds` | The cheater's original escape, still closed here. |
| `non_degenerate` | A declared exchange can manufacture runoff from nothing; the partition still has to be non-trivial. `non_degenerate` already counts `gwex` as supply in the runoff ratio. |

### The model that passes closure while doing no physics

`reference_noise_sink` is the exact bucket, conserving water internally, which
misreports its evaporation by a seeded ±30% and declares the difference as a
groundwater exchange. **It declares that it consumes the prescribed head, and
never reads it.** Its reported budget closes to **0.0000% of precipitation**. It
passes `closure`, `forcing_fidelity`, `state_bounds` and `non_degenerate`, and
it passes every one of the 27 probes currently in the suite. It is caught here
because its exchange is **identical with the head raised and with it lowered**:
the harness gives every variant of a seed the same model seed, its evaporation
noise is drawn from that seed, and nothing it computes reads `gwh`. Its
integrated response is exactly 0 mm both ways against a required 10 to 12 mm.

This is `reference_cheater` moved one column over. The cheater solves for
*storage* as whatever balances the budget and is caught by `state_bounds`,
because invented storage cannot stay physical. Until now nothing asked whether an
invented *flux* answers the driver it claims.

## Calibration

Across **all five gate seeds**, full ten-year record. The gate applies only to
the two models that declare the head; everything else is reported.

| Model | Declares `gwh` | Verdict | Δ raised (mm) | Δ lowered (mm) | ε (mm) | reversals, dry | head turns, dry | `Σ\|g\|/Σpr` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `reference_bucket` | no | PASS (negligible) | — | — | — | — | — | 0 |
| `flex_lumped` | no | PASS (negligible) | — | — | — | — | — | 0 |
| `flex_topo` | no | PASS (negligible) | — | — | — | — | — | 0 |
| `sacsma_snow17` | no | PASS (reported only) | 0 | 0 | 2e-7 | 0 | 15 | 0.003–0.004 |
| pycrest (dCREST) | no | PASS (reported only) | 0 | 0 | 2e-6 | 0 | 15 | 0.027–0.028 |
| `reference_driven_exchange` | **yes** | **PASS** | **+3650** | **−3197 to −3385** | 1.7e-5 | 15 | 15 | 0.195–0.233 |
| `reference_noise_sink` | **yes** | **FAIL** | **0** | **0** | 1.0e-5 to 1.2e-5 | 63–78 | 15 | 0.134–0.143 |

`reference_driven_exchange` is the positive control. Its exchange is
proportional to the difference between the prescribed head and a fixed
catchment reference head, so it answers a +0.5 m shift with exactly
`2 mm/m/day × 0.5 m × 3650 d = +3650 mm` and a −0.5 m shift with the mirror,
less what the soil could not give on the days it was empty. It reverses 360
times over the record and 15 times inside the rainless windows — the same
number of times as the head turns, which is a property of *this* control's
constant reference head and not a bound on head-driven exchange in general. A
model with a moving internal head may reverse more often and must still pass,
and under this gate it does.

**There is no calibrated threshold to report, because there is no physical
threshold.** The tolerance is `1e-8 · max(G_gross, 1 mm)`, of the order of
1e-5 mm for either declaring model, and exists so that a response of exactly
zero is distinguished from floating-point noise. The honest control clears it
by eight orders of magnitude; the adversarial control produces exactly zero.
A physical boundary with a conductance a thousand times smaller than the
control's would respond by about 3 mm and pass by the same margin it has here;
one a million times smaller would respond by 3 µm and still pass. What the
gate excludes is not a small response but no response, or one of the wrong
sign.

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

**It does not gate a model that declines the driver.** A residual sink that does
not declare `gwh` passes with its response and its reversals reported. That is
deliberate: the case cannot know what drives such a model's boundary, and a gate
would be a claim about a driver it never supplied.

**It does not check that the response is proportionate, only that it is there
and rightly signed.** A model that reads the head, ignores it for its budget and
adds a token head-proportional term would pass — which is also why no physical
minimum is asked for: a minimum would not stop the token term and would fail a
genuine boundary of small conductance. Every counterfactual probe in the suite
shares this limit, and the benchmark's answer to it is the same here: each probe
closes one door, and this one closes the door marked "declares the head and
does not read it".

**It does not bound the size of an exchange.** `Σ|gwex|/Σpr` and `|Σgwex|/Σpr`
are reported and never gated. Gross bidirectional movement is not bounded by
precipitation — the same water may cross the boundary repeatedly — and how much
external water a basin may legitimately depend on is not a question a
conservation probe can settle.

**A one-signed exchange is not thereby honest, and a reversing one is not
thereby invented.** `reference_leaky`'s vice is hiding its sink, not having one:
that same 15% loss declared honestly, without the head, passes with its zero
response reported. And a physical aquifer behind a fast boundary reverses its
flux as often as its own head moves; the earlier drafts of this probe would have
failed it, and this one does not.

## What this adds to the contract, and how it relates to #41

**`gwh` is a new optional forcing column.** The `/io` forcing contract in
`AGENTS.md` lists `time, pr, tas, pet` and, for the withdrawal probe, `abstr`.
This probe adds `gwh` in the same way `mass/human-abstraction` added `abstr`: an
optional prescribed driver that only this probe supplies, with its semantic
stated in `AGENTS.md`. The `gwex` *output* contract is unchanged. So the probe
is not, strictly, a test of the existing contract with nothing added; it keeps
the existing output channel and adds one prescribed input that gives that
channel something to be falsified against. `RunResult.model` follows from the
same fact: a criterion that applies only to models that opted in needs to see
the manifest.

**This revision moves #76 closer to #41 in mechanism.** Both perturb a
hydraulic driver and ask for a directionally correct exchange response. The
distinction is that #41 perturbs surface-water stage and is framed around gross
GW→SW and SW→GW fluxes that are not yet in `/io`, whereas #76 perturbs the
external head of the existing catchment-boundary `gwex` channel and requires
only that existing net flux. If the maintainers would rather have one generic
hydraulic-driver response criterion with two probes on it, the implementation
here factors that way without difficulty.

## Reproducing a failure

Every run records its seed. To re-run one:

```bash
ht run --model reference_noise_sink --probe mass/exchange-response --seed 515232987
ht run --model reference_driven_exchange --probe mass/exchange-response --seed 515232987
ht gate --probe mass/exchange-response
```

## Proposal

[#76](https://github.com/Flood-Lab/HydroTuring/issues/76). The prescribed
driver, the conditional gate, the reversing positive control and the move from
a reversal bound to a paired counterfactual all follow from the review there,
which supplied the MODFLOW counter-example and the `h_catchment = t + 0.1 sin(10t)`
counter-example to the second draft's claimed theorem.

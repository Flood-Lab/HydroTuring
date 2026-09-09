# Probes we want

The suite is deliberately small right now. This is the list of probes we want
next, so you can **claim one instead of inventing one**.

To claim: open a [probe proposal](../../issues/new?template=probe_proposal.yml)
naming the id below. We label it `accepted`, assign it to you, and it is
yours — nobody else will build it while your name is on it. Then fork the
repository and:

```bash
ht init-probe --list-templates           # the shapes available
ht init-probe --template <kind>          # get a template
ht init-probe --from probe-draft.yaml    # scaffold the directory
```

Each entry below names the template to start from. A templated probe scaffolds
into something that already passes `ht gate` against a placeholder case, so
you can watch it separate the reference models before you write any physics.

Contributors of merged probes are authors on the benchmark paper, as are
contributors of five accepted model proposals. Models are solicited too, and
this page does not list them: propose any model you think the benchmark
should have a verdict on, at
[model proposal](../../issues/new?template=model_submission.yml). You do not
need to be able to package it. See [CONTRIBUTING.md](CONTRIBUTING.md#credit).

Difficulty is about the physics and the discriminating case, not the code.
Every probe here is a few hundred lines at most.

---

## Cross-budget consistency

A model can close its water budget and close its energy budget while being
incoherent between them. The first probe that notices is merged; the second
is the one we most want next.

### `energy/latent-heat-et-consistency` &middot; **merged**
Latent heat must equal evapotranspiration times the latent heat of the phase
change it underwent, at every step, with a temperature-dependent &lambda;.
Filed under `energy` because the schema admits mass, energy and momentum.
Contributed by Changming Li (SCUT).

### `coupled/snowmelt-energy-water` &middot; hard &middot; **unclaimed**

Melt in the water budget must equal the energy consumed by melting divided by
the latent heat of fusion.

*Why it discriminates.* Same failure mode at the phase change, where it is
most consequential for runoff timing.

---

## Mass

### `mass/catchment-closure` &middot; **merged**
The reference implementation. Read it before writing your own.

### `mass/resolution-invariance` &middot; **merged**
The same month at the minute, the hour and the day: do the volumes agree?

### `mass/warming-response` &middot; **merged**
Same rain, air 3 °C warmer and cooler: does runoff move the right way, both ways?

### `mass/causality` &middot; **merged**
One storm added: nothing may change before it, something must after.

### `mass/dry-down` &middot; **merged**
Two years without rain: runoff can only fall, and only stored water can drain.

### `mass/steady-state` &middot; **merged**
Three years of the same day: everything settles and the budget balances.

### `mass/extreme-rain` &middot; **merged**
The largest storm scaled to ten times: runoff cannot fall, nor exceed the rain added.

### `mass/runoff-bounds` &middot; **merged**
Over ten years, is the runoff possible at all? The mass question a runoff-only model has to answer.

### `mass/area-invariance` &middot; **merged**
The same weather on a ten times larger catchment: every depth identical.

### `mass/response-nonnegativity` &middot; **merged**
An added storm may never lower the flow, on any day.

### `mass/antecedent-monotonicity` &middot; **merged**
The same storm after a dry month and a wet one: more runoff from the wet one, and no more than the extra water.

### `mass/phase-counterfactual` &middot; **merged**
The same water as rain instead of snow: timing moves, volumes do not.

### `mass/snowpack-mass-closure` &middot; starter &middot; **unclaimed**
Snowfall minus melt minus sublimation minus the change in SWE.
*Discriminates:* models that quietly lose water at the rain-snow transition,
a very common bug that a whole-catchment budget can absorb.

### `mass/multi-decadal-drift` &middot; starter &middot; **unclaimed**
Fifty years with no trend in the forcing. Total storage must not drift
secularly.
*Discriminates:* a leak too small to trip a ten-year 5 percent threshold but
large enough to be unphysical over a climate-relevant record. Good first
probe: the criterion already exists, the case is the contribution.

### `mass/routing-network-closure` &middot; standard &middot; **unclaimed**
A branching network. Mass must close reach by reach, not only basin-wide.
*Discriminates:* models that conserve globally while moving water between
reaches non-physically.

### `mass/human-abstraction` &middot; standard &middot; **unclaimed**
Irrigation withdrawal and return flow, which must both appear in the budget.
*Discriminates:* models that treat abstraction as an unaccounted sink.

---

## Generalisation

Conservation that holds only where a model was fitted is not conservation, it
is a coincidence of the training distribution. These four probes ask whether
the property survives a move — to another place, to another time, to a
different question. Each has a template, so the harness work is done and what
is left is the case.

### `mass/extreme-event-closure` &middot; hard &middot; **unclaimed**
`ht init-probe --template extrapolation-time`

Ordinary years, then conditions outside anything earlier in the record. The
budget must close over the anomalous stretch on its own terms, scored
separately so nine ordinary years cannot dilute it.

*Discriminates:* models that learned closure as a statistical regularity of
their training distribution rather than as a structural property. This is the
probe most likely to separate architecturally-constrained models from ones
that merely look conservative in-sample.

### `mass/ungauged-basin-closure` &middot; standard &middot; **unclaimed**
`ht init-probe --template extrapolation-space`

Every seed draws a catchment, and some draws sit outside the range models are
normally fitted over. Every seed must pass, so the verdict turns on the
corners.

*Discriminates:* models fitted to a gauged sample and deployed on an ungauged
one, which is the deployment case the field actually cares about. The work is
in defending where the hull boundary sits; push one attribute out at a time,
or you generate catchments no real place resembles and fail honest models.

### `mass/precipitation-counterfactual` &middot; standard &middot; **unclaimed**
`ht init-probe --template counterfactual`

The same seed twice, once wetter. The added water must appear in the
difference between the reported budgets, split across evaporation, runoff and
storage.

*Discriminates:* closure by construction, structurally rather than
circumstantially. A model that solves for a budget term as the residual closes
perfectly on every seed forever, and today only its storage bounds catch it. A
counterfactual asks where the extra water went, which construction cannot
answer.

### `mass/time-origin-invariance` &middot; starter &middot; **unclaimed**
`ht init-probe --template invariance`

The same weather under different dates. Nothing may move.

*Discriminates:* date features and trend terms that survived from training. The
cheapest probe in the suite and the hardest to tune towards, because there is
no tolerance worth arguing about: the two runs agree to floating point or they
do not. A good first contribution.

---

## Energy

### `energy/pet-consistency` &middot; **merged**
Evaporation follows demand when the model's own soil is wettest and water when it is driest; no energy flux needed.

### `energy/surface-energy-closure` &middot; starter &middot; **unclaimed**
Net radiation minus sensible minus latent minus ground heat flux, minus the
change in stored energy.
*Note:* the denominator is accumulated |Rn|, which crosses zero every night,
so this probe **must** set an absolute floor as well as the 5 percent rule.
The canonical energy probe and the natural first one.

### `energy/snowpack-cold-content` &middot; hard &middot; **unclaimed**
The full snowpack energy budget including cold content and phase change. Melt
must not occur while the pack is below freezing.
*Discriminates:* models that melt snow on a warm day regardless of whether
the pack has the energy to melt.

### `energy/radiation-consistency` &middot; standard &middot; **unclaimed**
Outgoing longwave must be consistent with the reported surface temperature
through Stefan-Boltzmann, given emissivity.
*Discriminates:* models that predict surface temperature and radiation with
separate heads that never have to agree.

---

## Momentum

### `momentum/routing-conservation` &middot; **merged**
The channel store is never negative and never holds more than its hydrograph can.

### `momentum/channel-routing-mass` &middot; starter &middot; **unclaimed**
Inflow minus outflow minus the change in channel storage, per reach.

### `momentum/stage-discharge-monotonic` &middot; standard &middot; **unclaimed**
Steady-flow rating must be monotonic. Where a loop rating appears, it must be
traversed in the physically correct direction, with the rising limb carrying
more discharge at a given stage than the falling limb.
*Discriminates:* models that fit a hydrograph while implying an impossible
relationship between depth and flow.

### `momentum/wave-celerity-bounds` &middot; hard &middot; **unclaimed**
Kinematic wave celerity must be positive and near the Manning expectation for
the reach geometry.
*Discriminates:* models that route a flood wave upstream, or at a speed the
channel cannot support.

### `momentum/froude-regime` &middot; standard &middot; **unclaimed**
Flow in a mild-sloped reach must stay subcritical.
*Discriminates:* spurious supercritical flow, which usually signals the model
is not solving anything resembling momentum.

---

## Beyond conservation

Counterfactual response and invariance have moved up into
[Generalisation](#generalisation): the harness runs paired cases now, and both
have templates. What remains out of scope for suite 0.1, listed so nobody
builds it twice:

- **Real-data track.** Internal closure under observed forcing. Note this
  tests something different from the synthetic track: observed budgets do not
  close, so observations can never be the reference.
- **Spatial permutation.** Reorder the reaches of a network, or the years of a
  record, and require the long-run totals to be unchanged. Weaker than it
  looks, because storage carries across the boundary: only the totals are
  invariant, not the series. Worth doing once the routing probes exist.
- **Cross-model agreement.** Not a conservation test at all, and a different
  kind of claim. Noted here only so it is clear it is deliberately absent.

---

## Adding something not on this list

Welcome, and please propose it first. The proposal form asks one question
that matters more than the rest: *how would a model pass your probe while
understanding no physics?* If you can answer that, you have a probe. If you
cannot, you may have a diagnostic rather than a test.

## Models we want

All of them, and this list deliberately does not exist. There is no roadmap
for models because there is no shortlist: any published rainfall-runoff
model, any LSTM, any foundation model making a hydrologic claim is in scope,
and the useful judgement is which ones the field would learn something from.
Open a [model proposal](../../issues/new?template=model_submission.yml) and
say why. The form's *packaging status* field decides whether the work lands
with you or with the maintainer; it does not affect your credit either way.

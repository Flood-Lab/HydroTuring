# Probes we want

The suite is deliberately small right now. This is the list of probes we want
next, so you can **claim one instead of inventing one**.

To claim: open a [probe proposal](../../issues/new?template=probe_proposal.yml)
naming the id below. We will label it `accepted` and it is yours. Then:

```bash
ht init-probe                            # get a template
ht init-probe --from probe-draft.yaml    # scaffold the directory
```

Contributors of merged probes are authors on the benchmark paper. See
[CONTRIBUTING.md](CONTRIBUTING.md#credit).

Difficulty is about the physics and the discriminating case, not the code.
Every probe here is a few hundred lines at most.

---

## Cross-budget consistency

**These are the two we most want.** A model can close its water budget and
close its energy budget while being incoherent between them, and nothing in
the suite currently notices. Closing that gap is the sharpest thing anyone
could contribute.

### `coupled/latent-heat-et-consistency` &middot; hard &middot; **unclaimed**

Latent heat flux must equal evapotranspiration times the latent heat of
vaporisation, at every step, with a temperature-dependent &lambda;.

*Why it discriminates.* A model with separate water and energy heads can
satisfy both budgets independently and still report an LE that implies a
different ET than the one it reported. No single-budget probe can see this.
Needs a reference model that closes both budgets separately but incoherently.

### `coupled/snowmelt-energy-water` &middot; hard &middot; **unclaimed**

Melt in the water budget must equal the energy consumed by melting divided by
the latent heat of fusion.

*Why it discriminates.* Same failure mode at the phase change, where it is
most consequential for runoff timing.

---

## Mass

### `mass/catchment-closure` &middot; **merged**
The reference implementation. Read it before writing your own.

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

### `mass/extreme-event-closure` &middot; hard &middot; **unclaimed**
A record-breaking event well outside anything in the generated record so far.
*Discriminates:* models that learned closure as a statistical regularity of
their training distribution rather than as a structural property. This is the
probe most likely to separate architecturally-constrained models from ones
that merely look conservative in-sample.

### `mass/human-abstraction` &middot; standard &middot; **unclaimed**
Irrigation withdrawal and return flow, which must both appear in the budget.
*Discriminates:* models that treat abstraction as an unaccounted sink.

---

## Energy

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

Out of scope for suite 0.1, listed so nobody builds them twice. The model
contract already declares `supports.perturbation`, so these are unblocked
whenever we decide to open the category.

- **Counterfactual response.** Double the precipitation and check that every
  budget term responds, rather than one term absorbing the whole change. This
  is the clean structural answer to closure-by-construction.
- **Symmetry and invariance.** Unit rescaling, time shift, spatial
  permutation. A model that changes its answer when you express rainfall in
  metres has not learned the physics.
- **Real-data track.** Internal closure under observed forcing. Note this
  tests something different from the synthetic track: observed budgets do not
  close, so observations can never be the reference.

---

## Adding something not on this list

Welcome, and please propose it first. The proposal form asks one question
that matters more than the rest: *how would a model pass your probe while
understanding no physics?* If you can answer that, you have a probe. If you
cannot, you may have a diagnostic rather than a test.

# momentum/wave-celerity-bounds

This probe implements accepted proposal #148. It asks a question deliberately
left open by `momentum/routing-lag-consistency`: does transient propagation
speed respond to hydraulic state, rather than merely staying causal and getting
slower in longer reaches?

## Why this is a separate momentum probe

The neighbouring momentum checks constrain different projections of the same
hydrograph, but none identifies the transient channel-wave speed itself.

| Check | What it constrains | What remains invisible |
| --- | --- | --- |
| `momentum/routing-conservation` | routed water is not created or destroyed | timing and propagation speed |
| `momentum/routing-lag-consistency` | gross rainfall-to-runoff lag grows on a geomorphic travel-time scale | runoff generation, storage and channel travel are mixed into one lag; a fixed positive celerity can pass |
| `momentum/stage-discharge-monotonic` | the reported rating is monotone and any hysteresis has the physical sign | how fast a transient moves down the reach |
| `momentum/uniform-flow-friction-consistency` | steady discharge and stage satisfy the declared Manning balance | the characteristic speed of a transient about that steady state |
| `momentum/routing-causality` (accepted proposal) | routed flow does not anticipate its driver | any positive but physically wrong propagation speed |
| `momentum/froude-regime` (accepted proposal) | the reported hydraulic state stays in an admissible regime | state-dependent wave speed inside that regime |

This sits next to, but is narrower than, existing hydrologic metamorphic and
sensitivity testing. Reichert et al. (2024,
https://doi.org/10.5194/hess-28-2505-2024) perturb hydrologic drivers and test
whether model responses remain physically plausible, and Dawson et al. (2014,
https://doi.org/10.2166/HYDRO.2013.222) use partial-derivative sensitivities to
judge the physical legitimacy of neural hydrologic models. #148 does not claim
those ideas as new. It specializes them to a governing hydraulic derivative:
a controlled black-box perturbation is used to identify a transient
characteristic speed, which is then compared quantitatively with the local
derivative of the declared steady hydraulic relation.

The physical gap is also one order deeper than a new state constraint. The
neighbouring steady probes test values of the hydraulic state or an equilibrium
constitutive relation. #148 tests the **local dynamic derivative** of that
relation: the characteristic speed seen by a small perturbation. A model can
sit on a Manning-consistent steady state and still use an unrelated transient
router. In that sense the probe asks whether the model's transient Jacobian is
consistent with the physics that defines its equilibrium state, not merely
whether the equilibrium itself looks plausible.

The identification step is the contribution as much as the final inequality.
The probe enters at the routing control volume directly: `q_in` is a prescribed
river inflow, not rainfall that a land model first has to turn into runoff.
That removes runoff-generation and soil-storage physics from the question.

There is still a reason to use a counterfactual length pair rather than time a
single reach from the forcing row. A black-box adapter or router can have a
common input clock, source-cell residence time, sub-step convention or other
latency that is not channel propagation. Write

```
t_response(L, Q) = t_common(Q) + t_route(L, Q).
```

The short and long variants have byte-identical `q_in` and identical static
attributes except `reach_length_m`. Their difference therefore removes that
common term,

```
Delta t(Q) = t_response(L_long, Q) - t_response(L_short, Q),
```

so `Delta x / Delta t(Q)` identifies the incremental propagation speed of the
additional reach using only externally visible inflow and discharge. Repeating
the same identification at low, medium and high base flow then asks a question
no existing probe asks: whether that local transient propagation operator
changes with hydraulic state as the declared reach physics requires.

A fixed-celerity router is the concrete blind spot. It can conserve mass, be
strictly causal, delay a longer reach more than a shorter one, remain
subcritical, and even have a Manning-consistent steady rating while still using
the same transient speed at every discharge. #148 is designed to reject that
specific construction.

## Hydraulic theory behind the estimator

This probe is a black-box system-identification experiment for the
one-dimensional flood-wave problem, not another empirical lag rule.

With no lateral inflow, continuity is

```
dA/dt + dQ/dx = 0.
```

In the kinematic limit the reach has an equilibrium rating `Q = Q(A)`. For a
small perturbation `a` about a steady state `A0`,

```
A = A0 + a
Q(A) = Q(A0) + (dQ/dA)|A0 a + O(a^2),
```

so the linearized continuity equation is

```
da/dt + c(A0) da/dx = 0,
c(A0) = (dQ/dA)|A0.
```

Thus the Kleitz--Seddon speed `dQ/dA` is the characteristic speed of the
small flood-wave perturbation in this controlled regime; it is not the mean
water velocity and it is not the gross rainfall-to-runoff lag. Mishra & Singh
(2001, https://doi.org/10.1080/02626660109492830) connect the Seddon formula to
the linearized Saint-Venant solution. Ponce & Simons
(1977, https://doi.org/10.1061/JYCEAJ.0004892) show that the full shallow-water
problem contains distinct gravity, dynamic and kinematic wave bands, which is
why this probe deliberately stays in the small, friction-dominated kinematic
regime.

For a wide rectangular Manning reach, `Q proportional to A^(5/3)`, hence

```
c = dQ/dA = (5/3) Q/A = (5/3) u.
```

The production criterion uses the exact rectangular hydraulic radius rather
than this wide-channel approximation.

### Why the response centroid is the measured clock

Keeping the next-order pressure effect gives the linear convection--diffusion
(or diffusion-wave) form

```
dq'/dt + c dq'/dx = D d2q'/dx2.
```

River-routing literature treats `c` and hydraulic diffusivity `D` as the two
physical propagation parameters; the Hayami solution is the classical
constant-parameter analytical case (Moussa 1996,
https://doi.org/10.1002/(SICI)1099-1085(199609)10:9%3C1209::AID-HYP380%3E3.0.CO;2-2).
A modern statement of the same diffusion-wave model gives
`c = dQ/dA` and, in the usual low-inertia approximation,
`D approximately Q/(2 B S)`
(https://doi.org/10.1029/2023WR034692).

For a reach of length `L`, the Hayami impulse kernel is

```
g_L(t) =
  L / (2 sqrt(pi D t^3))
  exp(-(L - c t)^2 / (4 D t)),    t > 0.
```

Its first two temporal moments are

```
E[T_L]   = L / c,
Var[T_L] = 2 D L / c^3.
```

This gives the paired experiment a stronger interpretation than "subtract two
lags". If any common upstream/input response is `f(t)`, the outlet response is the
convolution `f * g_L`. Temporal cumulants add under convolution, so

```
centroid(out_L) = centroid(f) + L/c.
```

Therefore

```
centroid(out_long) - centroid(out_short)
    = (L_long - L_short) / c,
```

and the unknown common timing cancels. Diffusion can broaden and
attenuate the hydrograph without moving this first-moment identity. This is why
the probe measures a response centroid rather than a peak index. The second
cumulant also gives a natural future diagnostic,

```
Var(out_long) - Var(out_short)
    = 2 D (L_long - L_short) / c^3,
```

so the same paired design can in principle separate translation (`c`) from
dispersion (`D`) without reading model internals. #148 keeps the verdict on
celerity only.

Temporal moments and cumulants are established tools in flood-routing theory,
not a new law introduced by this probe. Åkesson et al. (2015,
https://doi.org/10.1002/2014WR016279) derived central temporal moments from a
kinematic-diffusive wave description to study stage-dependent hydraulic
response in stream networks, and Romanowicz & Doroszkiewicz (2019,
https://doi.org/10.26491/MHWM/95023) review the use of impulse-response
cumulants for linearized Saint-Venant routing. The contribution here is the
controlled **counterfactual pairing**: evaluate the same black-box routing
model under the same prescribed river inflow and hydraulic state at two reach
lengths, then difference those moments so model-specific common timing cancels
before the hydraulic derivative is tested across states.

This interpretation also has an empirical analogue. Allen et al. (2018,
https://doi.org/10.1029/2018GL077914) estimated river-wave celerity from paired
upstream/downstream gauges by dividing reach distance by an observed hydrograph
lag. Meyer et al. (2019, https://doi.org/10.1080/02626667.2018.1557336) found
that celerity--discharge relations vary strongly across real rivers and can
reverse after floodplain activation. That is exactly why #148 restricts itself
to an in-bank prismatic reach instead of asserting monotonic celerity for
arbitrary natural-river states.

## Experiment

For each seed, the generator draws one mild prismatic rectangular reach
geometry and prescribes three positive river inflows: 8, 16 and 32 m3/s.
`pr` is zero throughout; the experiment is routing-only. Only two model runs
are needed per seed: a 4 km `short` reach and a 20 km `long` reach with
byte-identical `q_in`. Inside each run, low, medium and high
hydraulic states occupy successive eight-day blocks. Each block gets four days
to settle, then a six-hour +5% pulse and more than three days of response tail.

Both runs are read at the reach outlet, so the paired propagation distance is

```
Delta x = L_long - L_short.
```

The accepted proposal sketched a centre gauge, which would have used half this
distance. The production contract exposes outlet `dis`, so the implementation
uses the outlet consistently in both variants and therefore the full length
difference. This changes only the observation geometry, not the tested
state-dependent-celerity hypothesis. The production case also uses the public
routing input `q_in` directly instead of manufacturing the same inflow from
rainfall, removing a land-model confound identified during the physical-model
audit.

The criterion removes the pre-event discharge, takes the centroid of the
signed transient response about the settled base, and estimates

```
c_obs = Delta x / (t_long - t_short).
```

For each hydraulic state, differencing the paired reaches cancels any timing
component common to both variants and isolates the added routing distance.
A 72-hour local response window ends
before the next state transition, so the following plateau cannot pull the
previous pulse centroid downstream. Before taking the first moment, the
criterion also requires the final pulse-width block to contain at most 1% of
the captured response magnitude (absolute signed excess). This integrated guard cannot be evaded
by a response that happens to cross the baseline on the final sample while
remaining material just beforehand. A slower tail is reported as "response not
contained" rather than being truncated into a biased celerity. The generated
record is hourly so timing quantisation is smaller than the short/long
travel-time difference.

## Scope and identification assumptions

The paired estimator is intentionally narrower than a generic flood-wave test.
Its interpretation is valid only when changing `reach_length_m` changes the
propagation distance without changing the hydraulic experiment itself. The
criterion therefore enforces, rather than merely documents, the two conditions
on which the cancellation rests:

- short and long runs carry byte-identical forcing and identical static
  attributes except `reach_length_m`;
- over the 12-hour pre-pulse plateau, each outlet discharge agrees with the
  prescribed `q_in` operating point within one percent, varies by at most
  one percent of that inflow, and the two variants agree with each other.

A model/domain in which changing the declared reach length also changes the
settled base flow has not answered the paired experiment. That is why such a
configuration is not rescued by widening the celerity tolerance. The contract
can verify that a model consumes `reach_length_m`, but it cannot prove where
inside that model the value is used. If reach length is allowed to alter runoff
generation itself, the paired estimate is a net length-dependent response and
must not be interpreted as a pure channel celerity; the intended applicability
is to models that use this field as reach-routing geometry.

The physical comparison is also deliberately a kinematic-limit statement: a
mild prismatic rectangular reach, positive slope and discharge, no adverse
downstream boundary or backwater, and a small perturbation about a settled
state. Passing does not claim that `c = dQ/dA` describes tidal, strongly
backwatered, rapidly varied or supercritical flow. Those regimes require
different probes and boundary information.

## Independent expectation

For a rectangular section,

```
Q = (1/n) A R_h^(2/3) sqrt(S)
R_h = A / P
P = b + 2h
A = b h
```

and differentiating this declared rating gives

```
c_kin = dQ/dA
      = Q [5/(3A) - 4/(3 b P)].
```

The hydraulic operating point is the externally prescribed `q_in`, not a
discharge chosen by the model under test. Before evaluating the theory, the
criterion requires the model's 12-hour pre-pulse outlet plateau to match that
`q_in` within one percent and to vary by no more than one percent of it. It
then inverts the monotone Manning relation at the prescribed inflow and
evaluates the derivative above. A model therefore cannot alter its steady
discharge and make the theory target move with the alteration, and no
propagation time is reused from the reference model.

The settled-flow check above is an identification precondition, with a one-percent
relative allowance. Once it holds, the three celerity requirements are:

1. `c_obs > 0`;
2. `abs(c_obs - c_kin) / c_kin <= 0.05`;
3. `c_low < c_medium < c_high`.

The five-percent physical allowance is fixed before gate evaluation. Proposal
#148's feasibility audit found the original coarse Rusanov reference biased
low, with the bias shrinking under spatial refinement. A second implementation
audit exposed a subtler paired-design problem: giving both the 4 km and 20 km
variants the same cell count makes the short reach five times finer in space,
so reach length also changes numerical diffusion and CFL cost. The production
transient reference therefore fixes the physical mesh scale instead, using an
approximately 80 m cell size (50 cells over 4 km and 250 over 20 km). The
counterfactual now changes reach length without changing numerical resolution;
the five-percent physical tolerance is not widened to cover discretization.

The production generator keeps the synthetic section deliberately wide
(60--90 m for 8--32 m3/s) and uses a 4 km / 20 km pair. Because the reported
discharge is read at the reach outlet, the paired propagation distance is the
full 16 km length difference. The resulting short/long centroid delay remains
resolved at the one-hour output step even at high flow.

## Baselines and applicability

| Model | #148 status | Reason |
| --- | --- | --- |
| `reference_bucket` | N/A | no geometry-aware reach-wave process |
| `flex_lumped` | N/A | existing lag uses catchment channel-length fields and fixed 1 m/s celerity, not this reach contract |
| `flex_topo` | N/A | does not consume this reach geometry |
| `sacsma_snow17` | N/A | does not consume this reach geometry |
| `reference_saint_venant` | gate must-pass | independent finite-volume dynamic-wave solve |
| `reference_fixed_celerity` | gate must-fail | deliberately fixed 1 m/s propagation |
| `wflow_sbm` | submitted-model PASS | Deltares Wflow.jl native river kinematic wave, run through an independent Docker adapter path |

`reference_saint_venant` is the CI must-pass model. Its daily steady path is
unchanged; for sub-daily forcing it advances the same continuity and momentum
equations continuously through each output interval on the common approximately
80 m mesh scale. The criterion does not call the solver's fluxes or internal
wave speeds.

`reference_fixed_celerity` is the deliberately broken control. It transports
the prescribed discharge downstream at exactly 1 m/s for every hydraulic
state. It remains positive and length-dependent, but cannot satisfy the
low/medium/high state response.

The four general physical baselines are intentionally not forced through this
probe. `reference_bucket`, `flex_topo` and `sacsma_snow17` have no
reach-wave process, while `flex_lumped`'s existing geometry path is the
catchment-scale triangular lag used by the separate routing-lag probe. Extending
one of those models would invent the process rather than expose an existing
one. This is the process-missing exception documented in
`docs/writing-a-probe.md`: the exact reference remains the fast acceptance
gate, while a submitted physical model with the process supplies the
independent archived evidence.

### Independent submitted-model evidence: Wflow.jl

The first Wflow audit used its existing one-cell land/river/outlet
schematisation and rainfall-derived runoff. It was deliberately rejected rather
than tuned to pass: changing river geometry also changed the settled base
discharge, so the short/long pair no longer isolated propagation. That failure
identified a problem in the experiment interface, not a reason to widen the
five-percent celerity allowance.

The final evidence path uses interfaces that already exist on both sides of
the benchmark contract:

- HydroTuring's public `q_in` is a prescribed routing inflow.
- Wflow v1.0.4 exposes
  `river_water__external_inflow_volume_flow_rate` as a forcing of its native
  river kinematic wave.
- The adapter builds a two-cell river chain only when `q_in` is present.
  The first cell is a fixed 100 m source cell and receives `q_in`; the second
  cell is the tested 4 km or 20 km reach. The source cell is byte-identical
  between variants, so its residence time is a common term removed by the
  paired centroid difference.
- `width_m`, `slope`, `manning_n` and `reach_length_m` are written to
  Wflow's native river maps. The scored `dis` is Wflow's own downstream
  `q_av`, not total catchment runoff and not a discharge reconstructed by
  the criterion.

No Wflow parameter is fitted to this probe and the criterion is unchanged.
Wflow's kinematic wave uses its own `A = alpha Q^(3/5)` relation and a fixed
wetted-perimeter approximation rather than the criterion's exact rectangular
normal-depth inversion. Across the three deterministic gate geometries and
all three hydraulic states, the celerity implied by that native relation differs
from the criterion's exact rectangular `dQ/dA` by only -0.08% to +0.56%;
the production allowance remains 5%.

A clean GitHub Actions Docker run of adapter `1.0.4-ht.5` on all three gate
seeds passed #148 without changing the threshold:

```
low = 0.704 m/s < medium = 0.927 m/s < high = 1.225 m/s
```

This is the independent submitted-physical-model evidence for the probe.
The exact Saint-Venant reference and Wflow do not share the same numerical
scheme or routing implementation; they meet only at the externally visible
physics asserted by the criterion.

LISFLOOD was also explored as an independent candidate, but that audit used
adapter wiring that is not present in this repository and has not been
reproduced through the committed HydroTuring interface. Its numerical celerity
values are therefore not used as evidence for or against the probe. A future
LISFLOOD comparison should first expose the declared reach length through a
reproducible adapter path before interpreting any routing discrepancy as model
physics.

## Reproduction

```sh
ht validate
ht gate --probe momentum/wave-celerity-bounds
pytest -q
```

No generated forcing or model output is committed.

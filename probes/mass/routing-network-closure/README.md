Mass/routing-network-closure

Reach-by-Reach Mass Closure in a Branching Network

Status: implementation draft for mass/routing-network-closure (#108).

The proposal was accepted for PR, but the current HydroTuring /io contract
explicitly has no network dimension. This draft therefore keeps the forcing
side as close as possible to the existing contract and makes the network
requirement explicit rather than hiding it.

What it asserts

Two headwater sub-catchments feed reaches A and B. Those reaches join at a
zero-storage junction J and feed reach C.

headwater A -> Qin,A -> reach A -> Qout,A --\
                                             J -> Qin,C -> reach C -> Qout,C
headwater B -> Qin,B -> reach B -> Qout,B --/

For each reach i and daily interval t,

r_i,t = (Q_in,i,t - Q_out,i,t) * dt - (V_i,t+1 - V_i,t)

and for the zero-storage junction,

r_J,t = (Q_out,A,t + Q_out,B,t - Q_in,C,t) * dt

where Q_in and Q_out are interval-mean discharge in m3/s, V is
reach channel storage in m3, and dt = 86400 s.

The proposed reach score is

reach_error_i = sum_t(abs(r_i,t)) /
                sum_t(Q_in,i,t * dt)

and the proposed junction score is

junction_error = sum_t(abs(r_J,t)) /
                 sum_t((Q_out,A,t + Q_out,B,t) * dt)

Each reach and the junction must remain below 5%. Summing absolute interval
residuals prevents positive and negative errors at different times from
cancelling.

A denominator smaller than the equivalent of 1 mm of water over the
contributing area is treated as degenerate rather than as a pass. For
example, for a 60 km2 headwater reach the floor is 60,000 m3.

Case generation

generate.py creates:

365 daily spin-up intervals;

90 scored daily intervals;

one synthetic pr/tas/pet forcing series, matching the current forcing
contract rather than introducing pr_A/pr_B;

a fixed Y-shaped topology in static.json;

different contributing areas for headwaters A and B;

no additional lateral channel inflow, groundwater-channel exchange,
withdrawal, or other unreported routing source/sink.

During spin-up, intermittent precipitation events establish catchment and
channel states. One scored precipitation event occurs within the first
20 days, followed by at least 70 dry days. Temperature stays above freezing,
so the event is rain rather than snow, and PET stays moderate.

The same meteorological forcing is available to both headwater
sub-catchments. Their different contributing areas, and any model-native
differences in sub-catchment response, can produce distinct volumetric
headwater hydrographs. This avoids requiring two independent forcing-series
columns in the first implementation.

Groundwater and other catchment stores may contribute to headwater runoff
before it enters the routing network. Once water is in a routing reach,
however, no unreported groundwater-channel exchange or other extra source or
sink is allowed in this first case.

Contract feasibility

This is the part that does not exist in the current scalar /io contract.

The first implementation needs a minimal network representation, not a
general-purpose river-network API. For the fixed Y case it needs only:

static topology identifying reaches A, B, C and junction J;

contributing area for each reach, used for the denominator floor; and

model-reported reach-indexed inflow, outflow, and endpoint channel storage.

A concrete first schema could use columns such as:

q_in_A, q_out_A, channel_A
q_in_B, q_out_B, channel_B
q_in_C, q_out_C, channel_C

The exact upstream naming/representation should be settled in PR review.
A model that cannot expose these reach-scale quantities should be N/A for
this probe, not FAIL.

Gate plan

The gate needs more than an exact reference written alongside the criterion.

Proposed must_pass models:

reference_network_exact: an exact conservative three-reach reference,
useful for diagnosis;

flex_topo_network: a proposed extension of the existing physical
flex_topo baseline with a minimal conservative three-reach routing wrapper.
This provides an independent physical model that genuinely emits the
reach-indexed quantities required by the probe.

Proposed broken variants:

Broken model

Expected catch

reference_network_compensation

A loses water while B creates the same amount; whole-network closure is preserved but reach closure fails

reference_network_junction_loss

A and B close individually, but Qout,A + Qout,B != Qin,C; junction closure fails

reference_network_temporal_shift

Cumulative volume is preserved but water is shifted between intervals; sum(abs(r_t)) fails

The exact baseline names can be adjusted during PR review; the important point
is that at least one physical baseline must actually consume the network
representation and pass independently of the exact synthetic reference.

Scope

The probe tests local mass accounting and network connectivity. It does
not prescribe a routing equation, test wave speed, or judge hydrograph
accuracy beyond what is required to form the local conservation budgets.

# mass/routing-network-closure

Reach-by-reach mass closure in a branching network.

## What it asserts

Two headwater sub-catchments feed reaches A and B. These reaches join at a
zero-storage junction J and feed downstream reach C.

```text
headwater A -> Qin,A -> reach A -> Qout,A --\
                                             J -> Qin,C -> reach C -> Qout,C
headwater B -> Qin,B -> reach B -> Qout,B --/

For each reach i and daily interval t,
r_i,t = (Q_in,i,t - Q_out,i,t) * dt - (V_i,t+1 - V_i,t)

and for the zero-storage junction,
r_J,t = (Q_out,A,t + Q_out,B,t - Q_in,C,t) * dt

where Q_in and Q_out are interval-mean discharge in m3/s, V is
reach channel storage in m3, and dt = 86400 s.

The proposed scores are

reach_error_i = sum_t(abs(r_i,t)) / sum_t(Q_in,i,t * dt)

junction_error = sum_t(abs(r_J,t)) / sum_t((Q_out,A,t + Q_out,B,t) * dt)

Each reach and the junction must remain below 5%. Summing absolute interval
residuals prevents errors at different times from cancelling.

A denominator smaller than the equivalent of 1 mm of water over the
contributing area is treated as degenerate rather than as a pass.

## The case

generate.py creates:
365 daily spin-up intervals;
90 scored daily intervals;
one synthetic pr/tas/pet forcing series;
a fixed Y-shaped three-reach topology;
different contributing areas for headwaters A and B;
no additional lateral channel inflow, groundwater-channel exchange,
withdrawal, or other unreported routing source or sink.

During spin-up, intermittent precipitation events establish catchment and
channel states. During scoring, one precipitation event occurs within the
first 20 days, followed by at least 70 dry days. Temperature stays above
freezing and PET remains moderate.

The same meteorological forcing is supplied to both headwater
sub-catchments. Differences in contributing area and model-native
sub-catchment response can produce distinct headwater hydrographs.

Groundwater and other catchment stores may contribute to headwater runoff
before it enters the routing network, but no unreported exchange is allowed
within the routing reaches.

## Contract feasibility

The current scalar /io contract has no network dimension.

The first implementation therefore needs a minimal network representation for
the fixed Y-shaped case:

topology identifying reaches A, B, C and junction J;
contributing area for each reach;
model-reported reach inflow, outflow, and endpoint channel storage.

A possible first output schema is:

q_in_A, q_out_A, channel_A
q_in_B, q_out_B, channel_B
q_in_C, q_out_C, channel_C

The exact naming and representation can be settled during PR review. Models
that cannot expose these reach-scale quantities should be N/A for this probe,
not FAIL.

## Gate plan

Proposed must_pass models:

reference_network_exact: exact conservative three-reach reference;
flex_topo_network: proposed conservative three-reach extension of the
existing flex_topo physical baseline.

Proposed broken variants:

Model	Expected failure
reference_network_compensation	Local reach loss/gain while whole-network closure is preserved
reference_network_junction_loss	Qout,A + Qout,B != Qin,C
reference_network_temporal_shift	Cumulative volume is preserved but interval closure fails

At least one physical baseline should consume the network representation and
pass independently of the exact synthetic reference.

## Scope
[README.md](https://github.com/user-attachments/files/32397143/README.md)

This probe tests local mass accounting and network connectivity. It does not
prescribe a routing equation, test wave speed, or evaluate hydrograph accuracy
beyond what is required for local conservation.



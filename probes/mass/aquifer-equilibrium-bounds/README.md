# mass/aquifer-equilibrium-bounds

A water ledger can close while its hydraulic boundary or storage coefficient
is wrong. This probe asks a passive, fully confined aquifer to relax toward a
prescribed river head and checks that its storage never leaves the envelope
implied by that head and the supplied storativity.

## Physical statement

For positive transmissivity, uniform initial head h0, zero sources, a single
constant-head river connection of positive conductance, and impermeable other
boundaries, the diffusion maximum principle gives

    min(h0, H) <= h(x, t) <= max(h0, H).

Uniform confined storativity sigma converts head to storage change over the
same aquifer plan area:

    G(t) - G(0) = 1000 sigma (mean(h(t)) - h0)  [mm].

The control stays at H = h0. Raised and lowered variants change H only at the
first scored interval; common spinup and the model random seed are identical.
Every interval-end paired storage difference X = G_variant - G_control must
satisfy min(0, B) <= X <= max(0, B), where B = 1000 sigma (H - h0).
The maximum positive distance outside either bound is scored; a later recovery
cannot cancel an earlier excursion.

This necessary bound does not prescribe a recession shape, a response time,
or an exact spatial head field. A frozen aquifer can satisfy it. Such output
is explicitly reported as "no hydraulic response exercised", rather than
failed using an invented minimum hydraulic response.

## Case and control volume

Five fixed gate seeds draw area 0.5–2 km2, storativity 0.002–0.01 and head
changes 0.5–1.5 m. The daily record has two common equilibrium spinup days and
90 scored days, all with zero recharge. The minimum evaluation window is the
whole 90-day scored experiment. A positive conductance is drawn through the
lumped timescale A sigma / C = 2–10 days; a distributed solver need not
equilibrate on that lumped timescale to pass.

The initial head is 5 m above the aquifer top. The generator explicitly checks
that the entire head envelope remains above that top and that the highest
possible riverbed remains at least 1 m below the lowest admissible head.
This matters because the existing MODFLOW adapter sets riverbed elevation
to stage minus river_bottom_offset_m; the riverbed moves with the prescribed
stage. The offset is 4 m, so even the largest 1.5 m perturbation keeps the
connected RIV branch active. On that branch C(H - h) is independent of the
riverbed. Native cell heads must also verify the assumptions during physical
validation.

The aquifer is confined with constant density, geometry and storage properties.
There is no recharge, pumping, ET, compaction or other active aquifer boundary.
A declared zero gw_boundary is allowed; nonzero gross boundary activity makes
the case N/A (INCOMPATIBLE), even if its signed integral is zero. Missing
declared or non-finite boundary data are never silently replaced with zero.
Models lacking the required outputs are N/A (INCOMPLETE).

Storage is aquifer-area-equivalent mm. The MODFLOW adapter's specific-yield
initial inventory is a common additive offset; changes come from native STO
budget terms. The judge differences the reported inventory and never uses that
offset as a confined head estimate. The supplied confined storage coefficient,
not specific yield, sets B. River-aquifer exchange is gw_sw_exchange, an
internal transfer, and is never added to catchment-boundary gwex.

## Criterion and preliminary allowances

aquifer_equilibrium_bounds checks both head variants at every interval.
It also checks common-spinup inventory equality, neutral storage stationarity
and neutral gross exchange. All groundwater balance, component and state
checks run on every variant.

The initial bound allowance is max(0.05 |B|, 0.001 mm). The denominator is the
prescribed storage excursion, never the model's background storage or gross
exchange. The generator keeps |B| >= 1 mm. Initial and neutral guards initially
allow 0.001 mm; the additional-boundary activity floor is 1e-6 mm.
These numerical choices are pending the independent validation described below.
Derived arithmetic is checked for non-finite values as well as raw output.

## References and discrimination

The exact confined reservoir and its backward-Euler counterpart supply
Docker-free must-pass controls. An amplified-boundary control uses
H_eff = h0 + 1.4 (H - h0); an inflated-storativity control uses 1.4 sigma.
Both report their own conserved water ledger and correctly signed exchange,
yet can exceed the prescribed envelope. The gate must catch each only on
aquifer_equilibrium_bounds.

The standard four catchment references do not have river-aquifer exchange;
adding it would invent a process. This uses the process-specific exception
in docs/writing-a-probe.md: a domain reference in the gate and separately
archived evidence from the unchanged MODFLOW 6.7.0 adapter. Independent native
validation and audit-seed measurements will be recorded before submission.

## Relation to other probes

Accepted issue #159 is distinct from recharge ordering (#135): more recharge
can preserve storage/exchange ordering even when both paired runs relax
toward the same wrong equilibrium. In the amplified-head control, the
incremental recharge response still obeys d(delta G)/dt = added recharge
- k delta G with k > 0. The new head-derived bound tests an independently
specified hydraulic capacity.

Groundwater budget and component consistency test the ledger. Datum
invariance tests a coordinate transformation. A directional stage-response
test can accept a correctly directed response of excessive magnitude.
Those properties do not establish this storage envelope. Passing here does
not certify those other properties either.

## Sources

- Accepted proposal: https://github.com/Flood-Lab/HydroTuring/issues/159.
- MODFLOW 6.7.0 STO definitions:
  https://modflow6.readthedocs.io/en/6.7.0/_mf6io/gwf-sto.html.
- MODFLOW 6.7.0 RIV definitions:
  https://modflow6.readthedocs.io/en/6.7.0/_mf6io/gwf-riv.html.

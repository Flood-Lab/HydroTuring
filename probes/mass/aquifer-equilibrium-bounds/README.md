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

## Criterion and numerical allowances

aquifer_equilibrium_bounds checks both head variants at every interval.
It also checks common-spinup inventory equality, neutral storage stationarity
and neutral gross exchange. All groundwater balance, component and state
checks run on every variant.

The bound allowance is max(0.05 |B|, 0.001 mm). The denominator is the
prescribed storage excursion, never the model's background storage or gross
exchange. The generator keeps |B| >= 1 mm. Initial and neutral guards
allow 0.001 mm; the additional-boundary activity floor is 1e-6 mm.
These are provisional benchmark detection thresholds, not universal solver
accuracy estimates. They were retained through validation, without fitting
them to a model's response. The largest native discrepancy against the
requested storage law was 3.26e-6 mm, well below the absolute floor; the
40% fault controls exceed the bound allowance by at least a factor of 7.95.
The 5% band deliberately permits smaller constitutive errors and requires
maintainer judgment as a benchmark policy.
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
archived evidence from the unchanged MODFLOW 6.7.0 adapter.

## Measured validation

[The machine-readable audit](../../../results/aquifer-equilibrium-bounds/validation.json)
records all seeds, normalized criterion results, native margins and storage
traces. MODFLOW 6.7.0 and FloPy 3.11.0 ran in the adapter's pinned amd64 image;
its image digest is recorded in that report. The ordinary contract check and
all five gate-seed evaluations passed.

The native audit covers all three variants for five gate seeds plus independent
seeds 7 and 42 (21 runs). Six additional runs use the first gate seed with
IMS head tolerances tightened from 1e-8 to 1e-10 m or a 20 x 20 spatial grid
instead of 10 x 10. All retain the adapter's daily forcing and reporting.

| Native check, over all 27 runs | Measured result |
| --- | --- |
| Minimum cell head above aquifer top | 3.68186 m |
| Minimum cell head above the moving riverbed | 2.64142 m |
| Gross STO-SY release | 0 mm |
| Maximum head-envelope excess against serialized IC/RIV inputs | 2.85e-14 m |
| Maximum storage/head mapping error against serialized STO coefficient | 3.10e-10 mm |
| Maximum storage/head mapping error against requested inputs | 3.26e-6 mm |
| Maximum neutral storage change after common spinup | 1.62e-6 mm |
| Diagnostic versus ordinary result.csv, standard profile | 21 of 21 byte-identical |

FloPy serializes uniform IC/STO constants to eight decimal places and RIV
heads to eight-place scientific notation. The observer saves those native
inputs and verifies their expected rounding, rather than assuming lossless
serialization. The largest river-head rounding was 4.79e-7 m. Native heads
and the serialized confined storage coefficient independently reproduce the
reported storage change; this does not infer head or a missing flux from the
adapter's water-budget residual. Audit precision checks of 1e-7 m for native
head bounds and 1e-6 mm for the native storage mapping are separate from the
probe's provisional detection band. The native data satisfy much tighter
errors, as the table shows.

Both references pass every declared criterion on seeds 0–99 and the five
gate seeds (105 distinct seeds). Each negative control fails only
aquifer_equilibrium_bounds on all 105. Their smallest residual/allowance
ratios are 7.99677 for amplified head and 7.95696 for inflated storage,
against a failure threshold of 1. Tests also show both conservative faults
passing the existing recharge-ordering criterion on a common confined case
with an added recharge pulse. That is evidence of distinct logical content,
not a claim that either fault passes every existing probe.

Reproduce with a built MODFLOW image and optional host-side FloPy 3.11.0:

```bash
ht verify-adapter --model modflow6 --probe mass/aquifer-equilibrium-bounds
ht run --model modflow6 --probe mass/aquifer-equilibrium-bounds --gate-seeds
ht gate --probe mass/aquifer-equilibrium-bounds
python scripts/validate_aquifer_equilibrium.py --workdir "$HOME/aquifer-equilibrium-audit"
```

Use a fresh user-owned work directory; Docker Desktop needs a shared path.
The observer copies native files before the unchanged adapter removes its
temporary simulation, then extracts them outside the benchmark output
directory before normal contract validation. The raw files remain local;
the report is committed separately from the generated forcing.

The grid sensitivity confirms the bound and applicability assumptions, not
grid convergence or identical transient hydrographs. Only seven seeds were
audited with the distributed solver. This is an aggregate necessary bound:
local head errors can cancel in aggregate storage, and frozen responses can
pass. No minimum response or finite-time equilibrium claim is made.

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

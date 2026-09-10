# energy/surface-energy-closure

## What it asserts

For interval-mean fluxes over the same hourly interval, the surface residual is

```
r_i = Rn_i - H_i - LE_i - G_i                         [W m-2]
```

`Rn` (`rn`) is positive into the surface; sensible heat `H` (`hfss`) and
latent heat `LE` (`hfls`) are positive away from it; ground heat `G` (`hfg`)
is positive from the surface into the ground. Negative H and LE are allowed,
including at night. No rule requires nighttime evaporation to vanish.

For each **contiguous** solar daytime or nighttime block b, define

```
T_b = sum_i(dt_i)
E_b = sum_i(abs(r_i)  * dt_i) / T_b                  [W m-2]
D_b = sum_i(abs(Rn_i) * dt_i) / T_b                  [W m-2]
A_b = max(0.05 * D_b, 2 W m-2)
```

Every block must satisfy `E_b <= A_b`. At this probe's fixed hourly step the
duration-weighted means equal ordinary means. Absolute residuals are taken
**before** summing: opposite errors cannot cancel within a block, and a good
day cannot compensate for a bad night. Repeated labels denote separate
blocks, not one average over all daytime rows and another over all nights.
The criterion reports block duration, mean absolute residual, allowance,
failed block count and the worst block.

The 5% relative allowance and 2 W m-2 absolute floor carry forward the
existing energy-closure engineering conventions. They are starting values,
not observationally calibrated tolerances. The denominator uses absolute
net radiation; neither signed cancellation nor division by near-zero Rn is
involved. The floor provides a finite allowance under weak radiation.

### The budget boundary

The case is a homogeneous, snow-free bare surface with a **zero-capacity
skin**. It has no canopy, snow/ice melting or freezing, unreported phase-change
storage, or lateral energy transport. Evaporation or condensation is carried
by LE. G is at the **actual soil surface**. Subsurface heat storage lies
below that boundary and must not also be subtracted from this skin budget.

A model with G at a deeper soil boundary, or with finite skin/canopy heat
storage, needs a matching boundary or explicit storage accounting before
this residual can judge it. Output names and `PT1H` support alone cannot
establish this. The harness checks those declarations; adapter documentation
and review must establish the physical boundary. A failure should only be
interpreted as a conservation violation for a compatible boundary.

## The case

`generate.py` returns 384 hourly rows: two spinup days followed by 14 scored
days. The record starts at local-solar 06:00, leaving exactly 28 complete
12-hour scored blocks. A timestamp labels the start of `[time, time + 1 h)`.
Daytime is 06:00-18:00 and nighttime is 18:00-06:00 in this idealised solar
clock. `_regime` follows that clock, not the sign of Rn. The harness retains
annotations for scoring and strips them before sending forcing to a model.

Absorbed shortwave follows a 12-hour positive sine, integrated analytically
over each hour. Its daily amplitude is 420-540 W m-2 and hourly seeded cloud
attenuation is 0.45-1.0. Subtracting 25-40 W m-2 of daily longwave cooling
gives the imposed net radiation. Cloud attenuation and cooling are constant
within each interval, so Rn is an interval mean. Nighttime cooling and the
daylight shoulders include weak net radiation. Rn is a supplied boundary
flux, not a prediction based on surface radiative temperature.

Air stays warm (16-28 degC), canopy capacity is zero, and the water inputs
give the existing accounting reference modest evaporative demand and
occasional rainfall. Water fluxes remain rates in mm/day at every hour.
The reference starts with half of its 180 mm soil capacity. These inputs
support a warm bare-surface accounting case; no water or temperature states
are required from the model being tested.

`min_window_days: 14` prevents the default seven-day hourly evaluation
window from truncating the scored cycles. The generator uses the actual
`generate(seed)` contract; its row count is derived from the period, spinup
and hourly resolution, and the harness checks it against `probe.yaml`.

## Why this criterion

The existing `energy_closure` is unchanged. Its signed cumulative residual
can be zero while hourly residuals are large. Starting from an exact
accounting solution, `reference_diurnal_bias` changes H alone so the residual
is +20 W m-2 during the day and -20 W m-2 at night. The complete record
closes cumulatively, but every block fails `energy_closure_by_phase`. The
negative reference reads the current solar timestamp, never hidden phase
annotations, the evaluation horizon or future rows.

Only the new criterion is needed here: this probe tests that temporal
accounting property. The earlier energy probes continue to test latent-heat
identity and seasonal repartitioning. Unit regressions retain the old
criterion's passing result on cancelling errors, including alternating
hourly errors inside each phase, to show the specific difference.

### What can still pass

`reference_coupled` diagnoses H as `Rn - LE - G` and closes by construction.
It is an accounting positive control, not an independent thermal-inertia
model. A model that diagnoses a remainder, or makes equal and opposite
errors in simultaneous components, can still pass this probe. Passing does
not establish realistic turbulent exchange, component accuracy, radiative
temperature, soil heat storage or physical understanding. Heat
storage-temperature consistency is a separate question.

## Reproducing a failure

The directory was created with HydroTuring's `init-probe` scaffolding and
filled in for the scope accepted in
[proposal #15](https://github.com/Flood-Lab/HydroTuring/issues/15).

```bash
ht validate
ht gate --probe energy/surface-energy-closure
ht run --model reference_coupled --probe energy/surface-energy-closure --seed 0
ht run --model reference_diurnal_bias --probe energy/surface-energy-closure --seed 0
pytest -q tests/test_surface_energy_closure.py tests/test_surface_energy_probe.py
```

The first model must pass; the second must fail specifically on
`energy_closure_by_phase`. The tests also check weak-radiation allowances,
within-block cancellation, spinup exclusion, generator reproducibility,
complete scored phases, annotation stripping, and the minimum model window.
These checks establish the implemented accounting test, not a validation of
an AI or a comprehensive land-surface model.

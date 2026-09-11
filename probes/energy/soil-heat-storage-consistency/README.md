# energy/soil-heat-storage-consistency

## What it asserts

A model can close its surface energy budget while its reported soil temperature
never changes. This probe checks the energy retained in a declared soil layer:

```text
r = hfg - hfg_bottom - C_A * (T_end - T_start) / dt    [W m-2]
```

`hfg` is heat entering at the actual soil surface; `hfg_bottom` is heat leaving
at the layer's lower boundary. Both are downward-positive interval means.
`tsoil_layer` is the **mean temperature of that same layer at the interval
end**, in K. `C_A` is its constant areal heat capacity in J m-2 K-1; `dt` is
in seconds. Temperature is a diagnostic, never a water storage.

For each of the heating and recovery phases, require:

```text
mean(abs(r)) <= max(0.05 * mean(abs(hfg - hfg_bottom)), 1 W m-2)
```

Both phases must pass. Absolute residuals prevent signed cancellation; the
floor allows weak heat transfer. These are starting engineering tolerances,
not observational calibration. Phase averaging can still dilute brief errors.
Flux integration must match the temperature update; larger tolerances are
not a substitute for correcting an interval-mean versus endpoint mismatch.

## The case

The hourly record contains 24 hours of equilibrium spinup, 12 hours of heating
at positive net radiation, then 60 hours of recovery at zero net radiation.
The full three-day scored window is retained. Timestamps mark interval starts;
the last spinup row supplies the first scored starting temperature.

The surface is a bare, zero-capacity skin with no evaporation. Beneath it is
one fixed homogeneous layer, with no phase change, heat advection or internal
sources. Its lower boundary exchanges heat with a fixed-temperature reservoir
through a positive conductance. Air, layer and deep reservoir start at the
same temperature, so zero-radiation spinup is an exact equilibrium.

The seeded settings are:

| Setting | Range and units |
| --- | --- |
| Layer depth | 0.10-0.35 m |
| Volumetric heat capacity | 1.6-2.4 MJ m-3 K-1 |
| Sensible exchange coefficient | 10-15 W m-2 K-1 |
| Skin-to-layer conductance | 4-8 W m-2 K-1 |
| Bottom conductance | 0.5-1.5 W m-2 K-1 |
| Heating radiation | 120-240 W m-2 |
| Initial air, layer and deep temperature | 288.15-298.15 K |

`static.json` gives `soil_layer_depth_m`, `soil_heat_capacity_areal`,
`soil_temperature_initial`, `soil_deep_temperature`, `soil_top_conductance`,
`soil_bottom_conductance` and `soil_sensible_exchange`. The three conductances
use W m-2 K-1. These describe the test case; adapters must not fit heat capacity
or reconstruct a missing boundary flux from the residual being tested.

## Why these references

`reference_soil_heat` solves the layer's conductive heat equation analytically.
It integrates the native sensible and conductive flux equations over each
interval, and reports the endpoint layer temperature. No boundary flux is
inferred from the storage residual. The reference is a synthetic lumped-layer
model, not Noah or Noah-MP.

For other surface-energy probes that omit the thermal settings, the reference
uses a fixed 0.2 m layer with areal heat capacity 400,000 J m-2 K-1, sensible
exchange 12 W m-2 K-1, top conductance 6 W m-2 K-1 and bottom conductance
1 W m-2 K-1. Missing initial and deep temperatures use only the first forcing
row's air temperature. Supplied case settings always override these defaults;
the storage probe supplies every setting. Both broken controls share these
defaults. They configure the synthetic reference and do not fill missing
outputs from another model.

| Reference | Change to the correct solution | Expected storage verdict |
| --- | --- | --- |
| `reference_soil_heat` | None | PASS |
| `reference_frozen_soil` | Freeze only the reported temperature | FAIL |
| `reference_half_soil` | Halve only the reported temperature change | FAIL |

The broken controls retain the same surface and bottom heat fluxes. Their
surface energy budgets therefore still close; the new criterion catches the
disagreement between those fluxes and storage. Passing does not establish
accurate temperature or thermal inertia: equal top and bottom fluxes with
constant temperature, for example, satisfy this necessary budget condition.

The first version applies only when a model can represent the declared layer
and supply matching fluxes, heat capacity and temperature. It does not directly
score a multilayer model with evolving moisture or freeze-thaw storage against
this constant-capacity formula.

## Reproducing a failure

The optional [Noah-MP thermal-component check](../../../docs/noahmp-soil-heat-validation.md)
also evaluates native soil temperatures and conductive fluxes with this
criterion. It supplements these synthetic references; it is not a full
Noah-MP model adapter or a run of this probe's fixed-reservoir case.

```bash
ht gate --probe energy/soil-heat-storage-consistency
ht run --model reference_frozen_soil --probe energy/soil-heat-storage-consistency --seed 11
```

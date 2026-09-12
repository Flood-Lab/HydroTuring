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

The hourly record contains 24 hours of spinup without incoming shortwave,
12 hours of positive incoming shortwave, then 60 hours without shortwave.
The full three-day scored window is retained. Timestamps mark interval starts;
the last spinup row supplies the first scored starting temperature.

The case is warm and very dry, with a bare surface and no precipitation.
It specifies the depth, material heat capacity and initial temperature of
one homogeneous layer. Each model calculates its own surface exchange and
reports the heat flux through that layer's lower face. The deeper soil may
evolve; no fixed lower-temperature reservoir or particular exchange law is
required of evaluated models. No phase change, heat advection or internal
source is included in the layer balance.

Net radiation is a model calculation, not a prescribed forcing. All models
receive the same downward shortwave `rsds` and longwave `rlds` [W/m²], air
temperature `tas` [°C], wind speed `sfcWind` [m/s], specific humidity `huss`
[kg/kg], and surface pressure `ps` [Pa]. Radiation is constant within each
hour. Air is dry (`huss=0`), wind is 2 m/s and pressure is 101325 Pa.
Longwave is `sigma * T_initial**4`; the layer and air initially share that
temperature. This produces equilibrium for the ideal reference, while a
full model retains its native surface parameterizations.

The seeded settings are:

| Setting | Range and units |
| --- | --- |
| Layer depth | 0.10-0.35 m |
| Solid material volumetric heat capacity | 1.8-2.2 MJ m-3 K-1 |
| Heating-period incoming shortwave | 250-400 W m-2 |
| Initial air and soil temperature | 288.15-298.15 K |

`static.json` specifies layer depth, initial temperature, solid heat capacity,
porosity 0.464 and initial volumetric water content 0.005. The prescribed
constant areal capacity is calculated before any model runs:

```text
C_A = depth * [(1-porosity)*C_solid + theta*C_water + (porosity-theta)*C_air]
```

The test material uses `C_water=4188000` and effective pore-air coefficient
`C_air=1004.64`, in J m⁻³ K⁻¹. These explicit mixture settings match the
documented Noah-MP dry-soil parameterization; they are not values fitted from
its output. `soil_heat_capacity_areal` is the resulting fixed value used by
the criterion. Adapters configure the material and layer from these inputs,
never estimate capacity from fluxes or temperature changes.

A real model can retain tiny moisture changes from numerical surface updates.
In a dry-case validation, quantify their storage-rate effect using native
capacity and compare it with the existing allowance. Do not silently replace
the prescribed capacity with a time-varying or fitted value. Wet cases with
materially evolving heat capacity need a different physical assessment.

## Why these references

`reference_soil_heat` solves a massless radiative skin above a conductive
layer. It computes absorbed and emitted radiation, sensible exchange and
Fourier fluxes from their own equations, then integrates the temperature and
those fluxes with the same RK4 quadrature. No boundary flux is inferred from
the storage residual. Its fixed deep reservoir and exchange coefficients are
properties of this synthetic reference, not requirements on other models.

For other surface-energy probes that omit the thermal settings, the reference
uses a fixed 0.2 m layer with areal heat capacity 400,000 J m-2 K-1, sensible
exchange 12 W m-2 K-1, top conductance 6 W m-2 K-1 and bottom conductance
1 W m-2 K-1. Missing initial and deep temperatures use only the first forcing
row's air temperature. Supplied settings override those defaults in the older
net-radiation input format. Both broken controls share these
defaults. They configure the synthetic reference and do not fill missing
outputs from another model. The new shared-radiation case supplies the layer
and material; the reference's exchange coefficients remain model parameters.

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
and prescribed material, and supply matching fluxes and temperature. A
multilayer model can report the upper declared layer and its actual lower-face
flux. Materially changing heat capacity or freeze-thaw storage cannot be
ignored when applying this constant-capacity formula.

## Reproducing a failure

The optional [Noah-MP checks](../../../docs/noahmp-soil-heat-validation.md)
include a fixed-moisture thermal component, an exploratory moist HRLDAS run,
and complete HRLDAS runs on this generator's five actual cases. The generated
case validation uses the unchanged registered criterion and prescribed
capacity. Native capacity drift is quantified separately. This optional
physical validation does not register Noah-MP as a benchmark model or add it
to the mandatory reference gate.

```bash
ht gate --probe energy/soil-heat-storage-consistency
ht run --model reference_frozen_soil --probe energy/soil-heat-storage-consistency --seed 11
```

# mass/time-origin-invariance

A stationary catchment model given identical weather, static attributes,
initialization and random seed should give identical water fluxes and
storages when only the absolute calendar origin changes. A year-dependent
recession coefficient can violate this symmetry while closing the water
budget in both runs, so closure alone cannot detect the fault.

This probe implements [accepted proposal #10](https://github.com/Flood-Lab/HydroTuring/issues/10).

## Paired construction

Each seed generates 4,015 daily rows: 365 spinup steps followed by 3,650
scored steps under the full-window setting. `control` starts on 2000-01-01;
`shifted` starts on 1972-01-01. The 28-year shift preserves month, day,
weekday, day of year and leap-day positions throughout this specific
window. This is not a claim that every 28-year shift preserves a Gregorian
calendar, particularly across non-leap century boundaries.

The ten scored years give repeated wet and dry seasons and several leap
transitions over which an absolute-year dependency can appear. This is a
longer replay check than the three-year area-invariance generator; ten
years is a fixed test horizon, not a minimum required for the symmetry or
a claim of calibrated statistical power. Submitted models still use a
shorter event window as described below.

The non-time forcing columns and static attributes are exactly identical
between variants. Seasonal forcing uses actual calendar phase. Weather is
deliberately above freezing so a rainfall-response check is appropriate;
empty snow storage remains a valid outcome. Each run starts afresh using
the same model seed and initialization contract. Both variants use the
same scored row indices. The probe requests at least a 365-day scored
window so the annual runoff-ratio guard remains active for shorter model
evaluations. The 300-second budget is per variant, matching the other
paired probes and allowing for the longer submitted-model window.

## Criteria and tolerance

For each compared variable `x`, invariance requires:

```text
max_t |x_shifted(t) - x_control(t)|
-------------------------------- <= 1e-9
max(mean_t |x_control(t)|, 1e-12)
```

This uses the existing invariance criterion and its 1e-12 denominator
minimum. The relative limit is deliberately strict: all four physical
references replay exactly, including empty snow storage. No new absolute
tolerance is introduced. Required comparisons cover runoff
(`mrro`), evaporation (`evspsbl`), soil water (`mrso`), snow (`snw`) and
canopy water (`canopy`). Groundwater (`gw`), channel storage (`channel`)
and groundwater exchange (`gwex`) are compared when reported in both
variants. The criterion measures the maximum pointwise difference after
spinup, not just differences in integrated totals.

The existing `closure` criterion checks the control water budget:

```text
|sum(P + declared gwex - ET - runoff) - (S_end - S_start)| / sum(P) <= 0.05
```

Fluxes are integrated over the scored window; `S` includes required stores
and all additional recognized stores reported by the model. `S_start` is
the state immediately before that window. The denominator comes from the
input precipitation. A zero total fails closure rather than receiving an
artificial denominator floor. Its residual remains in the scorecard.
The shifted run is tested for invariance; it does not receive a separate
closure score. A calendar-independent leak in both runs still fails the
control budget. The reference-calendar test additionally checks that its
shifted budget closes, demonstrating a fault closure alone misses.

`non_degenerate` checks the control's runoff ratio, variation in runoff and
evaporation, and response to rainfall. Zero-runoff or constant-flux outputs
must not receive credit for invariance alone. A forcing-responsive fixed
partition such as `ET = 0.5 P` and `runoff = 0.5 P` can still pass these
criteria; the guards do not establish that a model represents hydrological
processes correctly.

## Discrimination

| Reference model | Expected result | Required failure |
| --- | --- | --- |
| `reference_bucket` | PASS | |
| `flex_lumped` | PASS | |
| `flex_topo` | PASS | |
| `sacsma_snow17` | PASS | |
| `reference_calendar` | FAIL | `invariance` |
| `reference_degenerate` | FAIL | `non_degenerate` |
| `reference_leaky` | FAIL | `closure` |

The calendar-dependent model is the key negative control: its water
budgets can close and its fluxes remain responsive, yet its paired outputs
change. The other two controls exercise the guards against trivial
invariance and an unbalanced budget.

## Scope and reproduction

```sh
ht validate
ht gate --probe mass/time-origin-invariance
pytest -q tests/test_time_origin_invariance.py tests/test_symmetry_outputs.py
```

This is a stationary synthetic experiment. It does not assert that a model
should ignore a changed physical driver such as land use or CO2, and it
does not evaluate such time-varying drivers. Adapters must provide the
required water-budget outputs and use repeatable initialization. The strict
tolerance targets deterministic replay with the same seed; stochastic
hardware variation is not calibrated here. Passing this probe is evidence
for this symmetry and budget check, not proof of general hydrological skill.

The generator adapts the project's `mass/area-invariance` weather and
catchment construction. No generated forcing or model output is committed.

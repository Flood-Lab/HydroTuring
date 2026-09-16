# Groundwater-surface-water exchange consistency

## What it asserts

This probe uses `gw_sw_exchange`, the net river-aquifer exchange, positive
into the aquifer. This is deliberately not `gwex`: `gwex` is a declared
exchange with the outside of the catchment (regional groundwater, an
inter-basin transfer) and `closure.py` counts it as a source in the
whole-catchment budget. River-aquifer exchange moves water between two
stores that are both inside a catchment model's control volume (`gw` and
`channel`), so a catchment model that also declared it as `gwex` would double
its water balance's opinion of that flow, and `mass/catchment-closure` would
misread an internal transfer as a source or sink crossing the boundary. This
probe also accepts two signed directional components:

- `sw_to_gw >= 0`: surface water losing to the aquifer (into the aquifer,
  same sign as `gw_sw_exchange`)
- `gw_to_sw <= 0`: aquifer losing to surface water (out of the aquifer)

They must satisfy `gw_to_sw + sw_to_gw = gw_sw_exchange` to reporting
precision: this is bookkeeping, not physics, and holds exactly for any model
that means it. The groundwater balance is:

```text
gw_recharge + gw_sw_exchange + gwex = change in groundwater storage (gw)
```

`gwex`, a source or sink crossing the control volume's boundary such as a
GHB or WEL package or a regional groundwater exchange, is added if the model
declares it and is otherwise zero: a model without that column is not
penalized for a boundary term it does not have, but one with a real
boundary term and no `gwex` column is scored on an incomplete budget and
fails. The only required store is `gw`. This is intentionally a groundwater
probe: canopy, soil, snow and catchment runoff are outside its control
volume, and a model is scored only on the groundwater state and exchange
fluxes it reports.

### Static parameters

| Key | Meaning | Units |
| --- | --- | --- |
| `area_km2` | catchment/model area the fluxes are converted over | km2 |
| `aquifer_specific_yield` | drainable porosity above the water table | dimensionless |
| `aquifer_storage_coefficient` | storage coefficient (storativity), not specific storage | dimensionless |
| `river_conductance_m2_per_day` | RIV package conductance between the river and the aquifer | m2/day |
| `river_bottom_offset_m` | river stage minus riverbed elevation | m |
| `aquifer_initial_head_m` | initial aquifer head, optional, defaults to 10.0 | m |

## The case

`generate.py` creates two years of daily forcing after a 90-day spinup
(`period_years: 2`, `spinup_days: 90`). It draws seeded daily groundwater
recharge (`gw_recharge`, mm/day, AR(1) noise on a mild seasonal cycle) and a
river stage (`sw_stage_m`, metres, an AR(1)-perturbed sine with a 90-day
period around a fixed reference elevation). The stage cycle drives the
aquifer through both gaining and losing conditions, so a correct model's
exchange genuinely reverses sign over the record. `case.min_window_days: 180`
widens a submitted model's default evaluation window to at least two stage
cycles, since a shorter window can land entirely on one side of the cycle
regardless of the model's physics.

## Criteria

| Criterion | Purpose |
| --- | --- |
| `groundwater_balance` | Ensures recharge, net exchange and groundwater storage agree. |
| `exchange_components` | Ensures the two signed directional terms sum to `gw_sw_exchange` and obey their signs. |
| `state_bounds` | Prevents negative aquifer storage. |
| `exchange_directions` | Requires both gaining and losing river exchange over the scored record. |

`reference_exchange_exact` is the physical `must_pass` baseline: a trusted,
Docker-free bookkeeping reference that closes its own groundwater balance
exactly and reports genuinely bidirectional exchange. CI's gate job runs
`ht gate`, which runs `reference_exchange_exact` as that `must_pass`, so it
is exercised on every run. `models/modflow6`, a real MODFLOW 6.7.0 binary
run through Docker, is a second, evidential baseline archived as a
submitted model rather than gated on, so the acceptance gate does not
depend on a 134 MB release download and a `flopy` build to score this
probe; it is not exercised by CI, and its own archive row is written
separately with
`ht run --model modflow6 --probe mass/gw-sw-exchange-consistency --gate-seeds --csv models/result.csv`.
`reference_bucket`, `flex_lumped`, `flex_topo` and `sacsma_snow17` do not
report `gw_sw_exchange`/`gw_to_sw`/`sw_to_gw`, so the probe is
`N/A (INCOMPLETE)` for them: a missing output is checked before an
unconsumed input, the same outcome as for any model without a declared
groundwater-exchange term.

The groundwater balance is checked two ways. Each step must be within a 1%
relative tolerance or a per-step floor, taken over the larger of the
recharge and exchange terms so the allowance is not set by whichever is
smaller. The floor is `max(2e-4 mm, 1e-5 * max|gw|)`: the fixed `2e-4 mm`
covers the recharge/exchange terms' own rounding, and the second term
covers the `gw` state's own rounding separately, scaled to the state's
magnitude rather than fixed, since that magnitude is the model's choice of
datum and not a probe constant. The record's cumulative residual must also
be within `max(1% of total recharge, 2x the per-step floor)`. Rounding a
stored state cancels out once the steps are summed; a per-step shortfall
that recurs in the same direction, such as an unmodeled leak out of the
aquifer, does not, and can otherwise hide under a per-step floor sized for
rounding. Honest output, at full precision or rounded to 6 significant
figures and at any datum, closes both checks. `exchange_components` is an
exact identity and uses a `1e-6` relative tolerance and a `1e-4 mm`
absolute floor: a reach that both gains and loses on the same step has
`gw_to_sw` and `sw_to_gw` each rounded to ordinary output precision before
they are summed, which the identity's own rounding does not absorb at a
tighter floor.

## Baselines

`reference_exchange_exact` and `reference_exchange_sign_error` share the
same lagged-head bookkeeping reference (its RIV-style exchange term is
computed the same way MODFLOW's RIV package would, but the head tracks a
smoothed river stage directly rather than integrating net flux through a
storage coefficient, so it is not a linear-reservoir aquifer), which closes
its own groundwater balance exactly and never reports negative storage.
`reference_exchange_sign_error` deliberately reports the aquifer-to-river
component (`gw_to_sw`) without its sign on every other step, so
`exchange_components` is the only criterion it trips.

## Local checks

From the repository root:

```bash
ht validate
ht gate --probe mass/gw-sw-exchange-consistency
```

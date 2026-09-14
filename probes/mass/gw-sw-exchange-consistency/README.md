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
gw_recharge + gw_sw_exchange = change in groundwater storage (gw)
```

The only required store is `gw`. This is intentionally a groundwater probe:
canopy, soil, snow and catchment runoff are outside its control volume, and a
model is scored only on the groundwater state and exchange fluxes it reports.

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
Docker-free linear-reservoir aquifer that closes its own groundwater balance
exactly and reports genuinely bidirectional exchange. `models/modflow6`, a
real MODFLOW 6.7.0 binary run through Docker, is a second, evidential
baseline archived as a submitted model rather than gated on, so the
acceptance gate does not depend on a 134 MB release download and a `flopy`
build to score this probe. Neither is exercised by CI's Docker-free harness
tests; `modflow6`'s own archive row is written separately with
`ht run --model modflow6 --probe mass/gw-sw-exchange-consistency`.
`reference_bucket`, `flex_lumped`, `flex_topo` and `sacsma_snow17` do not
consume `gw_recharge`/`sw_stage_m`, so the probe is `N/A (INCOMPATIBLE)` for
them, the same outcome as for any model without a declared
groundwater-exchange term.

The groundwater balance uses a 1% relative tolerance and a `2e-4 mm` absolute
per-step floor, taken over the larger of the recharge and exchange terms so
the allowance is not set by whichever is smaller. `exchange_components` is an
exact identity and uses a `1e-6` relative tolerance.

## Baselines

`reference_exchange_exact` and `reference_exchange_sign_error` share the
same lagged-head linear-reservoir aquifer (computed the same way MODFLOW's
RIV package would), which closes its own groundwater balance exactly and
never reports negative storage. `reference_exchange_sign_error` deliberately
reports the aquifer-to-river component (`gw_to_sw`) without its sign on
every other step, so `exchange_components` is the only criterion it trips.

## Local checks

From the repository root:

```bash
ht validate
ht gate --probe mass/gw-sw-exchange-consistency
```

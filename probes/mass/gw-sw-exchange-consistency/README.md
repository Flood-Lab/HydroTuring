# Groundwater-surface-water exchange consistency

## What it asserts

The HydroTuring output contract represents net groundwater exchange with the
signed `gwex` flux, positive into the aquifer (AGENTS.md). This groundwater-only
probe also accepts two signed directional components:

- `sw_to_gw >= 0`: surface water losing to the aquifer (into the aquifer,
  same sign as `gwex`)
- `gw_to_sw <= 0`: aquifer losing to surface water (out of the aquifer)

They must satisfy `gw_to_sw + sw_to_gw = gwex` to reporting precision: this is
bookkeeping, not physics, and holds exactly for any model that means it. The
groundwater balance is:

```text
gw_recharge + gwex = change in groundwater storage (gw)
```

`gwex` is a source in this budget, consistent with how `closure.py` treats it
in the whole-catchment probes.

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
| `exchange_components` | Ensures the two signed directional terms sum to `gwex` and obey their signs. |
| `state_bounds` | Prevents negative aquifer storage. |
| `exchange_directions` | Requires both gaining and losing river exchange over the scored record. |

`models/modflow6` is the executable baseline for this probe. It runs a real
MODFLOW 6.7.0 binary, reads the RIV and STO budget terms, converts cubic
metres per day to millimetres per day over the model area, and reports the
signed components and groundwater storage without fabricating them. It is
the only baseline that can be put to this probe today: `reference_bucket`,
`flex_lumped`, `flex_topo` and `sacsma_snow17` do not consume
`gw_recharge`/`sw_stage_m`, so the probe is `N/A (INCOMPATIBLE)` for them, the
same outcome as for any model without a declared groundwater-exchange term.

The groundwater balance uses a 1% relative tolerance and a `2e-4 mm` absolute
per-step floor, taken over the larger of the recharge and exchange terms so
the allowance is not set by whichever is smaller. `exchange_components` is an
exact identity and uses a `1e-6` relative tolerance.

## Baselines

`reference_exchange_sign_error` is a purpose-built reference: a linear-
reservoir aquifer with a genuine bidirectional river exchange (computed the
same way MODFLOW's RIV package would) that closes its own groundwater
balance exactly and never reports negative storage. It deliberately reports
the aquifer-to-river component (`gw_to_sw`) without its sign, so
`exchange_components` is the only criterion it trips.

## Local checks

From the repository root:

```bash
ht validate
ht gate --probe mass/gw-sw-exchange-consistency
```

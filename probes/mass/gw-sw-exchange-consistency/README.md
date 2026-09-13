# Groundwater-surface-water exchange consistency

## What it asserts

The current HydroTuring output contract represents net exchange with the signed
`gwex` flux, positive into the aquifer. This groundwater-only probe also accepts
two signed directional components:

- `gw_to_sw >= 0`: groundwater contribution to surface water
- `sw_to_gw <= 0`: surface-water loss to groundwater

They must satisfy `gw_to_sw + sw_to_gw = gwex`. The net term remains `gwex`,
so positive values are a net aquifer loss to the river and negative values are a
net river gain to the aquifer. The groundwater balance is:

```text
recharge - gwex = change in groundwater storage
```

The required stores include `gw`, `channel`, and the ordinary soil, snow and
canopy stores. Omitting groundwater or water in transit would leave the budget
open by exactly the amount held there.

## The case

`generate.py` creates five years of daily forcing after a 365-day spinup. It
uses seeded temperate weather with precipitation, air temperature and potential
evapotranspiration. It also creates private groundwater-head and surface-stage
annotations whose difference changes sign, preserving a genuinely bidirectional
case for a future gross-exchange criterion. Columns beginning with `_` are
removed by the harness before a model sees the forcing.

## Criteria

| Criterion | Purpose |
| --- | --- |
| `groundwater_balance` | Ensures recharge, net exchange and groundwater storage agree. |
| `exchange_components` | Ensures the two signed directional terms sum to `gwex` and obey their signs. |
| `state_bounds` | Prevents negative aquifer storage. |
| `exchange_directions` | Requires both gaining and losing river exchange over the scored record. |

`models/modflow6` is the executable baseline for this probe. It runs the local
MODFLOW 6 binary, reads the RIV budget and head output, converts cubic metres
per day to millimetres per day over the model area, and reports the signed
components without fabricating them.

The groundwater balance uses a 1% relative tolerance and a `2e-4 mm` absolute
per-step floor. The floor covers the small discrete-storage and binary-budget
rounding error of MODFLOW 6; it is still far below the generated daily recharge
and does not mask a hydrologic exchange term.

## Local checks

From the repository root:

```powershell
python -m py_compile probes\mass\gw-sw-exchange-consistency\generate.py
ht validate
```

`reference_bucket` is the first in-repository adapter updated to emit the two
components. Other models must expose genuine directional diagnostics before
they can be scored by this probe; they are not assigned fabricated zero columns.

This is a local addition only. No commit or push is performed by this work.
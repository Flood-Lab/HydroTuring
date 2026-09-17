# MODFLOW 6 groundwater exchange adapter

The adapter resolves MODFLOW 6 in this order:

1. `MODFLOW6_EXE`, if set;
2. on Linux, `/opt/modflow6/bin/mf6`, `/opt/modflow6/mf6`,
   `/usr/local/bin/mf6`, or `PATH`;
3. every executable named `mf6` (or `mf6.exe` on Windows) under
   `MODFLOW6_ROOT`, if set;
4. `mf6`/`mf6.exe` on `PATH`.

The supplied Dockerfile downloads the official MODFLOW 6.7.0 Linux asset
`mf6.7.0_linux.zip` and verifies SHA-256 before installing it. The release is
published at https://github.com/MODFLOW-ORG/modflow6/releases/tag/6.7.0.

The model is a transient 10 x 10 aquifer with one RIV cell; the grid's cell
size is derived from `area_km2` (not fixed), so the RIV/RCHA/STO flows are
always converted over the area the recharge actually fell on. The probe
supplies daily recharge and river stage, while `static.json` supplies the
aquifer properties (specific yield, storage coefficient), river conductance
and river-bottom offset. MODFLOW 6 supplies the RIV and STO budget terms; the
adapter converts them from m3/day to mm/day over the model area and
accumulates groundwater storage from the STO-SS and STO-SY release terms
rather than reconstructing it from mean head.

Two simplifications from a single RIV cell:

- `gw_to_sw` and `sw_to_gw` are never both nonzero on the same step, because
  one cell has one net flow direction per step. A case that reported "100 in,
  30 out, net 70" as "net 100" — the failure mode
  [Flood-Lab/HydroTuring#40](https://github.com/Flood-Lab/HydroTuring/issues/40)
  describes for a catchment with several gaining and losing reaches at once —
  cannot arise here. A second river cell at a different stage, with the
  adapter splitting each cell's `q` by sign before summing, would exercise
  that; this baseline does not.
- `river_bottom_offset_m` sets the riverbed at `sw_stage_m - river_bottom_offset_m`,
  so it moves with the stage rather than staying at a fixed elevation.
  MODFLOW's RIV package caps river-to-aquifer leakage at the conductance
  times this offset whenever the head drops below the riverbed, so
  `sw_to_gw` cannot exceed `river_conductance_m2_per_day * river_bottom_offset_m`
  converted to mm/day (about 0.15 mm/day at this case's defaults) no matter
  how low the head falls. A fixed riverbed elevation set below the lowest
  stage in the record would let that cap reflect an actual streambed instead
  of trailing the stage.

Sign mapping, per the probe's `gw_sw_exchange`-positive-into-the-aquifer
convention (AGENTS.md) — not `gwex`, which is a declared exchange with the
outside of the catchment rather than river-aquifer exchange within it:

- MODFLOW RIV `q > 0` is river into aquifer, reported as positive `sw_to_gw`.
- MODFLOW RIV `q < 0` is aquifer into river, reported as negative `gw_to_sw`.
- `gw_sw_exchange = gw_to_sw + sw_to_gw`.

The image is pinned to `linux/amd64` (native on x86-64, emulated on
arm64/Apple silicon): MODFLOW 6.7.0 ships no Linux arm64 build. `modflow6` is
archived as a second, evidential baseline alongside the trusted
`reference_exchange_exact`; it is not in the probe's `must_pass`, so the
acceptance gate does not need Docker or a network fetch of the MODFLOW 6
release to score `mass/gw-sw-exchange-consistency`.

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

The model is a transient 10 x 10 aquifer with one RIV cell. The probe supplies
daily recharge and river stage, while `static.json` supplies the aquifer
properties (specific yield, storage coefficient), river conductance and
river-bottom offset. MODFLOW 6 supplies the RIV and STO budget terms; the
adapter converts them from m3/day to mm/day over the model area and
accumulates groundwater storage from the STO-SS and STO-SY release terms
rather than reconstructing it from mean head.

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

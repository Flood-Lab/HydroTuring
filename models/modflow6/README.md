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

Sign mapping, per the repository's `gwex`-positive-into-the-aquifer
convention (AGENTS.md):

- MODFLOW RIV `q > 0` is river into aquifer, reported as positive `sw_to_gw`.
- MODFLOW RIV `q < 0` is aquifer into river, reported as negative `gw_to_sw`.
- `gwex = gw_to_sw + sw_to_gw`.

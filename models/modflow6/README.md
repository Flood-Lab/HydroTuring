# MODFLOW 6 groundwater exchange adapter

The adapter resolves MODFLOW 6 in this order:

1. `MODFLOW6_EXE`, if set;
2. on Windows, `D:\\hydrowang\\modflow\\bin\\mf6.exe`;
3. on Linux, `/opt/modflow6/bin/mf6`, `/opt/modflow6/mf6`,
   `/usr/local/bin/mf6`, or `PATH`;
4. every executable named `mf6` under `MODFLOW6_ROOT`, if set.

The supplied Dockerfile downloads the official MODFLOW 6.7.0 Linux asset
`mf6.7.0_linux.zip` and verifies SHA-256 before installing it. The release is
published at https://github.com/MODFLOW-ORG/modflow6/releases/tag/6.7.0.

The model is a transient 10 x 10 aquifer with one RIV cell. The probe supplies
daily recharge and river stage, while `static.json` supplies the aquifer
properties, river conductance and river-bottom offset. MODFLOW 6 supplies the
RIV budget flow and head field; the adapter converts the RIV flow from m3/day
to mm/day over the model area and reports groundwater storage from specific
yield and mean head.

Sign mapping:

- MODFLOW RIV `q > 0` is river into aquifer, reported as negative `sw_to_gw`.
- MODFLOW RIV `q < 0` is aquifer into river, reported as positive `gw_to_sw`.
- `gwex = gw_to_sw + sw_to_gw`.

The local subprocess run still uses the Windows installation when present; the
Dockerfile provides the Linux packaging path for a portable container run.

# sacsma_snow17

SAC-SMA (Burnash, Ferral & McGuire 1973; Anderson's 1979 code) with
Snow-17 (Anderson 1973; the 1995–2001 revisions by Duan and Koren) and a
gamma unit hydrograph, the US National Weather Service's operational
conceptual model. In the pool as a **physical reference**: every probe must
pass it, and a probe that fails it is examined before the model is.

## Provenance

A standard-library port (`sacsma_snow17.py`) of the legacy Fortran
packaged by [Upstream-Tech/SACSMA-SNOW17](https://github.com/Upstream-Tech/SACSMA-SNOW17)
at commit `8737b92`: `sacsma_source/sac/sac1.f`, `ex_sac1.f`, `duamel.f`
and `snow19/*.f`. The port keeps the Fortran's variable names and branch
order so it can be read against the original, and it was checked against
the f2py build of that Fortran (gfortran 15, `-fallow-argument-mismatch
-std=legacy`) on the closure probe's ten-year record:

- **SAC-SMA** reproduces the Fortran to 1e-4 mm per step in channel inflow
  and to 1e-6 mm in evaporation on identical input (float32 noise).
- **Snow-17** agrees with the Fortran in ten-year totals of rain-plus-melt
  to 0.02 % (8921.8 vs 8923.7 mm) but not step by step. Instrumenting the
  Fortran's `AESC19` shows it computing an areal cover of 0.99 on the
  depletion-curve branch while the carryover array holds `SBAESC = 0` on
  every one of 4015 steps: the compiled legacy code does not carry that
  variable across steps, and the port, which follows the source as written,
  does. The Fortran also leaves `tiny` in `AESC19` uninitialised (the port
  uses 1e-6) and never sets `SNOF` (zero in both).
- The **unit hydrograph** kernel agrees with `DUAMEL` (which takes an
  integer step of days; the port takes a real one so the hydrograph keeps
  its shape in days at the hourly step).

Left out of the port because they move no water or need observations: the
frozen-ground option (`IFRZE = 0` in the wrapper), the observation
updating routines (`UPDT19`, `ADJC19`, `AECO19`), the snow depth and
density diagnostics, and the user-specified melt-factor curve.

## Choices that make it a conservative reference

| Choice | Why |
| --- | --- |
| Snow-17 `SCF = 1.0` | the gauge-catch multiplier manufactures snow above one; a physical reference is fed what fell |
| SAC-SMA `SIDE = 0.02`, declared as `gwex` | deep baseflow leaves the catchment for good; declared as a negative exchange so the budget closes over what the model says it did. When the forcing carries `abstr`, the prescribed withdrawal is removed from the SAC stores and the day's runoff and added to `gwex` the same way |
| `RIVA = 0.02` in `evspsbl` | riparian evaporation from channel inflow is evaporation, as the Fortran counts it |
| `UZTWM + UZFWM + LZTWM` = the catchment's soil capacity | a physical model is told its catchment; the three tension/upper stores are rescaled from the default set keeping their ratios, lower-zone free water keeps its defaults |
| `PXTEMP`, `MBASE` = the catchment's snow threshold | same reason |
| `canopy` ≡ 0 | SAC-SMA has no interception store |

Default parameters are a mid-range set from the NWS calibration guidance
(Anderson 2002) and CAMELS-scale calibrations; they are printed in
`run.json`.

## What is reported

Catchment-average, weighting the pervious-area stores by `PAREA` and the
additional-impervious store by `ADIMP`, as the Fortran weights the runoff
components: `mrso` = UZTWC + UZFWC + LZTWC (+ ADIMC on its area), `gw` =
LZFSC + LZFPC, `snw` = Snow-17's total water equivalent (pack plus liquid
plus lagged excess plus storage), `channel` = channel inflow generated but
not yet released by the hydrograph, `gwex` = minus the non-channel
baseflow, plus any prescribed human withdrawal removed this step. Steps:
PT1D and PT1H, the whole hours Snow-17 is defined on.

## Result

**PASS, 17 of 17 probes passed**, with the gate seeds (`ht run --model
sacsma_snow17 --gate-seeds`).

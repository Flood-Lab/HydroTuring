# Noah-MP soil heat validation

These optional experiments supplement the generated probe's reference gate.
The component experiment checks `soil_heat_storage` against native Noah-MP
soil temperature calculations with fixed moisture. The full HRLDAS experiment
runs the complete column, including water updates, and inspects the native
heat equation and the fixed-capacity criterion's applicability. Neither
registers a Noah-MP model adapter or runs Noah-MP on the generated probe case.

## Run the fixed-moisture component

Install this checkout with `pip install -e .` and provide GNU Fortran. From the
repository root, choose a scratch directory for the downloaded sources and
outputs:

```bash
python scripts/noahmp_soil_heat/validate.py --workdir /tmp/noahmp-soil-heat --fc gfortran
```

The script downloads the required official modules from
[NCAR/noahmp v5.2.1, commit `17751dcd7a2442a3c57138f530f523daf08dc1c5`](https://github.com/NCAR/noahmp/tree/17751dcd7a2442a3c57138f530f523daf08dc1c5),
compiles in double precision, runs two time steps, and calls HydroTuring's
registered `soil_heat_storage` criterion. It writes two hourly CSV files and
`hydroturing_criterion_validation.json` in the scratch directory. No NetCDF,
MPI, forcing download, or GPU is needed. This is not a CI dependency.

To evaluate existing CSV files on a machine without a compiler or network:

```bash
python scripts/noahmp_soil_heat/validate.py --workdir /tmp/noahmp-soil-heat --skip-build
```

Only the custom driver and validation script are included here. Downloaded
Noah-MP code keeps its [UCAR license](https://github.com/NCAR/noahmp/blob/17751dcd7a2442a3c57138f530f523daf08dc1c5/LICENSE.txt)
and attribution in the scratch directory; HydroTuring's license does not replace
it. The Noah-MP modeling system was developed at the National Center for
Atmospheric Research (NCAR) with collaborations from university partners.
NCAR is sponsored by the United States National Science Foundation.

## Control volume and flux mapping

The driver calls the official `SoilThermalProperty`,
`SoilSnowThermalDiffusion`, `SoilSnowTemperatureSolver`, and
`MatrixSolverTriDiagonal` routines. Soil moisture stays fixed at 0.2, porosity
at 0.45, quartz fraction at 0.4, and solid heat capacity at 2 MJ m⁻³ K⁻¹.
Native thermal properties are calculated before stepping. There is no snow,
ice, water movement, phase change, or penetrating radiation.

The evaluated layer is the top 0.1 m of four native layers. Its heat capacity
per area comes from the native volumetric heat capacity times 0.1 m; its
temperature is `TemperatureSoilSnow(1)`. The other layers provide its lower
boundary interaction. This differs from the generated probe's fixed deep
temperature reservoir, so this experiment validates the budget criterion
under an additional physical boundary condition.

`HeatGroundTotMean` supplies the prescribed downward top flux. A diagnostic
output added to a copy of the native solver exposes its downward lower-face
Fourier flux, `2*k1*(T1-T2)/(-z2)`, after the implicit temperature update. The
temperature equations are unchanged. This flux comes from the native
conductivity and temperature gradient, not from a storage residual.

Hourly rows contain endpoint temperatures and averages of the fluxes used by
the native time integration. For the implicit solver, substituting the
hour-start flux changes that budget. The bottom-most `HeatFromSoilBot` is not
the lower flux for this top-layer control volume.

## Measured result

The complete download, compile, run, and criterion evaluation completed with
GNU Fortran 13.4.0 on 2026-09-12. Each run contains
24 h equilibrium at 285 K, 12 h with 80 W m⁻² top input, and 60 h with zero
top input. Both heating and recovery are scored separately. Native
`C_A = 193785.116 J m⁻² K⁻¹`; temperatures stayed between 285 and 296.133 K.

| Output supplied to the criterion | 300 s native step | 3600 s native step |
| --- | --- | --- |
| Native temperatures and matching fluxes | PASS | PASS |
| Temperature held at its initial value | FAIL, both phases | FAIL, both phases |
| Temperature departures halved | FAIL, both phases | FAIL, both phases |
| Hour-start lower flux substituted for its interval value | PASS | FAIL, heating |

Native phase mean absolute residuals were below `3.1e-12 W m⁻²`. In the 300 s
run, the frozen-temperature control had residuals of 49.94 W m⁻² during
heating and 7.41 W m⁻² during recovery, exceeding allowances of 2.50 and
1 W m⁻². These controls corrupt reported temperatures; they are not defects
observed in Noah-MP.

The deliberately misaligned flux produced 2.21 W m⁻² against 2.61 W m⁻²
allowed at 300 s, and 3.91 against 2.65 at 3600 s during heating. The script
records this diagnostic without requiring it to fail: the tolerance cannot
detect every timing error. It does require native output to pass and both
temperature controls to fail.

Conservation does not establish temperature accuracy. The two time steps
produced slightly different temperature responses while both conserved the
tested layer's energy.

## Run the complete HRLDAS column

The second experiment uses the official
[HRLDAS driver at `cd96df4`](https://github.com/NCAR/hrldas/tree/cd96df470220f7d7133cdbccd5f9c5355cf173e2)
and its Noah-MP submodule at the same `17751dc` revision as above. It executes
the complete land column, including surface exchange and water updates.
The meteorological forcing is synthetic; the bundled Bondville example
supplies the official input format and converter, not observations for this run.

On Linux, provide `git`, `make`, GNU Fortran, NetCDF C/Fortran and install this
checkout with `pip install -e '.[netcdf]'`. Use a fresh scratch directory:

```bash
python scripts/noahmp_soil_heat/run_full.py \
  --workdir /tmp/hrldas-soil-heat \
  --fc gfortran --netcdf-prefix /path/to/netcdf
```

The script downloads the pinned official source, builds a serial executable,
converts the point forcing to NetCDF, runs the original model, adds two
diagnostic CSV outputs, rebuilds and reruns the same case. It compares all
native NetCDF variable values between the two runs before evaluating the
budget. Source licenses remain in the downloaded repositories. The default
single precision is retained. The unrelated GRIB forcing converter is not
built; this experiment uses the official text-to-NetCDF converter.

The 96-hour case starts at 2020-06-20 06:00 UTC. It has four native soil
layers, bare ground, no rain or snow, initial soil and air temperature 295 K,
initial volumetric water content 0.20, wind 2 m/s and specific humidity 0.010.
After 24 h without shortwave radiation, incoming shortwave is 300 W/m² for
12 h, then zero for 60 h. Longwave stays at approximately 429.44 W/m².
Forcing, model, soil and output steps are all 300 s. End-time forcing records
align the pulse with the intended phases; hourly forcing would let the
driver interpolate across the phase boundaries.

`native_thermal_steps.csv` records the temperatures, heat capacity and
boundary fluxes used by each thermal solve. `native_column_steps.csv` records
the full model step, including water updates. The lower-face flux comes from
the same native Fourier relation as in the component experiment. Neither
flux is reconstructed from storage. The standard `SOILENERGY` output covers
the entire soil column with current heat capacities, so it cannot replace
this first-layer diagnostic.

To inspect saved CSV files without compiling or downloading again:

```bash
python scripts/noahmp_soil_heat/score_full.py \
  --workdir /tmp/hrldas-soil-heat/case/diagnostic
```

## Full-run result and applicability

The actual run completed 1,152 model steps and wrote 1,153 standard output
records, including the initial state. A fresh invocation of `run_full.py`
completed the entire download/build/run/evaluation sequence on `climet3` with
GNU Fortran 13.4.0, NetCDF Fortran 4.6.3 and NetCDF C 4.10.1. All 141 native
output variables were exactly equal before and after instrumentation.
Native `SWFORC` confirms exactly 144 heated steps with correctly aligned
phase boundaries. The evaluated top layer is 0.1 m.
Snow, soil ice, vegetation, penetrating shortwave and precipitation heat
advection were all zero. Temperature after the complete model step equalled
temperature immediately after the thermal solve.

Water content nevertheless changed through the native hydrology. Areal heat
capacity decreased from 190,986.531 to 185,744.797 J m⁻² K⁻¹, a 2.74% range.
The supplemental calculation therefore uses the native discrete equation,
`G_top - G_bottom = C_n * (T_n - T_(n-1)) / dt`, with each step's native
capacity before integrating to hourly means. This is not a check of total
soil-water enthalpy `Δ[C(T-T_ref)]`, or a formal fixed-capacity probe gate.

| Temperature supplied to the supplemental calculation | Heating mean absolute residual | Recovery mean absolute residual |
| --- | --- | --- |
| Native full-column output | 0.00158 W/m² | 0.00147 W/m² |
| Reported temperature frozen | 38.08 W/m² | 7.26 W/m² |
| Reported temperature departures halved | 19.04 W/m² | 3.63 W/m² |

The same numerical comparison allowance used by the probe is 1.90 W/m²
during heating and 1 W/m² during recovery. Both constructed reporting errors
exceed it in both phases; they are not observed defects in Noah-MP. Native
single-step residuals stay below 0.011 W/m². Top-layer temperatures range
from 294.50 to 303.33 K.

The script also replays the registered fixed-capacity criterion using the
initial capacity. It returns a numerical PASS for this case, but the report
marks it **inapplicable** because the capacity changed. Passing a tolerance
does not establish that a model meets the probe's physical assumptions.
The full run supplements the fixed-moisture component evidence without
changing the probe's equation, allowance or declared scope.

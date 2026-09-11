# Noah-MP soil heat component validation

This optional experiment checks `soil_heat_storage` against native Noah-MP
soil temperature calculations. It supplements the generated probe's reference
gate. It does not run the full Noah-MP/HRLDAS system or register a Noah-MP model
adapter in HydroTuring.

## Run it

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

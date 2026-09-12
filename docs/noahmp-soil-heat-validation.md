# Noah-MP soil heat validation

The complete HRLDAS/Noah-MP model was run on all five cases exported from the
current `energy/soil-heat-storage-consistency` generator. Its normal output
passed the unchanged `soil_heat_storage` criterion on every seed. Freezing
the reported temperature or halving its changes caused every seed to fail.

These optional runs validate the proposed probe with native model output.
They do not register Noah-MP as a HydroTuring benchmark model, provide a
model container, or add a dependency to CI or the mandatory reference gate.
The earlier component and moist-column experiments are retained below as
additional evidence.

## Current generated-case experiment

The model uses the official
[HRLDAS driver at `cd96df4`](https://github.com/NCAR/hrldas/tree/cd96df470220f7d7133cdbccd5f9c5355cf173e2)
and its
[Noah-MP v5.2.1 submodule at `17751dc`](https://github.com/NCAR/noahmp/tree/17751dcd7a2442a3c57138f530f523daf08dc1c5).
The full land column runs in its default single precision, including surface
exchange and water updates. The case is warm, bare and very dry, with no rain.
The model equations are unchanged.

`validate_generated.py export` calls the current harness's `build_case` and
`gate_seeds`. Each case has 24 hourly spinup intervals, 12 heating intervals,
and 60 recovery intervals. All models receive the same incoming shortwave and
longwave radiation and meteorology. The generator specifies layer thickness,
initial temperature, porosity, initial water content and material heat capacity.
The model computes its own surface exchange and lower-face heat flux.

This case specification revises the
[public proposal in issue #26](https://github.com/Flood-Lab/HydroTuring/issues/26):
incoming radiation replaces prescribed net radiation, and deeper soil may
evolve rather than follow a fixed temperature reservoir. The reference model's
exchange coefficients are no longer imposed on other models. These input and
boundary changes need maintainer review with the PR. The budget equation and
its 5% / 1 W m⁻² allowance have not changed.

### Input and output mapping

| Case quantity | Native HRLDAS/Noah-MP mapping |
| --- | --- |
| `rsds`, `rlds` [W m⁻²] | Incoming `SWDOWN`, `LWDOWN` |
| `tas` [°C] | `T2D` [K], after adding 273.15 |
| `sfcWind`, `huss`, `ps` | Wind along one horizontal axis, specific humidity and surface pressure |
| `pr` [mm day⁻¹] | `RAINRATE` [mm s⁻¹]; zero in these cases |
| Layer thickness | First `soil_thick_input` value and matching setup-file thickness and layer centre |
| Solid material heat capacity | Native table parameter `CSOIL_DATA` |
| Porosity and water content | STAS soil type 8, porosity 0.464; initial water content 0.005 |
| `hfg` [W m⁻²] | Native `HeatGroundTotMean`, downward into the first soil layer |
| `hfg_bottom` [W m⁻²] | Native lower-face Fourier flux, downward out of that layer |
| `tsoil_layer` [K] | First-layer temperature after the complete model step |

The native column has four layers. The first layer uses the case's sampled
thickness; deeper layer thicknesses are 0.3, 0.6 and 1.0 m. Every layer starts
at the case temperature and water content. Deeper soil follows the native
model, rather than the synthetic reference's prescribed reservoir.

Case timestamps label interval starts. The converter writes an initial forcing
record and one record at each interval end, carrying that interval's forcing.
The forcing, land-model, soil and output timesteps are all **3600 s**. There is
one native integration per case row, with no finer-step substitution or
resampling. Native `SWFORC` and `LWFORC` output confirms the pulse and all
96 interval values, to single-precision rounding. `_phase` is omitted from
the model's forcing CSV and restored from host-side metadata only for scoring.

Two read-only diagnostic outputs expose the native budget. The lower-face flux
is `2*k1*(T1-T2)/(-z2)` from the implicit temperature solve, using native
conductivity and the layer-centre gradient. It is not inferred from the
storage residual. The first layer's lower-face flux is distinct from the
bottom-of-column `HeatFromSoilBot`. Fluxes correspond to the native time
integration; temperature is the interval-end value.

### Prescribed capacity and its approximation error

The scorer uses the case's `soil_heat_capacity_areal`, calculated before the
model runs from layer depth and material properties. It never substitutes a
fitted or time-varying capacity. The mixture settings and formula are given in
the [probe description](../probes/energy/soil-heat-storage-consistency/README.md).

Native water updates still produce a very small capacity drift. Across all
five runs, the largest relative range was `1.294e-6` (0.0001294%).
The effect of replacing native capacity `C_n` with the prescribed `C_A`
was measured separately as:

```text
capacity approximation error = (C_A - C_n) * (T_end - T_start) / dt
```

Its maximum absolute value was `4.064e-5 W m⁻²`, far below the existing
1 W m⁻² floor. This is a quantified approximation in the dry test, not a
general allowance to ignore evolving moisture in wet soil. The largest native
ground latent heat flux was 0.002626 W m⁻². Snow, soil ice, vegetation,
penetrating radiation and precipitation heat advection were all zero, and
no later process changed the temperature returned by the thermal solve.
These checks and phase-specific capacity errors are included in the reports.

### Measured results

Each case completed 96 native steps and wrote 97 standard records, including
the initial state. For every seed, all **141 native output variables** were
exactly equal between the original executable and the executable with
diagnostic output. The runs used GNU Fortran 13.4.0, NetCDF Fortran 4.6.3 and
NetCDF C 4.10.1 on a Linux CPU.

| Gate seed | Largest normal phase mean absolute residual [W m⁻²] | Normal output | Frozen temperature | Halved temperature changes |
| --- | ---: | --- | --- | --- |
| 1239414052 | 0.000834 | PASS | FAIL | FAIL |
| 1746366166 | 0.000294 | PASS | FAIL | FAIL |
| 105834633 | 0.000623 | PASS | FAIL | FAIL |
| 612786747 | 0.000672 | PASS | FAIL | FAIL |
| 1119738861 | 0.000711 | PASS | FAIL | FAIL |

The controls alter only the reported temperature, retaining the native
fluxes. They are constructed reporting errors, not defects found in Noah-MP.
A failure in either phase fails a case. For seed 1239414052, the halved
temperature control fails during heating; its weaker recovery residual,
0.6973 W m⁻², remains inside the 1 W m⁻² floor. No tolerance was adjusted to
force every control phase to fail.

These results establish consistency between native fluxes and layer storage
under the specified test conditions. They do not establish temperature
accuracy or guarantee detection of every possible reporting error.

## Reproduce the generated-case runs

On Linux, provide `git`, `make`, GNU Fortran and NetCDF C/Fortran. Install
this checkout with `pip install -e '.[netcdf]'`. Use fresh scratch directories
for new exports and runs. From the repository root:

```bash
# Build the pinned native executables and official point-data converter.
# This also reproduces the additional moist-column experiment below.
python scripts/noahmp_soil_heat/run_full.py \
  --workdir /tmp/hrldas-soil-heat \
  --fc gfortran --netcdf-prefix /path/to/netcdf

# Export the current generator's five acceptance seeds.
python scripts/noahmp_soil_heat/validate_generated.py export \
  --cases /tmp/soil-heat-cases

# Run the original and instrumented full models on those exact inputs.
python scripts/noahmp_soil_heat/validate_generated.py run \
  --cases /tmp/soil-heat-cases \
  --workdir /tmp/soil-heat-native-runs \
  --native-case /tmp/hrldas-soil-heat/case

# Apply the current registered criterion and both temperature controls.
python scripts/noahmp_soil_heat/validate_generated.py score \
  --cases /tmp/soil-heat-cases \
  --workdir /tmp/soil-heat-native-runs
```

`run_full.py` downloads the pinned source, builds the official serial model
and text-to-NetCDF converter, then creates original and diagnostic executables.
The unrelated GRIB converter is not needed. Existing assets can be reused:
start with `export` and point `--native-case` at the existing case directory.

The generated-run directory contains per-seed native output, forcing,
`provenance.json`, `native_hydroturing_output.csv` and
`criterion_validation.json`. The combined criterion report is written at
the run directory's top level. Saved diagnostic CSVs and the exported inputs
can be scored elsewhere without compiling or downloading the model again.
The scorer exits with an error if a normal case does not pass or either
temperature control does not fail.

Downloaded model sources retain their
[UCAR license](https://github.com/NCAR/noahmp/blob/17751dcd7a2442a3c57138f530f523daf08dc1c5/LICENSE.txt)
and attribution; HydroTuring's license does not replace it. The Noah-MP
modeling system was developed at NCAR with university partners. NCAR is
sponsored by the United States National Science Foundation.

## Additional experiments

### Fixed-moisture thermal component

The earlier component experiment calls the official thermal-property,
diffusion and temperature-solver routines with a custom driver. It uses four
native layers, scores the top 0.1 m, and fixes moisture at 0.2, porosity at
0.45, quartz fraction at 0.4 and solid heat capacity at 2 MJ m⁻³ K⁻¹.
Its prescribed top flux is 80 W m⁻² for 12 hours, surrounded by 24 hours of
equilibrium and 60 hours of recovery. It does not run the full land column
or the current generated case.

```bash
python scripts/noahmp_soil_heat/validate.py \
  --workdir /tmp/noahmp-soil-heat --fc gfortran

# Re-score existing component CSVs without building again.
python scripts/noahmp_soil_heat/validate.py \
  --workdir /tmp/noahmp-soil-heat --skip-build
```

The double-precision component passed at both 300 s and 3600 s native steps,
with phase mean absolute residuals below `3.1e-12 W m⁻²`. Both temperature
controls failed in both phases. A deliberately substituted hour-start lower
flux passed at 300 s but failed during heating at 3600 s, illustrating that
the allowance does not detect every timing error. No NetCDF is needed for
this component experiment.

### Exploratory moist full column

The build command above also reproduces the original 96-hour full-column
experiment: initial moisture 0.20, air and soil temperature 295 K, specific
humidity 0.010, 300 W m⁻² incoming shortwave for 12 hours, and 300 s native
steps. It completed 1,152 steps; all 141 native output variables were identical
before and after instrumentation.

This moist run had a 2.74% heat-capacity range, so it was assessed separately
with the native stepwise equation `G_top-G_bottom=C_n*delta(T)/dt`. Normal
phase residuals were 0.00158 and 0.00147 W m⁻²; freezing and halving the
reported temperature changes caused both phases to fail the comparison.
This native thermal-equation check is not a total soil-water enthalpy budget,
nor the current generated-case fixed-capacity validation.

```bash
python scripts/noahmp_soil_heat/score_full.py \
  --workdir /tmp/hrldas-soil-heat/case/diagnostic
```

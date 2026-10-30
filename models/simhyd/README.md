# SIMHYD (MARRMoT m_18_simhyd_7p_3s) under HydroTuring

SIMHYD in the form the MARRMoT toolbox gives it, `m_18_simhyd_7p_3s` from
[wknoben/MARRMoT](https://github.com/wknoben/MARRMoT) v2.1.2 (`a95b925`). Proposed by Haoyang Su
as [Flood-Lab/HydroTuring#164](https://github.com/Flood-Lab/HydroTuring/issues/164).

**#164 was not labelled `accepted` when this was packaged. The maintainer asked for it anyway.**

MARRMoT's authors distinguish their formulations from the models that inspired them, so this is
not Chiew et al.'s original code, and a result here is about MARRMoT's formulation. It is also
not the differentiable SIMHYD decoder inside PD-HTHM (#120), a different model under a similar
name.

GPL-3.0; the image carries the licence at `/model/LICENSE`.

## What runs

MARRMoT's own MATLAB, unmodified, under GNU Octave. `simhyd_run.m` instantiates the class, sets
`delta_t`, `theta`, `S0`, `solver_opts` and `input_climate`, and calls `get_output()`, as
`workflow_example_1.m` does. `ht_adapter.py` is the units layer, standard library only. Nothing
is reimplemented.

Octave 6.4.0 with `optim` is the combination MARRMoT's README says it was tested with. Ubuntu
22.04 is the release that carries it, so the image is `ubuntu:22.04` pinned by digest, with
`octave=6.4.0-2` and `octave-optim=1.6.1-2build1` pinned by version (`octave-statistics 1.4.3-2`
and `octave-struct 1.0.17-3build1` came with them). The build checks the tag still points at the
commit, and asserts the interpreter version. The image is 219 MB.

## What is reported

| Column | Native | Note |
| --- | --- | --- |
| `pr` | the forcing | echoed from the CSV as text |
| `evspsbl` | `Ea` = `Ei` + `Et` | mm per step from MARRMoT, turned back into mm/day |
| `mrro` | `Q` = `Qt` = `SRUN` + `INT` + `BAS` | mm per step from MARRMoT, turned back into mm/day |
| `canopy` | `S1` | interception store, absolute mm |
| `mrso` | `S2` | soil store, absolute mm |
| `gw` | `S3` | groundwater store, absolute mm |

## What is not reported

- `snw`: the model has no snow store. All precipitation enters as rain and temperature is read
  by no equation. **Every closure probe requires `snw`, so the whole closure family is N/A
  (INCOMPLETE).** Writing a column of zeros would open them, and is what `flex_lumped` does, but
  a store the model lacks is not fabricated here.
- `channel`, `dis`: no routing. `Qt` leaves in the step it is generated.
- `gwex`: `S3` is closed inside the catchment, which answers the review's second question.
  `abstr` is not read, because the model has no withdrawal process.

These are named under `not_reported` in `run.json`.

## Parameters and initial stores

Fixed before any evaluation, as the proposal requires. No calibration exists for the synthetic
catchment. The proposal suggested calibrating on MARRMoT's example catchment and freezing the
result; that was not done, because a calibration to one real basin is no more applicable here
and would add a step the toolbox does not document as a default.

The values are the midpoints of the ranges the model file declares (`obj.parRanges`), which is
what the toolbox itself starts from: `MARRMoT_model.calibrate` falls back on
`par_ini = mean(obj.parRanges,2)`, and `workflow_example_4.m:105` sets the same. The two
capacities are taken from static.json, as `flex_lumped` and `sacsma_snow17` do.

| Parameter | Value | Range | Unit | Meaning | Source |
| --- | --- | --- | --- | --- | --- |
| `INSC` | `canopy_capacity_mm` | 0 to 5 | mm | interception capacity | **static.json**; midpoint 2.5 if absent |
| `COEFF` | 300 | 0 to 600 | mm/d | maximum infiltration loss | midpoint |
| `SQ` | 7.5 | 0 to 15 | - | infiltration loss exponent | midpoint |
| `SMSC` | `soil_capacity_mm` | 1 to 2000 | mm | soil moisture capacity | **static.json**; midpoint 1000.5 if absent |
| `SUB` | 0.5 | 0 to 1 | - | interflow constant | midpoint |
| `CRAK` | 0.5 | 0 to 1 | - | recharge constant | midpoint |
| `K` | 0.5 | 0 to 1 | 1/d | groundwater recession | midpoint |

The model file labels `COEFF` as mm; the equation uses it as a rate against the rain rate, so it
acts as mm/d.

Initial stores are zero: `input_s0 = zeros(numStores,1)` in `workflow_example_3.m:92` and
`workflow_example_4.m:57`. The second follows it with a 365-day warm-up (`warmup = 365`, line 117),
the spin-up most of the harness's daily probes supply.

Solver settings are `workflow_example_1.m`'s, which it notes are the defaults:
`resnorm_tolerance = 0.1`, `resnorm_maxiter = 6`. The solver's fifth and later restarts draw from
`randn`, so `rand` and `randn` are seeded from the request; repeat runs are bit-identical.

## The step

**The manifest declares PT1D**, as the proposal does. Every process is a rate in mm/d against
`delta_t` in days, so the adapter hands over the step as `delta_t` and rescales nothing. The
resolution probe also runs it at PT1H, and what that does is below.

## Evaluation

Gate seeds, all 37 probes: **9 PASS, 3 FAIL, 25 N/A** (24 INCOMPLETE, 1 INCOMPATIBLE).

Passing: `antecedent-monotonicity`, `area-invariance`, `causality`, `dry-down`, `extreme-rain`,
`phase-counterfactual`, `runoff-bounds`, `steady-state`, `warming-response`. Two of these say
little about this model: `phase-counterfactual` is trivially satisfied by a model that does not
distinguish snow from rain, and `area-invariance` by one with no area in it.

All three failures are in the model, in its equations at these parameters. None comes from the
adapter, and all three were checked on the harness's own cases against MARRMoT's internal fluxes.

### `energy/pet-consistency`: evaporation above demand

ET was 1.59, 1.56 and 1.53 of demand on the wettest fifth of soil states, against a limit of 1.0.

```matlab
flux_Ei = evap_1(S1,Ep,delta_t);          % min(S1/dt, Ep)
flux_Et = evap_2(10,S2,smsc,Ep,delta_t);  % min(10*S2/smsc, Ep, S2/dt)
```

Each is capped at `Ep` on its own, and neither is reduced by the other, so `Ea = Ei + Et` can
reach twice the demand. On the first gate seed ET exceeds PET on 1844 of 3650 days, on all of
them with both terms active, by up to exactly 2.0 times, and by 920 mm in total. Neither term
alone ever exceeds PET. No parameter tried brings the wet ratio under 1.4, so this is structural
in `m_18`. Whether Chiew et al.'s original reduces soil ET by interception was not checked.

### `mass/resolution-invariance`: runoff depends on the step

`mrro` moves by 20.3, 16.8 and 11.2% of the rain between PT1H and PT1D, against a limit of 10%.
Rain delivered in a shorter step arrives at a higher rate, and infiltration is capped by a rate:

```matlab
flux_INF = infiltration_1(coeff,sq,S2,smsc,flux_EXC);  % min(coeff*exp(-sq*S2/smsc), EXC)
```

On the first gate seed over the scored window, surface runoff `SRUN` is 142 mm hourly and 70 mm
daily, while total runoff is 154 and 111 mm: lower interflow and baseflow give part of it back.
Over the whole case, infiltration falls by 28% of the rain, and with it soil ET, interflow and
recharge.

Some of this is the physics the probe allows for, since within-day intensity is real. But `COEFF`
in SIMHYD is a daily parameter, and the proposal only claims PT1D, so the fair reading is that this
model does not transfer to a finer step at these parameters, not that its daily numbers are wrong.
It is sensitive to `SQ`: 11.5% at `SQ=1.5` and 6.9% at `SQ=0`, whole case, first seed. It is not
sensitive to the initial stores; see below.

### `mass/response-nonnegativity`: more rain, less flow later

Failed on one of three seeds: `mrro` falls 0.134 mm/day below the control, against an allowance of
0.12 mm/day, eight days after 120 mm was added. The sequence on seed 1713476937:

| | Control | With the added rain |
| --- | --- | --- |
| `S2` the day before the next storm, mm | 67.1 | 90.4 |
| Infiltration in that storm (45.7 mm), mm | 35.4 | 25.0 |
| Surface runoff that day, mm | 6.6 | 17.0 |
| Recharge that day, mm | 4.32 | 3.46 |
| `S3` at the end of that day, mm | 2.91 | 2.51 |
| Baseflow the day after, mm/day | 0.970 | **0.836** |

The wetter soil shuts infiltration through `exp(-sq*S2/smsc)`, which sends that day's water to
surface runoff (+9.2 mm/day) and starves the groundwater store, which at `K=0.5` is drained every
day, so the next day's baseflow is lower. No water is lost; it left earlier. This is structural,
but marginal and parameter-dependent: on this seed the dip disappears at `K=0.05`, at `SQ=1.5`, or
at the larger `SMSC`. The solver is not responsible: Newton-Raphson solved all 1095 steps of both
runs, with a worst residual norm of 7e-6.

## Water balance

Closed everywhere, to rounding. On the 11-year case: residual 1.7e-12 mm, and MARRMoT's own
`check_waterbalance` gives −5.7e-13 mm. Written to `run.json` on every run.

This holds by construction whatever the solver does: `MARRMoT_model.run` advances the stores from
the fluxes it evaluates at the solved state, so an unconverged step gets a wrong state, not a leak.
Over those 4015 steps Newton-Raphson solved 4003 and `fsolve` 12, `lsqnonlin` was never needed, and
the worst residual norm was 0.018. The stores stay non-negative; `S1`
reaches 2.09 mm over an `INSC` of 2 mm, because the interception threshold is a logistic smoother.
That would be scored by `state_bounds` on the closure probes if those could run.

## How much each choice moves the result

Over the 11-year record, spin-up included, for ET and runoff; on the worst nonnegativity seed for
the dip. The wet ratio is computed here the way the probe does, to within its windowing. "Step shift" is the
runoff moved between PT1H and PT1D, whole 40-day case, first seed, in % of the rain (18.2 here is
20.3 over the probe's 30-day window).

| Choice | Step shift | Wet ET / PET | Worst dip mm/day | Q / P |
| --- | --- | --- | --- | --- |
| **As packaged** | 18.2 | 1.60 | −0.134 | 0.377 |
| `SMSC`, `INSC` at the midpoints (1000.5, 2.5), not from the case | 21.6 | 1.40 | 0 | 0.229 |
| `SQ` = 1.5 | 11.5 | 1.45 | 0 | 0.307 |
| `SQ` = 0 | 6.9 | 1.45 | 0 | 0.307 |
| `COEFF` = 600 | 16.6 | 1.58 | −0.152 | 0.356 |
| `K` = 0.05 | 20.1 | 1.60 | 0 | 0.376 |
| `SUB` = 0.1, `CRAK` = 0.9 | 17.9 | 1.61 | −0.235 | 0.381 |

`SQ` 1.5 and 0 agree on the daily columns because a daily rain rate never reaches the capacity at
either, so `INF = EXC` throughout.

The initial stores matter only where the spin-up is short. Resolution probe, scored window, three
seeds, limit 10:

| S0 as a fraction of capacity | Seed 1 | Seed 2 | Seed 3 |
| --- | --- | --- | --- |
| **0 (packaged)** | 20.31 | 16.83 | 11.20 |
| 0.25 | 17.72 | 15.05 | 10.73 |
| 0.5 | 14.24 | 13.49 | 10.56 |
| 0.75 | 12.64 | 11.36 | 10.21 |

The verdict holds at every one. The 0 row reproduces the harness's numbers to the digits shown.
With a year of spin-up, half-full initial stores leave the wet ratio and the dip unchanged at the
digits shown, and move Q / P over the whole record from 0.377 to 0.391.

## Caveats

- Mid-range parameters are not a calibration. They are a documented, reproducible point in the
  model's space, and `SQ = 7.5` in particular is a strong response of infiltration to soil moisture.
  Two of the three failures move with parameters, and are reported as what this point does. The
  third, ET above demand, does not.
- Wall time including Octave's start: 0.27 s for 40 rows, 1.6 s for 960, 11.4 s for 4015.

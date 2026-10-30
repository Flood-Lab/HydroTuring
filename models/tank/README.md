# Tank Hydrologic Model

The four-tank conceptual model of Sugawara and Funiyuki (1956), in the equations of Paik et
al. (2005, doi:10.1002/hyp.5984), as implemented by Nazmul Ahasan in
[nzahasan/tank-model](https://github.com/nzahasan/tank-model) at tag `3.0.0b3`, commit
`aea105ca32d0957fee76fbcfe3a5b8993938c9a1`, MIT. Submission: issue #124.

The image vendors `tank_core/` from that commit, unmodified, and runs it through
`computation_helpers.compute_project` on one Subbasin. Nothing is reimplemented. Upstream's
CLI, plotting and calibration dependencies are not installed; `tank_core` needs numpy 2.5.1,
pandas 3.0.3 and scipy 1.18.0, the versions upstream's `uv.lock` resolves.

```
S[t+1, top] = max(S[t, top] + P[t+1] - E[t+1] - side[t, top] - bottom[t, top], 0)
S[t+1, i]   = max(S[t, i] + bottom[t, i-1] - side[t, i] - bottom[t, i], 0)
side[t, top] = S1 * max(S - H1, 0) + S2 * max(S - H2, 0);   side[t, i] = Si * max(S - Hi, 0)
bottom[t, i] = Bi * S
```

## Parameters

No calibration exists on a synthetic catchment. Every value is the midpoint of the range in
the "Sugawara" column of Table 5-3-2, which upstream ships as `docs/ParameterRange.pdf`. They
are per day, as Paik et al. state, and as `t0_*`, `t1_*`, `t2_*`, `t3_*` are named upstream.

| Parameter | Upstream name | Value | Range | Rescaled to the step |
| --- | --- | --- | --- | --- |
| S11, upper side outlet, top tank | `t0_soc_uo` | 0.3 /day | 0.1-0.5 | yes |
| S12, lower side outlet, top tank | `t0_soc_lo` | 0.3 /day | 0.1-0.5 | yes |
| B11, bottom outlet, top tank | `t0_boc` | 0.3 /day | 0.1-0.5 | yes |
| S21, second tank | `t1_soc` | 0.065 /day | 0.03-0.1 | yes |
| B21, second tank | `t1_boc` | 0.055 /day | 0.01-0.1 | yes |
| S31, third tank | `t2_soc` | 0.0075 /day | 0.005-0.01 | yes |
| B31, third tank | `t2_boc` | 0.0075 /day | 0.005-0.01 | yes |
| S41, fourth tank | `t3_soc` | 0.00275 /day | 0.0005-0.005 | yes |
| H11, upper outlet height | `t0_soh_uo` | 45 mm | 30-60 | no |
| H12, lower outlet height | `t0_soh_lo` | 20 mm | 10-30 | no |
| H21 | `t1_soh` | 25 mm | 0-50 | no |
| H31 | `t2_soh` | 15 mm | 0-30 | no |
| IST1, IST2, initial top two tanks | `t0_is`, `t1_is` | 25 mm | 0-50, "This study" column: the Sugawara column has none | no |
| IST3 | `t2_is` | 55 mm | 10-100 | no |
| IST4 | `t3_is` | 300 mm | 100-500 | no |

B41 is not a parameter: upstream's fourth tank has no bottom outlet, which is the conservative
end of the table's 0-0.005. These midpoints give the time constants Sugawara's chapter asks for
in the order it asks for them (upstream's `docs/207090972-Tank-Model.pdf`, section 6.1): a day
or so, about ten days, months, and about a year.

`soil_capacity_mm` is not used. The model has no capacity: a tank has outlets and no lid, so
nothing in it can be told the catchment's capacity. `area_km2` is upstream's `area` and enters
only the conversion to `dis`.

## Choices, and what each moves

The twelve scored probes were rerun with each choice changed, on the gate seeds, with no
container.

| Choice | Alternative | What moves |
| --- | --- | --- |
| Sugawara midpoints | Midpoints of `global_config.py`'s optimiser bounds, as upstream's `tests/tank_core_test.py` builds them | Nothing: all twelve pass |
| Sugawara midpoints | Low end of every range | `mass/steady-state` FAIL: `gw` still varies 1.3% of its level in the last year (limit 1%) and is not settling, its time constant being 5.5 years |
| Sugawara midpoints | High end of every range | `mass/extreme-rain` FAIL: 1.50 of the added rain returns. See the finding below |
| Coefficients rescaled, `1-(1-k)^dt` | Passed unscaled at PT1H | `mass/resolution-invariance` FAIL: `mrro` moves 99% of the rain instead of 3% |
| `mrso` = top + second tank, `gw` = third + fourth | `mrso` = top tank alone | Nothing: all twelve pass |
| States after the row's outflows | `tank_storage[t]` as upstream holds it, before them | Nothing scored; it would shift a budget by one row's runoff |
| `evspsbl` = demand less the reset | `evspsbl` = the demand | `energy/pet-consistency` FAIL, and the budget is open by 55% to 66% of rain |

**`evspsbl` is derived.** Upstream returns no evaporation. It subtracts the demand from the top
tank and resets a negative store to zero, so what actually left is the demand less what the
reset gave back, and `ht_adapter.py` computes that from upstream's own arrays. The reset fires
on about 2300 of the 4015 days of a ten-year record and takes back 5200 to 5700 mm of about
8300 mm of demand, so evaporation is 0.32 to 0.37 of potential over the record. All of it comes
from the top tank, which drains at 0.3 per day and more: rain moves into tanks evaporation
cannot reach. That low ratio is the model on these parameters, and calibration or the use of
actual evaporation as `E`, as Paik et al. do, would change it.

Row 0's `pr` and `pet` never enter the model, because upstream sets `S[0]` to the initial
storage. Row 0 is in every probe's spinup.

## Reported, and not reported

| Column | From |
| --- | --- |
| `pr` | echoed |
| `evspsbl` | demand less the top tank's reset, above |
| `mrro` | the sum of the four tanks' side outlets, as upstream's discharge is |
| `dis` | upstream's discharge, m3/s, on `area_km2` |
| `mrso` | top + second tank: surface and intermediate runoff stores |
| `gw` | third + fourth tank: sub-base and base-flow stores |

Not reported, named in `run.json`, and the reason 25 probes are N/A:

- `snw`: no snow store; precipitation enters the top tank as liquid whatever the temperature.
- `canopy`: no interception store.
- `channel`: one Subbasin has no reach. Upstream's Muskingum routing is optional and has no
  default `k` and `x`. It is not configured, because its `k` would have to be invented, and
  because a reach's travel time is far below a day: at PT1D `del_t/k + x > 1` trips upstream's
  own `WARNING-MUSKINGUM-01`.
- `stage`, `gwex`, and every energy term: the model has none. Withdrawals (`abstr`) are not read.

Nothing here reports a column of zeros to reach a probe. Every budget probe requires `snw` and
`canopy`, so `mass/catchment-closure` and its kin are N/A (INCOMPLETE).

## Steps

`[PT1D, PT1H]`. Upstream's step is a free parameter and enters only the conversion to m3/s, so
at any step its coefficients are fractions per row. PT1D is what the published values mean;
PT1H is reached by the rescaling above. Rows are never resampled.

## Results

Evaluated under Apptainer on gate seeds: 12 PASS, 0 FAIL, 25 N/A (23 INCOMPLETE, 2 INCOMPATIBLE).

| Probe | Result | Measured |
| --- | --- | --- |
| `energy/pet-consistency` | PASS | 1.00 of potential when wet, 0.08 when dry, 0.35 over the record |
| `mass/antecedent-monotonicity` | PASS | 10.0 mm more from the 60 mm storm, 0.17 of it |
| `mass/area-invariance` | PASS | identical |
| `mass/causality` | PASS | nothing moves before the storm; `mrro` returns 0.829 of the added rain |
| `mass/dry-down` | PASS | 282 mm over 730 dry days against a 625 mm bound |
| `mass/extreme-rain` | PASS | returns 0.99 of 997 mm added |
| `mass/phase-counterfactual` | PASS | 0.0%; the model never reads temperature |
| `mass/resolution-invariance` | PASS | worst 3.16% of rain, `mrso_end`, PT1D against PT1H |
| `mass/response-nonnegativity` | PASS | largest dip -2.1e-4 of the added 120 mm |
| `mass/runoff-bounds` | PASS | 0.71 of rain, inside [0.09, 1.09] |
| `mass/steady-state` | PASS | `gw` still settling, 36% as far in the last year as the one before |
| `mass/warming-response` | PASS | per unit of demand, warmer: `mrro` -0.219, `evspsbl` +0.229 |

`momentum/routing-lag-consistency` and `momentum/wave-celerity-bounds` are INCOMPATIBLE: no
channel lengths or `q_in` are consumed, because there is no reach.

**What was not tested.** Twelve passes do not say the model closes its budget, because no
scored probe differences its stores against its fluxes. Run by hand on the closure probe's five
gate seeds, the budget closes to about 1e-16 of rain with a largest row residual of about
1e-13 mm: `Q/P` 0.68 to 0.71, `E/PET` 0.32 to 0.37, `mrso` 0 to 113 mm, `gw` 204 to 423 mm.
That closure is close to an identity, because `evspsbl` is computed from the same store update,
and it would earn the closure criterion's `suspicious_exact` flag. What carries information is
in `run.json`: the reset never exceeded the demand and never fired in a lower tank (both totals
zero to 1e-12 mm), so the derived evaporation is never negative and no water appeared from
nowhere.

## Finding: the scheme can make water

Upstream's outlets are explicit and nothing limits a tank's outflow to its contents. It prints
`WARNING-TANK-02` when outflow exceeds storage, then carries on, and the reset to zero supplies
the difference. At these parameters the top tank's coefficients sum to 0.9 and it cannot happen
(`clamp_beyond_demand_mm` and `lower_tank_clamp_mm` in `run.json` are zero to rounding). At the
high end of the same published ranges they sum to 1.5, and `mass/extreme-rain` measures 1.50 of
added rain returned as runoff. Upstream's optimiser bounds, in `global_config.py`, permit that
sum, so a calibration can land there. This is about the scheme, and it is not in the verdict
above because the parameters here do not reach it.

## Reproducing

```bash
sbatch --chdir="$PWD" --job-name=ht-build-tank scripts/hpc/build_model.sbatch tank
sbatch --chdir="$PWD" --job-name=ht-tank scripts/hpc/run_model.sbatch tank
```

A ten-year case runs in about a second; `window_days: full`, so every probe sees its whole
record. The image is about 130 MB.

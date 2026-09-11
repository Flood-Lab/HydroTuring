# mass/precipitation-counterfactual

A model that reports one budget term as the residual of the others closes the water budget on every seed. This probe asks what construction cannot answer: when the same seed gets more or less rain, where did the difference go? It implements [accepted proposal #18](https://github.com/Flood-Lab/HydroTuring/issues/18). Every number below was measured with `ht run` and `ht gate` on the gate seeds (872466880, 1379418994, 1886371108) at upstream `273a5b2`.

## What it asserts

Over the scored window `W`, each perturbed variant `p` is compared with the control `c`, with `delta X = X_p - X_c`:

```text
P_v  = sum_W pr_v(t) dt          E_v = sum_W evspsbl_v(t) dt
Q_v  = sum_W mrro_v(t) dt        G_v = sum_W gwex_v(t) dt      (zero when not reported)
dS_v = S_v(end) - S_v(start-)    S   = every reported store
R_v  = P_v + G_v - E_v - Q_v - dS_v                                       [mm]

f_E = delta E / delta P     f_Q = delta Q / delta P
f_S = delta dS / delta P    f_X = -delta G / delta P    (signed export of the increment)

R_delta / delta P = 1 - (f_E + f_Q + f_S + f_X)                            [-]
```

`S_v(start-)` is the storage one row before the first scored row. The denominator is the precipitation added or removed, not the control total. For `drier20` it is negative, and each share is the fraction of the removed water that term gave up. These are incremental fractions, not elasticities; the measured elasticities (1.4 to 2.4) exceed one without violating anything.

| Criterion | Variants | Asserts |
| --- | --- | --- |
| `closure` | control | `abs(R_c) / P_c <= 0.05`; `gwex` counted as a declared source |
| `counterfactual_response` | control vs. each of wetter20, wetter10, drier20 | `f_E`, `f_Q` each `>= 0.03`; every share incl. `f_S` `<= 0.95`; `abs(f_E + f_Q + f_S - 1) <= 0.10`; any pair failing fails |
| `monotone_response` | all four, sorted by total `pr` | runoff does not fall between adjacent rungs; no rung adds more runoff than rain (1.05, slack 0.01); at least 0.1 of the rain added across the ladder runs off |
| `forcing_fidelity` | control | echoed `pr` equals the given `pr`, `rtol 1e-6` |
| `state_bounds` | control | `mrso` in `[0, soil_capacity_mm]`, `canopy` in `[0, canopy_capacity_mm]`, `snw >= 0` |

`counterfactual_response` leaves out `gwex`, so for `sacsma_snow17` the three-term sum is `1 - f_X`, with `f_X` measured at +0.0042 to +0.0049.

| Setting | Value | Origin |
| --- | --- | --- |
| `closure` threshold | 0.05 of `sum_pr` | suite-wide rule; loose against model numerics (baselines close to 1.7e-16; time-stepping errors are the subject of [R11]), tight against observed closure, where residuals of 5% to 25% are typical [R9, R10] |
| `min_share` / `sum_tolerance` | 0.03 / 0.10 | template defaults; no paper sets them, margins below |
| `max_share` | 0.95 | template default 0.90, raised here; see below |
| `monotone_response` 0.1 / 1.05 / 0.01 | as in `mass/extreme-rain` | criterion defaults; `0 <= dQ/dP <= 1` follows from Budyko [R6] and the budget |
| `forcing_fidelity` `rtol`, `state_bounds` capacities | 1e-6, `static.json` | numerical, physical |
| `min_window_days` | 3650 | conservative; 1825 days also passes every baseline, window table below |

**Margins.** At +20% on the worst seed, evaporation reaches 0.140 against 0.03, runoff 0.839 against 0.95, and the sum 0.994 against 0.10. Seeds differ by at most 0.02 in any share.

**Why 0.95.** The template's 0.90 leaves `sacsma_snow17` 0.061 of margin. In a wetter catchment near saturation the Budyko limit of the runoff share is one, so 0.90 would fail exact physics there. Here 0.95 changes no verdict.

**A physical bracket.** No paper sets the share bounds, so Budyko supplies one. With Fu's curve [R7] and the sensitivities of [R6], at PET/P 0.876 to 0.960 and curve parameters 2.0 to 3.0, the expected runoff share is 0.65 to 0.75. The expected evaporation share is 0.25 to 0.35; the baselines measure 0.59 to 0.84 and 0.14 to 0.37. The bounds therefore belong to this catchment. Measured elasticities lie inside the observed ranges for the United States, 1.0 to 2.5 [R8], and Australia, 2.0 to 3.5 [R12].

**The denominator cannot approach zero.** At +20% the added precipitation is 1589, 1634 and 1737 mm. Rounding to 1e-6 mm over 3650 rows leaves an arithmetic error of order 1e-2 mm, so no absolute floor is needed.

## The case

Each seed draws 4,015 daily rows: 365 of spinup and 3,650 scored. All random draws happen before the variant is applied. A variant multiplies the scored `pr` by 1.20, 1.10 or 0.80 and leaves `tas`, `pet`, wet-day occurrence, event order and the spinup untouched. Four variants is the schema's `maxItems`.

The control carries 794 to 869 mm a year on 91 to 95 wet days. Mean wet-day depth is 8.7 to 9.1 mm, and the largest daily total 70 to 97 mm. PET/P of 0.876 to 0.960 places it near the middle of the Budyko curve. A third of the precipitation falls below 0 °C, so snow must be reported in `snw`; `flex_lumped` and `flex_topo` have no snow module and treat it as rain. Temporal concentration is held fixed because it changes storage on its own [R3].

Shares (evaporation, runoff, storage) on the seed the harness reports as worst. Every cell passes the bounds:

| Amplitude | reference_bucket | flex_lumped | flex_topo | sacsma_snow17 |
| --- | --- | --- | --- | --- |
| −20% | 0.372, 0.593, 0.034 | 0.240, 0.745, 0.014 | 0.263, 0.727, 0.010 | 0.236, 0.725, 0.033 |
| −10% | 0.330, 0.660, 0.009 | 0.219, 0.767, 0.014 | 0.248, 0.742, 0.010 | 0.225, 0.737, 0.033 |
| −5% | 0.291, 0.674, 0.034 | 0.210, 0.776, 0.014 | 0.242, 0.744, 0.015 | 0.215, 0.748, 0.033 |
| +5% | 0.261, 0.712, 0.028 | 0.192, 0.786, 0.022 | 0.229, 0.757, 0.014 | 0.203, 0.760, 0.032 |
| +10% | 0.243, 0.730, 0.028 | 0.186, 0.793, 0.022 | 0.210, 0.775, 0.016 | 0.159, 0.795, 0.041 |
| +20% | 0.219, 0.754, 0.028 | 0.176, 0.811, 0.013 | 0.206, 0.784, 0.010 | 0.151, 0.803, 0.042 |
| +50% | 0.196, 0.777, 0.028 | 0.164, 0.820, 0.017 | 0.177, 0.814, 0.009 | 0.128, 0.827, 0.042 |

The runoff share rises smoothly with amplitude and model differences exceed seed differences. The runoff coefficient of `reference_bucket` goes from 0.41–0.43 to 0.47–0.49 at +20%, and `flex_topo` stays below 0.57, so ±20% keeps every baseline in its regime. The same ±20% grid in 10-point steps is used by [R1, section 3.2] and [R2, section 2.3.2].

The response is monotone but not symmetric: `reference_bucket` returns 0.593 of the removed water as runoff at −20% and 0.754 of the added water at +20%. [R5] find the same asymmetry on 353 observed catchments, so the probe asks only for monotonicity. Across the whole ladder the baselines return 0.69 to 0.82 of the added rain as runoff.

**The scoring window is the whole record.** A submitted daily model is scored by default on the largest 30-day flood event. A short window counts only the rain added inside it, while its runoff carries increment stored earlier. Verdicts with the shipped settings (`max_share` 0.95, all three pairs, gate seeds), naming the share that breaks:

| Window | reference_bucket | flex_lumped | flex_topo | sacsma_snow17 |
| --- | --- | --- | --- | --- |
| 30 days (default) | FAIL, runoff 2.870, storage −1.870 | FAIL, runoff 1.308 | FAIL, evaporation 0.013 | FAIL, runoff 1.053 |
| 365 days | FAIL, `drier20` runoff 0.961 | PASS | PASS | PASS |
| 730 days | PASS | PASS | PASS | FAIL, `wetter20` runoff 0.996 |
| 1825 days | PASS | PASS | PASS | PASS |
| 3650 days (full) | PASS | PASS | PASS | PASS |

The window follows the largest flood, so the sequence is not monotone in length: 365 days fails `reference_bucket` and 730 days fails `sacsma_snow17`. Sensitivity depends on aggregation time through storage [R4]. 1825 days passes all four baselines, so the full record is the conservative choice rather than the shortest window that works. `mass/phase-counterfactual` uses the same mechanism at 1095 days.

**Spinup.** One year is enough. Prepending 365 or 1460 days of independent weather moves the shares by at most 0.001, or 0.009 for `flex_topo`. Lengthening `SPINUP_DAYS` instead changes the scored weather and moves them by up to 0.05.

## Why these criteria

Closure alone is trivially satisfiable. Each other criterion closes off a specific cheat, measured here at +20%:

| Criterion | The cheat it closes off | Caught here |
| --- | --- | --- |
| `closure` | a silent sink | `reference_leaky`, residual 15.0% |
| `counterfactual_response` | ignoring the perturbation, giving it all to one term, or inventing water between two runs that each close | `reference_cheater` (E, Q, S = 0.000, 0.350, 0.650), `reference_degenerate` (1.000, 0.000, 0.000), `reference_leaky` (0.285, 0.542, 0.024; sum 0.850) |
| `monotone_response` | runoff that falls as rain rises, or exceeds it | `reference_degenerate`, 0.00 returned across the ladder |
| `forcing_fidelity` | rescaling the input to choose the denominator | |
| `state_bounds` | storage invented as the residual | `reference_cheater`, `mrso` reaches 1870 mm |

`reference_cheater` (`ET = 0.55 PET`, `Q = 0.35 P`, storage as the residual) closes exactly on every seed. `state_bounds` catches it only once its storage drifts out of range; `counterfactual_response` catches it regardless, because unchanged PET gives no evaporation response. All three broken models are pinned to `counterfactual_response` in `must_fail`.

**Limits.**

- `reference_in_sample` passes at every amplitude, with shares 0.219, 0.741, 0.028 and a sum of 0.977 to 0.987. Its loss needs daily rain above 55 mm, and +20% raises such days from 5–7 to 10–14 per record. The loss is 1.3% to 2.3% of the added water; tightening `sum_tolerance` to catch it would leave `sacsma_snow17` without margin.
- A fixed-ratio reporter whose ratios fall inside the bounds passes without dynamics.
- Equal and opposite errors inside one run cancel in a cumulative budget.
- The bounds are measured for this climate (PET/P 0.88–0.96, a third of precipitation as snow) and must be re-measured for another.

## Reproducing a failure

Every run records its seed. To re-run one:

```bash
ht run --model <name> --probe mass/precipitation-counterfactual --seed <n>
```

`ht gate --probe mass/precipitation-counterfactual` runs twelve ten-year daily runs per model, about 2 s per model with the subprocess runner. `pytest -q tests/test_precipitation_counterfactual.py tests/test_counterfactual_response.py` covers the generator and the criterion.

## References

- R1. Rasouli, Pomeroy & Whitfield (2022), *J. Hydrol.* 606, 127460. https://doi.org/10.1016/j.jhydrol.2022.127460
- R2. Li et al. (2024), *Hydrol. Earth Syst. Sci.* 28, 4521–4538. https://doi.org/10.5194/hess-28-4521-2024
- R3. Lesk & Mankin (2026), *Nature* 653, 425–432. https://doi.org/10.1038/s41586-026-10487-7
- R4. Zhang, Viglione & Blöschl (2022), *Water Resour. Res.* 58, e2021WR030601. https://doi.org/10.1029/2021WR030601
- R5. Tang et al. (2019), *J. Geophys. Res. Atmos.* 124, 11932–11943. https://doi.org/10.1029/2018JD030129
- R6. Roderick & Farquhar (2011), *Water Resour. Res.* 47, W00G07. https://doi.org/10.1029/2010WR009826
- R7. Zhang et al. (2004), *Water Resour. Res.* 40, W02502. https://doi.org/10.1029/2003WR002710
- R8. Sankarasubramanian, Vogel & Limbrunner (2001), *Water Resour. Res.* 37, 1771–1781. https://doi.org/10.1029/2000WR900330
- R9. Sahoo et al. (2011), *Remote Sens. Environ.* 115, 1850–1865. https://doi.org/10.1016/j.rse.2011.03.009
- R10. Lorenz et al. (2014), *J. Hydrometeorol.* 15, 2111–2139. https://doi.org/10.1175/JHM-D-13-0157.1
- R11. Kavetski & Clark (2010), *Water Resour. Res.* 46, W10511. https://doi.org/10.1029/2009WR008896
- R12. Chiew (2006), *Hydrol. Sci. J.* 51, 613–625. https://doi.org/10.1623/hysj.51.4.613

# Water balance integrity across declared catchment attributes

Shunan Zhou, Dalian University of Technology, Dalian, China. [ORCID](https://orcid.org/0009-0006-4072-6799), GitHub: Davidanan.

## Question and scope

Hold a model, its learned weights and its non-prescribed parameters fixed. Without recalibration against target-catchment discharge observations, does its reported budget remain closed across the catchment attributes supplied to it?

This is a synthetic daily rainfall–runoff experiment. The **experimental reference domain** is soil capacity 200–450 mm and canopy capacity 1–2 mm. Domain membership is defined by these intervals. It is not inferred from a model's training data. Capacity is the size of a store, not the quantity of water currently in it: capacities stay fixed within each independently initialised simulation.

## Native seed-per-catchment design

One seed generates one catchment and one weather record. All five native criteria see that run, exactly as on `mass/catchment-closure`. There are no paired variants, scoring wrappers, schema changes or default seed-stratification changes. Twelve fixed gate seeds cover all eight attribute categories below. Default evaluations use twelve fresh framework seeds; the categories and climates actually sampled are reported, without a guarantee that every category occurs in every random evaluation.

| `seed % 8` | Category | Soil capacity, mm | Canopy capacity, mm |
|---|---|---:|---:|
| 0 | reference | 320 | 2 |
| 1 | reference_medium | 230–280 | 2 |
| 2 | soil_small | 80–140 | 2 |
| 3 | soil_large | 500–700 | 2 |
| 4 | canopy_small | 320 | 0.5–0.8 |
| 5 | canopy_large | 320 | 2.2–2.5 |
| 6 | joint_small | 80–140 | 0.5–0.8 |
| 7 | joint_large | 500–700 | 2.2–2.5 |

Within a band, capacity is uniform. The experiment samples selected single-attribute changes and two joint corners, not every combination or a measured population frequency. It imposes no monotonicity or equality between catchments' hydrographs. Attribute and weather random streams are separate; the model receives capacities and forcing, not the generator seed or category label.

## Generated weather and duration

PT1D: 730 warm-up days followed by 3,650 evaluated days. Evaluation years are consecutive 365-day blocks, not calendar years. `min_window_days: 3650` prevents a submitted model's event-window setting from shortening this experiment. Two-year warm-up is not a claim of equilibrium for every model; actual storage change is always retained.

| `seed % 5` | Weather regime | Wet intercept / seasonal amplitude | Mean T / amplitude, °C | Mean PET / amplitude, mm/day | Latitude |
|---|---|---|---|---|---|
| 0 | temperate winter-wet | 0.28 / 0.10 | 15 / 9 | 2.2 / 1.3 | 38° |
| 1 | warm humid | 0.40 / 0.03 | 24 / 3 | 3.0 / 0.5 | 20° |
| 2 | summer rainfall | 0.30 / −0.20 | 23 / 5 | 3.0 / 1.0 | 25° |
| 3 | seasonally water-limited | 0.18 / 0.12 | 19 / 8 | 2.6 / 1.3 | 35° |
| 4 | cool maritime rain | 0.32 / 0.06 | 10 / 6 | 1.6 / 0.8 | 48° |

Wet-day probability is the intercept plus seasonal amplitude times cos(phase), plus 0.22 if the preceding day was wet. Wet-day amounts are drawn directly as `35 × Beta(0.9,3.2)` mm, without clipping or seed rejection. T uses uniform noise ±3°C and PET a factor uniform in [0.85,1.15]. Resulting bounds are P∈[0,35] mm/day, T∈[1,31]°C and PET∈[0.68,4.6] mm/day. The supplied conditions are rain-only; sustained ice, severe aridity, stronger storms and managed water systems are not represented by this generator. No model is exempted from any generated case on the basis of its familiar climate.

Area is 250 km². Bucket auxiliary defaults are baseflow coefficient 0.006/day, snow threshold 0°C and melt factor 3.2 mm/(°C day). These are declared attributes, not parameters fitted to an observed hydrograph.

The rainfall bound also stays inside a sufficient stability condition for the unchanged explicit FLEX soil partition: with `x=S/C` and beta=1.85, `1−x^beta ≤ beta(1−x)` for x∈[0,1]. Thus `p≤C/beta` prevents a capacity overshoot before subsequent drainage. At C=80 mm, the sufficient bound is 43.24 mm, above the generated 35 mm. This is a declared generator choice, not a repair to FLEX or a promise about stronger storms.

## Budget and native criteria

For the complete evaluation window W:

```text
R = sum_W[(P + G - ET - Q) * dt] - (S_end - S_before_W)
e = abs(R) / sum_W[P * dt]
```

P comes from supplied forcing; G is the signed external `gwex` if reported. ET=`evspsbl`, Q=`mrro`, dt=1 day. S contains required soil/snow/canopy and every reported groundwater/channel store, using the same catchment-area depth. Internal exchanges are not external inputs. A model must report every store it has, with its physical mapping documented in its adapter.

| Native criterion | Requirement |
|---|---|
| `closure` | Full-window e≤0.05; no pooling across seeds |
| `state_bounds` | Soil/canopy within supplied capacities; snow, groundwater and channel nonnegative where reported; native 1e-6 mm tolerance |
| `et_plausible` | ET nonnegative; full-window ET/PET≤1 with native numerical tolerances |
| `non_degenerate` | Runoff ratio against P+G in [0.02,0.98]; Q/ET CV≥0.1; seven-day rainfall/runoff correlation ≥0.05 |
| `forcing_fidelity` | Native maximum echoed-P deviation divided by mean supplied P≤1e-6 |

State bounds are checked over the native scored window; the budget uses the preceding measured inventory. The generated rainfall denominator is positive; the native criterion rejects a zero denominator. There is no custom threshold margin, clipping of model outputs or repair of the residual.

The native seven-day rainfall/runoff correlation screen is enabled (`min_response: 0.05`), consistent with catchment-closure and the rain-only guidance in the spatial-extrapolation template. It is a scoped anti-degeneracy screen, not a consequence of the mass-conservation identity. A conservative, sufficiently delayed response can fail this screen; that limitation does not justify disabling it only for this probe. Thresholds are fixed before evaluating models.

## Annual diagnostics and interpretation

`scripts/diagnose_ungauged_closure.py` reports ten annual residuals, each using the actual inventory immediately before that year. It regenerates the recorded seed, checks input identity and validates the stored output through the native protocol. It never modifies the standard verdict or archive.

Annual errors can cancel over the full record; the separate annual diagnostics expose this without changing the official full-window verdict.

## Positive and negative evidence

The positive references are the exact `reference_bucket` and three physical models: `flex_lumped`, `flex_topo` and `sacsma_snow17`. They run with unchanged numerical equations. Fitted models are evaluated with their native parameters; the probe does not require capacity-input consumption or override learned parameters. Capacities remain prescribed evaluation bounds, as in catchment-closure. Passing establishes the tested closure and bounds, not correct capacity use or catchment adaptation.

Seven attribute-dependent controls share the bucket implementation and are inactive inside the experimental reference domain. They read no case IDs, generator seeds, future rows or evaluation boundaries. The table lists the target check and additional failures observed across the fixed gate seeds.

| Fixture suffix | Out-of-domain fault | Target criterion | Additional failing criteria |
|---|---|---|---|
| loss | Export ET and Q at 0.8 of their values | closure | None |
| gain | Add 0.1P to Q without a debit | closure | None |
| capacity | Use fixed internal soil/canopy capacities of 320/2 mm | state_bounds | None |
| forcing | Simulate and echo 0.8P | forcing_fidelity | closure |
| et | ET=1.1PET with compensating Q | et_plausible | non_degenerate |
| negative | Reverse Q and compensate in ET | non_degenerate | et_plausible |
| frozen | Zero Q/ET with constant states | non_degenerate | closure |

The capacity control changes native dynamics, not reported storage. It preserves closure and is caught on the `soil_small` and `joint_small` gate cases. The `canopy_small` gate case is not caught because daily evaporation empties the fixed 2 mm canopy before the reported day-end state. Large-capacity cases test closure, not capacity consumption: an upper bound cannot detect an internal store that remains below it. Tests cover all eight attribute categories and five climates, verify the fixed-capacity dynamics, and check soil and canopy bounds separately: soil fails on every `soil_small`/`joint_small` seed in 0–39, and canopy fails on seeds 4, 20 and 28. The other six controls must fail their target check on every out-of-domain test case. Generic `reference_cheater` and `reference_degenerate` controls additionally test state bounds and non-degeneracy. All seven attribute-dependent controls pass the ordinary catchment-closure cross-check, whose capacities lie inside the reference domain.

## Reproduction

```bash
ht validate
ht gate --probe mass/ungauged-basin-closure
ht run --model flex_lumped --probe mass/ungauged-basin-closure --gate-seeds --workdir /tmp/ungauged --json /tmp/flex.json
python scripts/diagnose_ungauged_closure.py /tmp/ungauged/flex_lumped__mass__ungauged-basin-closure__1679270760 --seed 1679270760 --output /tmp/annual.json
pytest -q
```

Scientific-model demonstrations use the original adapters, serial execution and the standard per-invocation 60-second limit. The probe acceptance gate's trusted references run without Docker. Model runtime and probe overhead are distinct quantities.

## Sources and parameter semantics

- [HydroTuring spatial-extrapolation template](https://github.com/Flood-Lab/HydroTuring/blob/main/templates/probe.extrapolation-space.template.yaml): seed-per-catchment design, 730-day warm-up and broad capacity sampling.
- [van Oorschot et al. (2024)](https://doi.org/10.5194/hess-28-2313-2024): root-zone storage estimates reaching hundreds of millimetres, supporting magnitude rather than identity with every model's soil storage datum.
- [Zhang et al. (2006)](https://doi.org/10.5194/hess-10-65-2006): canopy interception storage magnitudes.
- [Zhong et al. (2022)](https://doi.org/10.5194/hess-26-5647-2022): leaf/canopy/land-area conventions; this probe uses catchment-area depths throughout.

These sources support physically plausible scales and definitions. They do not establish the frequency of the eight synthetic catchment categories or identify an AI training-data boundary.

## Review validation

Current evaluation results and per-model storage coverage are summarised below. SUMMA's general-adapter diagnostic is maintained separately from this probe; see the [SUMMA model card](../../../models/summa/README.md). Its archived FAIL is accepted by the maintainer and is not a water-leakage finding.

| Model | Native result |
|---|---|
| `reference_bucket` | PASS (OK) |
| `flex_lumped` | PASS (OK) |
| `flex_topo` | PASS (OK) |
| `sacsma_snow17` | PASS (OK) |
| `google_flood_forecast` | N/A (INCOMPLETE) |
| `wflow_sbm` | PASS (OK) |
| `summa` | FAIL (VIOLATION) |
| `cwatm` | PASS (OK) |
| `lisflood` | PASS (OK) |
| `dhbv2` | FAIL (VIOLATION) |

## Storage-bound coverage

| Model | Soil bound | Canopy bound | Other nonnegative stores exercised |
|---|---|---|---|
| `reference_bucket` | Nonzero storage | Nonzero storage | None |
| `flex_lumped` | Nonzero storage | Nonzero storage | gw, channel |
| `flex_topo` | Nonzero storage | Nonzero storage | gw, channel |
| `sacsma_snow17` | Nonzero storage | Zero throughout | gw, channel |
| `wflow_sbm` | Nonzero storage | Zero throughout | channel |
| `summa` | Nonzero storage | Nonzero storage | snw, gw, channel |
| `cwatm` | Nonzero storage | Nonzero storage | gw, channel |
| `lisflood` | Nonzero storage | Zero throughout | gw, channel |
| `dhbv2` | Nonzero storage | Zero throughout | gw, channel |

Coverage describes the 3,650 scored daily rows of all twelve cases. “Nonzero storage” means the bound is exercised on a nonzero reported store, not that it is reached or violated. Zero stores still undergo the formal check but do not exercise its upper bound; zero daily output does not imply interception is absent between outputs. Missing stores are not invented. Google is N/A for missing budget outputs and is excluded from this table.

δHBV passes cumulative closure (worst residual 0.1706% of precipitation) but fails the prescribed soil-storage bound. This is a boundary failure, not a water-leakage result.

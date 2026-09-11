# What the first two evaluations found

Read this before evaluating a third model or explaining a verdict: the
patterns repeat, and the numbers here are the baseline a new result is read
against. Both are archived in `models/result.csv`; the reports are in the
issues.

## google_flood_forecast (issue #1) — FAIL (INCOMPLETE), 2/7

Google Flood Hub's mean-embedding forecast LSTM via OpenHydroNet at commit
`828dfc5`, weights `google-floodhub-settings-55-epochs`. Discharge only,
declared as `mrro` + `dis`, so every budget probe is INCOMPLETE; the stress
probes still measure it.

| Probe | Result | Mechanism |
| --- | --- | --- |
| catchment-closure | INCOMPLETE | reports no ET or stores |
| resolution-invariance | 66% of rain between PT1H and PT1D | trained at one step, run at another |
| warming-response | −0.088 per unit of demand (needs 0.1) | marginal; temperature is a weak input |
| causality | PASS | after attributes were restricted to the first year |
| dry-down | non_degenerate: runoff CV 0.019 | a flat floor that never decays |
| steady-state | PASS | after the seasonal FAO-56 radiation mock was replaced |
| extreme-rain | 0.07 of a tenfold storm runs off | learned response saturates |

Adapter traps fixed along the way: whole-record attributes (causality), a
calendar-dependent radiation mock (false steady-state drift), CMAL sampling
noise (replaced by an exact mixture median). Inputs the probe does not
generate (HRES radiation, pressure) are mocked from the forcing with
Priestley–Taylor inverted from PET; 84 Caravan attributes at training mean
except the climate ones.

## dhbv2 (issue #2) — FAIL (VIOLATION), 3/7

δHBV 2.0 from `mhpi/dhbv2` v0.5.4 (`6d4bb68`), weights `hbv_2_ep100.pt`
from the MHPI S3 bucket, dmg 1.4.3 + hydrodl2 1.3.5, CPU torch 2.5.1. LSTM
writes `parBETA, parK0, parBETAET` daily; MLP writes 13 static HBV
parameters and 2 unit-hydrograph parameters; 3 components averaged. Reports
`pr, evspsbl, mrro, dis, mrso(SM), snw(SNOWPACK+MELTWATER), gw(SUZ+SLZ),
channel(UH store), canopy(0)`.

| Probe | Local (gate seeds) | CI (fresh seeds) | Mechanism |
| --- | --- | --- | --- |
| catchment-closure | residual 82.6% of rain; mrso 496–708 mm vs 320 cap; runoff ratio 1.15 | 64.2%; 555–746; 1.006 | regional-groundwater term; learned `parFC` |
| causality | PASS | PASS | bucket structure |
| dry-down | PASS | PASS | stores drain |
| extreme-rain | 1.07 of added rain returned | 1.05 | dynamic parameters release stored water (inferred, marginal) |
| resolution-invariance | 216% unscaled → with parameters scaled to the step (adapter .2): runoff 2.2%, ET 8.6%; failed the 5% limit, passes the 10% limit calibrated on the physical models (probe v2) | 482% (unscaled) | rate constants had no dt; what remains after scaling is the LSTM recurrence writing a lower `parBETAET` at the hourly step |
| steady-state | −1.08 mm/day unplaced | −1.04 | the same source term, constant |
| warming-response | PASS, −0.65/−0.70 per unit | PASS | ET follows PET and SM |

The headline mechanism, verified per step (residual constant to three
decimals):

    SLZ += parRT * clamp((Ac - parAC)/1000, -1, 1)   for Ac < 2500 km²

a learned inter-basin groundwater exchange (same idea as GR4J's X2), trained
on streamflow only, so nothing constrained it. On the 250 km² probe
catchment it is +1 to +2 mm/day. Sensitivity across attribute choices:
residual −48% to −76% of P in every configuration; runoff CV 0.07 → 1.1
depending only on the seasonality guess. Sensitivity to `uparea`:
(Q+AET)/P from 1.00 (7.8 km²) to 2.03 (2500 km²).

Packaging facts: shipped config does not fit weights (see SKILL.md);
`pip install /opt/dhbv2` fails under hatch-vcs in a one-commit checkout, so
the release wheel `dhbv2==0.5.4` is installed with `--no-deps` and the
checkout kept only for `LICENSE`; image 2.76 GB; one full 10-year record
runs in ~4 s inside the container, most of it loading 400 MB of MLP weights.

### Parameter scaling (adapter .2, 2026-09-05)

The maintainer's rule: a network that writes physical parameters in daily
units must have them rescaled to the case's step in the adapter, or the
probe measures the adapter's units. Per-day fractions → `1-(1-k)^dt`,
per-day amounts → `×dt`, unit hydrograph time constant `/dt` and length
`15/dt`; networks fed rates, physics fed depths. Verified: identity at
dt=1; hourly physics with the daily LSTM parameters repeated 24× matches
the daily volumes within 1%, so the residual 2.2%/8.6% is the recurrence.
AGENTS.md now states the rule.

## Harness changes the evaluations forced

- Time-based flood-event windows and per-variant time steps (issue #1).
- `resolution_invariance`, `response_sign`, `causality`, `dry_down`,
  `steady_state`, `monotone_response` criteria and their broken reference
  models (`reference_fixed_step`, `reference_anticipating`,
  `reference_climatology`, `reference_saturating`, `reference_restless`).
- Docker CPU request capped at the host count (CI has 4 cores).
- Extra store names `gw`, `channel`; closure over every reported store.
- Containers run as the invoking user (root without CAP_DAC_OVERRIDE could
  not write the runner-owned output directory on Linux).

## Physical reference models (2026-09-06)

`flex_lumped` and `flex_topo`, from chrimerss/HydrologicModels, are in
`must_pass` of every probe beside `reference_bucket`. Any probe PR is run
against the three first. Adapter deviations from the repository code are
listed in each model's README (ET `min(1,·)` limit; hillslope preferential
share moved once; wetland capillary rise with the wetland fraction;
capacities from static.json; step as a parameter). They calibrated one
tolerance: resolution-invariance 5% → 10% because FLEX-Topo moves 5.1%.
Trusted subprocess models are listed in `spec.TRUSTED_SUBPROCESS_MODELS`.

## SAC-SMA + Snow-17 (2026-09-07)

`sacsma_snow17`: a stdlib port of Upstream-Tech/SACSMA-SNOW17's Fortran
(`models/sacsma_snow17/sacsma_snow17.py`), fourth physical must-pass model;
passes all 14 probes. Validation recipe: f2py the Fortran in the scratch
venv with `PATH` including the venv (`meson`/`ninja` from pip) and
`--f77flags="-fallow-argument-mismatch -std=legacy"` for `ex_sac1`; drive
both step by step. SAC-SMA matches to 1e-4 mm; Snow-17 matches ten-year
totals to 0.02% but the compiled Fortran never carries `SBAESC` across
steps (instrumented `AESC19` computes 0.99, carryover holds 0), `tiny` is
uninitialised and `SNOF` never set. Snow-17 is defined on whole hours, so
the manifest lists PT1D and PT1H only. Conservative choices: SCF=1, SIDE
declared as negative `gwex`, capacities from static.json.

## Round-2 probes (2026-09-06)

Seven more probes, all gated on the three physical models: runoff-bounds
(the mass question a runoff-only model must answer), area-invariance,
response-nonnegativity, antecedent-monotonicity, phase-counterfactual,
energy/pet-consistency (regimes from the model's own soil range, not the
stated capacity), momentum/routing-conservation (`channel` <= max_lag ×
recent peak runoff). Six broken models were added for them. Lessons: a
broken model must be broken enough to cross the bound (overflowing at
+80% of rain, not ×1.3); `invariance` gained `optional` for stores a
model lacks; the response-nonnegativity storm has to be sharp (120 mm) to
expose a negative impulse lobe; FLEX-Topo's one-step-per-day beta
partition overshoots a small store and dipped 0.84 mm/day ten days after
an added storm, fixed by sub-stepping the partition (documented deviation).
`gwex` joined the contract: δHBV declares its regional groundwater term
and closes to 1e-8 (adapter .3).

## antecedent-monotonicity window (2026-09-11)

`antecedent_monotonicity` opened its 30-day window on the step after the
last antecedent rain, but the generator leaves ten rainless days before
the storm, so the window started ten days early: the wet run's recession
of antecedent water through those days was scored as storm runoff, and
`storm_mm` summed all rain in the window rather than the storm. Found on
wflow_sbm, whose `runoff_dry_mm`/`runoff_wet_mm` reproduced exactly when
summed from 9 July rather than the storm of 19 July. The criterion now
opens the window on the first rain after the last step on which the
variants differ, and takes the share against the storm itself (the
unbroken run of rain the window opens on), which is what probe.yaml and
the README always said. `storm_time` in the diagnostics shows where the
window opened; check it against the generator once when a criterion
locates an event from the forcing.

The shift took about a third of the apparent memory from the gate models
(bucket 0.108 → 0.072 of window rain on the first gate seed). SAC-SMA's
memory is genuinely small: its wet catchment holds 62 mm more soil water
at the storm but evaporates most of it, returning 1.4 mm in 30 days and
5.6 mm in a year. With the window moved and the old denominator it fell to
0.019, under the 0.02 bound; against the storm it is 0.024. Over 200 seeds
SAC-SMA falls below 0.02 on 4.0% before the fix, 11.5% with the window
moved alone, and 2.0% with the storm as the denominator (minimum 0.011);
the other physical models never go below 0.05 and the cheater is exactly 0.
Taking the share against the storm does loosen the bound on the gate
seeds, from 1.4–2.3 mm of extra runoff to 1.2 mm. Rain on the day after
the storm joins it, since the run ends at the first dry step; that happens
on 25 of the 200 seeds and changes none of these rates.
`min_share` stays 0.02; 0.01 would clear all 200 seeds. All five archived
models still pass; the lowest share on the gate seeds is dhbv2 0.21,
google_flood_forecast 0.10, flex_lumped 0.17, flex_topo 0.11, sacsma_snow17
0.024.

wflow_sbm merged while the fix was in review, with its row scored on the
old window. It still fails on the same two seeds, at 0.0013 and 0.0005 of
the storm instead of −0.0011 and −0.0010: its wet month's water has been
evaporated down to the rooting depth before the storm arrives, which no
window placement changes. Its README's sensitivity table was re-scored
from the kept per-setting runs, where the old criterion reproduces every
entry; with roots through 99% of the column the second seed now clears the
bound by a hair (0.020).

## Known open items

- The `probe` workflow's `container` job has failed on every recorded run
  since before these evaluations; its `gate` job passes. Not yet diagnosed.
- `reference_in_sample` and `reference_calendar` are baselines for the
  `regime_transfer` and `invariance` probes nobody has written yet.
- Upstream: tell MHPI that `dhbv_2.zip`'s `config.yaml` does not match its
  weights.
- Candidate next models: δHBV 2.0 MTS (hourly, same bucket, first hourly
  physics core for the resolution probe); NOAA-OWP NextGen LSTM
  (Apache-2.0, weights in repo, discharge-only so less new information).

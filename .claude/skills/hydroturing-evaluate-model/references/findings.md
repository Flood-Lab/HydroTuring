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

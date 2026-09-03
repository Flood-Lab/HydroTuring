# Writing a probe

A probe is a conservation law made executable. It defines a generated case, a
set of binary criteria, and the reference models it must be able to separate.

## Start from a template

```bash
ht init-probe --list-templates           # what shapes are available
ht init-probe                            # writes probe-draft.yaml
#   ... fill in the fields ...
ht init-probe --from probe-draft.yaml    # creates probes/<law>/<slug>/
```

`--template <kind>` starts from a probe of that shape instead of the plain
one, and brings a matching generator skeleton with it.

| Template | The question it asks |
| --- | --- |
| `default` | does the budget close over one generated case |
| `extrapolation-space` | does it still close on catchments outside the range models are fitted to |
| `extrapolation-time` | does it still close under conditions outside anything earlier in the record |
| `counterfactual` | add water to the same case and ask where it went |
| `invariance` | change something the physics does not depend on and require nothing to move |

Each of these ships filled-in baselines and passes `ht gate` as scaffolded, so
you can run the gate before you have written a line of physics and see it
separate the reference models. Then replace the placeholder case with yours.

The draft is heavily commented. Lines beginning `#!` are guidance and are
stripped from the generated `probe.yaml`; your own `#` comments are kept.

What you get already validates and already honours the length contract, so
your first `ht gate` fails on your physics rather than on scaffolding.

```
probes/<law>/<slug>/
  probe.yaml     spec, schema-validated, carries your authorship
  generate.py    generate(seed) -> (DataFrame, dict)
  README.md      the physics in prose
```

`probes/mass/catchment-closure/` is the reference implementation. Read it
before you start.

## The generator

```python
def generate(seed: int) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    ...
    return forcing, static
```

Requirements:

- deterministic given the seed, byte for byte
- a `time` column
- exactly the period plus the spinup, counted in rows at the probe's step:
  `period_years * 365 + spinup_days` rows for a daily probe, and for a
  sub-daily one `(period_days + spinup_days) / step`, where `period_days`
  may replace `period_years` in `case`
- fluxes as rates in mm per day whatever the step, so a minute of rain at
  one millimetre is a rate of 1440 mm/day
- no committed data files; CI rejects anything over 1 MB under `probes/`

Determinism is not a nicety. A benchmark whose failures cannot be reproduced
is unusable the first time a result is disputed.

## Criteria

Each entry in `criteria` names a registered criterion and its parameters.
Every one is binary.

| Criterion | Asserts | Scored over |
| --- | --- | --- |
| `closure` | the budget closes to within a share of the driving flux | one run |
| `state_bounds` | every reported storage stays physical | one run |
| `et_plausible` | ET is non-negative and bounded by potential ET | one run |
| `non_degenerate` | the partition and the response are non-trivial | one run |
| `forcing_fidelity` | the model reports back the forcing it was given | one run |
| `regime_transfer` | closure holds out of range as well as in range | labelled stretches |
| `counterfactual_response` | added water is partitioned, not absorbed | paired runs |
| `invariance` | a transform the physics ignores changes nothing | paired runs |

Picking a denominator for `closure` and `regime_transfer`:

| Denominator | Use for | Floor |
| --- | --- | --- |
| `sum_pr` | water budgets | not needed, precipitation is non-negative |
| `sum_abs_rn` | energy budgets | required, net radiation crosses zero nightly |
| `sum_inflow` | routing | not needed |

### Labelled stretches

`regime_transfer` scores parts of one record separately and compares them,
which is how a failure confined to a rare stretch is kept from being diluted
by an ordinary one. The generator marks each step with a `_regime` column.

Any column whose name begins with `_` is an annotation: the criteria see it and
the harness strips it before staging, so the model never receives it. That is
not a detail. A column called `_regime` whose values are `ordinary` and
`anomaly` would otherwise tell a model exactly which part of the record it is
being judged on.

### Paired runs

`case.variants` names the cases the generator builds for one seed:

```yaml
case:
  generator: generate.py
  variants: [control, wetter]
```

The model is run once per variant and paired criteria are handed all the
results. The generator takes the variant as well as the seed:

```python
def generate(seed: int, variant: str = "control") -> tuple[pd.DataFrame, dict]:
```

The first variant is the control. Every non-paired criterion is scored against
it alone, so adding a variant to a probe never silently changes what its
existing criteria measure. Declaring variants without a paired criterion, or a
paired criterion without variants, is rejected at load time rather than
becoming a silent no-op.

Draw everything that comes from the seed before you branch on the variant. Two
variants that differ in the weather as well as in the perturbation cannot
isolate the perturbation, and the comparison means nothing.

### Variants at different steps

A resolution transform is the same weather at two steps. Declare the step of
any variant that does not run at `case.timestep`:

```yaml
case:
  timestep: PT1M
  period_days: 30
  spinup_days: 10
  variants: [minute, hourly]
  timesteps:
    hourly: PT1H
```

The control keeps `case.timestep`. Each variant's row count follows its own
step, the harness hands the model the step in `request.json`, and criteria
integrate each run with its own step, so a paired criterion can compare a
minute record with its hourly aggregate. Aggregate, do not redraw: the coarse
variant must carry exactly the water of the fine one, and
`resolution_invariance` refuses a pair that does not. A model has to declare
every step the probe uses; one that runs at a single step is INCOMPATIBLE
with the probe, which is the honest verdict for it. See
`probes/mass/resolution-invariance` for the worked example.

## Baselines

```yaml
baselines:
  must_pass: [reference_bucket]
  must_fail:
    reference_leaky: closure
    reference_cheater: state_bounds
    reference_degenerate: non_degenerate
```

`must_pass` guards against tolerance drift: if a physically exact model ever
fails your probe, the probe is wrong. `must_fail` pins which criterion does
the catching, so a probe cannot appear to work while catching things for the
wrong reason.

The reference models available today:

| Model | What it does | Caught by |
| --- | --- | --- |
| `reference_bucket` | conserves water exactly by construction | nothing, it must pass |
| `reference_leaky` | hides a silent 15% sink | `closure` |
| `reference_cheater` | solves for storage as whatever balances the budget | `state_bounds` |
| `reference_degenerate` | evaporates all precipitation, produces no runoff | `non_degenerate`, `counterfactual_response` |
| `reference_in_sample` | exact in range, leaks outside it | `regime_transfer` |
| `reference_calendar` | recession drifts with the calendar year | `invariance` |
| `reference_streamflow_only` | reports too little to be checked | scored `INCOMPLETE` |

If your probe needs a broken model that does not exist yet, add it under
`models/` alongside the probe. A criterion with nothing that trips it is
untested.

## Before opening a PR

```bash
ht validate
ht gate --probe <law>/<slug>
pytest -q
```

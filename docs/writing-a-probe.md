# Writing a probe

A probe is a conservation law made executable. It defines a generated case, a
set of binary criteria, and the reference models it must be able to separate.

## Anatomy

```
probes/<law>/<slug>/
  probe.yaml     spec, schema-validated
  generate.py    generate(seed) -> (DataFrame, dict)
  README.md      the physics in prose
```

`probes/mass/catchment-closure/` is the reference implementation.

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
- exactly `period_years * 365 + spinup_days` rows
- no committed data files; CI rejects anything over 1 MB under `probes/`

Determinism is not a nicety. A benchmark whose failures cannot be reproduced
is unusable the first time a result is disputed.

## Criteria

Each entry in `criteria` names a registered criterion and its parameters.
Available today: `closure`, `state_bounds`, `et_plausible`, `non_degenerate`,
`forcing_fidelity`. Every one is binary.

Picking a denominator for `closure`:

| Denominator | Use for | Floor |
| --- | --- | --- |
| `sum_pr` | water budgets | not needed, precipitation is non-negative |
| `sum_abs_rn` | energy budgets | required, net radiation crosses zero nightly |
| `sum_inflow` | routing | not needed |

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

If your probe needs a broken model that does not exist yet, add it under
`models/` alongside the probe. A criterion with nothing that trips it is
untested.

## Before opening a PR

```bash
ht validate
ht gate --probe <law>/<slug>
pytest -q
```

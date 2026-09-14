<!-- Title: [PROBE: <law>] <short description> -->

Closes #<proposal issue>

## What this probe asserts

<!-- The residual equation and what a failure means physically. -->

## Discrimination

Every probe must separate the reference models. Fill in what `ht gate` reports:

| Reference model | Expected | Criterion that catches it |
| --- | --- | --- |
| `reference_bucket` | PASS | n/a |
| `reference_leaky` | FAIL | |
| `reference_cheater` | FAIL | |
| `reference_degenerate` | FAIL | |

## Checklist

- [ ] There is an `accepted` proposal issue and this PR closes it
- [ ] `authors` in `probe.yaml` names every author with `name`, `affiliation` and, where you have one, `orcid`, matching the proposal issue; CONTRIBUTORS.md, CITATION.cff and the paper's author list are built from it
- [ ] `ht validate` passes
- [ ] `ht gate --probe <id>` passes
- [ ] The generator is deterministic given a seed and commits no data
- [ ] The probe runs in under a minute on a two-core runner
- [ ] Tolerance and denominator are justified in `probe.yaml`
- [ ] I constructed a model that passes closure without doing physics, and a criterion catches it

# HydroTuring

**[flood-lab.github.io/HydroTuring](https://flood-lab.github.io/HydroTuring/)** &middot; English, Español, 中文

A benchmark that asks one question of any AI hydrologic model: **does it
conserve what physics says it must conserve?**

Not whether it fits a hydrograph. Whether its water budget closes, its energy
budget closes, and its routing conserves momentum. A model passes HydroTuring
only when every criterion of every probe passes.

```
$ ht run --model reference_bucket
reference_bucket v1.0.0  ->  PASS (OK)  [1/1 probes passed]
  PASS  mass/catchment-closure
        ok   closure            cumulative residual 0.0000% of sum_pr (limit 5.0%)
        ok   state_bounds       all storages stay physical
        ok   et_plausible       cumulative ET is 0.690 of potential ET (limit 1)
        ok   non_degenerate     partition and variability are non-trivial
        ok   forcing_fidelity   reported forcing matches the input
```

## The idea

A benchmark that only checks closure is trivially gamed. Three models in this
repository exist to prove it, and every probe must be able to catch all three
before it can be merged.

| Reference model | What it does | Caught by |
| --- | --- | --- |
| `reference_bucket` | conserves water exactly by construction | nothing, it must pass |
| `reference_leaky` | hides a silent 15% sink | `closure` |
| `reference_cheater` | solves for storage as whatever balances the budget | `state_bounds` |
| `reference_degenerate` | evaporates all precipitation, produces no runoff | `non_degenerate` |

`reference_cheater` is the one worth dwelling on. Its closure residual is
**exactly zero on every seed, forever**. Randomising the forcing cannot touch
it, because the cheat is in its internal wiring rather than its memory. It is
caught only because the model contract requires absolute storage states rather
than tendencies, so the storage it invents has to stay physical, and it does
not.

## How a case is generated

Forcing is generated fresh at run time from a recorded seed. Nothing is
committed as data, so there is nothing to memorise, and a reviewer reviews a
short deterministic script instead of a binary blob. Each probe runs several
seeds and all of them must pass, so no model gets through on a lucky draw.
Any run reproduces exactly with `ht run --seed <n>`.

## Verdicts

Binary, with the reason recorded separately, because these mean different
things:

- `VIOLATION` the model reported its budget and the budget did not close.
- `INCOMPLETE` the model never reported enough to be checked. Every
  streamflow-only model lands here today. It has not violated conservation;
  it has declined to be falsifiable.

The verdict is one bit. Everything under it stays quantitative, so a paper can
show that one model leaks 6% and another 40% long before anyone crosses the
line.

## Quick start

```bash
pip install -e '.[dev]'

ht init-probe                            # start a new probe
ht list                                  # probes and models
ht validate                              # schema-check everything
ht gate                                  # the probe acceptance gate
ht run --model reference_bucket          # evaluate one model
ht run --model my-model --json out.json  # machine-readable report
```

Without installing, `./ht` runs the CLI straight from `src/`.

## Submitting a model

Your model may be written in any language. It ships as a container plus a
thin adapter that reads `/io/request.json` and writes `/io/output/result.csv`.
The adapter is usually thirty lines. See [docs/adapting-a-model.md](docs/adapting-a-model.md),
and `AGENTS.md` if you are having a coding agent build the sandbox for you.

The container runs with no network and never sees the probe code, so a model
cannot read the tolerance it is being judged against.

## Contributing a probe

**[ROADMAP.md](ROADMAP.md) lists the probes we want.** Claim one rather than
inventing one. The two we most want are the cross-budget consistency probes:
a model can close its water budget and its energy budget while being
incoherent between them, and nothing in the suite currently notices.

Do not start from a blank page:

```bash
ht init-probe                            # writes probe-draft.yaml
#   ... fill in the fields ...
ht init-probe --from probe-draft.yaml    # creates probes/<law>/<slug>/
ht gate --probe <law>/<slug>             # prove it discriminates
```

**Contributors of merged probes are authors on the benchmark paper.** See
[CONTRIBUTING.md](CONTRIBUTING.md#credit) for the threshold and the process,
and [GOVERNANCE.md](GOVERNANCE.md) for how disagreements about tolerances get
settled.

## Status

Suite `0.1.0`, pre-release. One probe (mass), synthetic track only. Energy,
momentum and the real-data track are next. Scores are only comparable within
a suite version.

## Layout

```
src/hydroturing/   harness: protocol, runners, criteria, scoring
probes/<law>/<id>/ probe.yaml, generate.py
models/<name>/     model.yaml, Dockerfile, adapter
schemas/           JSON schemas for both spec files
```

<p align="center">
  <img src="res/HydroTuring.png" alt="HydroTuring Initiative" width="540">
</p>

<p align="center">
  <b><a href="https://flood-lab.github.io/HydroTuring/">flood-lab.github.io/HydroTuring</a></b>
  &middot; English, Español, 中文
</p>

A benchmark that asks one question of any AI hydrologic model: **does it
conserve what physics says it must conserve?**

Not whether it fits a hydrograph. Whether its water budget closes, its energy
budget closes, and its routing conserves momentum. A model passes HydroTuring
only when every criterion of every probe passes.

> **Write a probe, join the paper.** The suite is only as good as the physics
> people bring to it. Anyone may submit a probe — you do not need to be
> invited, affiliated, or known to us. **If your probe is merged, you are a
> co-author on the HydroTuring paper.** One merged probe is the whole
> threshold. Start at [ROADMAP.md](ROADMAP.md#probes-we-want) or propose your
> own.

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
| `reference_in_sample` | exact in range, leaks once the forcing leaves it | `regime_transfer` |
| `reference_calendar` | recession that drifts with the calendar year | `invariance` |

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

ht init-probe --list-templates            # probe shapes to start from
ht init-probe --template invariance      # start a new probe
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

**Please submit one.** A benchmark with four probes tests four things; the
reason this repository is open is that the physics worth testing is wider
than any one group knows. If you have spent time with a conservation law that
AI models get wrong, that law is a probe, and we would rather have it from you
than approximate it ourselves.

**Contributors of merged probes are co-authors on the benchmark paper.** The
threshold is one probe, merged and passing the acceptance gate. This is how
model intercomparison projects have always worked in this field: you
contribute an experiment, you are an author on the paper that reports it.
[CONTRIBUTING.md](CONTRIBUTING.md#credit) has the full terms — author order,
the right to decline, and what happens before anything is submitted.

**[ROADMAP.md](ROADMAP.md) lists the probes we want**, each with a difficulty
and a note on what it discriminates. Claiming one is easier than inventing
one, but inventing one is welcome too; propose it first so nobody builds it
twice. The two we most want are the cross-budget consistency probes: a model
can close its water budget and its energy budget while being incoherent
between them, and nothing in the suite currently notices.

Do not start from a blank page:

```bash
ht init-probe --list-templates           # the shapes available
ht init-probe --template <kind>          # writes probe-draft.yaml
#   ... fill in the fields ...
ht init-probe --from probe-draft.yaml    # creates probes/<law>/<slug>/
ht gate --probe <law>/<slug>             # prove it discriminates
```

The templates cover conservation over one case, extrapolation in space and in
time, counterfactual response, and invariance. Each scaffolds into a probe that
already passes the gate against a placeholder case, so you can watch it
separate the reference models before writing any physics, then replace the
case with yours. `docs/writing-a-probe.md` has the details.

The gate is the only bar that matters, and it is a technical one: your probe
must pass an exact physical model and catch the deliberately broken ones. See
[GOVERNANCE.md](GOVERNANCE.md) for how disagreements about tolerances get
settled.

## Status

Suite `0.1.0`, pre-release. One probe (mass), synthetic track only. Energy,
momentum and the real-data track are next. The harness runs paired cases and
scores labelled regimes, so the generalisation probes on the roadmap —
extrapolation in space and time, counterfactual response, invariance — are
unblocked and unclaimed. Scores are only comparable within a suite version.

## Layout

```
src/hydroturing/   harness: protocol, runners, criteria, scoring
probes/<law>/<id>/ probe.yaml, generate.py
models/<name>/     model.yaml, Dockerfile, adapter
schemas/           JSON schemas for both spec files
```

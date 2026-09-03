# Contributing to HydroTuring

Two very different kinds of contribution live in this repository. A **probe**
is a test of a conservation law. A **model** is something to be tested.

## Credit

Contributing to a benchmark is real scientific work, and this project treats
it that way.

**At the probe level.** Every probe carries its authors in `probe.yaml`. They
are named in every report that runs the probe and in each Zenodo release.
Nobody's contribution disappears into a commit log.

**On the benchmark paper.** Contributors of at least one merged probe are
invited to be co-authors on the HydroTuring paper. This follows how model
intercomparison projects have always worked in this field: you contribute an
experiment, you are an author on the paper that reports it.

The threshold is one probe that is merged and passes the acceptance gate.
Documentation improvements, bug fixes and model adaptations are genuinely
valuable and are credited in `CONTRIBUTORS.md`, but do not by themselves earn
authorship, because the paper is about the benchmark's design and findings.

Author order will be settled before submission and circulated to every
contributor for agreement. Anyone may decline authorship, and no manuscript
is submitted before all listed authors have seen and approved it.

If you are unsure whether what you have in mind clears the bar, ask in your
proposal issue. We would rather tell you up front than have you guess.

## Contributing a probe

### 0. Pick something from the roadmap

[ROADMAP.md](ROADMAP.md) lists the probes we want, with difficulty labels and
a note on what each one discriminates. Claiming one is easier than inventing
one, and it means two people do not build the same thing.

The two we most want are the cross-budget consistency probes, because a model
can close its water budget and its energy budget while being incoherent
between them, and nothing in the suite currently notices.

### 1. Propose it first

Open a [probe proposal](../../issues/new?template=probe_proposal.yml). The
form asks for the residual equation, the tolerance and its denominator, and
one question that matters more than the rest:

> How would an unphysical model pass this probe?

If you cannot construct such a model, the probe may not discriminate, and it
is better to find that out before you build it.

A maintainer labels the issue `accepted`. Then build.

### 2. Scaffold it

Do not start from a blank page. Get the template, fill in the fields, and let
the tool build the directory:

```bash
ht init-probe                            # writes probe-draft.yaml
#   ... fill it in ...
ht init-probe --from probe-draft.yaml    # creates probes/<law>/<slug>/
```

You get a `probe.yaml` that already validates, a `generate.py` skeleton that
already honours the length contract, and a `README.md` with the sections a
reviewer will look for. What is left is the physics, which is the part only
you can write.

```
probes/<law>/<slug>/
  probe.yaml      the spec: authors, what is required, the criteria
  generate.py     deterministic given a seed, returns (DataFrame, static dict)
  README.md       the physics, in prose
```

`probes/mass/catchment-closure/` is the worked example. Read it before you
start.

Two things to hold onto:

**Ship a generator, not a dataset.** Forcing is made fresh at run time so a
model cannot memorise the case. CI rejects any file over 1 MB under `probes/`.
Your generator must return byte-identical output for the same seed.

**Declare which criterion catches what.** `baselines.must_fail` maps each
broken reference model to the specific criterion that must trip. This is not
bookkeeping. It is what stops a probe from appearing to work while actually
catching things for the wrong reason.

### 3. The gate

```bash
ht validate
ht gate --probe <law>/<slug>
```

The gate runs your probe against every reference model and checks that the
physically exact one passes and each broken one fails on its declared
criterion. A probe that cannot separate them is measuring nothing, and a
tolerance tightened past what exact physics achieves breaks here rather than
silently mislabelling honest models later.

### 4. Submit

Title the PR `[PROBE: <law>] <description>` and fill in the template. Review
is two passes: one on the physics, one on the implementation.

## Contributing a model

See [AGENTS.md](AGENTS.md) for the contract and
[docs/adapting-a-model.md](docs/adapting-a-model.md) for the walkthrough. Your
model can be written in any language. It ships as a container with a small
adapter.

If you want the maintainers to package or evaluate an existing model, open a
[model submission](../../issues/new?template=model_submission.yml). Include
the model's code repository, exact version or commit, optional weights
repository, and whether it is AI-based, AI+physics or physics.

The form also asks whether there is a time window you want the test to run
over. Submitted models are scored on the largest flood event of the record,
a month for a daily model and a week for an hourly one unless you say
otherwise, so that a heavy model fits its time budget.

Submitting a model that fails is welcome and useful. `FAIL` with reason
`INCOMPLETE` is the current state of nearly every published rainfall-runoff
model, and recording that honestly is part of what the benchmark is for.

## Adding a criterion

New criteria go in `src/hydroturing/criteria/` and register with the
`@criterion` decorator. A criterion must be binary, must report the quantity
it measured, and should come with a reference model that it catches. Adding a
criterion without something that trips it means nothing tests the criterion.

## Running everything locally

```bash
pip install -e '.[dev]'
pytest -q          # the harness
ht validate        # the specs
ht gate            # the benchmark's own discrimination
```

No Docker needed for any of that: the reference models declare
`runner: subprocess`. Docker is only used for submitted models.

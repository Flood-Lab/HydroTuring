# Contributing to HydroTuring

Two very different kinds of contribution live in this repository, and **we are
asking for both**. A **probe** is a test of a conservation law. A **model** is
something to be tested. A benchmark with seven probes tests seven things, and
a benchmark that has scored a handful of models says very little about the
field. Neither gap closes without people outside this repository.

## How a contribution moves

The same three steps whichever kind you are bringing.

```
1. Open an issue        a probe proposal, or a model proposal
        |               a maintainer labels it `accepted` and assigns it
        v
2. Fork the repository  once it is accepted, not before
        |
        v
3. Open a pull request  from your fork, linking the issue it closes
```

**Propose before you build.** The issue is not a formality and it is not a
queue ticket. For a probe it is where we find out whether the thing
discriminates, which is cheaper to learn in a paragraph than in three hundred
lines. For a model it is where we find out whether it can be containerised at
all, and who is going to do it. Either way the answer arrives before you have
spent a weekend.

A pull request with no accepted issue behind it will be read, and may well be
good, but it starts the conversation at the most expensive point rather than
the cheapest one.

### Who the issue is assigned to

Whoever is actually going to do the work.

Usually that is you: you propose it, we label it `accepted`, you are assigned,
and it is yours. Nobody else will build it while your name is on it.

Sometimes it is not you, and that is fine and expected. Packaging a model can
need a GPU, a licence, weights you cannot redistribute, or a week you do not
have. **If you know of a model worth testing but cannot package it yourself,
propose it anyway.** The model proposal form asks directly, under *packaging
status*; say that help is needed and the issue is assigned to the maintainer
(@chrimerss) instead. The proposal still counts as yours.

A proposal you cannot execute is still a contribution. Knowing which models
are worth putting through the benchmark is a judgement about the field, and it
is not one the maintainer can make alone.

## Credit

Contributing to a benchmark is real scientific work, and this project treats
it that way.

**At the probe level.** Every probe carries its authors in `probe.yaml`. They
are named in every report that runs the probe and in each Zenodo release.
Nobody's contribution disappears into a commit log.

**On the benchmark paper.** The threshold for co-authorship is:

| Contribution | Counts as |
| --- | --- |
| one merged probe, passing the acceptance gate | co-authorship |
| five accepted model proposals | one probe, so co-authorship |

This follows how model intercomparison projects have always worked in this
field: you contribute an experiment, you are an author on the paper that
reports it.

Model proposals count whether or not you package the model yourself. A
proposal is *accepted* when the issue is labelled `accepted`; what is being
credited is the judgement of what belongs in the benchmark, and that judgement
is the same whoever builds the container afterwards. A proposal declined as
out of scope does not count, which is why the form asks you to make the case.

The two routes are additive: three model proposals plus one probe is well past
the line, and so is one probe on its own. Partial progress carries — nothing
resets between suite versions.

Documentation improvements, bug fixes, reviews and harness work are genuinely
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

A maintainer labels the issue `accepted` and assigns it to you. Then fork the
repository, and build.

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

Push to your fork, then open a pull request titled
`[PROBE: <law>] <description>` and fill in the template. Link the proposal
issue so it closes on merge. Review is two passes: one on the physics, one on
the implementation.

## Contributing a model

**We want more models, and we want the ones you think matter.** Every
published rainfall-runoff model, every LSTM, every foundation model with a
hydrologic claim is in scope. A benchmark is only as interesting as what has
been put through it.

### 0. Propose it first

Open a [model proposal](../../issues/new?template=model_submission.yml).
Include the model's code repository, the exact version or commit, the optional
weights repository, and whether it is AI-based, AI+physics or physics.

Two fields decide what happens next:

**Packaging status.** Say plainly whether the Dockerfile and adapter are
ready, whether you need help with the adapter, or whether you need help with
both. This is what assigns the issue: ready means it is yours to finish, help
needed means the maintainer takes it. It does not affect whether the proposal
is accepted, and it does not affect your credit.

**Evaluation window.** Submitted models are scored on the largest flood event
of the record rather than all ten years — a month for a daily model and a week
for an hourly one unless you say otherwise — so that a heavy model fits its
time budget.

A maintainer labels the issue `accepted` and assigns it. Then fork, and build.

### 1. Package it

See [AGENTS.md](AGENTS.md) for the contract and
[docs/adapting-a-model.md](docs/adapting-a-model.md) for the walkthrough. Your
model can be written in any language. It ships as a container with a small
adapter.

Declare what the model genuinely produces and nothing more. A model that
predicts discharge alone is scored `FAIL (INCOMPLETE)`, which is the honest
outcome and a completely different statement from `FAIL (VIOLATION)`.

### 2. Submit

Push to your fork, then open a pull request titled `[MODEL] <name>`, linking
the proposal issue.

**Submitting a model that fails is welcome and useful.** `FAIL` with reason
`INCOMPLETE` is the current state of nearly every published rainfall-runoff
model, and recording that honestly is part of what the benchmark is for.
Nobody is embarrassed by a result here; the point is to have one. Review is on
the contract — does the adapter honour `/io`, is `emits` honest, is the image
reproducible — never on the verdict.

## Adding a criterion

New criteria go in `src/hydroturing/criteria/` and register with the
`@criterion` decorator. A criterion must be binary, must report the quantity
it measured, and should come with a reference model that it catches. Adding a
criterion without something that trips it means nothing tests the criterion.

Criteria are harness work rather than science, so they do not by themselves
earn authorship. They usually arrive alongside a probe that needs them, which
does.

## Running everything locally

```bash
pip install -e '.[dev]'
pytest -q          # the harness
ht validate        # the specs
ht gate            # the benchmark's own discrimination
```

No Docker needed for any of that: the reference models declare
`runner: subprocess`. Docker is only used for submitted models.

# Contributing to HydroTuring

Two very different kinds of contribution live in this repository. A **probe**
is a test of a conservation law. A **model** is something to be tested.

## Contributing a probe

### 1. Propose it first

Open a [probe proposal](../../issues/new?template=probe_proposal.yml). The
form asks for the residual equation, the tolerance and its denominator, and
one question that matters more than the rest:

> How would an unphysical model pass this probe?

If you cannot construct such a model, the probe may not discriminate, and it
is better to find that out before you build it.

A maintainer labels the issue `accepted`. Then build.

### 2. Build it

```
probes/<law>/<slug>/
  probe.yaml      the spec: what is required, how it is generated, the criteria
  generate.py     deterministic given a seed, returns (DataFrame, static dict)
  README.md       the physics, in prose
```

`probes/mass/catchment-closure/` is the worked example.

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

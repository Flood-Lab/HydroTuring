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
reference_bucket v1.0.0  ->  ✅ PASS (OK)  [1/1 probes passed]
  ✅ PASS  mass/catchment-closure
        ✅  closure            cumulative residual 0.0000% of sum_pr (limit 5.0%)
        ✅  state_bounds       all storages stay physical
        ✅  et_plausible       cumulative ET is 0.664 of potential ET (limit 1)
        ✅  non_degenerate     partition and variability are non-trivial
        ✅  forcing_fidelity   reported forcing matches the input
```

<p align="center">
  <img src="res/hydroturing-infra.svg" alt="Any AI hydrologic model enters through the /io contract, HydroTuring probes it against mass, energy and momentum, and a verdict with its reason comes out" width="800">
</p>

## The idea

A benchmark that only checks closure is trivially gamed. A set of deliberately
broken models lives in this repository to prove it, and no probe is merged
until it passes an exactly conservative model and catches the broken ones it
names.

`reference_cheater` is the one worth dwelling on. Its closure residual is
**exactly zero on every seed, forever**. Randomising the forcing cannot touch
it, because the cheat is in its internal wiring rather than its memory. It is
caught only because the model contract requires absolute storage states rather
than tendencies, so the storage it invents has to stay physical, and it does
not.

## The probes

Seven, all of them mass so far. Each was merged only after the acceptance gate
saw it pass a model that conserves water exactly and fail a purpose-built
broken one on the named criterion. `ht list` prints them;
[ROADMAP.md](ROADMAP.md#probes-we-want) has the seventeen more we want, all
unclaimed.

| Probe | Law | What it asks | The broken model it catches |
| --- | --- | --- | --- |
| [`mass/catchment-closure`](probes/mass/catchment-closure) | mass | Does the water budget close over ten generated years? | `reference_leaky`, `reference_cheater`, `reference_degenerate` |
| [`mass/resolution-invariance`](probes/mass/resolution-invariance) | mass | The same month at the minute, the hour and the day: do the integrated volumes agree? | `reference_fixed_step`, `reference_degenerate` |
| [`mass/warming-response`](probes/mass/warming-response) | mass | The same rain with the air 3 °C warmer and 3 °C cooler: does runoff move the way physics says, in both directions? | `reference_degenerate`, `reference_streamflow_only` |
| [`mass/causality`](probes/mass/causality) | mass | One storm added mid-record: nothing may change before it, and runoff must answer after it. | `reference_anticipating` |
| [`mass/dry-down`](probes/mass/dry-down) | mass | Two years without rain: runoff can only fall, and no more may drain than the catchment held. | `reference_climatology` |
| [`mass/steady-state`](probes/mass/steady-state) | mass | Three years of the same day: does everything settle, and does the budget balance once it has? | `reference_restless` |
| [`mass/extreme-rain`](probes/mass/extreme-rain) | mass | The largest storm scaled up to ten times: runoff may not fall, nor exceed the rain that was added. | `reference_saturating` |

## Models

`models/` holds two kinds. A submitted model is there to be evaluated, and
every run of one is appended to [models/result.csv](models/result.csv). A
reference model is there to test the probes rather than to be tested: one that
is exactly conservative and must pass everything, and a set that is each broken
in one specific way, so that no criterion goes untested.

| Model | Kind | What it does | Standing |
| --- | --- | --- | --- |
| [`google_flood_forecast`](models/google_flood_forecast) | submitted | The mean-embedding forecast LSTM behind Google Flood Hub, at the published weights. Predicts discharge and nothing else. | **FAIL (INCOMPLETE)**, 2 of 7 probes passed |
| `reference_bucket` | exact | conserves water exactly by construction | must pass every probe |
| `reference_leaky` | broken | hides a silent 15% sink | caught by `closure` |
| `reference_cheater` | broken | solves for storage as whatever balances the budget | caught by `state_bounds` |
| `reference_degenerate` | broken | evaporates all precipitation, produces no runoff | caught by `non_degenerate`, `response_sign` |
| `reference_fixed_step` | broken | treats every row as a day whatever the step is | caught by `resolution_invariance` |
| `reference_anticipating` | broken | smooths runoff over a centred window, so three days of the future are in every value | caught by `causality` |
| `reference_climatology` | broken | emits the seasonal mean whatever falls, and keeps flowing without rain | caught by `dry_down` |
| `reference_saturating` | broken | caps its daily runoff, so an extreme storm adds rain and no runoff | caught by `monotone_response` |
| `reference_restless` | broken | a recession on an internal thirty-day clock, so it never settles | caught by `steady_state` |
| `reference_streamflow_only` | honest limit | reports discharge only, from a store that never reads the temperature | scored INCOMPLETE on budget probes; caught by `response_sign` |
| `reference_in_sample` | broken | exact in range, leaks once the forcing leaves it | waiting for a `regime_transfer` probe |
| `reference_calendar` | broken | a recession that drifts with the calendar year | waiting for an `invariance` probe |

## How a case is generated

Forcing is generated fresh at run time from a recorded seed. Nothing is
committed as data, so there is nothing to memorise, and a reviewer reviews a
short deterministic script instead of a binary blob. Each probe runs several
seeds and all of them must pass, so no model gets through on a lucky draw.
Any run reproduces exactly with `ht run --seed <n>`.

This is an **unseen-sample** guarantee, not a claim that the public generator's
distribution is secret. For evaluation against generators or data unavailable
during training, keep a second probe tree outside the repository and add it
with `--probe-root /secure/hidden-probes`. The model container receives only
opaque case metadata and the inputs needed for inference; the host retains the
probe identity, generator seed, annotations and scoring code. See
[docs/adapting-a-model.md](docs/adapting-a-model.md#private-evaluation-suites).

## Verdicts

Binary, with the reason recorded separately, because these mean different
things:

- `VIOLATION` the model reported its budget and the budget did not close.
- `INCOMPLETE` the model never reported enough to be checked. Every
  streamflow-only model lands here today. It has not violated conservation;
  it has declined to be falsifiable.
- `INCOMPATIBLE` the model and probe disagree on timestep, required forcing or
  paired-perturbation support, so running them would not be meaningful.
- `ERROR` the adapter or benchmark machinery failed. This is operational, not
  a scientific verdict, and makes model-evaluation CI fail.

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

Without installing, `./ht` runs the CLI straight from `src/`. Set `HT_ASCII=1`
for reports with the words and no marks, which is also what you get
automatically wherever the output stream cannot carry them.

## Submitting a model

Your model may be written in any language. It ships as a container plus a
thin adapter that reads `/io/request.json` and writes `/io/output/result.csv`.
The adapter is usually thirty lines. See [docs/adapting-a-model.md](docs/adapting-a-model.md),
and `AGENTS.md` if you are having a coding agent build the sandbox for you.

To ask the maintainers to package and evaluate a model, open a
[model submission](../../issues/new?template=model_submission.yml). The form
accepts AI-based, AI+physics and physics models, with separate links for the
model code and optional pretrained weights, and asks whether there is a time
window you want the test to run over.

A submitted model is scored on the largest flood event of the generated
record rather than on all ten years of it: by default a month for a daily
model and a week for an hourly one, with the full spinup in front, located
by the probe's own reference model. That is what keeps a model that takes
seconds per forecast inside the probe's time budget. `window_days` in
`model.yaml` changes it; `full` asks for the whole record.

The container runs with no network and never sees the probe code, so a model
cannot read the tolerance it is being judged against. Every evaluation is
appended to [models/result.csv](models/result.csv).

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

Suite `0.1.0`, pre-release. Seven probes, all mass, synthetic track only. Energy,
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

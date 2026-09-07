---
name: hydroturing-evaluate-model
description: End-to-end workflow for taking a hydrologic model into the HydroTuring conservation benchmark and reporting what it did — from a model-submission issue to a packaged models/<name>/ directory (model.yaml, Dockerfile, ht_adapter.py, README.md), a Docker evaluation on every merged probe, an archived row set in models/result.csv, the README/CONTRIBUTORS updates, the issue comment and the CI run. Use this whenever someone asks to evaluate, test, package, add, benchmark or "run the probes on" a model in this repository, mentions a [MODEL] issue, asks why a model passed or failed a probe, or wants to find the next model to test. Also read it before touching the harness (src/hydroturing), the probes, the Docker runner or CI for a model: it records the decisions and traps two evaluations have already hit, so a fresh session does not rediscover them.
---

# Evaluating a model under HydroTuring

HydroTuring asks one question of any hydrologic model: does it conserve what
physics says it must? A model goes in through a small `/io` contract, every
merged probe runs it on freshly generated weather, and one verdict comes out
with its reason and the numbers under it. This skill is the maintainer's
workflow for one model, written after the first two (Google Flood Hub's LSTM,
then δHBV 2.0). Read `references/findings.md` for what those two taught, and
`references/environment.md` for how this machine and the CI are set up.

## The shape of the job

1. **Pick or receive the model.** Requests arrive as GitHub issues from the
   `model_submission.yml` template (label `model-request`, title `[MODEL] …`).
   When asked to *find* a model, prefer one that changes what the benchmark can
   say: the first two were a discharge-only LSTM (INCOMPLETE on every budget
   probe by construction) and a differentiable HBV (the first that could be
   scored on closure). Check the repository for pinned revisions, published
   weights, licence, forcing needs and CPU inference before filing. File the
   issue with `gh issue create --label model-request --assignee <user>` and
   the template's `###` headings; issue #2 is a worked example of the level of
   detail that made the packaging go quickly.
2. **Branch** `model/<name>` from `main`. Work happens there; at the end `main`
   is fast-forwarded and pushed (the maintainer has approved that flow).
3. **Read the model's inference path, not its README.** Find how the network
   is constructed from its config and weights, how forcing and attributes are
   normalised and fed, and what the *full* output dictionary holds. BMI or
   operational wrappers usually expose one variable; the adapter must build the
   model underneath them to get stores and evaporation.
4. **Prototype outside Docker first.** A scratch venv with CPU torch and the
   model's package, a script that runs the probe's own `generate(seed)` record
   through the model and prints the budget (`P − Q − ET − ΔS`) and the runoff
   CV. This is where structural surprises show up in minutes rather than after
   a 15-minute image build.
5. **Write the package** in `models/<name>/`: `model.yaml`, `Dockerfile`,
   `ht_adapter.py`, `README.md`. Copy the structure of `models/dhbv2/`
   (stores reported) or `models/google_flood_forecast/` (discharge only).
6. **Build, verify, run, archive**:
   ```bash
   docker build --quiet -t hydroturing/<name>:<version> -f models/<name>/Dockerfile models/<name>
   ht verify-adapter --model <name> --csv models/result.csv
   ht run --model <name> --gate-seeds --json results/<name>/report.json --csv models/result.csv
   ```
   `results/` is untracked scratch; `models/result.csv` is the archive and is
   committed. Drop superseded rows from an adapter you changed in the same
   session before committing; the archive records evaluations, not attempts.
7. **Diagnose every failure before reporting it.** A number in the report is a
   claim about the model only if the adapter did not cause it. For each
   failing criterion, find the mechanism in the model's code (a term in the
   equations, a learned parameter at its bound, a step-size assumption) and
   check it with a targeted script. Both evaluations so far had one failure
   that was the adapter's fault and was fixed before the verdict was recorded.
8. **Report**: README models table row, `CONTRIBUTORS.md` models table,
   commit with `Closes #N`, fast-forward `main`, push, then `gh issue comment`
   with the verdict block and the mechanism behind each failure. The pages
   workflow regenerates the probe/model count badges on any `models/**` push.
9. **Optionally re-run on a clean runner.** The `model` workflow no longer
   fires on pushes (disabled 2026-09-05: too heavy for the 4-core runner,
   and the local archive is the record). `gh workflow run model --ref main
   -f model=<name>` runs one model on demand. Compare its scorecard with
   the local one: same pass/fail pattern is the expectation, different
   numbers are normal because CI draws fresh seeds.

## Adapter rules that were not obvious until they bit

The contract is in `AGENTS.md` and `docs/adapting-a-model.md`; these are the
places where following it took judgement.

- **Report every store the model has, invent none.** Closure is differenced
  over every known state column present (`mrso`, `snw`, `canopy`, `gw`,
  `channel`), the probe's required ones first. A groundwater box under
  `gw`, water inside a unit hydrograph under `channel` (cumulative unrouted
  minus routed flow). A store the model structurally lacks is reported as
  zero *and documented as such* in the README and `run.json` — that is a
  statement about the model, not a fabricated value. A flux the model does
  not compute is never fabricated; the model is declared INCOMPLETE instead.
- **Read structure from the checkpoint when the config disagrees.** δHBV's
  shipped `config.yaml` (4 components, no routing) did not fit its weights
  (3 components, unit hydrograph). Infer from tensor shapes, record both in
  `run.json`, and tell the upstream authors.
- **Derive only attributes with an unambiguous definition.** Mean annual P
  and PET, their ratio, mean temperature, snowfall mass fraction, area:
  yes. Seasonality indices or snow-cover fractions whose training-set
  definition you cannot establish: training mean. A guessed definition that
  lands outside the training range moved δHBV's learned parameters to their
  caps and flattened its runoff — a result about the guess, not the model.
  Always run the attribute sensitivity table (all derived / each group at the
  mean / everything at the mean) and put it in the model README.
- **Climatology from the first year only.** Any statistic over the whole
  record leaks the future into earlier rows and fails `causality`.
- **Never resample.** Feed rows at the step in `request.json`; the adapter is
  the units layer (rate × dt in, flux ÷ dt out). Step dependence is what the
  resolution probe measures.
- **Determinism.** Seed the model from `request.json`'s `seed`; paired
  variants share it. Prefer an exact statistic (a mixture median by quantile
  search) over sampling so stress probes see the model, not its sampler.
- **Mock only from the forcing, with no calendar dependence**, when a model
  needs inputs the probe does not generate (radiation, pressure). Label them
  in `run.json`. Seasonal mocks masqueraded as model drift on the
  steady-state probe once.
- **Weights are baked at build time**; the container has no network. Pin a
  commit and a weights URL as `ARG`s. If a package takes its version from
  git tags, install the release wheel of the same tag rather than fighting
  hatch-vcs in a one-commit checkout.
- **Licence.** Non-commercial model licences (PSU) permit this evaluation;
  copy `LICENSE` into the image and never push the image to a registry.

## Harness facts a model evaluation depends on

- Verdict = worst reason across probes: `OK < VIOLATION < INCOMPATIBLE <
  INCOMPLETE < ERROR`. INCOMPLETE is decided by `ModelManifest.missing_for`
  against `requires` in `probe.yaml`; INCOMPATIBLE only on single-step probes.
- Submitted models are scored on a flood-event window (`window_days` in
  `model.yaml`; default 30 d daily / 7 d hourly; `full` for the whole record;
  probes may set `min_window_days`). Reference models always see the full
  record because the gate is defined on it.
- Containers run `--network none --read-only --cap-drop ALL`, as the invoking
  user (`--user uid:gid`, `HOME=/tmp`), CPUs capped at the host's count.
  Both of those were CI failures first.
- `ht validate`, `ht gate`, `ht list`, `pytest -q` must all be green before a
  push; `ht run` exit 1 is a scientific FAIL (fine), exit 2 is a harness
  ERROR (not fine).
- Every probe must pass four physical models (`reference_bucket`,
  `flex_lumped`, `flex_topo`, `sacsma_snow17`) and fail its named broken one; a probe PR is
  gated on the physical models first. The seven merged probes and what each
  catches are in the README tables; `docs/writing-a-probe.md` is the probe
  author's guide.

## Writing it up

Commit messages explain why in prose (see `git log`); no bullet lists, no
"Closes" except when a push should close an issue. The README's models table
row carries the verdict and the one-line mechanism. The issue comment carries
the verdict block verbatim plus the mechanism behind each failure and the
packaging caveats. When the finding is about the model's design (a learned
source term, a missing dt), say so plainly and say what would have hidden it.

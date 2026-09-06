# Local machine and CI

## This Mac

- No torch in the system Python. The harness runs from a scratch venv:
  `PYTHONPATH=src <venv>/bin/python -m hydroturing.cli <cmd>` (or `./ht`
  after `pip install -e '.[dev]'`). `HT_ASCII=1` for plain-text reports.
- Model prototyping: a second venv with CPU torch
  (`pip install --index-url https://download.pytorch.org/whl/cpu torch==2.5.1`)
  plus the model's package. Keep it in the session scratchpad, not the repo.
- Docker Desktop is available; a submitted model runs through it at roughly
  4–20 s per case. `ht verify-adapter` rebuilds the image (cached) on every
  call, so its reported time includes that; the probe's `max_runtime_s`
  applies to the container run only.
- Weights downloads: Google's from GitHub at build time; MHPI's from a
  public S3 URL (`--no-sign-request`, or plain HTTPS), 374 MB for δHBV daily.
- The scratch venv's tests: `PYTHONPATH=src <venv>/bin/python -m pytest -q
  -p no:cacheprovider tests/` (~20 s, 101 tests as of 2026-09-05).

## Git flow the maintainer approved

Work on `model/<name>`; `git checkout main && git merge --ff-only model/<name>
&& git push origin main`; back to the branch. Commit messages: prose, the
why, a `Co-Authored-By` trailer as instructed by the session. `Closes #N`
in the message closes the submission issue on push.

## CI (GitHub Actions, ubuntu-latest, 4 cores)

- `model` workflow: manual only since 2026-09-05,
  `gh workflow run model --ref main -f model=<name>` (blank = every model); 45-minute limit;
  builds the image without secrets, `ht verify-adapter`, then `ht run` with
  fresh seeds; exit 1 (FAIL) is green, exit 2 (ERROR) is red; scorecard in
  the log and `results/` as an artifact. δHBV on the full record took ~10
  minutes; Google on 30-day windows could not finish within 45 minutes on
  4 cores before the CPU cap fix and has not been re-run there since.
- `pages` workflow: deploys `site/` on `site/**`, `probes/**`, `models/**`
  pushes and generates `site/badges/{probes,models}.json` with
  `scripts/badges.py` for the README badges.
- `probe` workflow: `gate` job passes; `container` job has a pre-existing
  failure (unrelated to models).

## Where things live

```
models/<name>/            model.yaml, Dockerfile, ht_adapter.py, README.md
models/result.csv         the archive (committed); one row per probe per run
results/<name>/           local scratch reports (untracked)
probes/<law>/<id>/        probe.yaml, generate.py, README.md
src/hydroturing/          harness; criteria/ has one file per criterion family
tests/                    pytest; test_states.py covers extra reported stores
scripts/badges.py         probe/model counts for the site badges
.github/ISSUE_TEMPLATE/model_submission.yml   the [MODEL] form
```

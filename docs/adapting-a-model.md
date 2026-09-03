# Adapting a model

Walkthrough for wrapping an existing hydrologic model so HydroTuring can run
it. See `AGENTS.md` for the contract stated compactly, which is also what to
hand a coding agent.

## 1. Decide what your model honestly emits

This is the only decision that requires judgement. List the variables the
model genuinely produces. Do not pad the list.

A model that predicts only discharge declares:

```yaml
emits:
  fluxes: [mrro, dis]
  states: []
```

and is scored `FAIL (INCOMPLETE)`. That is the correct outcome, and it is a
different statement from `FAIL (VIOLATION)`. Inventing an evapotranspiration
column to escape `INCOMPLETE` converts an honest limitation into a false
claim, and the budget will not close anyway.

## 2. Write the adapter

Read `/io/request.json`, read the forcing, call your model, write the table.
The plumbing is about thirty lines and does not depend on your model.

```python
request = json.loads(Path(sys.argv[-1]).read_text())
io_dir  = Path(sys.argv[-1]).parent
forcing = list(csv.DictReader(open(io_dir / request["input"]["forcing"])))
static  = json.loads((io_dir / request["input"]["static"]).read_text())

rows = my_model.run(forcing, static)          # <- the only model-specific line

with open(io_dir / request["output"]["table"], "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=COLUMNS)
    writer.writeheader(); writer.writerows(rows)
```

`models/reference_bucket/ht_adapter.py` is the full worked version, standard
library only.

## 3. Write the Dockerfile

Install the model and everything it needs. There is no network at run time, so
weights, lookup tables and parameter files must be baked into the image.

```dockerfile
FROM python:3.11-slim
RUN pip install --no-cache-dir my-hydro-model==1.2.3
COPY ht_adapter.py /model/ht_adapter.py
WORKDIR /model
ENTRYPOINT []
```

The entrypoint stays empty; `model.yaml` supplies the argv.

## 4. Verify the contract before the physics

```bash
ht verify-adapter --model my-model
```

This always invokes the adapter, even when the model cannot emit enough
variables for the selected scientific probe. It asks for every output declared
in `model.yaml` and checks the row count, exact time axis, declared columns,
finite values and output-size limit. Get it green first. A residual computed
from a malformed table tells you nothing.

## The evaluation window

A submitted model is scored on a flood event, not on the full generated
record. The harness generates the whole record, asks the probe's reference
model where the largest flood of the requested length is, and hands the
submitted model that stretch with the full spinup in front of it. The
default is 30 days for a daily model and 7 days for an hourly one, which is
long enough to see the event and short enough that a model taking seconds
per forecast fits the probe's time budget. Every criterion is scored on the
window, and the report says which dates were scored for each seed.

The submission form asks whether you want a particular window. Whatever was
agreed goes in the manifest, and the command line can override it for a
one-off:

```yaml
window_days: 30        # any positive number of days, or "full"
```

```bash
ht run --model my-model --window 90
ht run --model my-model --window full
```

Two consequences for an adapter. `n_steps` in the request is the length to
emit, spinup included, and it is not the length of the full record. And a
model that needs a long history behind every prediction, as a sequence model
does, gets the probe's spinup for that purpose and nothing more: the first
rows of the record have less history behind them than the last, and the
adapter has to produce a finite value for them anyway.

Reference models always see the full record, because the acceptance gate is
defined on it. The one criterion that is climatological by construction, the
runoff ratio inside `non_degenerate`, is reported but not judged on a window
shorter than a year, since a melt flood returns more water than fell on it
that month and that is physics rather than degeneracy.

## Private evaluation suites

The probes committed to this repository are public development tests. They
generate fresh cases, so exact rows cannot be memorised, but their generating
distribution is intentionally reviewable and is not secret.

A trusted evaluator can add an uncommitted probe tree at run time:

```bash
ht validate --probe-root /secure/hidden-probes
ht run --model my-model --probe-root /secure/hidden-probes --json report.json
```

The external tree uses the same `probes/<law>/<id>/` layout. Its generator and
criteria execute only on the host. A submitted container receives read-only
forcing and static inputs, an opaque case identifier, a model-specific random
seed, and a writable output directory; it never receives the probe path,
generator seed, annotations, criterion names or tolerances. Freeze the model
image before running this private suite if the result is meant to demonstrate
generalisation beyond training.

## 5. Run it

```bash
ht run --model my-model --markdown
ht run --model my-model --json results/my-model/report.json
ht run --model my-model --gate-seeds --csv models/result.csv
```

`models/result.csv` is the archive of every evaluation: one dated row per
probe, plus one for the adapter contract check when `verify-adapter` is
given the same `--csv`. For a model that reports only discharge, the
contract row is the only line saying it was actually built and run, because
its scientific verdict is INCOMPLETE before the container is started.

## Common failures

**"result has N rows, expected M"** — the spinup was trimmed. Emit one row per
input row; the harness slices the window itself.

**`forcing_fidelity` fails** — the reported `pr` does not match the input,
usually a unit conversion. Forcing is mm/day.

**`state_bounds` fails but `closure` passes** — the model is computing a
storage term as the budget residual. This is the pattern the benchmark exists
to detect, and it is worth checking whether it is deliberate.

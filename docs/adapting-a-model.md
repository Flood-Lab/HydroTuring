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

This runs one seed and checks the shape of the result: right number of rows,
requested columns present, no non-finite values. Get it green first. A
residual computed from a malformed table tells you nothing.

## 5. Run it

```bash
ht run --model my-model --markdown
ht run --model my-model --json results/my-model/report.json
```

## Common failures

**"result has N rows, expected M"** — the spinup was trimmed. Emit one row per
input row; the harness slices the window itself.

**`forcing_fidelity` fails** — the reported `pr` does not match the input,
usually a unit conversion. Forcing is mm/day.

**`state_bounds` fails but `closure` passes** — the model is computing a
storage term as the budget residual. This is the pattern the benchmark exists
to detect, and it is worth checking whether it is deliberate.

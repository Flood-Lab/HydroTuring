# Building a HydroTuring sandbox

Instructions for a coding agent asked to make an existing hydrologic model
runnable under HydroTuring. The work is always the same three files, and it
does not require understanding the model's internals.

## What you are building

```
models/<model-name>/
  model.yaml      declares what the model emits and needs
  Dockerfile      installs the model and the adapter
  ht_adapter.py   translates between the /io contract and the model's own API
```

## The contract

The harness exposes a small `/io` filesystem and runs your entrypoint once:

```
<entrypoint> --request /io/request.json
```

Paths inside `request.json` are relative to the request file's directory.

```
/io/request.json          read-only: opaque case id, model seed, timestep, n_steps, outputs
/io/input/forcing.csv     read-only: columns time, pr, tas, pet (mm/day, degC, mm/day)
/io/input/static.json     read-only: catchment attributes
/io/output/result.csv     write: one row per forcing row, spinup included
/io/output/run.json       write: {"status": "ok"}
```

The model seed is deterministic for reproducible stochastic inference, but it
is not the generator seed recorded in the host-side report. Exit 0 on success.
There is no network. Do not attempt to download weights or data at run time;
bake them into the image.

## The evaluation window

A submitted model is not run over the full generated record. The harness
cuts the record down to the largest flood event, located by the probe's own
reference model, with the full spinup in front of it: by default 30 days for
a daily model and 7 days for an hourly one, so a heavy model fits the probe's
time budget. `n_steps` in `request.json` is the length the adapter must
emit, spinup included; it is 395 rows for a daily model on the default
window, not the 4015 of the full record. Ask for a different window, or the
whole record, in `model.yaml`:

```yaml
window_days: 30      # or any positive number of days, or "full"
```

A probe may set a minimum window when its expectation only holds over a
longer stretch; the warming-response probe asks for a year. Reference
models always see the full record, because the acceptance gate is defined
on it.

## The timestep

`request.json` names the step the case runs at, as an ISO 8601 duration:
`PT1D`, `PT1H`, `PT15M`, `PT5M` or `PT1M`. Fluxes are rates in mm per day at
every step, so the depth moved in one row is the rate times the step
length, and a one-millimetre burst in one minute arrives as 1440 mm/day.
`model.yaml` declares the step the model was built at, and may list others
it also runs at, native step first:

```yaml
timestep: PT1D             # or a list: [PT1H, PT1M]
```

A resolution probe serves the same weather at several steps and runs every
model at its own step and at the next finer one, whatever the manifest
lists. It then measures how far the integrated runoff moved, in percent of
the rain. A model whose arithmetic assumes its native step moves a lot.

**Feed the model the rows at the step you were given.** Do not aggregate
hourly rows to days inside the adapter, or split days into hours, to reach
the step the model prefers. The adapter is the model's units layer, not a
resampler: it turns rates into the depths the model wants for the step it
was handed, and resampling would hide exactly the dependence on the step
that the probe exists to measure. On a single-step probe the manifest's
step is respected: a daily closure probe given an hourly-only model reports
INCOMPATIBLE, because closure is what it measures and it cannot measure it
on a model it cannot feed.

## Variable names and units

| Name | Meaning | Units |
| --- | --- | --- |
| `pr` | precipitation, echoed back from the forcing | mm/day |
| `evspsbl` | evapotranspiration | mm/day |
| `mrro` | total runoff | mm/day |
| `dis` | river discharge | m3/s |
| `mrso` | soil water storage | mm |
| `snw` | snow water equivalent | mm |
| `canopy` | canopy interception storage | mm |

## Three rules that are easy to get wrong

1. **States are absolute, never tendencies.** Report the storage itself at
   each step. The harness differences it. If you report `dS/dt` you have
   removed the quantity the benchmark checks, and the model will be scored as
   if it invented its storage.

2. **Echo `pr` exactly as given.** The harness compares your reported
   precipitation against the forcing it handed you. Rescaling it, even for
   unit reasons, fails `forcing_fidelity`.

3. **Emit one row per input row, spinup included.** Do not trim the spinup
   yourself. The harness knows where the window starts and slices it.

## Emit only what the model actually produces

Do not invent values to fill columns. A model that genuinely produces only
discharge should declare exactly that:

```yaml
emits:
  fluxes: [mrro, dis]
  states: []
```

It will be scored `FAIL` with reason `INCOMPLETE`, which is the honest
outcome. Fabricating an `evspsbl` column to avoid `INCOMPLETE` produces
`VIOLATION` instead, which is worse and is also dishonest.

## Verify before you submit

```bash
ht verify-adapter --model <model-name>   # contract only, no physics
ht run --model <model-name>              # the actual evaluation
```

`verify-adapter` runs a single seed, on the same window the evaluation will
use, and checks the shape of what came back. Get that green before looking
at any residual. Both commands take `--csv models/result.csv` to append what
they found to the archive, and `--window DAYS|full` to override the
manifest.

## A minimal adapter

`models/reference_bucket/ht_adapter.py` is a complete worked example in the
standard library alone. Copy its structure. The model-specific part is the
`simulate` function; everything around it is contract plumbing that does not
change.

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

The harness mounts a directory at `/io` and runs your entrypoint once:

```
<entrypoint> --request /io/request.json
```

Paths inside `request.json` are relative to the request file's directory.

```
/io/request.json          read: case id, seed, timestep, n_steps, requested variables
/io/input/forcing.csv     read: columns time, pr, tas, pet (mm/day, degC, mm/day)
/io/input/static.json     read: catchment attributes
/io/output/result.csv     write: one row per forcing row, spinup included
/io/output/run.json       write: {"status": "ok"}
```

Exit 0 on success. There is no network. Do not attempt to download weights or
data at run time; bake them into the image.

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

`verify-adapter` runs a single seed and checks the shape of what came back.
Get that green before looking at any residual.

## A minimal adapter

`models/reference_bucket/ht_adapter.py` is a complete worked example in the
standard library alone. Copy its structure. The model-specific part is the
`simulate` function; everything around it is contract plumbing that does not
change.

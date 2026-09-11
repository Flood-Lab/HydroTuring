# Building a HydroTuring sandbox

Instructions for a coding agent asked to make an existing hydrologic model
runnable under HydroTuring. A maintainer's agent evaluating a submission
end to end should start from `.claude/skills/hydroturing-evaluate-model/`,
which records the whole workflow and what the first evaluations taught. The work is always the same three files, and it
does not require understanding the model's internals.

Before any of it: contributions here start with an issue, not a pull request.
If there is no accepted [model proposal](https://github.com/Flood-Lab/HydroTuring/issues/new?template=model_submission.yml)
for this model, open one and wait for it to be labelled `accepted` and
assigned. Then fork, build, and open the pull request from the fork. See
[CONTRIBUTING.md](CONTRIBUTING.md).

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
/io/input/forcing.csv     read-only: columns time, pr, tas, pet (mm/day, degC, mm/day),
                          and, when a probe prescribes a human withdrawal, abstr (mm/day, net)
/io/input/static.json     read-only: catchment attributes
/io/output/result.csv     write: one row per forcing row, spinup included
/io/output/run.json       write: {"status": "ok"}
```

The model seed is deterministic for reproducible stochastic inference, but it
is not the generator seed recorded in the host-side report. Exit 0 on success.
There is no network. Do not attempt to download weights or data at run time;
bake them into the image.

Some probes prescribe a human withdrawal in the forcing as an `abstr` column
(mm/day, net of return flow). Honouring it means removing that water from the
stores and fluxes the model reports, and declaring whatever was actually
removed as a negative `gwex`. A model that never reads the column reports a
budget that closes on its own yet misses the withdrawal; the probe scores that
as a failure, not as INCOMPLETE, because such a model reports everything the
criterion needs.

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

**Put learned parameters in the units of the step.** A model whose network
writes physical parameters (a recession coefficient per day, a percolation
cap in mm per day, a unit hydrograph in days) has to have those parameters
rescaled to the step the case runs at before its physics uses them, or the
physics drains its stores as many times too fast as there are steps in a
day, and the probe measures the adapter's units rather than the model. A
per-day fraction k becomes 1 - (1 - k)^dt, a per-day amount becomes
amount × dt, a hydrograph keeps its shape in days. The networks themselves
are fed what they were trained on, rates in mm per day; the physics is fed
depths per row. `models/dhbv2/ht_adapter.py` does this.

**Feed the model the rows at the step you were given.** Do not aggregate
hourly rows to days inside the adapter, or split days into hours, to reach
the step the model prefers. The adapter is the model's units layer, not a
resampler: it turns rates into the depths the model wants for the step it
was handed, and resampling would hide exactly the dependence on the step
that the probe exists to measure. On a single-step probe the manifest's
step is respected: a daily closure probe given an hourly-only model reports
INCOMPATIBLE, because closure is what it measures and it cannot measure it
on a model it cannot feed.

## Nothing from the future

The stress probes add one storm to a record and check that nothing the
model reports changes before it. Any path from later rows to earlier
output fails that: a bidirectional pass, a centred filter, attention over
the whole record, or statistics taken over the record you were handed,
such as normalising inputs by their mean or deriving catchment attributes
from all of it. If the model needs a climatology, take it from the first
year of the record and nothing later. Paired variants of one seed carry
the same `seed` in the request, so a stochastic model draws the same
numbers in both runs; keep it that way and do not reseed from the clock.

## Variable names and units

| Name | Meaning | Units |
| --- | --- | --- |
| `pr` | precipitation, echoed back from the forcing | mm/day |
| `evspsbl` | evapotranspiration | mm/day |
| `mrro` | total runoff | mm/day |
| `dis` | river discharge | m3/s |
| `gwex` | a declared exchange with the outside: regional groundwater, inter-basin transfer; positive into the catchment | mm/day |
| `sbl` | the sublimating share of `evspsbl`: a component of it, never an addition; report it if the model knows which kilograms left as ice | mm/day |
| `hfls` | latent heat flux, positive away from the surface | W/m2 |
| `hfss` | sensible heat flux, positive away from the surface | W/m2 |
| `hfg` | ground heat flux at the actual soil surface, positive into the ground; a flux taken below the surface must be corrected for heat storage above that depth | W/m2 |
| `mrso` | soil water storage | mm |
| `snw` | snow water equivalent | mm |
| `canopy` | canopy interception storage | mm |
| `gw` | groundwater storage below the soil column | mm |
| `channel` | water generated as runoff but not yet released by the model's routing | mm |

For `hfg`, an adapter mapping a plate-depth or deeper-boundary flux must use
`G_surface = G_depth + (E_above_end - E_above_start) / dt`, with downward
fluxes positive, `E_above` in J/m2 and `dt` in seconds. The storage and fluxes
must cover the same area, layer and time interval; this relation assumes no
other energy sources or sinks in that layer. Use the model's actual heat
storage, not a value inferred from the surface-budget residual. Document the
mapping and any unavailable terms. Once `hfg` is mapped to the surface, do
not subtract subsurface heat storage again in the surface budget.

A probe names the stores it requires. Report every store the model
actually has, including ones the probe did not name: a groundwater zone
under `gw`, water in transit through a unit hydrograph or a channel under
`channel`. The budget is closed over all of them. Leaving one out opens the
budget by exactly the water it holds, and a model with a groundwater store
that reports only soil water will be measured as leaking. A store the model
does not have at all is not fabricated to fill a column; say so in the
adapter's README and in `run.json`.

A model with an explicit exchange with the outside, a regional groundwater
term or an inter-basin transfer, reports it as the flux `gwex`, positive
into the catchment. Declared, it is a source in the budget and the budget
can close; hidden, it is the residual, and the closure probe reports the
model as leaking or making water by exactly that term. The runoff ratio is
taken against rain plus the declared source. `models/dhbv2` does this for
HBV 2.0's regional groundwater term.

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

## When a probe or a model merges: what has to move with it

This section is for whoever lands the pull request, human or agent. The
harness discovers a probe from its directory and a model from its manifest,
but six other places describe the suite in prose and none of them are
generated. `tests/test_docs_in_sync.py` fails when any of them is behind,
and the probe workflow runs it on every pull request, so a merge that skips
this list is caught before it lands rather than noticed a week later.

For a merged probe:

1. `README.md`: a row in the probes table, and a row in the reference-models
   table for every reference model the probe adds. If the probe requires a
   variable the physical models do not report, their standing changes from
   "PASS, N of N" to "N of N+1, INCOMPLETE on ..." and the sentence above the
   table that says what they must pass changes with it.
2. `CONTRIBUTORS.md`: a row in the probes table naming the author with
   their affiliation, as `Name (Institution)`. A merged probe earns
   co-authorship, so this row is the record of that. The affiliation comes
   from `authors` in the probe's `probe.yaml`; if it is missing there, ask
   the author before merging rather than after, because the test that
   checks it (`test_probe_authors_carry_an_affiliation`) fails the build.
3. `CITATION.cff`: an `authors` entry for every probe author not already
   listed, with `given-names`, `family-names`, `affiliation` and `orcid` as
   a full `https://orcid.org/` URL. This is the software's author list and
   the paper's starting point. Contact emails stay out of the repository;
   ask for one privately at merge time.
4. `ROADMAP.md`: the entry moves from unclaimed to `**merged**`, under the
   id the probe actually took, with the author named. Rewrite the paragraph
   above it if it counted the unclaimed entries.
4. `site/index.html`, three times, once per language block: the row in
   `probes.rows` moves from the wanted block to the merged block under its
   real id; every `N / M` in `models.rows` takes the new probe count; the
   flowchart under "How it works" gains a labelled entry in the right pillar,
   with an `infra.probe.<key>` translation in each language, measured against
   the pillar's width. The badge count is generated and needs nothing.
5. `models/result.csv`: one row per evaluated model on the new probe, written
   by `ht run --model <name> --probe <id> --gate-seeds --csv models/result.csv`
   rather than by hand. INCOMPLETE is a verdict and is archived like any
   other.
6. This file, if the probe introduced a variable: a row in the table above.
   `spec.py` accepts the name the moment it is in `FLUX_VARS`; nothing tells
   an adapter author it exists except this table.

For a merged model: the README models table, the CONTRIBUTORS models table,
the site's `models.rows` in all three languages, and the archive rows from
`ht run` on every probe. The pages workflow redeploys the site on any push
that touches `site/`, `probes/` or `models/`.

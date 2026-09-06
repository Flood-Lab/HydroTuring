<!-- Title: [MODEL] <name> -->

Closes #<proposal issue>

## The model

<!-- What it is, which revision, and a link to the code and any weights. -->

| | |
| --- | --- |
| Category | <!-- AI-based / AI+physics / physics --> |
| Code revision | |
| Weights | |
| Licence (code, weights) | |

## What it honestly emits

```yaml
emits:
  fluxes: []
  states: []
```

<!-- If this list is short, say so plainly. A model that reports discharge
     alone is scored FAIL (INCOMPLETE), which is an honest result and one
     worth recording. Do not pad the list to escape it. -->

## Result

<!-- Paste `ht run --model <name> --markdown`. A FAIL is a legitimate
     outcome and is not a reason to hold the PR back. -->

## Checklist

- [ ] There is an `accepted` proposal issue and this PR closes it
- [ ] `ht verify-adapter --model <name>` passes
- [ ] `ht run --model <name>` completes, whatever the verdict
- [ ] `emits` lists only what the model genuinely produces
- [ ] States are absolute storages, not tendencies
- [ ] The forcing is echoed back exactly as given, with no unit conversion
- [ ] One row per input row, spinup included
- [ ] Everything the model needs is baked into the image; it runs with no network
- [ ] The licence permits distributing it this way, and `model.yaml` names it

## Conflicts of interest

<!-- If you are an author of this model, say so. It is common, it is not a
     problem, and it means a second reviewer is assigned. See GOVERNANCE.md. -->

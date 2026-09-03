# {id}

## What it asserts

<!-- The residual equation, with units. -->

```
R = ...        [units]
```

Accept when the cumulative residual is under **5 percent of <driver>**. Say
what the denominator is and whether it can approach zero. If it can, say what
absolute floor you set and why.

## The case

`generate.py` produces, from a seed:

<!-- What the forcing looks like, with rough annual totals. If this is a
     catchment probe, say where it sits on the Budyko curve, because a
     reviewer will check. -->

<!-- Say how much spinup precedes the scored window. -->

## Why these criteria

Closure alone is trivially satisfiable. State what each additional criterion
is closing off.

| Criterion | The cheat it closes off |
| --- | --- |
| `closure` | |
| `state_bounds` | |
| `non_degenerate` | |

<!-- The question a reviewer will ask: describe a model that passes closure
     while doing no physics, and name the criterion that catches it. -->

## Reproducing a failure

Every run records its seed. To re-run one:

```bash
ht run --model <name> --probe {id} --seed <n>
```

# mass/catchment-closure

## What it asserts

Over a ten year record of generated daily weather, a lumped catchment's water
budget must close:

```
R = pr - evspsbl - mrro - d(mrso + snw + canopy)/dt        [mm/day]
```

Accept when the cumulative residual is under **5 percent of total
precipitation**. Precipitation is strictly non-negative, so the denominator
never approaches zero and no absolute floor is needed.

The reference budget is exact by construction. The model is judged against the
forcing it was handed, never against observations, so the large closure errors
in published hydrologic datasets do not enter this probe at all.

## The case

`generate.py` produces, from a seed:

- seasonal precipitation occurrence with gamma-distributed wet-day depths,
  about 830 mm per year
- an annual temperature cycle with a persistent AR(1) weather anomaly, giving
  a real seasonal snowpack, roughly a fifth of days below freezing
- temperature-driven potential ET with a daylength factor, about 765 mm per
  year

Aridity index is near 0.92, so the catchment sits in the middle of the Budyko
curve where both water-limited and energy-limited behaviour appear. The
reference model partitions at ET/P = 0.59, which is where Budyko puts it.

One year of spinup precedes the scored window. Storage change is measured from
the state at the last spinup step, not the first scored step.

## Why five criteria and not one

Closure alone is trivially satisfiable. Each additional criterion exists
because a specific cheat defeats the ones before it.

| Criterion | The cheat it closes off |
| --- | --- |
| `closure` | a silent sink, the honest failure mode |
| `state_bounds` | solving for storage as whatever balances the budget |
| `et_plausible` | evaporating all precipitation so nothing runs off |
| `non_degenerate` | any remaining trivial partition, or ignoring the forcing |
| `forcing_fidelity` | closing a budget against a privately rescaled driver |

`state_bounds` is the one that matters most. A model that computes storage as
the budget residual has a closure error of exactly zero on every seed, so
neither a tighter tolerance nor fresher random forcing will ever catch it. Its
invented storage is not a storage, and over ten years it drifts far outside
what a 320 mm soil column can hold. `reference_cheater` demonstrates this and
is required to fail here.

## Reproducing a failure

Every run records its seed. To re-run one:

```bash
ht run --model <name> --probe mass/catchment-closure --seed <n>
```

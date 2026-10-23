# mass/interception-response

A larger canopy intercepts more water. That water has to go somewhere: less
runoff, more evaporation, or a reported store. The probe runs the same
three-year weather twice, changing only `static.canopy_capacity_mm` from 2 mm
to 10 mm, and scores where the difference went.

## The residual

Write `Delta X = X(capacity_high) - X(capacity_low)` over the scored window,
with fluxes integrated to depths in mm:

```
R = Delta evspsbl + Delta mrro + Delta gwex
    + Delta (mrso + snw + canopy + gw + channel)        [mm]
```

A model that conserves water has `R = 0`. The capacity change is
`Delta cap = 8 mm`, and the four checks are:

1. `max(0, Delta mrro) / Delta cap <= 0.05` — a larger canopy may not raise
   runoff;
2. `max(0, -Delta evspsbl) / Delta cap <= 0.05` — it may not lower
   evaporation;
3. `abs(R) / Delta cap <= 0.05` — the paired budget must account for the
   difference, using every store the model reports;
4. `Delta evspsbl / Delta cap >= 0.5` — the response must exist and must be an
   evaporation response.

The absolute residual floor is 1 mm, so the effective limit is
`max(0.05, 1 mm / Delta cap) = 12.5%` at `Delta cap = 8 mm`. The floor is
there so a numerical residual on a multi-year record does not bind; the
measured exact bucket closes to better than 0.01 mm.

Criterion 4 is what makes this more than a closure check. A model that
declares `canopy_capacity_mm` and ignores it produces two identical runs:
criteria 1-3 pass trivially and 4 fails. The floor is on evaporation rather
than on the canopy-store peak because a model that routes interception loss
directly to evaporation is physically sound and should not be failed for not
exposing a separate store.

## Why 2 -> 10 mm

10 mm is a closed-canopy capacity and is representable by every model that
currently declares the static; LISFLOOD's `lai_for_canopy_capacity`
conversion tops out near 11.7 mm, and 20 mm raises before the native run. The
change is still a factor of five. Native CWatM 1.11, run outside Docker at
2 -> 20 mm on an earlier dry-down seed, moved +63.32 mm of ET and -62.84 mm of
runoff with a 0.784 mm paired residual (reported by @cehw in #161), which is
independent evidence that the design asks a real question. Its row on this
probe is generated through the Docker adapter, not by hand.

## Measured on the gate seed

Seed `1135680359`, 1,095 scored days (365 spinup days in front of them):

| model | Delta sum evspsbl | Delta sum mrro | paired residual |
| --- | ---: | ---: | ---: |
| `reference_bucket` | +82.18 mm | -104.91 mm | 0.00 mm |
| `flex_lumped` | +144.39 mm | -151.57 mm | 0.00 mm |
| `flex_topo` | +180.58 mm | -180.64 mm | 0.00 mm |

The response magnitudes differ by more than a factor of two, which is why
criterion 4 is a floor and not a magnitude test. All three pass on every gate
seed.

## Baselines

`must_pass`: `reference_bucket`, `flex_lumped`, `flex_topo`.

`must_fail`, one control per criterion:

| control | fault | caught by |
| --- | --- | --- |
| `reference_interception_blind` | declares the capacity and ignores it | criterion 4 (zero response) |
| `reference_interception_wetter` | larger capacity raises runoff | criterion 1 |
| `reference_interception_drier` | larger capacity lowers evaporation | criterion 2 |
| `reference_interception_leaky` | returns part of the difference as evaporation and drops the rest | criterion 3 |

Each is listed under `baselines.must_fail_only`, so the gate requires it to
trip `interception_response` alone rather than a contract accident.

`cwatm` and `summa` also declare `canopy_capacity_mm` and report `canopy`;
they are candidates for `must_pass` once their Docker rows pass the gate.
`sacsma_snow17` and `dhbv2` report `canopy` as an identically zero state but
do not declare the static, so they are `N/A (INCOMPATIBLE)` and their
standing does not move.

## What this probe does not catch

- A model that does not declare `canopy_capacity_mm` is not asked; the static
  is the opt-in.
- A model whose conduction of the parameter is physically wrong but still
  raises evaporation and closes the budget can pass. The probe scores
  allocation, not the internal interception formulation.
- The case fixes 2 -> 10 mm. A model that cannot represent 10 mm is out of
  scope until a second pair is agreed; the first gate does not silently score
  it at a different range.

## Reproducing

```bash
ht run --model reference_bucket --probe mass/interception-response
ht run --model reference_interception_leaky --probe mass/interception-response
ht gate --probe mass/interception-response
```

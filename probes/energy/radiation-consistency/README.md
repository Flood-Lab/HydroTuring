# energy/radiation-consistency

## What it checks

For a uniform, opaque, snow-free gray surface with fixed emissivity, total
upward longwave must equal surface emission plus reflected downward longwave:

```
r_i = rlus_i - [eps * sigma * ts_i^4 + (1 - eps) * rlds_i]
|r_i| <= max(0.005 * |rlus_i|, 0.5 W m-2)
```

`ts` is skin temperature in K; `rlus` and `rlds` are upward and downward
longwave in W m-2. `sigma` is the Stefan-Boltzmann constant. The criterion
reads `eps` from the supplied `static.json`, never from model output.

Every scored step must pass. Opposite errors cannot cancel, and interval
means are out of scope: mean(T^4) is not mean(T)^4. The relative tolerance
uses the **reported** upward flux. `floor_steps` counts steps where the
absolute floor governs; none occur in the reference cases.

The 0.5 percent tolerance and 0.5 W m-2 floor follow
[proposal #22](https://github.com/Flood-Lab/HydroTuring/issues/22).
They are consistency thresholds, not observational uncertainty. For example,
at 290 K, `eps = 0.98` and `rlds = 300 W m-2`, reflected longwave is 6 W m-2,
about 1.5 percent of total upward flux. Omitting it passes a 5 percent bound
but fails this one.

## Case and controls

The generator supplies 768 hourly rows starting at local-solar 06:00:
two spinup days and 30 scored days. `min_window_days: 30` preserves all
scored day/night cycles for submitted models, including dry hours outside
the default seven-day flood window.

`tas`, `rlds`, `ts` and `rlus` are instantaneous. `pr`, `pet` and `rn` are
means over the following hour and drive the reference water/energy budgets.
Net radiation is prescribed and is **not reconciled with the longwave
components**.

Air temperature spans 16-30 degC with a daily weather offset and a 10 K
diurnal range. Hourly cloud transmission (0.45-1.0) both dims the sun and
raises daily clear-sky emissivity (0.74-0.82) towards one. Downward longwave
follows from that sky emissivity and air temperature; the scan below spans
297-456 W m-2. Surface emissivity is drawn once per seed in [0.95, 0.99],
rounded to four decimals, and held fixed.

| Model | Construction | Expected result |
| --- | --- | --- |
| `reference_radiative` | `Ts = Ta + H / g_H`, with `g_H = 20 W m-2 K-1`; longwave uses this Ts, supplied eps and reflected sky | PASS |
| `reference_air_emitter` | Reports the same Ts but emits at air temperature; reflection unchanged | FAIL on `radiative_identity` |
| `reference_no_reflection` | Reports surface emission alone as total upward longwave | FAIL on `radiative_identity` |

The positive control retains `reference_coupled`'s water and energy columns
exactly. It treats interval-mean sensible heat as holding at the sampling
instant: an explicit approximation, not a solved instantaneous energy
balance. Its radiation identity holds by construction; the diagnosed
temperature is approximate. Each negative control changes the radiation
mode; tests verify that only `rlus` differs from the positive control.

## Limits

A wrong Ts paired with longwave computed from that Ts can pass. This checks
output consistency, not temperature accuracy, energy balance or thermal
inertia. Adapters must expose native model outputs, not manufacture `rlus`
from `ts`. Missing either output yields INCOMPLETE, as for every archived
physical and submitted model on this probe.

Around 290 K under a 300 W m-2 sky, temperature differences below about
0.4 K can pass; resolution is coarser at warmer temperatures. Unit tests
also document a linearised Stefan-Boltzmann law: around 288 K, excursions
of 5 K pass while 10 K fails, with a neglected second-order term near
2.8 W m-2. That optional control is not part of the gate.

## Reproduce

Run from the repository root after installing HydroTuring:

```bash
ht validate
ht gate --probe energy/radiation-consistency
ht run --model reference_radiative --probe energy/radiation-consistency --seed 0
ht run --model reference_no_reflection --probe energy/radiation-consistency --seed 0
pytest -q tests/test_radiation_consistency.py tests/test_radiation_probe.py tests/test_radiation_references.py
python3 scripts/radiation_margins.py
```

The first model must pass; the second must fail on `radiative_identity`.
The gate makes 15 adapter invocations across five seeds. The scan adds
seeds 0-19, seed 36 (eps 0.9891, highest among seeds 0-49), and the five
gate cases forced to each emissivity endpoint: 36 cases in total.

Ratios below are absolute residual divided by tolerance. The positive
column gives the largest ratio in the group. Negative columns give the
**minimum across cases of each case's maximum ratio**; the design target
is at least 1.3 for both negative controls, with the criterion threshold
unchanged at 1.

| Group | Cases | Positive maximum | Air emitter | No reflection |
| --- | --- | --- | --- | --- |
| gate | 5 | 3e-14 | 45.7 | 3.10 |
| additional | 20 | 3e-14 | 44.5 | 3.12 |
| drawn eps 0.9891 | 1 | 3e-14 | 49.9 | 2.15 |
| forced eps 0.99 | 5 | 3e-14 | 46.1 | 1.95 |
| forced eps 0.95 | 5 | 3e-14 | 44.7 | 10.2 |

The air emitter violates 679-690 of 720 scored steps per case: every night
step and all daytime steps except some 06:00 and 17:00 rows near Ts = Ta.
No reflection violates all 720 steps. At eps 0.99 its ratio is
`2 * rlds / rlus`: per-case maxima are 1.95-1.98, while the minimum over
**all individual steps** is 1.20. Skin temperature stays within 287-320 K
and never freezes. The script also checks unchanged water/energy columns
and that negative controls differ only in `rlus`; an optional output path
writes the per-case results as JSON.

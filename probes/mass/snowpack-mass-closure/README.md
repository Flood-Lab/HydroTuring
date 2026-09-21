# mass/snowpack-mass-closure

## What it asserts

The probe closes the water balance over the snowpack control volume. For each
contiguous regime stretch \(k\), the residual is

$$
R_k =
\sum_{t \in k}
\left(pr_t - snm_t - sbl_t\right)\Delta t
-
\left(snw_{k,\mathrm{end}} - snw_{k,\mathrm{start}}\right).
$$

Here, `pr` is precipitation supplied by the forcing, `snm` is liquid water
delivered by the model's snow module to the underlying hydrologic system,
`sbl` is reported snow sublimation, and `snw` is total snowpack water storage.
The probe uses a daily timestep, so \(\Delta t = 1\) day.

If a model does not report `sbl`, sublimation is taken as zero. Any real
unreported sublimation therefore appears as a snowpack mass residual rather
than being inferred by the benchmark.

The generator labels six contiguous `_regime` stretches:

`accumulation -> storage -> melt -> accumulation -> storage -> melt`.

`_regime` is a host-side annotation used by the criteria and is not supplied
to the model as a forcing variable.

Closure is evaluated separately for each contiguous stretch. The absolute
stretch residuals are then accumulated so that water lost in one stage cannot
be cancelled by water created in another:

$$
\epsilon =
\frac{\sum_k |R_k|}
     {\sum_t pr_t \Delta t}
\le 0.05.
$$

The denominator is total precipitation over the full experiment. The synthetic
case is designed to provide precipitation during both accumulation stages so
that the snowpack response and closure checks are meaningfully exercised.

## The case

`generate.py` creates two seasonally aligned synthetic snow cycles from a
recorded random seed. Each cycle contains three 120-day stages.

During **accumulation**, temperature remains well below the prescribed
snow/rain threshold and precipitation events build the snowpack.

During **storage**, temperature remains below freezing and precipitation is
zero. Any decrease in `snw` must be accounted for by liquid water leaving the
snow module as `snm` and/or by reported sublimation as `sbl`.

During **melt**, temperature remains well above freezing and precipitation is
zero. The stored snow is expected to be strongly depleted, primarily through
liquid water delivered to the underlying hydrologic system.

Storm occurrence, storm magnitude, and temperature variability change with
the seed, while the stage structure is fixed. The same seed always produces
the same forcing.

A single cold, precipitation-free spinup day precedes the scored record. The
first scored accumulation stage begins on 1999-09-01. Each cycle occupies approximately one year, so repeated accumulation and melt stages occur at comparable times of year.

The synthetic case sets `canopy_capacity_mm = 0` for models that consume this
case attribute, so the repository reference bucket has no canopy interception
in this experiment. The closure formulation assumes that forcing precipitation
enters the tested snow-system control volume without an upstream interception
store; models with upstream snow interception require an extended control
volume and are outside the current probe scope unless explicitly supported.

## Why these criteria

| Criterion | The cheat it closes off |
| --- | --- |
| `closure` | Water disappears from or is created inside the snowpack, including water that bypasses snow storage while the whole-catchment budget can still close. |
| `state_bounds` | A model obtains apparent closure by allowing physically impossible negative snow storage. |
| `snowpack_response` | A trivial model keeps `snw = 0` and passes precipitation directly through the snow module, satisfying closure without representing snow accumulation and melt. |

`state_bounds` requires

$$
snw \ge 0.
$$

For each accumulation-plus-storage stretch, `snowpack_response` requires

$$
\frac{snw_{\mathrm{peak}}}
     {P_{\mathrm{accumulation+storage}}}
\ge 0.80.
$$

At the end of the following melt stage, it requires

$$
\frac{snw_{\mathrm{melt,end}}}
     {snw_{\mathrm{peak}}}
\le 0.1.
$$

The first condition requires a material snowpack to form during the cold
stages. The second requires at least 90% of that peak snow storage to be
depleted by the end of the following melt stage.

The 0.80 floor sits in a window with a hard edge on each side. Above it,
`sacsma_snow17` reached a minimum peak fraction of 0.9084 over 200 seeds, and
`reference_bucket` is at 1.0 by construction. Below it, `reference_snow_bypass`
retains `1 - BYPASS_FRACTION` of every snowfall exactly, so its peak fraction
is the constant 0.85 and it must stay above the floor: that model is declared
to fail on `closure`, and a floor above 0.85 would make it fail
`snowpack_response` as well, so it would no longer isolate the criterion it
exists to demonstrate.

What the floor does not close off is a model that reports a large share of
cold-stage precipitation as `snm` and stores the rest. Such a model conserves
mass — the water it passes through is declared, so the snowpack budget closes
— and the probe will pass it as long as the retained share clears 0.80.
Delivering a fifth of midwinter snowfall to the ground as liquid water is
physically doubtful, but it is not a conservation violation, and this probe
measures conservation. The energy side of that question belongs to
`energy/snowpack-cold-content`.

`reference_snow_bypass` conserves total catchment water but sends 15% of
snowfall directly to the soil without reporting that bypass as snow-module
outflow. It therefore preserves the whole-catchment water balance while
violating the internal snowpack balance, and must fail `closure`.

`reference_snowless` reports `snm = pr` and keeps `snw = 0`. Its snowpack
budget therefore closes algebraically, but it never develops a snowpack and
must fail `snowpack_response`.

## Reproducing a failure

Every run records its seed. To reproduce a `reference_snow_bypass` failure:

```bash
ht run --model reference_snow_bypass --probe mass/snowpack-mass-closure --seed <n>
```
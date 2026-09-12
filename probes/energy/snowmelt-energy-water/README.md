# energy/snowmelt-energy-water

## What it asserts

Melting ice is not free. Every kilogram of snow converted to liquid consumes
the latent heat of fusion, `lambda_f = 3.337e5 J kg-1`.

A model can close its water budget and close its energy budget while never
charging one for what the other spent. The common construction is a
**degree-day melt head** — melt is a function of air temperature — bolted to an
**energy head** that closes `Rn = H + LE + G` by itself. Each ledger balances.
The pack melts because the air is warm, whether or not the energy to melt it
existed. No single-budget criterion can see this, because within either budget
the numbers are consistent.

Over the melt block the surface residual is therefore **not** required to
vanish. It is required to be the melt energy:

```
M       = -(snw_end - snw_before) - sum_i(sbl_i * dt_i)      [mm == kg m-2]
R_avail = mean_i(Rn_i - LE_i - H_i - G_i)                    [W m-2]
R_demand = lambda_f * M / (N * dt * 86400)                   [W m-2]

D = abs(R_avail - R_demand) <= max(0.05 * R_demand, 2 W m-2)
```

`snw_before` is the pack one step **before** the block opens, taken through the
shared `storage_at` helper, so the change is never measured from the block's
own first row.

## The budget boundary

The residual is written at the **snow surface**, and every term must be
reported there.

`Rn` (`rn`) is positive into that surface. `LE` (`hfls`), `H` (`hfss`) and `G`
(`hfg`) are positive away from it, and any of the three may be negative. Over a
melting pack `H` is normally negative: the air is warmer than a surface held at
0 °C, so the turbulent flux is energy arriving, and it melts ice.

`G` is the conduction **into the pack** across that surface, not into the soil
beneath it. The distinction is not cosmetic. A deep pack insulates the ground,
so soil heat flux and snow-surface conduction are different quantities, and a
model reporting the first where the probe reads the second moves the residual
by whatever the pack is absorbing. An adapter whose model reports `G` at a
different boundary must map it to the snow surface the same way
`energy/surface-energy-closure` requires, using the model's own heat storage in
the layer between, never the budget residual.

The catchment is snow-covered throughout the scored block, so this is also the
catchment's surface: while a pack is present the soil and the canopy are
underneath it and do not exchange latent heat with the atmosphere through it.
`LE` over the block is the pack's own sublimation.

## Why the case is built the way it is

Two properties of the scored stage make the assertion possible, and both are
built by the generator rather than assumed of the model.

### No precipitation falls during it

The `/io` contract carries **neither a melt flux nor a snowfall flux**. Melt
must therefore be inferred from the pack the model reports, and that inference
is only exact if nothing is being added to the pack at the same time. With
`pr = 0` throughout the block, snowfall is zero whatever rain–snow threshold a
model uses internally, and the pack's loss is melt and sublimation and nothing
else.

This is what lets the probe exist without extending `FLUX_VARS`. It is also a
real limitation: the probe says nothing about rain-on-snow, which is
operationally the most consequential melt case. See *Scope* below.

### The pack is ripe before it opens

The full snowpack energy budget carries a cold-content term:

```
Rn - LE - H - G  =  lambda_f * M  +  dU/dt
```

`U` is the internal energy of a sub-freezing pack. Energy spent warming it
toward 0 °C does no melting, and `U` is not a contract variable — no criterion
can observe it. So the case spends that cold content in an **unscored ripening
stage** and opens the scored block only once the pack is isothermal at 0 °C and
`dU/dt ~ 0`. What remains is fusion and nothing else.

That the stage is *sufficient* is a property of the case, not an assumption.
For a pack of `m` kg m-2 sitting at `T` below freezing,

```
U = c_ice * m * abs(T)            c_ice = 2100 J kg-1 K-1
```

and the ripening stage must deliver more than that. The generator bounds the
accumulation temperature to -10..-2 °C precisely because `U` scales with it,
and holds 40 days of rising radiation against it.

`generate.py`'s `__main__` computes that margin over a 348-seed sweep and
exits non-zero below 2.0. The worst seed delivers **2.89x** the cold content
of the pack it has to ripen. Both sides of the ratio are deliberately taken
at their conservative end: the cold content uses the *coldest* accumulation
day rather than the mean, since a pack integrates its season; and the energy
delivered has the pack's own sublimation subtracted, because latent heat
leaving the surface cannot also ripen the pack. Computing either the generous
way inflates the margin — done both ways, the same case reports 2.89x
honestly and 4.6x flatteringly.

A model carrying cold content that is *slower* to ripen than this margin allows
would meet the block still cold, spend part of the residual on `dU/dt`, and
fail. That is the probe's main false-positive risk and it is why the margin is
stated as a number rather than asserted.

Not assuming a ripe pack is
[`energy/snowpack-cold-content`](../../../ROADMAP.md), a separate and harder
probe that this one does not attempt.

## The case

One snow year at `PT1D`, in four stages labelled in a `_regime` column. The
harness strips every underscore-prefixed column before staging, so the model is
handed the weather and not the answer: it has to tell the stages apart from the
forcing itself.

| Stage | Days | Precipitation | Daily mean T | Rn | Scored |
| --- | ---: | --- | --- | --- | :---: |
| `spinup` | 30 | none | ~ -8 °C | low | no |
| `accumulation` | 110 | snowfall | -10..-2 °C | low | no |
| `ripening` | 40 | none | rising through 0 °C | rising to 55 W m-2 | no |
| `melt` | 60 | none | > +3 °C | 55 to 70-100 W m-2 | **yes** |

The cold stages sit clearly below, and the melt stage clearly above, any
plausible rain–snow threshold, so the phase split never depends on a model's
own PXTEMP and `snow_threshold_degC` is carried in `static` as the other probes
do. The melt stage staying well above freezing also suppresses nightly
refreeze, which would reopen the cold-content term the ripening stage closed.

Seeds vary the pack depth, the temperatures, the melt-season radiation and an
AR(1) weather anomaly, all drawn before any stage branch. Across seeds the pack
peaks at roughly 460–840 mm of SWE. A pack that size melting out entirely over
the 60-day block would demand **29–54 W m-2**; what a given model actually
demands depends on how much of it that model melts, and is lower whenever the
pack outlives the block.

## Tolerance

The 5 percent engineering rule against the melt energy itself, with the
2 W m-2 absolute floor carried from `energy/surface-energy-closure` and for the
same reason: the denominator goes to zero for a model that reports no pack.
These are the suite's engineering conventions, not observationally calibrated
values.

`1 mm day-1` of melt is `lambda_f / 86400 = 3.86 W m-2`, which is the scale to
keep in mind when reading the criterion's messages.

**Which of the two rules binds.** The crossover is at a demand of 40 W m-2:
below it five percent is under 2 W m-2 and the floor is the operative bound,
above it the relative rule is. Both happen here. Across the 24 to 43 W m-2 the
reference models actually produce, the floor binds on the lighter melts and
the relative rule on the heavier ones, so the effective band runs **+/-5 to
+/-8 percent** rather than a flat 5. Both rules being exercised seems better
than tuning the case until only one is ever reached, but say so if the probe
should sit wholly in one regime.

### Degeneracy guards

A model reporting no pack at all would satisfy the identity with two zeroes, so
the criterion enforces, as shares rather than depths so they hold for any seed:

- peak `snw` at least **half** the cumulative snowfall of the accumulation
  stage — the pack must appear;
- `M > 0` — it must have melted rather than merely sublimated away.

Failing a guard is a `FAIL` with the measured quantity, not an error. The case
supplied the snow; reporting no pack is a physics answer, and a wrong one.

**How much of the pack a model melts is deliberately not scored.** An earlier
draft required it to melt out, which fails a model that reports a physically
consistent budget and simply melts slowly - a calibration choice scored as a
conservation violation, and scored *ahead* of the identity the probe exists to
test. A model that melts a tenth of its pack and pays exactly the fusion energy
for that tenth is coupled correctly, which is the only thing this criterion is
entitled to judge. The share of the peak left is reported as a diagnostic so a
reviewer can see how far the melt got.

## What it catches

| Model | Water budget | Melt reported | Energy available | Verdict |
| --- | --- | ---: | ---: | --- |
| `reference_snow_energy` | closes 0.0000% | 379.5 mm | 24.43 vs 24.43 W m-2 | PASS |
| `reference_coupled` | closes 0.0000% | 667.6 mm | 0.00 vs 42.97 W m-2 | FAIL |
| `reference_two_head` | closes 0.0000% | 628.7 mm | 0.00 vs 40.47 W m-2 | FAIL |

Both failing models close their water budget to floating point and still pay
**nothing** for well over half a metre of fusion. That is the gap the probe
exists to find, and it is why a water-budget probe cannot find it.

`reference_coupled` melts a degree-day depth
(`models/reference_coupled/ht_adapter.py:117`) and then solves
`sensible = rn - ground - latent`, which drives its residual to exactly zero.
`reference_two_head` replaces the latent flux with a climatological evaporative
fraction and closes the same way, so it fails from the other direction with the
same residual.

Both therefore fail through one mechanism — an energy budget that closes on
itself and so has nothing left for fusion — rather than two independent ones.
That is the failure this probe was proposed to catch, and no reference model
available today pays the *wrong* amount rather than none: one that did would
sharpen the `must_fail` set, and adding it is worth doing when a model that
melts on a partial energy budget exists to imitate.

`reference_snow_energy` is the reference this probe adds: the same bucket with
melt bought from the pack's own surface budget,
`M = max(0, Rn - H - LE - G) / lambda_f`, capped by the ice that is there. Its
turbulent flux is a bulk form over a surface held at 0 °C,
`H = K_H * (0 - T_air)` with `K_H = 3 W m-2 K-1`, which is the melting-surface
assumption and applies only while the pack is melting. On a cold pack, where
the available energy is negative and no ice melts, the surface is below
freezing, the assumption does not hold, and `H` is the closing term — the same
statement `reference_coupled` makes, that the surface has no other place to put
the energy.

While a pack survives a step the soil and canopy do not evaporate through it,
so the latent flux is the pack's sublimation alone and the flux that decided
the melt is the flux that is reported. Over a melt season that identity holds
to under 1 W m-2, the residue being the single step on which the pack melts out
and the ice cap binds, where the surplus warms ground that has just been
uncovered. Letting soil evaporation run under a metre of snow, as an earlier
draft did, left the model deciding melt on one energy budget and reporting
another, 26 W m-2 apart.

### Baselines

`must_pass` names `reference_snow_energy` alone. That is structural, not an
exception: `must_pass` can only name models emitting what the probe requires,
and `reference_bucket`, `flex_lumped`, `flex_topo` and `sacsma_snow17` emit no
energy fluxes at all, so they are `N/A (INCOMPLETE)` here. Every merged energy
probe carries a single-model `must_pass` for the same reason.

`closure` is declared ahead of `melt_energy` as a precondition: a model that
leaks water reports a pack loss that was never melt. Nothing is pinned to it in
`baselines`, so it is never reported as this probe's result.

## Scope

What this probe does **not** test, stated so that a passing verdict is not read
for more than it is worth:

1. **Cold content.** The block is engineered to have none. A model that melts a
   sub-freezing pack is caught by `energy/snowpack-cold-content`, not here.
2. **Rain-on-snow.** Melt is inferred from `snw` and `sbl`, which requires a dry
   window. Closure during rain-on-snow is out of reach until the contract
   carries a melt flux — see
   [#17](https://github.com/Flood-Lab/HydroTuring/issues/17).
3. **Sublimation a model does not report.** `sbl` is in `requires`, so a model
   that omits it is `N/A (INCOMPLETE)` rather than failed. Unreported
   sublimation would inflate `M` while its energy cost already sits in `hfls`,
   and a reporting gap must not be scored as a physics failure. Of the models
   evaluated today `summa` misses this probe by that one variable alone.
4. **Melt timing within the block.** The identity is integrated over the block,
   so a model melting the right total on the wrong days passes. Catching that
   needs a per-step form and a tolerance that survives a pack melting out
   mid-block.
5. **A model that had the energy and declined to melt.** With no melt-out
   requirement, such a model passes if it reports a surface residual as small
   as its melt. Over a ripe pack under strong radiation that means putting the
   energy into `hfss` while the air is warmer than a surface at 0 C, which is
   not physical - but catching it needs the surface temperature, and the
   contract does not carry one.

## Provenance

Synthetic, generated at run time by `generate.py` from a recorded seed. Nothing
is committed as data, so there is nothing to memorise. Proposal and agreed
scope: <https://github.com/Flood-Lab/HydroTuring/issues/66>.

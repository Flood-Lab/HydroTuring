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

Over the scored block the surface residual is therefore **not** required to
vanish. It is required to equal what the pack absorbed:

```
M        = -(ice_end - ice_before) - sum_i(sbl_i * dt_i)     ice = snw - lwsnl
dU       = -(csnow_end - csnow_before)                       [J m-2]
R_avail  = mean_i(Rn_i - LE_i - H_i - G_i)                   [W m-2]
R_demand = (lambda_f * M + dU) / (N * dt * 86400)            [W m-2]

D = abs(R_avail - R_demand) <= max(0.05 * abs(R_demand), 2 W m-2)
```

with `c_ice = 2100 J kg-1 K-1` and `T0 = 273.15 K`.

A failure means the model moved water across a phase change, or changed the
temperature of its own snowpack, on energy its surface budget never accounted
for — or charged for a phase change that never appeared in the pack.

### Why the ice is `snw - lwsnl` and not `snw`

**`snw` is the pack's total water in this suite's own adapters.** Snow-17
reports `WE + LIQW + lagged excess + storage`
(`models/sacsma_snow17/ht_adapter.py:28`), and SUMMA reports `scalarSWE`, "ice
plus liquid in the pack" (`models/summa/ht_adapter.py:124`).

So a fall in `snw` is **net water leaving the pack**, not evidence of a phase
change, and differencing it misclassifies in both directions:

- **melt retained as liquid** moves no `snw` at all, yet fusion has happened
  and the energy for it was spent — an honest model would be failed for
  spending it;
- **drainage of water that melted days earlier** moves `snw` downward with no
  fusion happening now — an honest model would be failed for not spending
  energy it had no reason to spend.

Requiring `lwsnl`, the liquid held in the pore space, makes the ice observable.
Refreezing is not scored, and is not claimed to be — see *Tolerance* below for
why an integrated balance cannot see it.

### Why cold content is a term and not an assumption

A pack below freezing spends energy warming towards zero that does no melting.
That energy is real, it is in the surface residual, and no case can rule it out
for a model whose snow physics differs from the reference's — Snow-17 carries
exactly such a store in `NEGHS`.

An earlier version of this probe ripened the pack in an unscored stage so the
residual would be fusion alone. That made ripeness a **precondition no
criterion could observe**: the sweep that established it used the reference
model's own turbulent-exchange and radiation parameters, so it validated that
implementation rather than every compatible one. Warming and fusion stayed
indistinguishable in principle.

Requiring `csnow` and carrying `dU` in the balance removes the assumption. The
scored block now deliberately **opens on a cold pack** and warms through zero
inside the scored stretch, so the term has something to measure.

### What the case still rests on

One property, and the criterion enforces it rather than trusting it: **the
scored block carries no precipitation.** The contract has neither a melt flux
nor a snowfall flux, so ice arriving during the block could not be told from
ice melting. `melt_energy` refuses a block that is not dry.

## The control volume

The budget is over the **whole snowpack as a control volume**:

- the **snow surface**, where `rn` arrives and `hfls` and `hfss` leave;
- the **snow-soil interface**, where `hfg` leaves into the ground.

`hfg` keeps the definition it has everywhere else in the suite — ground heat
flux at the actual soil surface, positive into the ground. **This probe does
not relocate it** and asks nobody for a conductive flux at the top of a pack.

`Rn` is positive into the surface; `LE`, `H` and `G` are positive away from the
volume, and any of the three may be negative. Over a melting pack `H` is
normally negative: the air is warmer than a surface held at 0 °C, so the
turbulent flux is energy arriving.

`ice_before` and the cold content it carries are taken one step **before** the
block opens — from `state0` where the block opens the window, and from the
previous row otherwise — so nothing is measured from the block's own first
row. The criterion does this itself rather than through `storage_at`, because
that helper sums water storages and the ice here is a difference of two
reported columns.

## The case

One snow year at `PT1D`, in three stages labelled in a `_regime` column. The
harness strips every underscore-prefixed column before staging, so the model is
handed the weather and not the answer.

| Stage | Days | Precipitation | Daily mean T | Rn | Scored |
| --- | ---: | --- | --- | --- | :---: |
| `spinup` | 30 | none | -24..-12 °C | low | no |
| `accumulation` | 150 | snowfall | -24..-12 °C | low | no |
| `melt` | 60 | none | opens below -14 °C, ends above +4 °C | rising from winter values to a 70-100 W m-2 peak | **yes** |

The accumulation stage sits clearly below any plausible rain-snow threshold, so
the phase split never depends on a model's own PXTEMP, and
`snow_threshold_degC` is carried in `static`.

Both ends of the scored block are held rather than left to the draw, and how
cold the opening has to be is not a matter of taste. The ratio of the two terms
is fixed by the constants alone:

```
cold content / fusion  =  c_ice * dT / lambda_f  =  0.0063 * dT
```

independent of how big the pack is. A pack 6 K below freezing carries a cold
content worth 3.8 percent of its own fusion, which the 5 percent relative rule
cannot see; at 8 K it is 5.0 percent and at 18 K it is 11.3 percent. The
opening twenty days are therefore held below **-14 °C**, and the block is 60
days rather than 100 so the term also clears the 2 W m-2 floor in absolute
terms. The closing twenty days are held above +4 °C so the fusion term is
exercised too. A case that tested only half the balance would prove half as
much, and a case that opened at -6 °C would carry the cold-content term without
ever being able to catch anything with it.

Seeds vary the pack depth, the temperatures, the melt-season radiation and an
AR(1) weather anomaly, all drawn before any stage branch.

## Tolerance

The 5 percent engineering rule against the energy demanded, with the 2 W m-2
absolute floor carried from `energy/surface-energy-closure` and for the same
reason: the denominator goes to zero for a model that reports no pack. These
are the suite's engineering conventions, not observationally calibrated values.
The relative part is taken on `abs(R_demand)` defensively, but note that
refreezing is not scored here and no credit for it exists: a block that gains
ice over its whole length fails the `M > 0` guard, and refreezing that happens
inside the block and re-melts leaves both endpoints unchanged and so is
invisible to an integrated balance. Catching refreeze needs a per-step form,
and is not claimed.

`1 mm day-1` of melt is `lambda_f / 86400 = 3.86 W m-2`, which is the scale to
keep in mind when reading the criterion's messages.

**Which of the two rules binds.** The crossover is at a demand of 40 W m-2:
below it five percent is under 2 W m-2 and the floor is the operative bound,
above it the relative rule is. Both occur here. That is the floor doing its job
rather than a tolerance that got away, and it is how the floor already behaves
in the merged suite: `energy/surface-energy-closure` records that "each night's
|Rn| is 25-40 W m-2, so every night block uses the 2 W m-2 floor". A floor that
never binds is not a floor, and the case is deliberately not tuned to keep the
demand above the crossover, because that would remove the protection on exactly
the models that need it.

### Guards

`lwsnl` and `csnow` are the model's own statements about its own pack, and both
sit on one side of the balance. A single boundary row could therefore decide a
verdict, so the criterion bounds what those rows can say. Seven checks in all.
Four are arithmetic consistency checks on a single row:

- `snw`, `lwsnl` and `csnow` must be **non-negative**. A model that pays for
  600 mm and melts 300 could otherwise report `lwsnl = -300` on the step before
  the block and have the difference read as ice that left.
- `lwsnl` may not exceed `snw`. The liquid is held inside the pack.
- `csnow` must be **zero where there is no ice**. Cold content is the energy
  needed to bring ice to 0 °C, and there is none to bring.
- `csnow` may not exceed what the reported ice would hold **10 K below the
  coldest air in the record** — applied only to the rows the identity reads,
  the step before the block and its last step. It bounds a model declaring its
  pack cold on either of them; with the specific-cold-content check on, that
  check reaches the last-row case first, and the cap is what stands when that
  one is relaxed or when the opening row is the row inflated. Applied to *every*
  row it would assert a temperature bound on states the criterion never
  evaluates, and would refuse an honest Snow-17-type pack: Snow-17 bounds its
  deficit as a mass fraction, `NEGHS <= 0.33 * WE`, worth 52 K, which a thin
  early-season pack legitimately reaches. The claim is not that a pack can
  never be colder than the air, but that the reported state must stay
  physically plausible where this identity reads it.
- **no liquid on the row where the block opens.** `ice_before` is a self-report
  on a single row that no other term constrains, and declaring part of it
  liquid shrinks the melt attributed to the model. The share guard below cannot
  oppose that, because it is evaluated on the same row and moves with it. This
  case holds the air at or below −12 °C through accumulation and −14 °C into
  the block, so an honest model opens frozen — a property of this case, not of
  snowpacks, which is why the tolerance is a parameter.
- **specific cold content may not rise across the block**:
  `csnow_after <= ice_after * csnow_before / ice_before`. Inflating the opening
  cold content to buy room raises the demand by at least as much. Tied to the
  same precondition: it may rise legitimately over a block that ends *colder*
  than it began, so a probe whose block does not end warm turns this off. This
  check closes a **construction, not a route** — it sees the flip only when
  `csnow` is copied from the block's last row; scaled to the declared ice
  instead, the excess is round-off and only the opening-liquid check stands.
- **`snw` may not gain water**: `Δsnw + sbl·dt <= 0.5 mm`, read on **every
  step** of the block rather than on a single row, because the construction it
  closes is a one-step excursion that returns. `snw` alone, not `snw + canopy`, because the
  identity differences `snw - lwsnl` and summing the canopy in would let the
  same dip hide there. One consequence, since this case carries no canopy water
  through the block: a probe whose case did would see canopy unloading faster
  than the tolerance refused, and would need a different bound. This bounds gains in one store, so it is **narrower
  than `closure`** rather than a per-step form of it: `closure` is cumulative
  and spans every store, which is exactly why a one-step excursion that
  returns nets to nothing there. Deposition a model reports in `sbl` is netted
  out, so an honest pack that gains water by saying so passes at any
  tolerance; what the tolerance covers is *unreported* gain.

  The tolerance is coupled to the floor and the two should be changed together.
  A model that lowers the opening pack and closes the offset gradually at `g`
  mm per step understates the demand by `lambda_f · g / (dt · 86400)` W m-2,
  independent of block length. At a daily step, keeping that inside the
  2 W m-2 floor requires `g <= floor · 86400 / lambda_f = 0.518` mm, which is
  why the tolerance is 0.5 and not the 1 mm that honest deposition alone would
  justify.

  That figure is the **least** the tolerance lets through, not a ceiling. Each
  step's gain is net of the drainage the model reports, so a pack reported not
  to drain while its runoff carries the water away gets more room than this.
  That is a multi-row fiction rather than a boundary-row one; it is named in
  *Scope* rather than guarded, and closing it needs a pack outflow flux, which
  [#17](https://github.com/Flood-Lab/HydroTuring/issues/17) will add as
  `snmsl`.

Two are about the pack being real at all, as shares rather than depths so they
hold for any seed:

- the ice standing **at the step before the block opens** must be at least
  **half** the accumulation stage's cumulative snowfall — the pack must appear.
  Measured there rather than at the peak anywhere in the window, because a peak
  lets a model that melts and drains before the block still satisfy it, and
  lets one row flip it;
- `M > 0` — it must have melted rather than merely sublimated away.

Each of these is reported as a **contract** error rather than a physics one,
because that is what it is.

None fires on an honest model. Over 505 seeds and the three reference models
no guard trips; on the gate seeds the opening row carries exactly
0 mm of liquid, and the largest excess over the specific-cold-content
inequality is +2.2×10⁻⁸ J m⁻² over a 20-seed sweep, which the 1 J m⁻²
tolerance clears by eight
orders of magnitude.

Together they close a construction that the share guard alone does not. Taking
`reference_degree_day`'s own output and changing **one row** — the step before
the block, setting `lwsnl = snw - (ice_end + block sublimation + 1 mm)` with
`csnow` copied from the block's last row — leaves the declared ice above half
the snowfall, so every earlier guard passes, while the melt attributed to the
model falls to 1 mm and a model that paid nothing reports
`melt of 1.0 mm demands 0.06 W m-2`. It passes on **four of the five gate
seeds** without the two checks above, and on none of them with either one.

**How much of the pack a model melts is deliberately not scored.** Melting part
of it and paying for exactly that part is correct coupling; failing it for
melting slowly would score a calibration choice as a conservation violation.
The exact model leaves 44–66 % of its opening ice. The share left is reported
as a diagnostic.

## What it catches

On **gate seed 1415560089**, one recorded case rather than a per-model worst:

| Model | Melt | Fusion | Cold content | Demanded | Available | Verdict |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `reference_snow_energy` | 291.8 mm | 18.78 | 7.76 | 26.54 | 26.53 | PASS |
| `reference_degree_day` | 300.3 mm | 19.33 | 2.39 | 21.72 | **0.00** | FAIL |
| `reference_warming_free` | 291.8 mm | 18.78 | 7.76 | 26.54 | **18.78** | FAIL |

All figures W m-2 except melt. Two different failures, not one repeated:

`reference_degree_day` melts on air temperature and closes its energy budget
around its evaporation alone, so **nothing** pays for 300 mm of fusion. Its own
cold content is smaller than the exact model's because melting on temperature
takes ice out faster, and the cold content that ice was holding leaves with it
— an earlier draft of this model removed the ice and left the deficit behind,
which drove its reported pack to -30 °C while the air stood at +10 °C. It is
the construction the probe was proposed to catch.

`reference_warming_free` melts on the energy actually available — so its melt
matches `reference_snow_energy` exactly, as it must, since they share a water
side and differ only in what the energy head reports — but it charges the
fusion alone and lets the pack climb from sub-freezing towards zero for
nothing. **It fails by 7.76 W m-2, which is precisely the cold-content term** — nearly
four times the 2.00 W m-2 allowed, and by 3.6 to 4.2 times across the five gate
seeds. `reference_degree_day` misses by 10.9 to 18.0 times over the same seeds.
A probe that assumed a ripe pack instead of measuring one would not notice it.

`reference_snow_energy` is the reference this probe adds: an energy-balance
pack that buys melt with `max(0, Rn - H - LE - G) / lambda_f`, holds liquid in
its pore space up to 5 percent of the ice, carries its cold content as an
energy deficit, and refreezes liquid when the budget turns negative. Every
joule is accounted at the phase change it drove, so the identity holds by
construction, as `reference_bucket` closes the water budget by construction.

### Baselines

`must_pass` names `reference_snow_energy` alone. That is structural:
`must_pass` can only name models emitting what the probe requires, and the four
physical baselines emit no energy fluxes at all, so they are `N/A (INCOMPLETE)`.
Every merged energy probe carries a single-model `must_pass` for the same reason.

`closure` is declared ahead of `melt_energy` as a precondition: a model that
leaks water reports a pack change that was never a phase change. Nothing is
pinned to it, so it is never reported as this probe's result.

## Scope

What this probe does **not** test:

1. **Rain-on-snow.** The scored block must be dry, because ice arriving cannot
   be told from ice melting without a snowfall flux. Operationally the most
   consequential melt case, and out of reach until the contract carries one —
   see [#17](https://github.com/Flood-Lab/HydroTuring/issues/17).
2. **Melt timing within the block.** The balance is integrated, so a model
   melting the right total on the wrong days passes.
3. **A model that had the energy and declined to melt**, if it reports a
   surface residual as small as its melt. Over a ripe pack under strong
   radiation that means putting the energy into `hfss` while the air is warmer
   than a surface at 0 °C, which is not physical, but catching it needs a
   surface temperature the contract does not carry at this boundary.
4. **Sublimation a model does not report.** `sbl` is in `requires`, so a model
   omitting it is `N/A (INCOMPLETE)` rather than failed: a reporting gap must
   not be scored as a physics failure.
5. **Everything reachable from the two boundary rows.** Every construction
   found in review has landed on one of them, and it is worth saying why rather
   than listing the tricks. The identity reads the pack on exactly two rows —
   the step before the block and its last step — and differences them. Neither
   is scored by anything else: the block's interior is covered by the identity,
   and the record's other rows by `closure` in the aggregate. So anything those
   two rows say about `snw`, `lwsnl` or `csnow` moves the melt or the warming
   attributed to the model, and nothing else notices. Seven checks now bound
   what they can say, four of them arithmetic on a single row. A reader looking for the next construction should look
   there first.

   The ones already known, none of which rescues a `must_fail` reference
   because none of them over-pays:

   - an **under**-payer lowering the opening `snw` or `ice` — the flip and the
     dip, both now closed — or raising the last row's, which is **bounded
     rather than closed**: the gain check still allows that row to rise by the
     step's own reported decrease plus the 0.5 mm tolerance;
   - an **under**-payer declaring the last row's retained liquid as ice, which
     raises the melt it appears to have done;
   - an **over**-payer declaring its remaining pack liquid on the last row,
     raising the melt it appears to have done, with the floor absorbing the
     rest;
   - an **over**-payer raising `snw` on the opening row, or lowering it on the
     last row, which are the same move from the other side;
   - an **over**-payer raising the opening `csnow` to the cap while reporting
     zero on the last row, which the specific-cold-content check allows in
     that direction.

6. **A model that reports an internally consistent fiction.** This is the limit
   of the method rather than a gap to be closed. The criterion reads the
   model's own states, so it can bound *inconsistency between* what a model
   reports — liquid appearing where the block opens, specific cold content
   rising, a pack gaining water in a dry block — and it cannot reach a model
   whose whole reported series is coherent and wrong. A pack reported low
   through the entire accumulation stage is consistent, and no criterion
   reading only self-reports will catch it. Specifically still open:

   - `reference_warming_free` passes if it reports `csnow = 0` — on the row
     before the block alone, not only throughout, since that is the row the
     cold-content term is differenced from.
   - A pack reported not to drain while its runoff carries the water away gains
     more room under the pack-gain tolerance than the `lambda_f · g / (dt ·
     86400)` figure suggests, because each step's gain is net of the reported
     drainage. Closing that needs a pack outflow flux, which
     [#17](https://github.com/Flood-Lab/HydroTuring/issues/17) will add as
     `snmsl`.


   Closing any of these needs a quantity the case supplies rather than the
   model, which for a snowpack means an observable the contract does not carry.

## Provenance

Synthetic, generated at run time by `generate.py` from a recorded seed. Nothing
is committed as data. Proposal and agreed scope:
<https://github.com/Flood-Lab/HydroTuring/issues/66>.

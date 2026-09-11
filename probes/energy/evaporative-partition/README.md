# energy/evaporative-partition

Two runs of one seed. Net radiation, air temperature and evaporative demand
are byte-identical between them. One of them loses 120 days of summer rain.

Because the radiation did not change,

```
sum_W (dLE + dH + dG) = sum_W dRn = 0
```

is an identity. The latent heat a drying surface gives up has nowhere to go
but the sensible and ground fluxes. There is no parameterisation in it, no
atmospheric feedback to model, and no tolerance worth arguing about.

The sibling probe, `energy/latent-heat-et-consistency`, asks whether a
model's two budgets agree on one number. This one asks whether they agree
under a change neither of them has seen.

## What this case is, and what it is not

**It is an idealised seasonal drydown. It is not a flash drought**, and
nothing here claims otherwise.

The flash drought literature settled its definition on the *rate of
intensification* rather than on duration (Otkin et al. 2018), and the field
has since organised around the triad of rapid onset, drought development and
impacts (Christian et al. 2024). Onset criteria sit on the order of one to six
pentads, and duration caps exist precisely to keep these events distinct from
seasonal drought. A 120-day rainfall-free window is on the wrong side of that
line.

The window is long because the identity needs a signal large enough to score,
and the case takes a long time to produce one. That is measured, not assumed:

| Days into the window | Share of the reference model's total shift |
| --- | --- |
| 1 – 60 | under 1 percent |
| by day 82 – 85 | 10 percent |
| by day 105 | 50 percent |
| days 76 – 120 | about 98 percent |

Stable across seeds. The reason is worth stating because it is not obvious:
the control run keeps receiving rain, so the two runs do not diverge until the
drought run has spent its soil buffer. Ninety days was the first attempt and
produced falls as small as 8 percent of the window's net radiation on some
seeds — a perturbation, not a regime transition.

So the case is a seasonal drydown that *contains* a rapid intensification, and
the `_intensification` column marks days 76 to 120 where that intensification
happens. The criterion reports the shift over that stretch alongside the
four-month total, which turns a total into a rate:

```
ok  partition_shift  ... dhfls -1266, dhfss +1266, dhfg +0, sum -2.3e-13, dEF -0.333;
                     intensification stretch (45 steps) dhfls -1236, dhfss +1236, dEF -0.464
```

That stretch is reported and never gated. Its placement comes from the
reference model, which is how this suite already locates evaluation windows,
and a threshold whose position depends on the reference model is a threshold
worth not having.

The repartitioning being tested is the one that governs flash drought onset.
The case is simply slower than a flash drought, and says so.

## What it catches that the sibling probe does not

`reference_ground_dodge` **passes every criterion of the sibling probe.** Its
latent heat agrees with its evaporation at every step, its water budget closes
to floating point, and so does its energy budget. It fails here:

```
FAIL  partition_shift   hfg absorbed 1.00 of the shift (limit 0.20); the energy
                        went into the soil instead of the air
                        [dhfls -1266, dhfss +0, dhfg +1266 (W m-2 day),
                         dEF -0.112; intensification stretch dEF -0.204]
```

Its sensible heat flux is a fixed share of net radiation that never reads the
soil, so when the soil dries the ground flux takes the whole shift. Over four
months that is a claim about soil heat capacity that no soil has. Bounding the
ground share is what closes the dodge, and it is why the criterion asserts
three things rather than only the sum.

`reference_two_head` fails earlier, on the response gate: its latent flux does
not move at all, because its energy side has no idea its water side ran out of
water.

```
FAIL  partition_shift   hfls barely moved: +0 against a required fall of 422.4
                        (10% of the 4224 the window received)
```

## How a model would pass this while understanding no physics

It would not. Passing requires moving a specific quantity of energy from
latent to sensible in response to a water deficit that never occurred, under
radiation that did not change. There is no regularity in any training
distribution that supplies it, because the counterfactual half of the pair
does not exist in any observed record. That is also why this probe cannot be
ported to the real-data track: the case is not a measurement, it is a
comparison between a measurement and something that did not happen.

There is a consequence beyond this benchmark. The transfer from latent to
sensible heat as soil dries is the mechanism by which the land surface
amplifies heat extremes (Zhou et al. 2021; Dong et al. 2020 shows that the
evaporation stress operator this probe interrogates is what separates land
surface models from one another). A model that fails `partition_shift` is not
merely inconsistent. It is structurally incapable of representing a hot
drought, however well it fits a hydrograph. Failing this is a statement about
what the model can never be used for.

## Related work

Two lineages meet here, and the probe is new at their intersection rather than
in either one.

**Conservation in learned emulators.** Beucler et al. (2021) enforce analytic
constraints in neural networks emulating convection, architecturally or
through the loss, and quantify the "physical inconsistency" that remains
across mass, momentum, radiation and energy budgets. That is the nearest prior
art for the sibling probe. What is new here is not the idea that a learned
model should conserve; it is that the assertion is made about the *response to
a counterfactual* rather than about a same-instant residual, which is a
property no penalty term evaluated on a validation set can reach.

**Metamorphic testing.** Asserting how an output must change when an input is
changed in a known way is the software-testing idea this criterion belongs to,
and the suite already contains two of them: `mass/causality` and
`mass/extreme-rain` are metamorphic relations on the water budget. This is the
first one on the energy partition.

## Why the ground-share bound is a judgement and stays

`max_ground_share = 0.20` is the one number in either probe that is a physical
judgement rather than an identity. Over a four-month window the ground heat
store is a small residual of the surface energy budget: G/Rn is a bounded
fraction that land surface schemes parameterise explicitly and that varies
with cover rather than with soil moisture (Kustas & Daughtry 1990; Yang et al.
1999; Liang et al. 1999).

The number is open to argument. Dropping the bound is not: without it
`reference_ground_dodge` passes, and it is the model that most resembles what
an unconstrained energy head would learn.

## The case, in short

Three years at a daily step with a year of spinup, 3 seeds, 120-day window
from 1 June of the second scored year, soil capacity 120 mm — smaller than the
sibling probe's 320 mm, deliberately, because a bucket that cannot empty
cannot become water-limited. `min_shift: 0.10` sits below a measured range of
29 to 43 percent, a factor of three in hand rather than a threshold fitted to
a result. The window is sized against the moisture-limited / energy-limited
framework of Haghighi et al. (2018).

Holding radiation fixed removes a confounder, not the phenomenon. Real
droughts arrive with coupled warming, which would make the shift larger rather
than smaller; an honest model is not made to look worse by the simplification.

`_perturbed` and `_intensification` are stripped before the model sees the
forcing, so the model is never told which stretch it is being judged on.

## References

- Beucler, T., Pritchard, M., Rasp, S., Ott, J., Baldi, P., & Gentine, P. (2021).
  Enforcing analytic constraints in neural networks emulating physical systems.
  *Physical Review Letters*, 126, 098302. https://doi.org/10.1103/PhysRevLett.126.098302
- Christian, J. I., Hobbins, M., Hoell, A., Otkin, J. A., Ford, T. W., Cravens, A. E.,
  Powlen, K. A., Wang, H., & Mishra, V. (2024). Flash drought: A state of the science
  review. *WIREs Water*, 11(3), e1714. https://doi.org/10.1002/wat2.1714
- Dong, J., Dirmeyer, P. A., Lei, F., Anderson, M. C., Holmes, T. R. H., Hain, C.,
  & Crow, W. T. (2020). Soil evaporation stress determines soil moisture-evapotranspiration
  coupling strength in land surface modeling. *Geophysical Research Letters*, 47(21).
  https://doi.org/10.1029/2020GL090391
- Haghighi, E., Short Gianotti, D. J., Akbar, R., Salvucci, G. D., & Entekhabi, D.
  (2018). Soil and atmospheric controls on the land surface energy balance.
  *Water Resources Research*, 54(3), 1831-1851. https://doi.org/10.1002/2017WR021729
- Kustas, W. P., & Daughtry, C. S. T. (1990). Estimation of the soil heat flux/net
  radiation ratio from spectral data. *Agricultural and Forest Meteorology*.
  https://doi.org/10.1016/0168-1923(90)90033-3
- Liang, X., Wood, E. F., & Lettenmaier, D. P. (1999). Modeling ground heat flux in land
  surface parameterization schemes. *Journal of Geophysical Research: Atmospheres*.
  https://doi.org/10.1029/98JD02307
- Otkin, J. A., Svoboda, M., Hunt, E. D., Ford, T. W., Anderson, M. C., Hain, C.,
  & Basara, J. B. (2018). Flash droughts: A review and assessment of the challenges
  imposed by rapid-onset droughts in the United States. *Bulletin of the American
  Meteorological Society*, 99(5), 911-919. https://doi.org/10.1175/BAMS-D-17-0149.1
- Yang, Z.-L., Dai, Y., Dickinson, R. E., & Shuttleworth, W. J. (1999). Sensitivity of
  ground heat flux to vegetation cover fraction and leaf area index. *Journal of
  Geophysical Research: Atmospheres*. https://doi.org/10.1029/1999JD900230
- Zhou, S., Williams, A. P., Lintner, B. R., Berg, A. M., Zhang, Y., Keenan, T. F.,
  Cook, B. I., Hagemann, S., Seneviratne, S. I., & Gentine, P. (2021). Soil
  moisture-atmosphere feedbacks mitigate declining water availability in drylands.
  *Nature Climate Change*, 11(1), 38-44. https://doi.org/10.1038/s41558-020-00945-z

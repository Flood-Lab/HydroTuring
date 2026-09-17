# momentum/stage-discharge-monotonic

A staff gauge is the oldest instrument in hydrology and it imposes two
constraints that no water budget does. This probe asks a model that reports
the stage of its channel to satisfy both.

The first is monotonicity. Stage has to rise with the flow it is drawn
against: read the reach at a higher discharge and the gauge reads a higher
level, while a gauge that ratchets up on a running maximum and relaxes only a
fraction of a percent a day steps down against its own history instead. The
share of the span a binned rating may dip is 8%, because a stage read off a
store is hysteretic by construction and dips for that reason alone — the 2%
the criterion first used was only survivable while the rating spanned tens of
metres, which no reach does.

The second is the loop. A flood wave steepens as it arrives, so at the same
water in transit the rising limb sits *lower* than the falling limb. The sign
of that loop is physical rather than conventional: a rating read backwards,
high while the flow is still climbing and low once it is leaving, inverts it,
and the criterion fails that direction and only that one. The limbs are paired
on the channel store — the quantity the hysteresis is about; a reach that
holds nothing in transit is read on its discharge instead.

A gauge that is a function of the discharge alone — `f(dis)` and nothing else
— has no loop to show, and the criterion reads it as the single-valued rating
it is rather than failing it on the store-axis noise the pairing produces.
`flex_lumped` and `sacsma_snow17` are exactly that gauge: they solve Manning's
depth from their own outflow, so their rating is single-valued by construction.

Four years of temperate weather with five multi-day storms, in a single reach
whose width, slope and roughness are handed to every model through the static
file. `stage` is a diagnostic the probe asks the model to report in metres. It
is deliberately not one of the storages the suite differences its budget over:
a stage is a reading, not a volume, and summing it into `reported_states`
would corrupt every closure test in the suite.

| Criterion | Asserts |
| --- | --- |
| `rating_monotonic` | stage does not fall against the running maximum as the abscissa rises |
| `rating_loop` | at matching water in transit, the rising limb sits below the falling limb |
| `non_degenerate` | the stage varies with the weather |

Must-fail: `reference_rating_drift`, whose stage is a decaying running maximum
of its own discharge; `reference_rating_inverted`, which solves its gauge from
the floodplain's release instead of the channel's and so loops backwards;
and `reference_flat_stage`, whose gauge reports a constant.

The probe requires `dis`, so the rating is drawn against discharge and not
against a storage: a model that reports only its channel store is read on
that store instead, and the criterion says so in its reason rather than
claiming a discharge check it did not make. `reference_rating` is in
`must_pass` — with two separated time constants it is the only baseline that
traces a correctly-signed loop, and it is what keeps the loop criterion's
"the right way" branch live rather than dead.

## Limitations

A gauge that is hydrograph-shaped but decoupled from the model's own discharge
is not caught by either criterion. A stage taken from another seed's record, or
one lagged by a fortnight, can pass while being physically meaningless, because
the probe only checks that the gauge moves with the flow it claims to describe.

Where `reference_rating`'s loop clears the 2 mm size floor (90 of 130 seeds in
the sweep) it measures 2.0–5.8 mm, so it sits just above the floor; on the
other 40 seeds it takes the size escape, which passes either way.

## Where the implementation departs from the design settled in #53

Stated so they are accepted knowingly rather than discovered in review:

- the limbs are paired on the channel store by nearest falling step, rather
  than on discharge bins with a turning-point dead-band;
- twelve bins rather than eight;
- the size floor is an absolute 2 mm — a staff gauge's reading error — rather
  than 5% of `dH` plus 0.02 m. It is absolute because the reference gauges are
  different instruments whose spans differ by two orders of magnitude between
  seeds, while the loop a gauge can resolve does not;
- there is no `dH >= 0.5 m` degenerate case, so a flat stage passes both rating
  criteria and only `non_degenerate` catches it;
- one fixed reach geometry, handed to every model, rather than one drawn per
  seed;
- two linear reservoirs rather than a Muskingum prism plus wedge;
- a decaying peak memory (`0.997` a day) rather than a random walk.

The two deviations that were load-bearing, and why, sit in the criterion
docstrings: the loop and the monotonicity share are opposite constraints on
the same rating, and the reference gauges had to be built at a physical scale
before both could hold at once.

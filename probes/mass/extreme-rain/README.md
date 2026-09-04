# mass/extreme-rain

**More rain cannot mean less runoff, nor more runoff than was added.**

## The physics

Take a record and scale its largest storm by two, five and ten. Ten times
the largest day of a temperate record is close to a metre of rain in a
day, outside anything a training record holds. Three things stay exact
along that ladder. Integrated runoff cannot fall when more rain is added.
No rung can add more runoff than it added rain, because runoff has to be
rain or stored water. And across the whole ladder most of the added rain
must run off: evaporation is bounded by demand and storage by capacity,
and a storm several times larger than either has nowhere else to go. On
every rung, what runs off cannot exceed what fell plus what the catchment
could have held.

## What it measures

Four variants of the same weather, with one day scaled. The criterion sorts
them by rain, checks each rung for a fall in runoff or a rise beyond the
rain added, checks the whole ladder for at least a tenth of the added rain
returning as runoff, and checks every variant against the absolute bound.

What fails it: `reference_saturating`, the reference bucket with its daily
runoff capped at 25 mm, which is what a learned response looks like at the
edge of its data. Scaling the storm from one to ten adds 997 mm of rain
and not a millimetre of its runoff. The in-sample leaker, which loses just
over half of any rain above what it has seen, passes: it still returns 0.45
of the added rain on every rung, monotone and bounded. Extrapolation that
loses water is caught by the closure probes; this one catches responses
that flatten, turn down, or amplify.

## Results on the gate seeds

The exact bucket returns 1.00 of the rain added on every rung: on the day
of the largest storm its soil is already full, so every extra millimetre
runs off. The capped bucket returns 0.00.

## Criteria

| Criterion | What it asserts |
| --- | --- |
| `monotone_response` | runoff does not fall along the ladder, no rung adds more runoff than rain, at least a tenth of the rain added at the top runs off, and no variant runs off more than fell plus the storage bound |
| `non_degenerate` | runoff varies with the weather |

## Baselines

- `must_pass: reference_bucket`
- `must_fail: reference_saturating` on `monotone_response`

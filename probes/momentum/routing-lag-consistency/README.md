# momentum/routing-lag-consistency

A larger drainage network should not move a comparable storm response faster
merely because it is larger. This probe gives byte-identical daily weather and
the same isolated design storm to four synthetic catchments, then asks whether
the runoff-peak delay is on a Snyder travel-time scale and whether that delay
grows across a Hack-derived geometry ladder.

This implements [accepted proposal #8](https://github.com/Flood-Lab/HydroTuring/issues/8).
The expected lag is calculated only from public fields in `static.json`; it is
never supplied to the model.

## Synthetic experiment

Every seed produces 365 daily spinup rows followed by a 90-day scored record.
The spinup is stochastic, temperate and rain dominated. The scored record is
dry except for one 48--52 mm storm occupying one daily row selected at scored
day 20 +/- 2. This gives at least ten dry pre-event days for a baseline and
more than 30 post-event days for the routed response.

Proposal #8 used a three-hour storm as an example. This implementation adapts
it to one `PT1D` row so the probe can exercise daily model adapters. A one-step
pulse gives one unambiguous rainfall row. Row timestamps mark interval starts,
so the continuous centroid of this 24-hour rectangular pulse is 12 hours after
that row starts; the criterion measures from this interval centre.

The four variants have byte-identical `time`, `pr`, `tas`, `pet` and hidden
`_event_pr` columns. Area is prescribed, main-channel length follows Hack's
length-area relation, and centroid channel length is half the main-channel
length:

```text
A_mi2 = A_km2 / 2.589988110336
L_mi  = 1.4 * A_mi2^0.6
L_km  = 1.609344 * L_mi
Lc_km = 0.5 * L_km
```

| Variant | Area (km2) | L (km) | Lc (km) | Raw Snyder lag (h) | Corrected lag (days) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `small` | 30 | 9.797 | 4.898 | 9.582 | 0.631 |
| `medium` | 300 | 39.001 | 19.500 | 21.951 | 1.123 |
| `large` | 3,000 | 155.264 | 77.632 | 50.287 | 2.250 |
| `xlarge` | 10,000 | 319.741 | 159.870 | 77.570 | 3.335 |

The relation and units follow Hack's drainage-basin analysis in
[USGS Professional Paper 294-B](https://pubs.usgs.gov/publication/pp294B).
Hack's relation is a regional empirical fit rather than a universal channel
law, and the upper end of this ladder, including the 10,000 km2 catchment,
extends beyond the roughly 375 mi2 upper end of Hack's original sample. These
lengths define an idealised synthetic geometry ladder; they are not estimates
for four representative real basins.

The harness removes `_event_pr` before staging `forcing.csv`; the criterion
retains it only to locate the added storm. The probe requires a model to
consume `pr` and the three public geometry fields `area_km2`,
`main_channel_length_km` and `centroid_channel_length_km`. Standard synthetic
catchment attributes are also present, but do not change among variants. A
runoff-only model is therefore eligible only when its adapter truthfully
declares and consumes those geometry inputs.

## Reference travel time

For lengths in kilometres, the probe uses a Snyder coefficient of 4.0 and the
SI conversion factor 0.75:

```text
t_raw = 0.75 * 4.0 * (L_km * Lc_km)^0.3 hours
```

`Ct` is a regional calibration coefficient, not a universal constant. The
value 4.0 is a declared synthetic-experiment choice. The
[HEC-HMS Snyder documentation](https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/transform/snyder-unit-hydrograph-model)
describes 1.8--2.2 as typical and reports calibrated values from about 0.4 to
8.0, so applying this benchmark to a real basin would require local
calibration.

Snyder's standard excess-rainfall duration is `t_raw / 5.5`. Because this
probe's storm lasts 24 hours, the non-standard-duration correction gives:

```text
D_standard = t_raw / 5.5
t_expected = t_raw - (D_standard - 24) / 4
t_reference_peak_from_storm_start = t_expected + 24 / 2
```

The resulting expected lags are approximately 0.631, 1.123, 2.250 and 3.335
days from the excess-rainfall centroid to the runoff peak. The added 12 hours
only converts that origin to the storm-start origin used to index the reference
kernel; it is not additional physical travel time. The Snyder relation comes from
[Snyder (1938), "Synthetic Unit-Graphs"](https://doi.org/10.1029/TR019i001p00447).

Snyder's unit hydrograph is defined for **excess rainfall**, whereas the probe
observes prescribed gross rainfall (`pr`) and modelled total runoff (`mrro`).
Gross-rainfall-to-total-runoff lag is therefore an observable end-to-end proxy:
it includes any interception, infiltration, storage, runoff generation and
routing delays inside the model. The Snyder lag is an idealised comparison
scale here, not a direct prediction of the model's effective-rainfall response.
The probe cannot isolate channel routing from those other processes.

## What is measured

The precipitation time is the volume centroid of the hidden one-day event:

```text
t_P = sum(interval_centre(t) * event_pr) / sum(event_pr)
```

For each variant, the criterion estimates pre-event runoff from the five rows
immediately before the event and finds the total-runoff peak from event onset
to the end of the record:

```text
t_observed = t_peak(mrro) - t_P
```

The positive response above the pre-event median must integrate to at least
0.01 of the event-rain depth. This prevents a zero or effectively flat series
from acquiring a meaningless lag through `argmax`.

`lag_time_bounds` requires every variant to satisfy:

```text
max(0, 0.5 * t_expected - 0.5 day)
    <= t_observed <=
2.0 * t_expected + 0.5 day
```

The factor-of-two interval allows broad empirical variation. The half-day term
accounts for quantisation of a daily peak. `scaling_monotonicity` sorts the runs
by `area_km2`, allows at most a half-day apparent reversal between adjacent
daily measurements, and requires at least two days of total lag growth from
the smallest to the largest catchment.

| Criterion | Asserts |
| --- | --- |
| `lag_time_bounds` | every observed rainfall-to-runoff peak lag lies within its broad Snyder bounds |
| `scaling_monotonicity` | lag does not materially reverse and grows observably across the full area ladder |

## Discrimination and limits

| Reference model | Expected result | Required failure |
| --- | --- | --- |
| `reference_snyder_router` | PASS | |
| `flex_lumped` | PASS | |
| `reference_instant_router` | FAIL | `lag_time_bounds` |
| `reference_inverse_router` | FAIL | `scaling_monotonicity` |

The positive reference converts a fixed share of rain to runoff and routes it
with a causal, conservative triangular unit hydrograph. Snyder lag is measured
from the excess-rainfall centroid, so the triangle's continuous mode relative
to the storm start is the corrected lag plus half the 24-hour event duration.
Integrating that triangle over daily bins produces a kernel whose weights sum
to one. A `Ct = 2` regression gives 0, 1, 1 and 2 days on all three gate seeds
and passes the unchanged criteria. The instantaneous control returns runoff in
the rainfall row and therefore has zero lag. The inverse control deliberately
applies discrete lags of 1, 2, 1 and 3 days in ascending area order. Each lies
inside its individual broad Snyder bound, but the one-day reversal between the
second and third catchments fails the scaling criterion.

`flex_lumped` is the registered physical control. It retains its nonlinear
interception, soil-moisture partition, fast and slow reservoirs, then maps the
centroid-to-outlet channel distance to its native triangular `Tlag` through
the independent constant-celerity relation `travel time = Lc / c0`, with
`c0 = 1.0 m s-1` fixed before gate evaluation. The choice follows
[Beven (2020, Appendix equations A18--A20 and Figure A3)](https://doi.org/10.5194/hess-24-2655-2020),
which distinguishes mean water velocity from kinematic-wave celerity, derives
`c = dQ/dA`, and uses this value for an upland-channel example. The underlying
generalised kinematic routing method is described by
[Beven (1979)](https://doi.org/10.1029/WR015i005p01238). This value also lies
within the approximately 0.8--1.6 m s-1 regional range independently measured by
[Le Mesnil et al. (2021)](https://doi.org/10.5194/hess-25-1259-2021).
It is a deliberately simple synthetic prior rather than a universal river
constant. Applying it across the probe's synthetic catchment-size ladder is a
first-order benchmark, not a site-calibrated hydraulic relation or an exact
peak-time prediction. It shares no coefficient, length exponent or duration
correction with the criterion's Snyder relation.

FLEX's legacy triangle is indexed from the start of a generated-runoff row.
To assign the channel travel time to this native row-based routing kernel using
the source interval centre, the adapter uses
`Tlag = 2 * (travel time + dt / 2)`: the half-step aligns time origins, while
the factor two converts the target triangle mode to its full base. The final
timing still includes the model's interception, soil, fast-reservoir and
slow-reservoir effects. The routing kernel conserves generated runoff, while
the total runoff peak is evaluated independently against the probe's
duration-corrected Snyder range.

The fixed-celerity channel travel times are 0.057, 0.226, 0.899 and 1.850 days
from the small through xlarge variants. On all three registered gate seeds,
the complete FLEX response measures 0, 1, 2 and 3 days at the daily output
step; the small travel time is shorter than one output interval. The
corresponding criterion expectations are 0.631, 1.123, 2.250 and 3.335 days.
Every measured lag lies in its broad bound, the sequence never reverses, and
its three-day span exceeds the declared two-day minimum. The reproducible
archive row is in `models/result.csv`.

Existing models that do not declare consumption of both channel-length fields
are formally N/A (INCOMPATIBLE) on this probe. This avoids scoring a geometry
response against a model that never received the geometry as an input it can
use.

A model can pass by reproducing the end-to-end timing pattern without solving
a momentum equation. A calibrated model may fail because its runoff-generation
or storage response differs from this synthetic experiment. Passing supports
only the two stated lag claims; it does not establish a momentum balance,
identify channel routing in isolation, or demonstrate forecast accuracy in
real basins.

## Reproduction

From the repository root:

```sh
ht validate
ht gate --probe momentum/routing-lag-consistency
pytest -q
```

No generated forcing or outputs are committed. Repeating a seed reproduces the
complete case; changing the seed moves the event within its five-day placement
window and changes its depth and spinup weather.

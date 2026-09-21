# mass/spinup-cycle-invariance

## Question

Does a model reach the same seasonal state cycle after repeated identical
forcing, or does its answer still depend on how many cycles preceded it?

The generator draws one 365-day rain-dominated weather sequence from a seed
and maps it onto a fixed 3652-day Gregorian record. February 29 repeats
February 28's value, so every January 1 and month/day is aligned between
cycles even when a model reads the calendar. A preceding 365-day calendar
year is emitted as real spin-up, giving closure a state row before the scored
record. The `short`, `plus3`, and `long` variants start in 2001, 1998, and
1997, respectively. Their period offsets are five, eight, and nine cycles;
including the preceding spin-up year, six, nine, and ten complete prior cycles
are present before they all score the same non-leap calendar year, 2007. Their visible dates therefore differ,
but their static metadata and their scored `pr`, `tas`, and `pet` values are
identical. The host-only phase labels are stripped before any adapter runs.

The probe asks whether all three selected cycles agree. It does not require a
model to settle under a constant climate, as `mass/steady-state` does. The
forcing remains seasonal and variable; the expected attractor is a periodic
annual orbit, not a constant point.

## Physical statement

For a deterministic catchment driven by periodic forcing

$$
X(t+365\ \mathrm{d})=X(t),
$$

an initialized physical model should converge to a repeatable seasonal orbit.
Once it has reached that orbit, giving it three or four additional copies of
the same forcing cannot change the following copy:

$$
Q_N(t)\approx Q_{N+3}(t)\approx Q_{N+4}(t),
$$

and likewise for actual evapotranspiration and every reported physical store.

The probe first establishes the conditioning empirically: the exact bucket,
FLEX-Lumped, FLEX-Topo, and SAC-SMA/Snow-17 must all agree between the sixth
and tenth copies for every gate seed. Their pass is evidence that the five-cycle
period offset (six complete prior cycles including spin-up) is enough for this
generated catchment, not a claim that five or six years is a universal spin-up
length.

## Criterion

For each required variable and any store all three runs elect to report, compare
every pair of the 365-day evaluation series point by point. The criterion also compares
the sum of every reported physical storage column, because several small
changes can exceed the individual-store floor when combined. For a flux $Y$ and
pair $(i,j)$ use

$$
d_{Y,i,j} =
\frac{\max_t |Y_i(t)-Y_j(t)|}
{\max(\mathrm{mean}_t |Y_i(t)|, 0.05\ \mathrm{mm\ day^{-1}})}.
$$

Storage floors are calibrated per reported variable for this case: 1.0 mm
for `mrso`, 1.0 mm for the empty `snw` channel, and 0.01 mm for `canopy`;
`state_floor_mm: 1.0` remains the fallback for other stores and for the
aggregate. The largest variable-level departure must remain below 5 percent.
For the independent seed `20260912`, the scored-year means are about 118.5 mm
for `mrso`, 0 for `snw`, and 0.0565 mm for `canopy`. Thus the soil floor is
less than 1% of the soil signal, the snow floor is an explicit absolute
conditioning bound for an empty store, and the canopy floor leaves the 5%
relative rule active instead of allowing a 1 mm floor to dominate a 0.0565 mm
signal.

For `total_reported_storage`, the scale is the mean absolute sum of all stores
reported by the model. In these cases soil moisture (`mrso`) is the largest
store, so it naturally dominates that aggregate scale; the 1 mm fallback
still protects the score
when every reported store is nearly empty.

The criterion checks the complete visible driver schema and values for exact
equality before comparing outputs; the absolute `time` labels may differ so
that every variant can score the same calendar year. It rejects an optional
output that appears in only some variants, because the runs would then make
different reporting claims. The ordinary closure and state-bound preconditions
are evaluated over the complete post-spinup record in all variants, with the row immediately
before that record used as the initial state. This covers the full fixed-length
record; the `evaluation` label is reserved for the paired comparison of the
three phase-aligned 365-day cycles.

## What it catches

`reference_restless` is the exact conservative bucket with an internal
30-day clock that changes its recession coefficient. It closes its water
budget and keeps each reported store within bounds, but its clock has a
different phase among the three history lengths: period offsets five, eight,
and nine, or six, nine, and ten prior cycles including spin-up. It fails only
`spinup_cycle_invariance`.

This is the intended hidden-state failure: an unmodelled state keeps evolving
under a forcing sequence whose physical seasonal state has already repeated.
The probe does not identify the state or require a neural model to expose it;
it detects its observable consequence.

## What it does not claim

A pass does not prove that a model has the right groundwater, routing, or
energy mechanisms. A physical system can have a long transient, multiple
stable regimes, or a slow external driver; then the probe's periodic-attractor
assumption is not applicable until its case is redesigned or the model has
made that external driver explicit.

A FAIL means that the model did not reproduce the same evaluation cycle after
different amounts of identical prior history. The result identifies a failure
of cycle invariance under the prescribed spin-up, but it does not identify the
mechanism. A hidden state and a physical store that has not yet settled can
produce the same observable signature.

Longer spin-up can help distinguish these cases. In the packaged CWatM
diagnostic, a slow-groundwater configuration converged after 40 years of
repeated forcing. This is evidence for that configuration, not a universal
40-year requirement for every model or generated case.

Similarly, a model with no reported evaporation or storage is incomplete for
this mass probe rather than nonphysical. The required outputs are the minimum
needed to combine the cycle-invariance claim with closure and state-bound
preconditions.

The comparison is deliberately finite. A hidden clock with a period longer than the tested offsets can still evade a finite test. The coprime three- and four-cycle offsets eliminate the short integer-period alias demonstrated during review (periods 2, 3, 4, 6, and 8), but they do not prove convergence for every possible hidden process. Extend the case if a model exposes evidence of a longer calendar lock.

The deliberately permissive negative control `reference_cheater` also escapes
all criteria on four of 200 gate seeds (43, 103, 164 and 186). Those seeds are
reported as an explicit limitation of the discriminator; the other seeds and
the named `reference_restless` control still establish that the probe catches
the intended failure mode.

## Relationship to existing probes

- `mass/steady-state` holds forcing constant and asks whether all variables
  become constant. This probe keeps realistic seasonal variation and asks
  whether that variation repeats.
- `mass/time-origin-invariance` changes the absolute calendar while holding a
  finite forcing sequence fixed. This probe holds the forcing cycle fixed and
  changes the amount of prior identical history.
- `mass/catchment-closure` detects a net water-budget error in one run. A
  non-convergent hidden state can conserve water perfectly, so closure alone
  does not detect it.
- `mass/multi-decadal-drift` follows a single long repeated-forcing run for
  storage that slowly escapes physical bounds. This
  probe compares three phase-aligned seasonal cycles after different spin-up
  lengths, so it also catches a bounded internal clock that produces no
  storage drift.

## Long-spin-up diagnostic

The default case uses period offsets of five, eight, and nine cycles (six, nine,
and ten complete prior cycles including spin-up) because it must remain small
enough for the gate. To check the slow-store concern raised for CWatM, an
isolated copy with `recessionCoeff=0.001` was run with 40 years of repeated
history before the ten-year scored record. It differed by about 10% at the
prescribed spin-up, then passed after 40 years; the packaged default CWatM also
passed its 40-year diagnostic at 0.02% (the `gw` store). These are empirical
results for the two packaged configurations, not a universal 40-year
requirement for every model or generated case.

A 100-year version was also attempted in the same Docker environment. It reached
the project runner's 600-second per-model time budget before producing a result.
That timeout is a computational limit, not a FAIL or evidence of a physical
problem. The completed 40-year run is evidence that this packaged CWatM setup
had converged well before 100 years; a 100-year archive would require a faster
execution path or a separately approved runtime budget.

## Acceptance gate

Each seed runs the model three times (one per history length). The four
physical models must pass. `reference_restless` must fail specifically
on `spinup_cycle_invariance`; `reference_leaky`, `reference_cheater`, and
`reference_degenerate` must still fail `closure`, `state_bounds`, and
`non_degenerate` respectively. `tests/test_spinup_cycle_invariance.py`
repeats these checks on an independent seed and verifies that adapters receive
neither phase labels nor variant-specific case metadata, while their three
visible evaluation-year forcings remain aligned. The acceptance threshold and
per-variable floors are calibrated against the four gate references, whose
departures are 0.00% in the calibrated run. The
independent-seed diagnostic is intentionally conservative: `reference_restless`
has a 1.36 departure on `mrro`, while its `evspsbl` and `mrso` departures are
about 4.88% and 4.81%, respectively. The archived
LISFLOOD submission was 0.47%; it is reported as a model observation rather
than folded into the gate-reference range. The 5% engineering rule therefore
retains a conservative margin while the absolute floors keep near-zero
variables well-conditioned.

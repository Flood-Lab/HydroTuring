# mass/spinup-cycle-invariance

## Question

Does a model reach the same seasonal state cycle after repeated identical
forcing, or does its answer still depend on how many cycles preceded it?

The generator draws one 365-day rain-dominated weather sequence from a seed
and repeats its numerical `pr`, `tas`, and `pet` values exactly for nine
cycles. It runs the same model twice with the same visible forcing and the
same model seed. Host-only labels select the evaluation cycle after five
repetitions in `short` and after eight repetitions in `long`; the labels are
stripped before either adapter runs.

The probe asks whether the two selected cycles agree. It does not require a
model to settle under a constant climate, as `mass/steady-state` does. The
forcing remains seasonal and variable; the expected attractor is a periodic
annual orbit, not a constant point.

## Physical statement

For a deterministic catchment driven by periodic forcing

$$
X(t+365\ \mathrm{d})=X(t),
$$

an initialized physical model should converge to a repeatable seasonal orbit.
Once it has reached that orbit, giving it three additional copies of the same
forcing cannot change the following copy:

$$
Q_N(t)\approx Q_{N+3}(t),
$$

and likewise for actual evapotranspiration and every reported physical store.

The probe first establishes the conditioning empirically: the exact bucket,
FLEX-Lumped, FLEX-Topo, and SAC-SMA/Snow-17 must all agree between the sixth
and ninth copies for every gate seed. Their pass is evidence that five cycles
are enough for this generated catchment, not a claim that five years is a
universal spin-up length.

## Criterion

For each required variable and any store both runs elect to report, compare
the two 365-day evaluation series point by point. For a flux $Y$ use

$$
d_Y =
\frac{\max_t |Y_{N+3}(t)-Y_N(t)|}
{\max(\operatorname{mean}_t |Y_N(t)|, 0.05\ \mathrm{mm\ day^{-1}})}.
$$

For a storage use the same expression with a 1 mm floor. The largest
variable-level departure must remain below 5 percent. Floors make empty snow
or channel stores well-defined rather than granting them an accidental exact
score.

The criterion checks that the selected `pr`, `tas`, and `pet` cycles are
identical before comparing outputs. It rejects an optional output that appears
in only one variant, because the two runs would then make different reporting
claims.

## What it catches

`reference_restless` is the exact conservative bucket with an internal
30-day clock that changes its recession coefficient. It closes its water
budget and keeps each reported store within bounds, but its clock has a
different phase after five and eight annual cycles. It fails only
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

Similarly, a model with no reported evaporation or storage is incomplete for
this mass probe rather than nonphysical. The required outputs are the minimum
needed to combine the cycle-invariance claim with closure and state-bound
preconditions.

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
- The open `mass/multi-decadal-drift` proposal follows a single long
  repeated-forcing run for storage that slowly escapes physical bounds. This
  probe compares two phase-aligned seasonal cycles after different spin-up
  lengths, so it also catches a bounded internal clock that produces no
  storage drift.

## Acceptance gate

The four physical models must pass. `reference_restless` must fail specifically
on `spinup_cycle_invariance`; `reference_leaky`, `reference_cheater`, and
`reference_degenerate` must still fail `closure`, `state_bounds`, and
`non_degenerate` respectively. `tests/test_spinup_cycle_invariance.py` repeats
these checks on an independent seed and verifies that adapters receive neither
phase labels nor distinct case metadata.

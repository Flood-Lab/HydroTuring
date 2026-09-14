# mass/human-abstraction

**Law:** mass conservation · **Track:** synthetic · **Step:** daily

> The same weather twice — once natural, once with a prescribed net
> irrigation withdrawal. The abstracted water must leave the reported budget,
> and the difference between the two runs must account for exactly the
> prescribed volume.

## The question

Irrigation withdrawal is the oldest unaccounted sink in hydrology: water a
model was told left the catchment, that it never removed from any store or
flux it reports. Every single-run criterion in the suite is blind to the
omission, because nothing leaks — the water was never taken. The only way to
catch it is the paired question: run the same seed with and without the human
term and ask where the abstracted water went.

## The case

`generate.py` draws ten years of daily weather once from the recorded seed,
plus one spinup year. The `natural` variant carries the abstraction column as
zeros; the `irrigated` variant carries the schedule. Precipitation,
temperature and PET are byte-identical between the two.

The schedule is a smooth growing-season bump (May–September, sin² shape,
peak 0.5 mm/day net), giving ~38 mm/yr catchment-mean against ~830 mm/yr of
rain — a modest irrigated fraction, not an irrigation district. The amplitude
is bracketed from both sides: a blind model must fail by a wide margin (at
~380 mm over the record the residual is twenty times the tolerance), and the
withdrawal must stay honourable (at 2.5 mm/day peak the soil column empties
every summer and even an exact model under-removes 27-30% of the prescription;
at 0.5 mm/day the honest model's residual is ~0.0% on the three scored seeds,
and stays below ~1% on arbitrary seeds — far under the tolerance). The schedule
is a NET withdrawal (mm/day, gross abstraction minus return flow); that the
budget term is net rather than gross follows Perry (2007, doi:10.1002/ird.323)
and Döll & Siebert (2002, doi:10.1029/2001WR000355). No gross-to-net fraction
is applied: the generator's schedule is already net.

The human term is a visible forcing column `abstr` (mm/day), the prescribed
NET withdrawal. It is prescribed by the case, not reported by the model: an
external driver the model must honour, in the same sense that `pr` is. A model
honours it by removing the water from its stores and fluxes; a model that
routes the removal outside its reported budget declares it as a (negative)
`gwex`. Note that `human_abstraction` scores the difference between the runs
and does not include `gwex`, and nothing on this probe scores the irrigated
run's own `gwex`; a declared exchange that responds to the drier stores
therefore lands in the residual. That is harmless here (SAC-SMA's deep loss,
SIDE, is the whole of its 0.15-0.29%), but a model with a large
storage-dependent exchange would meet it first.

## The criteria

| criterion | what it asserts |
|---|---|
| `closure` | the control run closes its own budget to 5% of the rain — the same rule as `mass/catchment-closure` |
| `human_abstraction` | the paired budget difference plus the prescribed abstraction integrates to zero, to 5% of the **abstracted** volume (floor 5 mm) |
| `state_bounds` | every reported storage stays physical |

The denominator of `human_abstraction` is the abstracted volume, not total
precipitation — ten years of unabstracted rain would otherwise dilute the
signal by two orders of magnitude.

`closure` is a single-run criterion, so it is scored on the control (natural)
variant only. A PAIRED closure is deliberately not added: against `sum_pr` the
whole withdrawal is 4.3 to 4.6 percent of the window's rain on the gate seeds,
under the 5 percent limit, so a model that removed the water and declared none
of it would still close — the check could never fail for the reason it would
exist. A leak of that size in both runs already shows in the control's closure,
and a leak confined to the irrigated run moves the difference that
`human_abstraction` scores.

## What it catches

- **Abstraction ignored** — the model never reads the driver, both runs are
  identical, and the residual is the whole abstracted volume
  (`reference_abstraction_blind`).
- **Unreported sink** — the model removes the water from nothing it reports;
  each run may close while the difference does not move by the prescribed
  amount.
- **Inconsistent removal** — the right total on one seed, the wrong total on
  the next; every seed must pass.

A fixed-fraction abstraction — a set share of the withdrawal always taken
from runoff, whatever the state — would pass. That is accepted and stated
plainly: such a model still conserves and still accounts for the human term.
Where abstraction is drawn from is a partition question, and
`mass/antecedent-monotonicity` already probes partition state-dependence.

## Baselines

| model | expected | why |
|---|---|---|
| `reference_abstraction_blind` | FAIL `human_abstraction` | the exact bucket blind to the driver; the residual is the whole prescribed volume |
| `reference_leaky` | FAIL `closure` | its silent 15% sink shows in the control run's own closure |

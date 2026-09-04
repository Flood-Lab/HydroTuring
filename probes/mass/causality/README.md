# mass/causality

**A catchment cannot respond to rain that has not fallen yet.**

## The physics

Runoff on a given day is a function of the weather up to that day. There is
no mechanism by which tomorrow's storm reaches today's stream. So if two
records are identical up to some day and differ only from that day on, a
model's output must be identical up to that day too, to floating point,
and must differ afterwards, because a storm has to run off somewhere.

## What it measures

The `pulse` variant is the control's weather with one 40 mm storm added on
a single day in the second scored year. The criterion finds that day from
the forcing, compares every reported variable between the runs on all the
days before it, and requires the difference to be zero to a relative
tolerance of one part in a million. After the day, the runoff must answer
by at least 2 percent of the added rain, so a model cannot pass by
ignoring the rain altogether.

What fails it is any path from the future to the present:

- a bidirectional sequence model, or attention over the whole record;
- a centred filter, which is what `reference_anticipating` does: the
  reference bucket reporting its runoff smoothed over seven days, three of
  them ahead;
- statistics taken over the record the model is handed: normalising inputs
  by their mean, or deriving catchment attributes from the whole record.
  The Google Flood Hub adapter derives its climate attributes from the
  first year only for exactly this reason.

A stochastic model draws the same random numbers in both runs, because the
harness hands every variant of one seed the same model seed. Without that,
the comparison would see the sampling noise and call it a response.

## Results on the gate seeds

The exact bucket differs by exactly nothing before the storm and returns
0.47 of the added rain afterwards. `reference_anticipating` differs by 16
percent of its mean runoff before the storm, on the three days its window
reaches ahead. Only runoff is required, so every runoff model is asked.

## Criteria

| Criterion | What it asserts |
| --- | --- |
| `causality` | nothing reported changes before the added storm; runoff changes after it |
| `non_degenerate` | runoff varies with the weather, so a constant cannot pass |

## Baselines

- `must_pass: reference_bucket`
- `must_fail: reference_anticipating` on `causality`

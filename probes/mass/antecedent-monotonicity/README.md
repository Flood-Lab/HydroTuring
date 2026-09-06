# mass/antecedent-monotonicity

Memory is the property neural hydrology is proudest of, and no probe had
asked for its sign. This one does. The same storm falls after a dry month
and after a wet one; over the month that follows, the wetter catchment must
run off more of it, by at least two percent of the storm, and by no more
than the extra water it was given.

The lower bound is what a catchment does: fuller stores pass a larger
share of what arrives, whatever the mechanism, saturation excess, a beta
partition, a nonlinear reservoir. The upper bound is mass: the difference
between the two hydrographs over the window cannot exceed the water that
distinguishes the two cases.

| Criterion | Asserts |
| --- | --- |
| `antecedent_monotonicity` | runoff over the 30 days after the storm is larger on the wet catchment by at least 0.02 of the storm and at most the 120 mm of antecedent rain |
| `non_degenerate` | runoff varies with the weather |

Must-fail: `reference_cheater`, which runs off a fixed share of each day's
rain and so answers both storms identically.

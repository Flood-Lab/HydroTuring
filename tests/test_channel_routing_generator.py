"""The recessions this probe scores have to be the ones the storms leave.

`momentum/channel-routing-mass` asks what the channel store does on a step where
nothing enters the reach, and it asks it of a record built from storms each
followed by a long dry spell. The generator draws stochastic weather first and
used to add the storms on top without ever zeroing the spells between them: the
record kept raining at `p_wet` ~ 0.25 throughout, the `_recessions` mask the
module documents was never called, and the steps the criterion scored were
chance runs of four dry days — `1095 * 0.75**4`, about 346 of them — rather
than the recessions the design describes.

These tests pin the property rather than the storm table, so lengthening a storm
or moving an event cannot quietly take the dry spells away again.
"""

from __future__ import annotations

import numpy as np
import pytest

from hydroturing import registry
from hydroturing.harness import load_generator
from hydroturing.seeds import gate_seeds

PROBE = "momentum/channel-routing-mass"

# The storm table spaces the events 215-240 days apart, so a spell that ever
# falls below this means the generator stopped zeroing the rain over them.
MIN_SPELL_DAYS = 200
# What `recession_drainage` counts on a record whose spells are real: the
# recessions the four events leave, not a handful of chance dry days.
MIN_SCORED_RECESSION_STEPS = 500


@pytest.fixture(scope="module")
def generator():
    return load_generator(registry.find_probe(PROBE))


@pytest.fixture(scope="module")
def seeds():
    return gate_seeds(PROBE, 3)


def _dry_run_lengths(pr: np.ndarray) -> list[int]:
    """Lengths of the maximal runs of rainless steps, in order."""
    runs, current = [], 0
    for value in pr:
        if value == 0.0:
            current += 1
        elif current:
            runs.append(current)
            current = 0
    if current:
        runs.append(current)
    return runs


def test_each_post_storm_spell_is_rainless(generator, seeds):
    """The rain is zero over every spell, not merely sparse.

    This is the property the review found missing: the storms were added on top
    of weather that kept raining, so the spells carried rain of their own.
    """
    for seed in seeds:
        frame, _ = generator.generate(seed)
        pr = frame["pr"].to_numpy()
        mask = generator._recessions(len(pr))
        assert int(mask.sum()) >= MIN_SPELL_DAYS
        assert (pr[mask] == 0.0).all(), "a post-storm spell still carries rain"


def test_the_dry_spells_are_long_enough_to_be_the_generators_own(generator, seeds):
    """An unbroken rainless run the weather could not have produced.

    Chance dryness at `p_wet` of about 0.25 gives runs of two or three days and
    almost never a hundred, so the length of the longest run separates a record
    the generator shaped from one it merely drew.
    """
    for seed in seeds:
        frame, _ = generator.generate(seed)
        longest = max(_dry_run_lengths(frame["pr"].to_numpy()), default=0)
        assert longest >= MIN_SPELL_DAYS, (
            f"the longest rainless run is {longest} days; the probe scores what a "
            "storm leaves behind, so the spells between storms have to be dry"
        )


def test_every_storm_still_arrives(generator, seeds):
    """Zeroing the spells must not touch the events themselves."""
    for seed in seeds:
        frame, _ = generator.generate(seed)
        pr = frame["pr"].to_numpy()
        for start, days, rate in generator.STORMS:
            lo = generator.SPINUP_DAYS + start
            total = float(pr[lo:lo + days].sum())
            assert total >= rate * days, (
                f"the storm at step {lo} carries {total:.1f} mm, less than the "
                f"{rate * days:.1f} mm its own spell places"
            )


def test_the_criterion_scores_a_long_recession_not_a_handful_of_days(generator, seeds):
    """What `recession_drainage` counts on the record the generator builds.

    The mask mirrors `src/hydroturing/criteria/drainage.py`: a step is in
    recession when it and the three before it were rainless (`settle_days`), and
    the first thirty steps are excluded as startup. Before the fix this counted
    about 346 steps out of 1460 — exactly what chance gives.
    """
    settle_days, settle_steps, dry_threshold = 3, 30, 0.05
    for seed in seeds:
        frame, _ = generator.generate(seed)
        pr = frame["pr"].to_numpy()
        dry = pr <= dry_threshold
        recession = np.copy(dry)
        for k in range(1, settle_days + 1):
            recession[k:] &= dry[:-k]
        recession[:settle_steps] = False
        scored = int(recession.sum())
        assert scored >= MIN_SCORED_RECESSION_STEPS, (
            f"only {scored} steps are in recession on seed {seed}; the record is "
            "supposed to hand the criterion the spells between the storms"
        )

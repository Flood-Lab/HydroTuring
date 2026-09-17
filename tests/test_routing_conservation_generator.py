"""The record this probe runs on has to be the one its bound needs.

`momentum/routing-conservation` allows a channel at most `max_lag_days` times
the largest runoff of the preceding month, and it reads that peak over a
rolling 30-day window. On a record that never goes dry the window always holds
a storm, the ceiling stays at a storm-time level, and a kernel that keeps a
fraction of a percent of every day's runoff hides inside it — a residue of a
few millimetres against an allowance of hundreds.

Four storms, each followed by a long rainless spell, are what change that. Over
a spell the runoff falls to the baseflow and then to nothing, the peak decays
out of the window, and the ceiling comes down with it. That is when a residue
of a few millimetres stops being invisible.

These tests pin the weather's properties and the separation they buy, rather
than the storm table: moving an event or lengthening a storm must not quietly
take the dry spells away again.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.harness import build_case, load_generator
from hydroturing.protocol import RunResult
from hydroturing.runner import get_runner
from hydroturing.seeds import gate_seeds

PROBE = "momentum/routing-conservation"

# The bound reads its peak over `lookback_days`, so a spell shorter than that
# never lets the ceiling fall. The storm table spaces the events at least 200
# days apart, and the spell after the last one runs to the end of the record.
MIN_SPELL_DAYS = 200
# A storm has to be a storm: the generator adds 30 to 48 mm/day on top of the
# weather, and the weather's own daily rain is a gamma draw with mean ~9 mm.
MIN_STORM_MM = 30.0


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe(PROBE)


@pytest.fixture(scope="module")
def generator(probe):
    return load_generator(probe)


@pytest.fixture(scope="module")
def seeds(probe):
    return [int(s) for s in gate_seeds(PROBE, 3)]


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


def _params(probe, name: str) -> dict:
    return dict(next(c.params for c in probe.criteria if c.name == name))


@pytest.mark.parametrize("which", [0, 1, 2])
def test_every_storm_is_followed_by_a_spell_longer_than_the_bound_s_window(
    generator, seeds, which
):
    forcing, _static = generator.generate(seeds[which])
    runs = _dry_run_lengths(forcing["pr"].to_numpy(dtype=float))
    assert max(runs) >= MIN_SPELL_DAYS, (
        f"the longest rainless spell is {max(runs)} days; the bound's window is "
        "30, so a short one leaves the ceiling at a storm-time level and a small "
        "leak inside it"
    )
    assert len(runs) >= 4, f"only {len(runs)} rainless runs in the record"


@pytest.mark.parametrize("which", [0, 1, 2])
def test_the_spells_follow_storms_rather_than_dry_weather(generator, seeds, which):
    """A spell is only a recession if something filled the store first."""
    forcing, _static = generator.generate(seeds[which])
    assert forcing["pr"].max() >= MIN_STORM_MM, (
        f"the wettest day is {forcing['pr'].max():.2f} mm, which is the weather "
        "rather than one of the named storms"
    )


def test_the_generator_is_deterministic(generator, seeds):
    first, static = generator.generate(seeds[0])
    second, static_again = generator.generate(seeds[0])
    assert first.equals(second)
    assert static == static_again


def test_a_kernel_summing_to_0999_is_caught_by_this_weather(probe):
    """The case the weather changed for.

    A kernel summing to 0.999 releases all but a tenth of a percent of each
    day's runoff, so it leaves a residue of `0.001 x the runoff so far` behind
    it. That residue is built here on the exact bucket's own runoff, so it is
    the fault and nothing else: on a record whose peak never leaves the
    window the allowance would be hundreds of millimetres and the residue would
    pass, and on this one it fails.
    """
    model = registry.find_model("reference_bucket")
    case = build_case(probe, int(gate_seeds(PROBE, 3)[0]))
    run = get_runner(model).run(model, probe, case, Path(tempfile.mkdtemp()))

    table = run.table.copy()
    residue = 0.001 * table["mrro"].cumsum()
    table["channel"] = table["channel"] + residue
    leaking = RunResult(case=run.case, table=table, meta=run.meta, wall_seconds=0.0)

    params = _params(probe, "routing_conservation")
    caught = get("routing_conservation")(leaking, probe, params)
    exact = get("routing_conservation")(run, probe, params)

    # An exact router passes on the same record, so the failure is the residue.
    assert exact.status == PASS
    assert caught.status == FAIL

    # And it is caught while holding far less than a storm-time ceiling would
    # allow, which is the whole point: `max_lag_days` of the record's own peak
    # runoff is orders of magnitude above the residue, so a criterion reading
    # the peak alone would pass it.
    peak = table["mrro"].to_numpy(dtype=float)
    storm_time_ceiling = params["max_lag_days"] * (1.0 + params["tolerance"]) * peak.max()
    assert residue.max() < 0.02 * storm_time_ceiling, (
        f"the residue reaches {residue.max():.3f} mm against a storm-time "
        f"ceiling of {storm_time_ceiling:.1f} mm; this test is only meaningful "
        "while the residue is far below it"
    )

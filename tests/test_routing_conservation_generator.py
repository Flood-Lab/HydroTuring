"""The weather this probe runs on, and the separation it buys.

`momentum/routing-conservation` allows a channel at most `max_lag_days` times
the largest runoff of the preceding month, read over a rolling 30-day window.
On a record that never goes dry the window always holds a storm, the ceiling
stays at a storm-time level, and a kernel that keeps a fraction of a percent of
every day's runoff hides inside it — a residue of under a millimetre against an
allowance of hundreds.

Four storms, each followed by a long rainless spell, are what change that. Over
a spell the runoff falls to the baseflow and then to nothing, the peak decays
out of the window, and the ceiling comes down with it. That is when the residue
stops being invisible.

These tests pin the weather's properties, the separation they buy, and the size
of the absolute allowance, rather than the storm table: moving an event or
lengthening a storm must not quietly take the dry spells away, and neither must
deleting the storm injection, which the background weather would otherwise cover
for.
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
# A storm is not the weather: the table adds 30 to 48 mm/day over five to seven
# days, so a window that carries less than half of that has lost the injection.
MIN_STORM_SHARE = 0.5
# Outside the spells the record is temperate, rain-dominated weather. The band
# is wide on purpose: it pins the character the description claims, not the seed.
MIN_ANNUAL_RAIN_MM = 500.0
MAX_ANNUAL_RAIN_MM = 1300.0
# What `wflow_sbm` keeps above the hydrograph bound through the droughts, in mm.
# The probe README gives the seeds, the command and the calculation.
WFLOW_NEED_MM = 0.0026


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe(PROBE)


@pytest.fixture(scope="module")
def generator(probe):
    return load_generator(probe)


@pytest.fixture(scope="module")
def seeds(probe):
    return [int(s) for s in gate_seeds(PROBE, 3)]


@pytest.fixture(scope="module")
def params(probe):
    return dict(next(c.params for c in probe.criteria if c.name == "routing_conservation"))


class _Router:
    """A three-day kernel summing to `total`, run on a model's own runoff.

    The retained share is left in the store and the released share is reported,
    both derived from the same generation, so the fault is a router that keeps a
    tenth of a percent of every day's runoff rather than a store that gains
    water from nowhere.
    """

    WEIGHTS = np.array([0.5, 1.0 / 3.0, 1.0 / 6.0])

    def __init__(self, generated: np.ndarray, total: float):
        weights = self.WEIGHTS / self.WEIGHTS.sum() * total
        released = np.zeros_like(generated)
        for k, weight in enumerate(weights):
            released[k:] += weight * generated[: len(generated) - k]
        self.released = released
        self.store = np.cumsum(generated) - np.cumsum(released)


def _bucket_run(probe, seed):
    model = registry.find_model("reference_bucket")
    case = build_case(probe, int(seed))
    return get_runner(model).run(model, probe, case, Path(tempfile.mkdtemp()))


def _altered(run, **columns) -> RunResult:
    table = run.table.copy()
    for name, values in columns.items():
        table[name] = values
    return RunResult(case=run.case, table=table, meta=run.meta, wall_seconds=0.0)


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


def test_every_storm_is_followed_by_a_spell_longer_than_the_bound_s_window(
    generator, seeds
):
    for seed in seeds:
        forcing, _static = generator.generate(seed)
        runs = _dry_run_lengths(forcing["pr"].to_numpy(dtype=float))
        assert max(runs) >= MIN_SPELL_DAYS, (
            f"seed {seed}: the longest rainless spell is {max(runs)} days; the "
            "bound's window is 30, so a short one leaves the ceiling at a "
            "storm-time level and a small leak inside it"
        )
        assert len(runs) >= 4, f"seed {seed}: only {len(runs)} rainless runs"


def test_each_storm_supplies_its_rain_and_its_spell(generator, seeds):
    """The storms are pinned event by event, not by the record's wettest day.

    The background weather reaches 30 mm/day on its own on some seeds, so a
    record-wide maximum cannot tell the injection from the weather — deleting
    the storms outright would still clear it, because the dry mask is applied
    whether or not the rain was added.
    """
    spinup = generator.SPINUP_DAYS
    table = generator.STORMS
    for seed in seeds:
        forcing, _static = generator.generate(seed)
        pr = forcing["pr"].to_numpy(dtype=float)
        for index, (start, days, rate) in enumerate(table):
            lo = spinup + start
            hi = min(lo + days, len(pr))
            supplied = float(pr[lo:hi].sum())
            assert supplied >= MIN_STORM_SHARE * days * rate, (
                f"seed {seed}, storm {index}: {supplied:.1f} mm over {days} days, "
                f"less than half of the {days * rate:.0f} mm the table adds"
            )
            stop = spinup + table[index + 1][0] if index + 1 < len(table) else len(pr)
            spell = pr[hi:stop]
            assert len(spell) > 0, f"seed {seed}, storm {index}: no spell follows it"
            assert np.all(spell == 0.0), (
                f"seed {seed}, storm {index}: {int((spell > 0).sum())} of the "
                f"{len(spell)} steps after it carry rain"
            )


def test_the_record_outside_the_spells_stays_temperate_rainfall(generator, seeds):
    """The description says about 880 mm/yr, and the dry spells do not pin that:
    a mutant that doubles the wet-day frequency leaves the spells just as dry,
    so the total is checked here to keep the record the one the probe describes.
    """
    spinup = generator.SPINUP_DAYS
    for seed in seeds:
        forcing, _static = generator.generate(seed)
        pr = forcing["pr"].to_numpy(dtype=float)
        annual = float(pr.sum()) / ((len(pr) - spinup) / 365.0)
        assert MIN_ANNUAL_RAIN_MM <= annual <= MAX_ANNUAL_RAIN_MM, (
            f"seed {seed}: {annual:.0f} mm/yr over the record, outside the "
            f"[{MIN_ANNUAL_RAIN_MM:.0f}, {MAX_ANNUAL_RAIN_MM:.0f}] band"
        )


def test_the_generator_is_deterministic(generator, seeds):
    first, static = generator.generate(seeds[0])
    second, static_again = generator.generate(seeds[0])
    assert first.equals(second)
    assert static == static_again


def test_a_kernel_summing_to_0999_is_caught_on_every_gate_seed(probe, params, seeds):
    """The case the weather changed for, run as the fault it names.

    A kernel summing to 0.999 releases all but a tenth of a percent of each
    day's runoff. The released series and the store it leaves behind are built
    here from the exact bucket's own runoff, so what is scored is a router that
    retains a tenth of a percent rather than a store that invents water. It
    fails on every gate seed, and the exact kernel passes on the same record.
    """
    bound = get("routing_conservation")
    ratios = []
    for seed in seeds:
        run = _bucket_run(probe, seed)
        generated = run.table["mrro"].to_numpy(dtype=float)

        leaking = _Router(generated, 0.999)
        caught = bound(
            _altered(run, mrro=leaking.released, channel=leaking.store), probe, params
        )
        assert caught.status == FAIL, f"seed {seed}: the 0.999 kernel passed"
        ratios.append(caught.value)

        exact = _Router(generated, 1.0)
        assert bound(
            _altered(run, mrro=exact.released, channel=exact.store), probe, params
        ).status == PASS, f"seed {seed}: the exact kernel failed"

    # The figure the README quotes has to be the one that binds: every seed has
    # to pass, so it is the smallest ratio across them, not the largest.
    assert min(ratios) > 4.0, f"ratios {ratios} leave less than 4x of margin"
    assert max(ratios) < 5.0, f"ratios {ratios} are larger than the README says"


def test_the_residue_is_far_inside_what_a_peak_reading_bound_would_allow(probe, params, seeds):
    """Why the dry spells are the mechanism rather than a detail.

    The residue the 0.999 kernel leaves is under a millimetre, while
    `max_lag_days` of the record's own peak runoff is hundreds. A bound reading
    its peak from a storm would pass the same leak, so what fails it here is the
    ceiling coming down during the drains.
    """
    run = _bucket_run(probe, seeds[0])
    generated = run.table["mrro"].to_numpy(dtype=float)
    # The store carries the lag's water in transit as well as the retention, so
    # the residue is the difference between the two kernels' stores: what the
    # 0.999 kernel holds that an exact one does not.
    residue = _Router(generated, 0.999).store - _Router(generated, 1.0).store
    storm_time_ceiling = (
        params["max_lag_days"] * (1.0 + params["tolerance"]) * generated.max()
    )
    assert residue.max() < 0.02 * storm_time_ceiling, (
        f"the residue reaches {residue.max():.3f} mm against a storm-time "
        f"ceiling of {storm_time_ceiling:.1f} mm; this test is only meaningful "
        "while the residue is far below it"
    )


def test_the_allowance_is_what_lets_an_honest_residue_through(probe, params, seeds):
    """The parameter has an effect, and this is the size of it.

    A store sitting `WFLOW_NEED_MM` above the hydrograph bound at one step is
    what `wflow_sbm` keeps through the droughts. With the default 1e-6 mm it
    fails; with this probe's allowance it passes. Nothing else pins that the
    parameter is read at all.
    """
    run = _bucket_run(probe, seeds[0])
    runoff = run.table["mrro"].to_numpy(dtype=float)
    peak = np.array([runoff[max(0, i - 30): i + 1].max() for i in range(len(runoff))])
    proportional = params["max_lag_days"] * peak * (1.0 + params["tolerance"])
    channel = np.zeros_like(runoff)
    middle = len(channel) // 2
    channel[middle] = proportional[middle] + WFLOW_NEED_MM

    bound = get("routing_conservation")
    altered = _altered(run, channel=channel)
    assert bound(altered, probe, {**params, "min_allowance_mm": 1e-6}).status == FAIL
    assert bound(altered, probe, params).status == PASS


def test_one_row_cannot_defeat_the_bound(probe, params, seeds):
    """The allowance comes from the runoff, not from the store.

    A criterion whose scale is read off the store itself is at the mercy of a
    single row: one infinite value makes its floor infinite and the criterion
    passes. This bound is `max_lag_days` of the recent runoff, so neither an
    infinite store nor a finite spike can raise it — both have to fail. The row
    goes inside the scored window: a row in the spinup is not scored, though the
    last spinup row does reach every criterion as its `state0`
    (`criteria/base.py`), so a non-finite value there is visible to the criteria
    that read it as an initial state, even though this one works on the window's
    own steps.
    """
    run = _bucket_run(probe, seeds[0])
    channel = run.table["channel"].to_numpy(dtype=float).copy()
    scored = len(channel) // 2
    bound = get("routing_conservation")

    spiked = channel.copy()
    spiked[scored] = 1e6
    assert bound(_altered(run, channel=spiked), probe, params).status == FAIL

    infinite = channel.copy()
    infinite[scored] = np.inf
    assert bound(_altered(run, channel=infinite), probe, params).status == FAIL


def test_a_non_finite_or_unknown_parameter_is_rejected(probe, params, seeds):
    """A parameter that would disable the bound is refused, not read.

    `min_allowance_mm: .inf` passes every model and `.nan` compares false
    against everything, so both would look satisfied. A misspelled name is
    refused too rather than silently dropped: an author writing
    `min_allowance` means to change the bound, and would not learn that nothing
    happened.
    """
    run = _bucket_run(probe, seeds[0])
    bound = get("routing_conservation")
    for value in (np.inf, np.nan, -1.0):
        with pytest.raises(ValueError, match="min_allowance_mm"):
            bound(run, probe, {**params, "min_allowance_mm": value})
    with pytest.raises(ValueError, match="max_lag_days"):
        bound(run, probe, {**params, "max_lag_days": 0.0})
    with pytest.raises(ValueError, match="unknown"):
        bound(run, probe, {**params, "min_allowance": 0.05})

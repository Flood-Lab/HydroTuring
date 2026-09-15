"""Counterexamples for the recession-drainage criterion, built by hand.

`recession_drainage` is a one-sided statement about a channel store: on steps
where the forcing has been rainless, the store must not rise. The acceptance
gate drives it with the physical models (which drain) and with
`reference_stuck_router` (which fills), but a gate only exercises the paths its
baselines happen to take. These tests build the frames directly so the branches
the gate never reaches — the relative floor, the `settle_days` window, the
startup exclusion and the degenerate record — are pinned by construction.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.protocol import Case, RunResult


def build(channel, pr, spinup=0, timestep="PT1D"):
    """A RunResult carrying the store and the forcing the criterion reads.

    `channel` is the water in transit and `pr` the rain, both per step. Nothing
    downstream reads the timestamps beyond their count.
    """
    channel = np.asarray(channel, dtype=float)
    pr = np.asarray(pr, dtype=float)
    n = len(channel)
    time = pd.date_range("2001-01-01", periods=n, freq="D")
    case = Case(
        probe_id="momentum/channel-routing-mass", seed=1,
        forcing=pd.DataFrame({"time": time, "pr": pr}),
        static={}, spinup_steps=spinup, timestep=timestep,
    )
    return RunResult(case, pd.DataFrame({"time": time, "channel": channel}), {}, 0.0)


def storm_then_dry(n, storm_at, storm_days=3, rate=40.0):
    """Rain on the named steps only: one storm, then a long recession."""
    pr = np.zeros(n)
    pr[storm_at:storm_at + storm_days] = rate
    return pr


def test_a_store_that_drains_on_recessions_passes():
    """The honest case: it fills in the storm and only falls afterwards."""
    n = 200
    pr = storm_then_dry(n, 10)
    channel = np.zeros(n)
    channel[13:] = np.linspace(50.0, 0.0, n - 13)
    result = get("recession_drainage")(build(channel, pr), None, {})
    assert result.status == PASS
    assert result.diagnostics["rising_steps"] == 0
    assert "drains on recession steps" in result.message


def test_a_store_that_fills_from_nothing_on_recessions_fails():
    """The leaky-router case: it rises on steps where nothing enters the reach."""
    n = 200
    pr = storm_then_dry(n, 10)
    channel = np.zeros(n)
    channel[13:] = np.linspace(1.0, 100.0, n - 13)
    result = get("recession_drainage")(build(channel, pr), None, {})
    assert result.status == FAIL
    assert result.diagnostics["rising_steps"] > 0
    assert "cannot rise" in result.message


def test_the_measure_is_a_share_of_recession_steps_not_a_worst_case():
    """One isolated rise in a settled record must not decide it.

    A discrete unit hydrograph renormalising over a partial history moves the
    store by a fraction of a percent of its scale on a single step. Scoring the
    worst step would fail the model for arithmetic; the share does not.
    """
    n = 400
    pr = storm_then_dry(n, 10)
    channel = np.zeros(n)
    channel[13:] = np.linspace(80.0, 20.0, n - 13)
    channel[200] += 0.5  # one step up, far above the floor, but a single step
    result = get("recession_drainage")(build(channel, pr), None, {})
    assert result.status == PASS
    assert result.diagnostics["rising_steps"] == 1


def test_a_rise_below_the_relative_floor_is_not_counted():
    """Floating-point dust is not a rise."""
    n = 300
    pr = storm_then_dry(n, 10)
    channel = np.full(n, 100.0)
    channel[13:] += 1e-9 * np.sin(np.arange(n - 13))
    result = get("recession_drainage")(build(channel, pr), None, {})
    assert result.status == PASS
    assert result.diagnostics["rising_steps"] == 0


def test_the_storm_tail_is_not_scored_until_the_reach_settles():
    """A step only counts as recession after `settle_days` dry steps.

    The fastest surface response of an event is still arriving a day or two
    after the rain stops, and a store is allowed to move then.
    """
    n = 200
    pr = storm_then_dry(n, 50)
    channel = np.zeros(n)
    channel[53:56] = [5.0, 10.0, 15.0]  # rises while the tail is still draining
    channel[56:] = np.linspace(15.0, 0.0, n - 56)
    result = get("recession_drainage")(build(channel, pr), None, {})
    assert result.status == PASS
    assert result.diagnostics["rising_steps"] == 0


def test_the_startup_transient_of_the_window_is_excluded():
    """A store filling from zero at the start is arithmetic, not a violation.

    Every recession step here rises through the excluded startup and only then
    falls, so counting the startup would fail a record whose scored stretch is
    flawless.
    """
    n = 100
    pr = np.zeros(n)  # every step is dry
    channel = np.concatenate([
        np.linspace(0.0, 100.0, 30),
        np.linspace(100.0, 0.0, n - 30),
    ])
    result = get("recession_drainage")(build(channel, pr), None, {})
    assert result.status == PASS
    assert result.diagnostics["rising_steps"] == 0


def test_a_record_with_no_recession_steps_is_a_failure_not_a_pass():
    """A case the criterion cannot score is not a case it passed."""
    n = 200
    pr = np.full(n, 5.0)  # it never stops raining
    channel = np.linspace(0.0, 50.0, n)
    result = get("recession_drainage")(build(channel, pr), None, {})
    assert result.status == FAIL
    assert "no recession steps" in result.message
    assert result.diagnostics["recession_steps"] == 0


def test_the_threshold_is_the_declared_share():
    """`max_rising_fraction` decides, not the raw count."""
    n = 400
    pr = storm_then_dry(n, 10)
    channel = np.zeros(n)
    # Half the recession rises: past a 0.25 limit, inside a 0.75 one.
    channel[13:] = np.where(np.arange(n - 13) % 2 == 0, 50.0, 49.0)
    strict = get("recession_drainage")(build(channel, pr), None,
                                       {"max_rising_fraction": 0.25})
    loose = get("recession_drainage")(build(channel, pr), None,
                                      {"max_rising_fraction": 0.75})
    assert strict.status == FAIL
    assert loose.status == PASS


def test_the_store_column_must_be_present():
    """A missing store is a failure, not a silent pass."""
    n = 100
    run = build(np.linspace(0.0, 10.0, n), np.zeros(n))
    del run.table["channel"]
    assert get("recession_drainage")(run, None, {}).status == FAIL


def test_the_forcing_must_carry_rain():
    """Without `pr` the recession cannot be located, so it must not pass."""
    n = 100
    time = pd.date_range("2001-01-01", periods=n, freq="D")
    case = Case(
        probe_id="momentum/channel-routing-mass", seed=1,
        forcing=pd.DataFrame({"time": time}), static={},
        spinup_steps=0, timestep="PT1D",
    )
    table = pd.DataFrame({"time": time, "channel": np.linspace(10.0, 0.0, n)})
    run = RunResult(case, table, {}, 0.0)
    assert get("recession_drainage")(run, None, {}).status == FAIL

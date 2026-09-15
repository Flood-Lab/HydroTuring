"""Counterexamples for the rating criteria, on ratings built by hand.

`rating_monotonic` and `rating_loop` read a `stage` diagnostic against a
discharge or a store, and the shape they are looking for — a single-valued
rating on one hand, a correctly-signed hysteresis loop on the other — is
difficult to produce from a reference model alone. The acceptance gate runs
them on the four must-pass baselines and three broken gauges, but a gate only
ever drives the paths its baselines happen to take: a branch no baseline
reaches can be wrong and stay green.

These tests build the frames directly, so every branch of `rating_loop` is
exercised by construction: a correctly-signed loop (which must PASS *and*
report the loop it found), an inverted loop (which must FAIL), a
single-valued rating with no loop to read, and a constant gauge with nothing
to read at all.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.criteria.rating import _rating_scatter
from hydroturing.protocol import Case, RunResult


def build(q, stage, channel=None, spinup=0, timestep="PT1D", discharge=True):
    """A RunResult carrying just the columns the rating criteria read.

    `q` is the discharge (`dis`), `stage` the gauge reading and `channel` the
    store the limbs are paired on. Times are synthetic; nothing downstream of
    the criteria reads them beyond the count. `discharge=False` drops the `dis`
    column, which is how a model that reports only a store is represented.
    """
    q = np.asarray(q, dtype=float)
    stage = np.asarray(stage, dtype=float)
    n = len(q)
    table = {"time": pd.date_range("2001-01-01", periods=n, freq="D"),
             "stage": stage}
    if discharge:
        table["dis"] = q
    if channel is not None:
        table["channel"] = np.asarray(channel, dtype=float)
    case = Case(
        probe_id="momentum/stage-discharge-monotonic", seed=1,
        forcing=pd.DataFrame({"time": table["time"]}),
        static={}, spinup_steps=spinup, timestep=timestep,
    )
    return RunResult(case, pd.DataFrame(table), {}, 0.0)


def loop_run(sign=+1, separation=0.5, periods=6, n=600):
    """A rating that traces a hysteresis loop of the requested sign.

    The discharge sweeps a sinusoid, so the record has both limbs. The stage
    is a function of discharge *plus* a term that depends on whether the flow
    is rising or falling: with `sign=+1` the gauge sits higher on the fall
    (the physical direction), with `sign=-1` higher on the rise (inverted).

    The store is set to the discharge multiplied by a slowly-varying factor,
    so pairing on it is a genuine pairing rather than a relabelling of the
    discharge — the case the probe's `abscissa: store` is meant to read.
    """
    t = np.linspace(0.0, periods * 2 * np.pi, n)
    q = 1.0 + 0.5 * np.sin(t)
    dq = np.gradient(q)
    rising = dq > 0
    # A stage that is a function of the flow, with the loop on top of it.
    stage = 2.0 * q + sign * separation * np.where(rising, -0.5, +0.5)
    channel = q * (1.0 + 0.3 * np.sin(0.25 * t))
    return build(q, stage, channel=channel)


def test_a_correctly_signed_loop_passes_and_reports_the_loop_it_found():
    """The branch that asserts "the loop runs the right way" must be reachable.

    With the loop clipped at zero before the size test, every correctly-signed
    loop reads as zero and takes the no-loop escape, so the criterion can say
    "no loop" or "wrong loop" but never "a loop, the right way". This pins the
    reachable branch.
    """
    result = get("rating_loop")(loop_run(sign=+1), None, {"bins": 12, "tolerance": 0.0})
    assert result.status == PASS
    # It reports a real loop, not the zero the escape path reports.
    assert result.diagnostics["loop_size_m"] > 0.1
    assert result.diagnostics["loop_median_m"] < 0.0
    assert result.diagnostics["single_valued"] is False
    assert "below the falling limb" in result.message


def test_an_inverted_loop_fails():
    result = get("rating_loop")(loop_run(sign=-1), None, {"bins": 12, "tolerance": 0.0})
    assert result.status == FAIL
    assert result.diagnostics["loop_m"] > 0.1
    assert "wrong way" in result.message


def test_a_single_valued_rating_has_no_loop_to_read():
    """A stage that is an exact function of the abscissa has no hysteresis."""
    t = np.linspace(0.0, 12 * np.pi, 600)
    q = 1.0 + 0.5 * np.sin(t)
    result = get("rating_loop")(build(q, 3.0 * q), None, {"bins": 12, "tolerance": 0.0})
    assert result.status == PASS
    assert result.diagnostics["single_valued"] is True
    assert "single-valued" in result.message


def test_the_direction_of_a_loop_is_not_read_from_a_store_that_never_moves():
    """A constant store is no abscissa, so the limbs are read on `dis` instead."""
    t = np.linspace(0.0, 12 * np.pi, 600)
    q = 1.0 + 0.5 * np.sin(t)
    dq = np.gradient(q)
    stage = 2.0 * q + 0.5 * np.where(dq > 0, -0.5, +0.5)
    result = get("rating_loop")(build(q, stage, channel=np.zeros_like(q)), None,
                                {"bins": 12, "tolerance": 0.0, "abscissa": "store"})
    assert result.status == PASS
    assert result.diagnostics["abscissa"] == "dis"


def test_stage_missing_is_a_failure_not_a_pass():
    run = build([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
    del run.table["stage"]
    assert get("rating_loop")(run, None, {}).status == FAIL
    assert get("rating_monotonic")(run, None, {}).status == FAIL


def test_rating_monotonic_catches_a_stage_that_steps_down():
    q = np.linspace(0.1, 5.0, 200)
    stage = 2.0 * q
    stage[100:] -= 3.0  # one downward step that a running max sees
    result = get("rating_monotonic")(build(q, stage), None, {"bins": 12, "tolerance": 0.0})
    assert result.status == FAIL
    assert result.diagnostics["deficit_m"] > 0.0


def test_rating_monotonic_names_the_store_when_the_rating_is_read_off_it():
    """The archived reason must not claim a discharge check that was not made."""
    q = np.linspace(0.1, 5.0, 200)
    run = build(q, 2.0 * q, channel=np.linspace(0.2, 9.0, 200), discharge=False)
    result = get("rating_monotonic")(run, None, {"bins": 12, "tolerance": 0.0})
    assert result.status == PASS
    assert result.diagnostics["discharge_source"] == "channel"
    assert "channel store" in result.message
    assert "with discharge" not in result.message


def test_an_explicit_tolerance_is_honoured_even_when_tighter_than_the_floor():
    """`rating_loop`'s `min_loop_m` is a true override; `tolerance` matches it."""
    # A rating that rises, peaks and then falls away: the running-maximum
    # deficit over the binned medians is a visible share of the span. The two
    # tolerances are chosen to straddle it, so the outcome turns on which one
    # the criterion actually used. Both are far below the 2%-of-span floor of
    # the upper case, which is what makes this a test of the override.
    q = np.linspace(0.1, 5.0, 400)
    stage = np.where(q < 3.0, 10.0 * q, 30.0 - 2.0 * (q - 3.0))
    loose = get("rating_monotonic")(build(q, stage), None, {"bins": 12, "tolerance": 8.5})
    tight = get("rating_monotonic")(build(q, stage), None, {"bins": 12, "tolerance": 0.01})
    assert loose.status == PASS
    assert tight.status == FAIL
    # The tight run kept its declared tolerance rather than the floor.
    assert tight.threshold == pytest.approx(0.01)
    assert loose.diagnostics["deficit_m"] > tight.threshold


def test_a_tie_in_the_abscissa_does_not_mask_the_rest_of_the_rating():
    """A short-circuit on ties would report any rating as single-valued."""
    x = np.array([0.0, 0.0, 1.0, 2.0, 3.0])
    h = np.array([0.0, 0.0, 9.0, 1.0, 9.0])  # no single-valued rating fits these
    assert _rating_scatter(h, x) >= 1.0


def test_a_tie_with_equal_stage_still_measures_the_separated_pairs():
    x = np.array([0.0, 0.0, 1.0, 2.0])
    h = np.array([4.0, 4.0, 5.0, 1.0])
    assert _rating_scatter(h, x) > 0.0


def test_an_exact_function_has_scatter_at_rounding_scale():
    x = np.linspace(0.0, 10.0, 50)
    assert _rating_scatter(2.0 * x + 1.0, x) < 1e-9


def test_a_stage_that_is_a_function_of_the_discharge_has_no_loop():
    """A gauge computed from the flow is single-valued in the flow, whatever
    the store does.

    The probe pairs the limbs on the store, and a physical model that computes
    its stage from its own discharge — `flex_lumped`, `sacsma_snow17` — shows
    centimetres of apparent loop on that axis. It is pairing noise: the store
    lags the flow, so two steps at one store carry different discharges. The
    criterion has to read the rating as a function of the discharge it is
    drawn against, or it fails an honest model on a noise floor no size test
    can separate from a real loop.
    """
    t = np.linspace(0.0, 12 * np.pi, 600)
    q = 1.0 + 0.5 * np.sin(t)
    stage = 3.0 * q  # an exact function of the discharge
    # A store that lags the flow, as a routing store does.
    store = q * (1.0 + 0.4 * np.sin(0.25 * t + 1.0))
    result = get("rating_loop")(build(q, stage, channel=store), None,
                                {"bins": 12, "tolerance": 0.0, "abscissa": "store"})
    assert result.status == PASS
    assert result.diagnostics["single_valued"] is True
    assert result.diagnostics["single_valued_in"] == "dis"
    assert "single-valued" in result.message


def test_a_small_but_consistent_inversion_still_fails():
    """A loop that is small in absolute terms but inverted on every bin is the
    physics, not noise.

    The size floor exists to excuse a separation the reach cannot resolve, and
    the previous escape also excused a loop that was *consistently* on the
    wrong side but under the floor — which is how `reference_rating_inverted`
    came to trip nothing. With the floor set below a real instrument's
    resolution, a consistent inversion has to fail whatever its width.
    """
    t = np.linspace(0.0, 12 * np.pi, 600)
    q = 1.0 + 0.5 * np.sin(t)
    dq = np.gradient(q)
    # Inverted and small: a fortieth of the rating's span, not half of it.
    stage = 2.0 * q + 0.02 * np.where(dq > 0, +1.0, -1.0)
    store = q * (1.0 + 0.3 * np.sin(0.25 * t))
    result = get("rating_loop")(build(q, stage, channel=store), None,
                                {"bins": 12, "tolerance": 0.0, "abscissa": "store",
                                 "min_loop_m": 1e-4})
    assert result.status == FAIL
    assert "wrong way" in result.message
    assert result.diagnostics["loop_m"] > 0.0

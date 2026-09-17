"""Unit tests for exchange_components: the sign arm's two halves and the
residual arm must each catch something on their own, and its `abs_tol` is
pinned to the value probe.yaml actually declares.

reference_exchange_sign_error trips both the sign arm and the residual arm
on the same steps (a component whose magnitude is right but whose sign is
wrong is also a component that no longer sums to the declared net), so
deleting either arm from the criterion would leave the gate green. Swapping
gw_to_sw and sw_to_gw's values also trips both halves of the sign arm at
once (gw_to_sw turns positive and sw_to_gw turns negative together), so
that alone would let either half be deleted without a test noticing. These
tests isolate each half of the sign arm, and the residual arm, with a
synthetic table the sign-error reference does not produce.

`abs_tol` used to be 1e-6 and no test pinned it: reverting it to 1e-6 left
every test in this file green. A model with two reaches active on the same
step -- one gaining from the aquifer, one losing to it, both nonzero every
day -- has `gw_to_sw` and `sw_to_gw` each rounded to ordinary precision
(here, 6 decimal places) before they are summed into `gw_sw_exchange`.
Rounding each component separately from rounding their sum can disagree by
up to 1e-6 (half a unit in the last decimal place, from each of the two
independently-rounded components), which an honest model has no way to
avoid and which `rel_tol` alone (a relative test) does not cover when the
components themselves are small. `abs_tol: 1.0e-4` clears that; `1.0e-6`
does not. The tests below build that two-reach case at the probe's own
gate seeds and check both values, so raising `abs_tol` back to 1e-6 (or
lowering it from 1e-4) is caught here rather than only in a maintainer's
manual review.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.criteria.groundwater import exchange_components
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult
from hydroturing.seeds import gate_seeds


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/gw-sw-exchange-consistency")


@pytest.fixture(scope="module")
def params(probe):
    return dict(next(c.params for c in probe.criteria if c.name == "exchange_components"))


def _run(probe, table: pd.DataFrame) -> RunResult:
    case = build_case(probe, gate_seed(probe))
    return RunResult(case=case, table=table, meta={}, wall_seconds=0.0)


def gate_seed(probe) -> int:
    from hydroturing.seeds import gate_seeds

    return gate_seeds(probe.id, 1)[0]


def _base_table(probe) -> pd.DataFrame:
    case = build_case(probe, gate_seed(probe))
    n = len(case.forcing)
    # A genuinely bidirectional, correctly-signed exchange: enough steps
    # each way to clear exchange_directions' minimum, and a gw state that
    # need not itself balance for this criterion, which only checks the
    # components against the declared net.
    half = n // 2
    gw_to_sw = np.where(np.arange(n) < half, -0.2, 0.0)
    sw_to_gw = np.where(np.arange(n) < half, 0.0, 0.3)
    return pd.DataFrame(
        {
            "time": case.forcing["time"],
            "gw_sw_exchange": gw_to_sw + sw_to_gw,
            "gw_to_sw": gw_to_sw,
            "sw_to_gw": sw_to_gw,
            "gw": np.full(n, 2000.0),
        }
    )


def test_correctly_signed_and_summing_table_passes(probe, params):
    table = _base_table(probe)
    run = _run(probe, table)
    result = exchange_components(run, probe, params)
    assert result.status == PASS, result.message


def test_sign_arm_alone_is_caught(probe, params):
    """Swap gw_to_sw and sw_to_gw's values (not their column names): the
    magnitudes still sum to the declared net exactly, since addition does
    not care which column held which value, but gw_to_sw is now positive
    and sw_to_gw negative, violating the sign convention with zero residual."""
    table = _base_table(probe)
    swapped = table.copy()
    swapped["gw_to_sw"] = table["sw_to_gw"]
    swapped["sw_to_gw"] = table["gw_to_sw"]
    run = _run(probe, swapped)
    result = exchange_components(run, probe, params)
    assert result.status == FAIL
    assert result.diagnostics["max_abs_residual_mm"] < 1e-9, (
        "the residual arm should be exactly satisfied here; only the sign "
        "arm should be responsible for the failure"
    )


def test_sign_arm_gw_to_sw_half_alone_is_caught(probe, params):
    """Only gw_to_sw turns positive; sw_to_gw keeps its legal sign. The net
    is adjusted to match exactly, so the residual arm stays exactly
    satisfied and only the gw_to_sw half of the sign arm can fail this."""
    table = _base_table(probe)
    broken = table.copy()
    broken["gw_to_sw"] = table["gw_to_sw"].abs()
    broken["gw_sw_exchange"] = broken["gw_to_sw"] + broken["sw_to_gw"]
    run = _run(probe, broken)
    result = exchange_components(run, probe, params)
    assert result.status == FAIL
    assert result.diagnostics["max_abs_residual_mm"] < 1e-9, (
        "the residual arm should be exactly satisfied here; only the "
        "gw_to_sw half of the sign arm should be responsible"
    )


def test_sign_arm_sw_to_gw_half_alone_is_caught(probe, params):
    """Only sw_to_gw turns negative; gw_to_sw keeps its legal sign. The net
    is adjusted to match exactly, so the residual arm stays exactly
    satisfied and only the sw_to_gw half of the sign arm can fail this."""
    table = _base_table(probe)
    broken = table.copy()
    broken["sw_to_gw"] = -table["sw_to_gw"].abs()
    broken["gw_sw_exchange"] = broken["gw_to_sw"] + broken["sw_to_gw"]
    run = _run(probe, broken)
    result = exchange_components(run, probe, params)
    assert result.status == FAIL
    assert result.diagnostics["max_abs_residual_mm"] < 1e-9, (
        "the residual arm should be exactly satisfied here; only the "
        "sw_to_gw half of the sign arm should be responsible"
    )


def test_residual_arm_alone_is_caught(probe, params):
    """Keep every sign correct, but report a net that omits sw_to_gw: the
    components no longer sum to the declared net, with no sign violation."""
    table = _base_table(probe)
    broken = table.copy()
    broken["gw_sw_exchange"] = table["gw_to_sw"]
    run = _run(probe, broken)
    result = exchange_components(run, probe, params)
    assert result.status == FAIL
    assert result.diagnostics["minimum_gw_to_sw"] <= 0.0
    assert result.diagnostics["maximum_sw_to_gw"] >= 0.0
    assert result.diagnostics["max_abs_residual_mm"] > 0.0


def _two_reach_table(probe, seed: int) -> pd.DataFrame:
    """A synthetic two-reach model: on every step, one reach is gaining
    from the aquifer (gw_to_sw < 0) and a second reach is losing to it
    (sw_to_gw > 0) at the same time, at the O(1e-1) to O(1e0) mm magnitudes
    the probe README describes. Both signed components are written at
    ordinary 6-decimal precision, and so is their sum (as an adapter would
    report gw_sw_exchange, a column of its own rather than a live
    recomputation of gw_to_sw + sw_to_gw). Rounding each component
    separately from rounding their exact sum disagrees by up to 1e-6 on a
    real fraction of steps -- not a contrived one-off -- which is the
    rounding gap `abs_tol` exists to absorb.

    The two frequencies (37 and 53 days) are coprime-ish and deliberately
    different from each other and from the case's own 90-day spinup or
    180-day window, so the two reaches' gaining/losing cycles are never in
    lockstep and both stay strictly nonzero (never crossing zero at the
    same time as each other) across the whole record.
    """
    case = build_case(probe, seed)
    n = len(case.forcing)
    t = np.arange(n)
    gw_to_sw = -(0.15 + 0.35 * (1 + np.sin(2 * np.pi * t / 37 + seed % 7)) / 2)
    sw_to_gw = 0.20 + 0.80 * (1 + np.cos(2 * np.pi * t / 53 + seed % 5)) / 2
    net_exact = gw_to_sw + sw_to_gw
    return pd.DataFrame(
        {
            "time": case.forcing["time"],
            "gw_sw_exchange": np.round(net_exact, 6),
            "gw_to_sw": np.round(gw_to_sw, 6),
            "sw_to_gw": np.round(sw_to_gw, 6),
        }
    )


@pytest.mark.parametrize("seed", gate_seeds("mass/gw-sw-exchange-consistency", 5))
def test_two_reach_rounding_passes_at_the_probes_declared_abs_tol(probe, params, seed):
    """At probe.yaml's actual `abs_tol` (1e-4 as of this test), the two-reach
    case's rounding gap must not fail an honest model on any of the five
    gate seeds."""
    assert params["abs_tol"] == pytest.approx(1.0e-4), (
        "this test is pinned to probe.yaml's declared abs_tol; update it "
        "if that value ever changes"
    )
    case = build_case(probe, seed)
    table = _two_reach_table(probe, seed)
    run = RunResult(case=case, table=table, meta={}, wall_seconds=0.0)
    result = exchange_components(run, probe, params)
    assert result.status == PASS, result.message


@pytest.mark.parametrize("seed", gate_seeds("mass/gw-sw-exchange-consistency", 5))
def test_two_reach_rounding_fails_if_abs_tol_reverted_to_one_e_minus_six(probe, params, seed):
    """The same two-reach case, at the old `abs_tol: 1.0e-6`, must fail on
    every one of the five gate seeds: this is the regression that motivated
    raising `abs_tol` to 1e-4, and it must stay caught if the value is ever
    lowered again."""
    case = build_case(probe, seed)
    table = _two_reach_table(probe, seed)
    run = RunResult(case=case, table=table, meta={}, wall_seconds=0.0)
    tightened = dict(params, abs_tol=1.0e-6)
    result = exchange_components(run, probe, tightened)
    assert result.status == FAIL, (
        f"seed {seed}: expected the two-reach rounding gap to fail at "
        f"abs_tol=1e-6, got {result.status} ({result.message})"
    )


def test_two_reach_case_still_catches_the_sign_error_reference(probe, params):
    """Raising `abs_tol` to cover honest rounding must not also let a real
    sign fault through: apply reference_exchange_sign_error's fault (flip
    the sign of gw_to_sw on alternating steps) onto the two-reach table and
    confirm it still fails at the probe's declared abs_tol."""
    case = build_case(probe, gate_seed(probe))
    table = _two_reach_table(probe, gate_seed(probe))
    broken = table.copy()
    n = len(broken)
    flipped = np.where(np.arange(n) % 2 == 0, -table["gw_to_sw"], table["gw_to_sw"])
    broken["gw_to_sw"] = flipped
    run = RunResult(case=case, table=broken, meta={}, wall_seconds=0.0)
    result = exchange_components(run, probe, params)
    assert result.status == FAIL, result.message


def test_two_reach_case_still_catches_a_net_that_omits_a_component(probe, params):
    """Raising `abs_tol` must not let a net that drops one of the two
    signed components through either."""
    case = build_case(probe, gate_seed(probe))
    table = _two_reach_table(probe, gate_seed(probe))
    broken = table.copy()
    broken["gw_sw_exchange"] = table["gw_to_sw"]
    run = RunResult(case=case, table=broken, meta={}, wall_seconds=0.0)
    result = exchange_components(run, probe, params)
    assert result.status == FAIL, result.message

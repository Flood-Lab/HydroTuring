"""Unit tests for exchange_components: the sign arm and the residual arm
must each catch something on their own.

reference_exchange_sign_error trips both arms on the same steps (a
component whose magnitude is right but whose sign is wrong is also a
component that no longer sums to the declared net), so deleting either arm
from the criterion would leave the gate green. These tests isolate each
arm with a synthetic table the sign-error reference does not produce:
swapped-sign columns that still sum correctly (sign arm only), and a net
column that omits one component while every sign stays correct (residual
arm only).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.base import FAIL, PASS, make_window
from hydroturing.criteria.exchange import exchange_components
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/gw-sw-exchange-consistency")


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


def test_correctly_signed_and_summing_table_passes(probe):
    table = _base_table(probe)
    run = _run(probe, table)
    result = exchange_components(run, probe, {})
    assert result.status == PASS, result.message


def test_sign_arm_alone_is_caught(probe):
    """Swap gw_to_sw and sw_to_gw's values (not their column names): the
    magnitudes still sum to the declared net exactly, since addition does
    not care which column held which value, but gw_to_sw is now positive
    and sw_to_gw negative, violating the sign convention with zero residual."""
    table = _base_table(probe)
    swapped = table.copy()
    swapped["gw_to_sw"] = table["sw_to_gw"]
    swapped["sw_to_gw"] = table["gw_to_sw"]
    run = _run(probe, swapped)
    result = exchange_components(run, probe, {})
    assert result.status == FAIL
    assert result.diagnostics["max_abs_residual_mm"] < 1e-9, (
        "the residual arm should be exactly satisfied here; only the sign "
        "arm should be responsible for the failure"
    )


def test_residual_arm_alone_is_caught(probe):
    """Keep every sign correct, but report a net that omits sw_to_gw: the
    components no longer sum to the declared net, with no sign violation."""
    table = _base_table(probe)
    broken = table.copy()
    broken["gw_sw_exchange"] = table["gw_to_sw"]
    run = _run(probe, broken)
    result = exchange_components(run, probe, {})
    assert result.status == FAIL
    assert result.diagnostics["minimum_gw_to_sw"] <= 0.0
    assert result.diagnostics["maximum_sw_to_gw"] >= 0.0
    assert result.diagnostics["max_abs_residual_mm"] > 0.0

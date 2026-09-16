"""Unit tests for exchange_components: the sign arm's two halves and the
residual arm must each catch something on their own.

reference_exchange_sign_error trips both the sign arm and the residual arm
on the same steps (a component whose magnitude is right but whose sign is
wrong is also a component that no longer sums to the declared net), so
deleting either arm from the criterion would leave the gate green. Swapping
gw_to_sw and sw_to_gw's values also trips both halves of the sign arm at
once (gw_to_sw turns positive and sw_to_gw turns negative together), so
that alone would let either half be deleted without a test noticing. These
tests isolate each half of the sign arm, and the residual arm, with a
synthetic table the sign-error reference does not produce.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.criteria.exchange import exchange_components
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult


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

"""Unit tests for exchange_directions, using the parameters
probes/mass/gw-sw-exchange-consistency/probe.yaml declares.

The criterion asks that, over the scored record, the aquifer both loses water
to the river (`gw_to_sw` below zero) and gains water from it (`sw_to_gw` above
zero), each by at least `minimum_gross_mm`. `reference_exchange_sign_error`
does not trip it and `reference_exchange_exact` clears it by a wide margin, so
without these tests the gate would notice neither a criterion that accepted a
one-sided exchange nor one that counted a wrongly signed component.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.criteria.groundwater import exchange_directions
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult
from hydroturing.seeds import gate_seeds

PROBE_ID = "mass/gw-sw-exchange-consistency"


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe(PROBE_ID)


@pytest.fixture(scope="module")
def params(probe):
    return dict(next(c.params for c in probe.criteria if c.name == "exchange_directions"))


@pytest.fixture(scope="module")
def case(probe):
    return build_case(probe, gate_seeds(probe.id, 1)[0])


def _score(probe, params, case, gw_to_sw: np.ndarray, sw_to_gw: np.ndarray):
    table = pd.DataFrame(
        {
            "time": case.forcing["time"],
            "gw_sw_exchange": gw_to_sw + sw_to_gw,
            "gw_to_sw": gw_to_sw,
            "sw_to_gw": sw_to_gw,
            "gw": np.full(len(case.forcing), 2000.0),
        }
    )
    run = RunResult(case=case, table=table, meta={}, wall_seconds=0.0)
    return exchange_directions(run, probe, params)


def _every_other_step(n: int, value: float, offset: int) -> np.ndarray:
    return np.where(np.arange(n) % 2 == offset, value, 0.0)


def test_the_minimum_gross_volume_is_the_reviewed_value(params):
    assert params["minimum_gross_mm"] == pytest.approx(0.1)


def test_an_aquifer_that_both_loses_and_gains_passes(probe, params, case):
    n = len(case.forcing)
    result = _score(
        probe, params, case, _every_other_step(n, -0.2, 0), _every_other_step(n, 0.3, 1)
    )
    assert result.status == PASS, result.message


@pytest.mark.parametrize("missing", ["gw_to_sw", "sw_to_gw"])
def test_a_one_sided_exchange_fails(probe, params, case, missing):
    n = len(case.forcing)
    gw_to_sw = _every_other_step(n, -0.2, 0)
    sw_to_gw = _every_other_step(n, 0.3, 1)
    if missing == "gw_to_sw":
        gw_to_sw = np.zeros(n)
    else:
        sw_to_gw = np.zeros(n)
    result = _score(probe, params, case, gw_to_sw, sw_to_gw)
    assert result.status == FAIL, result.message


def test_a_wrongly_signed_gw_to_sw_does_not_count_as_the_aquifer_losing_water(
    probe, params, case
):
    """A positive `gw_to_sw` breaks the sign convention. It must not be
    counted as the aquifer losing water to the river."""
    n = len(case.forcing)
    result = _score(
        probe, params, case, _every_other_step(n, 0.2, 0), _every_other_step(n, 0.3, 1)
    )
    assert result.status == FAIL, result.message
    assert result.diagnostics["gw_to_sw_mm"] == 0.0


def test_a_wrongly_signed_sw_to_gw_does_not_count_as_the_aquifer_gaining_water(
    probe, params, case
):
    """The mirror case: a negative `sw_to_gw` must not be counted as the
    aquifer gaining water from the river."""
    n = len(case.forcing)
    result = _score(
        probe, params, case, _every_other_step(n, -0.2, 0), _every_other_step(n, -0.3, 1)
    )
    assert result.status == FAIL, result.message
    assert result.diagnostics["sw_to_gw_mm"] == 0.0


@pytest.mark.parametrize(("total", "expected"), [(0.09, FAIL), (0.11, PASS)])
def test_the_minimum_is_a_total_over_the_scored_record(probe, params, case, total, expected):
    """The aquifer-to-river volume is summed over the scored window only: a
    large loss during spinup does not count, and a scored total just under
    or over `minimum_gross_mm` decides the result."""
    n = len(case.forcing)
    start = case.spinup_steps
    gw_to_sw = np.zeros(n)
    gw_to_sw[:start] = -1.0
    gw_to_sw[start:] = -total / (n - start)
    sw_to_gw = np.full(n, 0.3)
    result = _score(probe, params, case, gw_to_sw, sw_to_gw)
    assert result.status == expected, result.message
    assert result.diagnostics["gw_to_sw_mm"] == pytest.approx(total)

"""Unit tests for groundwater_balance's storage-scaled per-step floor.

The floor has two terms: a fixed `abs_tol_mm` covering the recharge and
exchange terms' own rounding, and a term scaled to `max|gw|` covering the
storage state's own rounding at whatever datum the model chose. `gw` is an
absolute storage, not a recharge or exchange volume, so a model can raise
the second term arbitrarily by reporting `gw` at a larger offset, and
without a cap that hides a real, fixed-size leak or an unexplained
one-step storage jump behind a floor sized for a datum the probe never
asked for. These tests build an exact balance directly (bypassing a model
adapter) and perturb only `gw`, so any resulting PASS or FAIL is
attributable to the floor alone.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from hydroturing import registry
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.criteria.groundwater import groundwater_balance
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult
from hydroturing.seeds import gate_seeds


def _probe():
    return registry.find_probe("mass/gw-sw-exchange-consistency")


def _gate_seed(probe) -> int:
    return gate_seeds(probe.id, 1)[0]


def _params(probe):
    return dict(next(c.params for c in probe.criteria if c.name == "groundwater_balance"))


def _exact_table(probe, gw_datum: float = 0.0):
    """A hand-built table, the same length as `case.forcing` (spinup rows
    included, matching what a real adapter reports), that balances exactly
    at any chosen `gw` datum. `make_window` reads `state0` from the row at
    `spinup_steps - 1`, so `gw` must accumulate from row 0, not from the
    start of the scored window, or state0 would not match the window's own
    running sum and manufacture a spurious residual.
    """
    case = build_case(probe, _gate_seed(probe))
    n = len(case.forcing)
    recharge = case.forcing["gw_recharge"].to_numpy(dtype=float)
    # An arbitrary, genuinely bidirectional exchange (sign only matters for
    # exchange_directions, which this test does not exercise).
    exchange = 0.05 * np.sin(np.arange(n) * 2 * np.pi / 90)
    storage_change = recharge + exchange
    gw = gw_datum + np.cumsum(storage_change)
    table = pd.DataFrame(
        {
            "time": case.forcing["time"],
            "gw_sw_exchange": exchange,
            "gw": gw,
        }
    )
    return table, case


def _run(case, table: pd.DataFrame) -> RunResult:
    return RunResult(case=case, table=table, meta={}, wall_seconds=0.0)


def test_exact_balance_passes_at_any_datum():
    probe = _probe()
    params = _params(probe)
    for datum in (0.0, 1.0e5, 1.0e6):
        table, case = _exact_table(probe, gw_datum=datum)
        result = groundwater_balance(_run(case, table), probe, params)
        assert result.status == PASS, f"datum {datum}: {result.message}"


def test_uncapped_floor_would_hide_a_fixed_leak_but_the_cap_catches_it():
    """A steady leak out of the aquifer (storage grows less than recharge
    plus exchange says it should) must fail regardless of the `gw` datum.
    Before the cap, shifting the datum up scaled the floor past the leak's
    size and hid it; the cap keeps the floor bounded so the leak still
    trips groundwater_balance."""
    probe = _probe()
    params = _params(probe)
    table, case = _exact_table(probe, gw_datum=1.0e6)
    leaky = table.copy()
    n = len(leaky)
    daily_leak = 0.0199
    leaky["gw"] = leaky["gw"].to_numpy(dtype=float) - np.arange(n) * daily_leak
    result = groundwater_balance(_run(case, leaky), probe, params)
    assert result.status == FAIL, result.message


def test_uncapped_floor_would_hide_a_one_step_jump_but_the_cap_catches_it():
    """An unexplained one-step storage jump must fail regardless of the
    `gw` datum, for the same reason as the steady leak above."""
    probe = _probe()
    params = _params(probe)
    table, case = _exact_table(probe, gw_datum=1.0e6)
    jumped = table.copy()
    values = jumped["gw"].to_numpy(dtype=float).copy()
    mid = len(values) // 2
    values[mid:] -= 5.0
    jumped["gw"] = values
    result = groundwater_balance(_run(case, jumped), probe, params)
    assert result.status == FAIL, result.message


def test_honest_output_survives_rounding_at_eight_significant_figures():
    """The probe's README asks for `gw` at full precision or at least 8
    significant figures; confirm that precision clears the capped floor
    even under a large datum shift."""
    probe = _probe()
    params = _params(probe)
    table, case = _exact_table(probe, gw_datum=1.0e5)
    rounded = table.copy()
    rounded["gw"] = [float(f"{v:.8g}") for v in rounded["gw"]]
    result = groundwater_balance(_run(case, rounded), probe, params)
    assert result.status == PASS, result.message


def test_six_significant_figures_can_fail_under_a_large_datum_shift():
    """Below the README's precision floor, rounding noise at a large datum
    can itself exceed the capped floor: this is expected, and is why the
    probe asks for higher precision rather than leaving the floor
    unbounded."""
    probe = _probe()
    params = _params(probe)
    table, case = _exact_table(probe, gw_datum=1.0e5)
    rounded = table.copy()
    rounded["gw"] = [float(f"{v:.6g}") for v in rounded["gw"]]
    result = groundwater_balance(_run(case, rounded), probe, params)
    assert result.status == FAIL, result.message

"""Unit tests for groundwater_balance, using the probe's own declared
parameters (probes/mass/gw-sw-exchange-consistency/probe.yaml) throughout
rather than parameters invented in the test.

groundwater_balance runs two independent checks over the same residual
(recharge + exchange + declared sources - storage change):

  * a per-step check: |residual[i]| <= max(rel_tol * denom[i], step_floor)
  * a cumulative check: |sum(residual)| <= max(rel_tol * total_recharge,
    2 * step_floor)

Before this file, no test exercised either check against the probe's real
rel_tol/abs_tol_mm, and no test proved that either check is actually load-
bearing: deleting the cumulative check, or silently ignoring `sources`,
left pytest and `ht gate` green. That is the gap this file closes. Each
test below is built so that, if the mechanism it names were removed or
disabled, that specific test (not some other one, by coincidence) would
flip from FAIL to PASS or PASS to FAIL.

The five cases asked for:
  1. exact output passes;
  2. a leak sized at 10% of daily recharge fails (naturally on both checks,
     since a 10% relative leak also exceeds the 1% per-step rel_tol; a
     second, smaller leak is added that fails ONLY the cumulative check, so
     that check is pinned in isolation too);
  3. an alternating +/-0.5 mm per-step error fails only the per-step check
     (it cancels to zero cumulatively, so the cumulative check alone would
     wave it through);
  4. a large `gw` datum shift hides none of the faults above;
  5. `sources` scoping (the `gw_boundary` fix) behaves as documented: a
     boundary term is credited only when declared under the name the probe
     actually asks for, not under `gwex`.

Also covered:
  * `value` is the worse of the two residuals as a share of its own
    allowance, against a threshold of 1, so a run fails exactly when its
    value exceeds 1;
  * the probe's own rel_tol and abs_tol_mm are the reviewed values;
  * the storage-scaled floor lets honest output through at the precision the
    probe README asks for, scales with the reported store, and is capped at
    0.05 mm;
  * a non-zero exchange enters the budget with its sign.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.criteria.groundwater import groundwater_balance
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult
from hydroturing.seeds import gate_seeds

PROBE_ID = "mass/gw-sw-exchange-consistency"


def _probe():
    return registry.find_probe(PROBE_ID)


def _gate_seed(probe) -> int:
    # The gate's own seed, not an arbitrary one, so the recharge series
    # these tests reason about is the same series `ht gate` itself uses.
    return gate_seeds(probe.id, 1)[0]


def _params(probe):
    # Read rel_tol / abs_tol_mm / sources straight from probe.yaml. A test
    # that hard-codes its own tolerance numbers can drift out of step with
    # the probe and stop meaning anything; reading the real params is what
    # "pins" these tests to the probe's actual acceptance criteria.
    return dict(next(c.params for c in probe.criteria if c.name == "groundwater_balance"))


def _case_and_recharge(probe):
    """The gate seed's case and its full recharge series, spinup included.

    `make_window` scores from `case.spinup_steps` on but reads `state0`
    from the row just before it, so every array here spans the FULL
    forcing record (spinup + window), not just the scored window -- a
    table shorter than that misaligns state0 and manufactures a spurious
    residual that has nothing to do with the fault under test.
    """
    case = build_case(probe, _gate_seed(probe))
    recharge = case.forcing["gw_recharge"].to_numpy(dtype=float)
    return case, recharge


def _table(case, gw: np.ndarray, exchange: np.ndarray | None = None) -> pd.DataFrame:
    n = len(case.forcing)
    return pd.DataFrame(
        {
            "time": case.forcing["time"],
            "gw_sw_exchange": exchange if exchange is not None else np.zeros(n),
            "gw": gw,
        }
    )


def _run(case, table: pd.DataFrame) -> RunResult:
    return RunResult(case=case, table=table, meta={}, wall_seconds=0.0)


def _window_slice(case, array: np.ndarray) -> np.ndarray:
    """The scored window's own slice of a full-length array, matching what
    `make_window` hands the criterion (spinup rows dropped)."""
    return array[case.spinup_steps :]


def _tolerances(probe, params, window_recharge: np.ndarray):
    """The exact per-step and cumulative allowances `groundwater_balance`
    computes for this gate seed, from the probe's own rel_tol/abs_tol_mm.
    `gw` is held at a small, unshifted datum here (its cumsum starting at
    0), so the storage-scaled term of the per-step floor never dominates
    over abs_tol_mm within these helper calculations; the datum-shift tests
    below re-derive the floor at their own, much larger datum instead of
    reusing these numbers.
    """
    rel_tol = float(params["rel_tol"])
    abs_tol = float(params["abs_tol_mm"])
    total_recharge = float(window_recharge.sum())
    per_step_allowed = np.maximum(rel_tol * window_recharge, abs_tol)
    cumulative_allowed = max(rel_tol * total_recharge, 2.0 * abs_tol)
    return rel_tol, abs_tol, total_recharge, per_step_allowed, cumulative_allowed


def _tiny_recharge_case(probe, scale: float = 1.0e-6):
    """The gate seed's case with its recharge scaled down, so that
    `rel_tol * total_recharge` becomes far smaller than `2 * abs_tol_mm`
    and the cumulative floor term is the one that actually binds. At this
    probe's own gate-seed scale, `rel_tol * total_recharge` always
    dominates (see `_tolerances`'s docstring), so this scaled-down case is
    the only way to exercise the floor-dominated branch of
    `cumulative_allowed` at all.
    """
    case, _ = _case_and_recharge(probe)
    scaled_forcing = case.forcing.copy()
    scaled_forcing["gw_recharge"] = scaled_forcing["gw_recharge"] * scale
    scaled_case = dataclasses.replace(case, forcing=scaled_forcing)
    return scaled_case, scaled_forcing["gw_recharge"].to_numpy(dtype=float)


# ---------------------------------------------------------------------------
# 1. Exact output passes.
# ---------------------------------------------------------------------------


def test_exact_output_passes():
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    gw = np.cumsum(recharge)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == PASS, result.message


# ---------------------------------------------------------------------------
# 2. A leak sized at 10% of daily recharge fails, and a smaller leak that
#    clears the per-step floor on every single day still fails on the
#    cumulative check alone -- proving that check is load-bearing by
#    itself, not just riding along with the per-step check.
# ---------------------------------------------------------------------------


def test_a_leak_of_ten_percent_of_daily_recharge_fails():
    """A leak equal to 10% of each day's recharge (the aquifer keeps only
    90% of what recharge plus exchange say it should). At this probe's
    rel_tol = 1%, a 10% relative shortfall is ten times the per-step
    allowance on every day it occurs, so this leak fails the per-step
    check as well as the cumulative one -- both are real findings here,
    and the test asserts both diagnostics to say so rather than picking
    only the passing story.
    """
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    window_recharge = _window_slice(case, recharge)
    _, _, _, per_step_allowed, cumulative_allowed = _tolerances(probe, params, window_recharge)

    leak = 0.10 * recharge
    gw = np.cumsum(recharge - leak)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == FAIL, result.message

    window_leak = _window_slice(case, leak)
    assert np.any(window_leak > per_step_allowed), (
        "test setup error: the 10% leak was expected to also exceed the "
        "per-step floor on at least one day"
    )
    assert window_leak.sum() > cumulative_allowed, (
        "test setup error: the 10% leak was expected to exceed the "
        "cumulative allowance too"
    )


def test_a_smaller_systematic_leak_fails_the_cumulative_check_alone():
    """A same-signed leak sized to 99% of each day's own per-step
    allowance clears the per-step check on every individual day (0.99x is
    always <= 1x, by construction), yet the daily allowances themselves
    sum to more than the record's cumulative allowance -- so a leak riding
    just under the per-step floor every day still drifts the budget past
    what the cumulative check permits. This is the case that isolates the
    cumulative check: were it deleted (or its threshold loosened to match
    the per-step one), this leak would report PASS instead of FAIL, while
    every other test in this file would be unaffected.
    """
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    window_recharge = _window_slice(case, recharge)
    _, _, _, per_step_allowed, cumulative_allowed = _tolerances(probe, params, window_recharge)

    # Only the scored window's per-step allowance is well-defined here (the
    # spinup rows are never checked by make_window), so the leak is applied
    # window-relative and zero during spinup.
    leak_window = 0.99 * per_step_allowed
    assert np.all(leak_window <= per_step_allowed)  # every day clears its own floor
    assert leak_window.sum() > cumulative_allowed  # but the record does not

    n = len(case.forcing)
    leak_full = np.zeros(n)
    leak_full[case.spinup_steps :] = leak_window
    gw = np.cumsum(recharge - leak_full)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == FAIL, result.message
    # The message and diagnostics must name the cumulative mechanism, not
    # the per-step one, confirming which check actually caught this.
    assert "cumulative" in result.message
    assert result.diagnostics["cumulative_residual_mm"] > result.diagnostics["cumulative_allowed_mm"]
    assert result.diagnostics["max_abs_step_residual_mm"] <= per_step_allowed.max()


# ---------------------------------------------------------------------------
# 3. An alternating +/-0.5 mm per-step error fails only the per-step check.
# ---------------------------------------------------------------------------


def test_an_alternating_half_millimetre_error_fails_the_per_step_check_alone():
    """A +0.5 mm / -0.5 mm error that alternates sign every step sums to
    (near) zero over the record -- the cumulative check alone would not
    notice it -- but 0.5 mm exceeds this probe's per-step allowance
    (max(1% of daily recharge, 2e-4 mm)) on every single day, so the
    per-step check must catch it on its own. This is the mirror image of
    the cumulative-only leak above: it isolates the per-step check instead.
    """
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    window_recharge = _window_slice(case, recharge)
    _, _, _, per_step_allowed, cumulative_allowed = _tolerances(probe, params, window_recharge)

    n = len(case.forcing)
    alt_full = np.where(np.arange(n) % 2 == 0, 0.5, -0.5)
    alt_window = _window_slice(case, alt_full)
    assert np.all(np.abs(alt_window) > per_step_allowed)  # exceeds every day's own floor
    assert abs(alt_window.sum()) <= cumulative_allowed  # yet clears the cumulative one

    # storage_change = recharge - alt, so residual = recharge - storage_change = alt.
    gw = np.cumsum(recharge - alt_full)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == FAIL, result.message
    assert "per step" in result.message
    assert result.diagnostics["max_abs_step_residual_mm"] > per_step_allowed.max()
    assert abs(result.diagnostics["cumulative_residual_mm"]) <= result.diagnostics["cumulative_allowed_mm"]


def test_value_exceeds_the_threshold_exactly_when_a_check_fails():
    """`value` is the larger of the worst per-step residual / allowance and
    the cumulative residual / cumulative allowance, against a threshold of
    1. Both failure modes must show a value above 1, and a pass a value of
    at most 1: the alternating error fails the per-step check while its
    cumulative residual is ~0, and the leak riding under the per-step floor
    fails the cumulative check while every step is within its allowance.
    """
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    window_recharge = _window_slice(case, recharge)
    _, _, _, per_step_allowed, _ = _tolerances(probe, params, window_recharge)
    n = len(case.forcing)

    exact = groundwater_balance(_run(case, _table(case, np.cumsum(recharge))), probe, params)
    assert exact.status == PASS
    assert exact.threshold == 1.0
    assert exact.value <= 1.0

    alt_full = np.where(np.arange(n) % 2 == 0, 0.5, -0.5)
    alternating = groundwater_balance(
        _run(case, _table(case, np.cumsum(recharge - alt_full))), probe, params
    )
    assert alternating.status == FAIL
    assert alternating.value > 1.0
    assert alternating.diagnostics["cumulative_ratio"] <= 1.0 < alternating.diagnostics["step_ratio"]

    leak_full = np.zeros(n)
    leak_full[case.spinup_steps :] = 0.99 * per_step_allowed
    drifting = groundwater_balance(
        _run(case, _table(case, np.cumsum(recharge - leak_full))), probe, params
    )
    assert drifting.status == FAIL
    assert drifting.value > 1.0
    assert drifting.diagnostics["step_ratio"] <= 1.0 < drifting.diagnostics["cumulative_ratio"]


def test_a_pass_on_the_floor_dominated_allowance_still_reports_value_at_most_one():
    """With recharge scaled down, `2 * step_floor` is the cumulative
    allowance that binds. A systematic leak at 90% of that allowance passes,
    although it is far more than rel_tol of the (tiny) recharge. The run's
    value must still be at most 1, and the share of recharge is kept as a
    diagnostic.
    """
    probe = _probe()
    params = _params(probe)
    rel_tol = float(params["rel_tol"])
    abs_tol = float(params["abs_tol_mm"])

    case, recharge = _tiny_recharge_case(probe)
    window_recharge = _window_slice(case, recharge)
    total_recharge = float(window_recharge.sum())
    cumulative_allowed = max(rel_tol * total_recharge, 2.0 * abs_tol)
    assert cumulative_allowed == pytest.approx(2.0 * abs_tol), (
        "test setup error: expected the floor term to dominate "
        "cumulative_allowed at this scaled-down recharge"
    )

    n = len(case.forcing)
    # A systematic per-step leak, small enough to clear the per-step floor
    # on every day, but summing to 90% of cumulative_allowed -- comfortably
    # within it, so the run PASSES -- while sitting far above rel_tol as a
    # share of this tiny recharge total.
    target_cumulative_leak = 0.9 * cumulative_allowed
    leak_per_day = target_cumulative_leak / len(window_recharge)
    assert leak_per_day <= abs_tol  # clears the per-step floor everywhere
    expected_value = target_cumulative_leak / total_recharge
    assert expected_value > rel_tol, (
        "test setup error: expected this leak's value to exceed rel_tol"
    )

    leak_full = np.zeros(n)
    leak_full[case.spinup_steps :] = leak_per_day
    gw = np.cumsum(recharge - leak_full)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)

    assert result.status == PASS, (
        f"test setup error: expected this leak to pass on the floor-"
        f"dominated cumulative allowance; got {result.message}"
    )
    assert result.diagnostics["cumulative_residual_share_of_recharge"] == pytest.approx(expected_value)
    assert result.diagnostics["cumulative_ratio"] == pytest.approx(0.9)
    assert result.threshold == 1.0
    assert result.value <= 1.0


# ---------------------------------------------------------------------------
# 4. A shifted `gw` datum hides none of the faults above.
# ---------------------------------------------------------------------------


LARGE_DATUM = 1.0e6


def test_exact_output_still_passes_under_a_large_datum_shift():
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    gw = LARGE_DATUM + np.cumsum(recharge)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == PASS, result.message


def test_the_cumulative_only_leak_still_fails_under_a_large_datum_shift():
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    window_recharge = _window_slice(case, recharge)
    _, _, _, per_step_allowed, _ = _tolerances(probe, params, window_recharge)

    n = len(case.forcing)
    leak_full = np.zeros(n)
    leak_full[case.spinup_steps :] = 0.99 * per_step_allowed
    gw = LARGE_DATUM + np.cumsum(recharge - leak_full)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == FAIL, result.message
    assert "cumulative" in result.message


def test_the_alternating_error_still_fails_under_a_large_datum_shift():
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    n = len(case.forcing)
    alt_full = np.where(np.arange(n) % 2 == 0, 0.5, -0.5)
    gw = LARGE_DATUM + np.cumsum(recharge - alt_full)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == FAIL, result.message
    assert "per step" in result.message


def test_a_fixed_leak_still_fails_under_a_large_datum_shift():
    """A plain, constant-rate leak (not scaled to recharge or to the
    per-step floor) is the simplest possible fault, and the one the
    storage-floor cap fix was written for: before that cap, shifting `gw`
    to a large enough datum scaled the floor past this leak's size and hid
    it. It must still fail here."""
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    n = len(case.forcing)
    daily_leak = 0.02  # mm/day, comparable to this seed's mean daily recharge
    gw = LARGE_DATUM + np.cumsum(recharge) - np.arange(n) * daily_leak
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == FAIL, result.message


# ---------------------------------------------------------------------------
# 5. `sources` scoping: a boundary term is credited only when declared
#    under the flux this probe actually names (`gw_boundary`), not under
#    the whole-catchment `gwex`, matching probe.yaml and the README.
# ---------------------------------------------------------------------------


def test_an_aquifer_outflow_declared_as_gw_boundary_passes():
    """A model with a real GHB/WEL-style term acting on the aquifer must
    declare it under the flux named in probe.yaml's `sources` (currently
    `gw_boundary`) to have it credited; declared, the budget closes
    exactly."""
    probe = _probe()
    params = _params(probe)
    assert params["sources"] == ["gw_boundary"], (
        "this test is pinned to the probe's declared sources; update it "
        "if probe.yaml's groundwater_balance.sources ever changes"
    )
    case, recharge = _case_and_recharge(probe)
    n = len(case.forcing)
    daily_outflow = 0.05
    gw = np.cumsum(recharge) - np.arange(n) * daily_outflow
    table = _table(case, gw)
    table["gw_boundary"] = -daily_outflow
    result = groundwater_balance(_run(case, table), probe, params)
    assert result.status == PASS, result.message


def test_the_same_outflow_left_undeclared_fails():
    """The same aquifer outflow with no `gw_boundary` column is an
    unexplained residual, as it should be."""
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    n = len(case.forcing)
    daily_outflow = 0.05
    gw = np.cumsum(recharge) - np.arange(n) * daily_outflow
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == FAIL, result.message


def test_a_withdrawal_taken_from_the_channel_and_reported_as_gwex_still_passes():
    """A model that takes a withdrawal from the channel, not the aquifer,
    and correctly reports it as `gwex` (a whole-catchment boundary term)
    rather than `gw_boundary`, must still pass groundwater_balance: `gw`
    genuinely did not change, and `gwex` is not one of this criterion's
    sources, so it is neither credited nor needed. Before `sources`
    defaulted to `[]` and probe.yaml named `gw_boundary` instead of
    `gwex`, this scenario mis-scored: a model correct about which control
    volume its withdrawal left was penalized for a residual that was
    never there.
    """
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    gw = np.cumsum(recharge)  # gw itself is untouched by the channel withdrawal
    table = _table(case, gw)
    table["gwex"] = -0.05
    result = groundwater_balance(_run(case, table), probe, params)
    assert result.status == PASS, result.message


def test_declaring_the_outflow_as_gwex_instead_of_gw_boundary_still_fails():
    """The mirror image of the previous two tests: an outflow that DID
    leave the aquifer, but is declared under `gwex` rather than the
    `gw_boundary` name probe.yaml's `sources` actually names, must not be
    credited -- `gwex` is deliberately not one of this criterion's
    sources, so a model that mislabels an aquifer-boundary term as `gwex`
    is scored on an incomplete budget and fails, exactly as it would if it
    had not declared the term at all."""
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    n = len(case.forcing)
    daily_outflow = 0.05
    gw = np.cumsum(recharge) - np.arange(n) * daily_outflow
    table = _table(case, gw)
    table["gwex"] = -daily_outflow
    result = groundwater_balance(_run(case, table), probe, params)
    assert result.status == FAIL, result.message


# ---------------------------------------------------------------------------
# 6. The reviewed tolerances, the storage-scaled floor and its cap, and the
#    exchange term's sign.
# ---------------------------------------------------------------------------


def test_the_probes_tolerances_are_the_reviewed_values():
    """Every test above reads rel_tol and abs_tol_mm from probe.yaml, so
    loosening them there would loosen those tests with them. Pin the values
    that were reviewed."""
    params = _params(_probe())
    assert params["rel_tol"] == pytest.approx(0.01)
    assert params["abs_tol_mm"] == pytest.approx(2.0e-4)


def _rounded(values: np.ndarray, fmt: str) -> np.ndarray:
    return np.array([float(fmt % v) for v in values])


@pytest.mark.parametrize(
    ("datum", "fmt"),
    [(2000.0, "%.6g"), (1.0e5, "%.8g"), (5.0e5, "%.8g")],
)
def test_honest_output_passes_at_the_precision_the_readme_allows(datum, fmt):
    """6 significant figures below 10,000 mm and 8 below 1,000,000 mm: the
    storage-scaled term of the floor absorbs that rounding. Without it (a
    floor of abs_tol_mm alone), or with the cap set far below 0.05 mm, these
    fail."""
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    gw = _rounded(datum + np.cumsum(recharge), fmt)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == PASS, result.message


@pytest.mark.parametrize(("datum", "fmt"), [(1.0e4, "%.6g"), (1.0e6, "%.8g")])
def test_output_coarser_than_the_capped_floor_fails(datum, fmt):
    """Once `gw` is written in steps coarser than the 0.05 mm cap, honest
    output fails on rounding alone, which is why the README asks for full
    precision."""
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    gw = _rounded(datum + np.cumsum(recharge), fmt)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == FAIL, result.message


def _alternating_storage_error(n: int, amplitude: float) -> np.ndarray:
    """A storage error of `amplitude` on every other row. The storage change,
    and so the residual, alternates between +amplitude and -amplitude, and
    the record's residual sums to at most one amplitude."""
    return amplitude * (np.arange(n) % 2)


@pytest.mark.parametrize(("amplitude", "expected"), [(0.045, PASS), (0.055, FAIL)])
def test_the_storage_floor_is_capped_at_five_hundredths_of_a_millimetre(amplitude, expected):
    """At a 1e6 mm datum the uncapped term would be 10 mm. With the cap, a
    0.045 mm per-step error passes and a 0.055 mm one fails."""
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    n = len(case.forcing)
    gw = LARGE_DATUM + np.cumsum(recharge) + _alternating_storage_error(n, amplitude)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.diagnostics["step_floor_mm"] == pytest.approx(0.05)
    assert result.status == expected, result.message


@pytest.mark.parametrize(("datum", "expected"), [(2000.0, PASS), (0.0, FAIL)])
def test_the_storage_floor_scales_with_the_reported_store(datum, expected):
    """A 0.015 mm per-step error sits inside the 0.02 mm floor of a store
    near 2,000 mm, but not inside the 2e-4 mm floor of a store that never
    exceeds about 16 mm."""
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    n = len(case.forcing)
    gw = datum + np.cumsum(recharge) + _alternating_storage_error(n, 0.015)
    result = groundwater_balance(_run(case, _table(case, gw)), probe, params)
    assert result.status == expected, result.message


def test_the_exchange_term_enters_the_budget_with_its_sign():
    """Every other test here reports no exchange. A bidirectional exchange
    that the store follows closes the budget; the same table with the
    exchange's sign flipped does not."""
    probe = _probe()
    params = _params(probe)
    case, recharge = _case_and_recharge(probe)
    n = len(case.forcing)
    exchange = 0.2 * np.sin(2.0 * np.pi * np.arange(n) / 90.0)
    gw = np.cumsum(recharge + exchange)
    honest = groundwater_balance(_run(case, _table(case, gw, exchange)), probe, params)
    assert honest.status == PASS, honest.message
    flipped = groundwater_balance(_run(case, _table(case, gw, -exchange)), probe, params)
    assert flipped.status == FAIL, flipped.message

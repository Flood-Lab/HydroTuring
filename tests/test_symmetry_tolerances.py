"""The invariance denominator must remain meaningful when a variable is dry.

The explicit floor gives a probe an absolute tolerance near zero without
weakening its relative tolerance for larger signals. Existing probes retain
their original denominator when they do not opt in.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.symmetry import invariance
from hydroturing.protocol import Case, RunResult


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/area-invariance")


def _runs(control, transformed, *, spinup_steps=0):
    runs = {}
    for variant, values in (("control", control), ("transformed", transformed)):
        table = pd.DataFrame(values)
        forcing = pd.DataFrame({"time": pd.date_range("2000-01-01", periods=len(table))})
        case = Case("test/invariance", 0, forcing, {}, spinup_steps)
        runs[variant] = RunResult(case, table, {}, 0.0)
    return runs


def test_near_zero_output_uses_declared_native_unit_floor(probe):
    control = np.array([0.0, 1e-15, -1e-15])
    runs = _runs({"mrro": control}, {"mrro": control + 5e-10})
    params = {"unchanged": ["mrro"], "rtol": 1e-9}

    original = invariance(runs, probe, params)
    explicit_default = invariance(runs, probe, {**params, "normalization_floor": 1e-12})
    floored = invariance(runs, probe, {**params, "normalization_floor": 1.0})

    assert not original.passed
    assert original.value == pytest.approx(500.0)
    assert original == explicit_default
    assert floored.passed
    assert floored.value == pytest.approx(5e-10)
    assert floored.diagnostics["normalization_floor"] == 1.0


@pytest.mark.parametrize(
    "delta, passes",
    [
        (0.0, True),
        (np.nextafter(1e-9, 0.0), True),
        (1e-9, True),
        (np.nextafter(1e-9, np.inf), False),
    ],
)
def test_floor_tolerance_includes_boundary_but_rejects_next_float(probe, delta, passes):
    runs = _runs({"mrro": [0.0, 0.0]}, {"mrro": [0.0, delta]})
    result = invariance(
        runs, probe, {"unchanged": ["mrro"], "rtol": 1e-9, "normalization_floor": 1.0}
    )
    assert result.passed is passes
    assert result.value == delta
    assert result.threshold == 1e-9


def test_signal_above_floor_uses_mean_magnitude_not_peak_or_signed_mean(probe):
    runs = _runs({"mrro": [-2.0, 6.0]}, {"mrro": [-2.0, 6.5]})
    result = invariance(
        runs, probe, {"unchanged": ["mrro"], "rtol": 0.125, "normalization_floor": 1.0}
    )
    assert result.passed
    assert result.value == 0.125  # max departure 0.5 / mean magnitude 4


def test_scaled_expectation_and_optional_variables_receive_same_floor(probe):
    runs = _runs(
        {"dis": [0.125, 0.25], "snw": [0.0, 0.0]},
        {"dis": [0.5, 1.0], "snw": [0.0, 5e-10]},
    )
    result = invariance(
        runs, probe,
        {
            "scaled": {"dis": 4.0}, "optional": ["snw", "canopy"],
            "rtol": 1e-9, "normalization_floor": 1.0,
        },
    )
    assert result.passed
    assert result.diagnostics["deviations"] == {"dis": 0.0, "snw": 5e-10}
    assert result.diagnostics["worst_variable"] == "snw"


def test_scaled_signal_above_floor_keeps_expected_magnitude_denominator(probe):
    runs = _runs({"dis": [1.0, 2.0]}, {"dis": [4.0, 8.75]})
    result = invariance(
        runs, probe, {"scaled": {"dis": 4.0}, "rtol": 0.125, "normalization_floor": 1.0}
    )
    assert result.passed
    assert result.value == 0.125  # max departure 0.75 / mean expected magnitude 6


def test_floor_and_departure_are_computed_only_after_spinup(probe):
    runs = _runs(
        {"mrro": [1e6, 0.0, 0.0]},
        {"mrro": [-1e6, 0.0, 5e-10]},
        spinup_steps=1,
    )
    result = invariance(
        runs, probe, {"unchanged": ["mrro"], "rtol": 1e-9, "normalization_floor": 1.0}
    )
    assert result.passed
    assert result.value == 5e-10


@pytest.mark.parametrize("floor", [0, -1, np.nan, np.inf, -np.inf, None, True, "bad", []])
def test_invalid_normalization_floor_is_rejected(probe, floor):
    runs = _runs({"mrro": [0.0]}, {"mrro": [0.0]})
    with pytest.raises(ValueError, match="normalization_floor must be a positive finite number"):
        invariance(runs, probe, {"unchanged": ["mrro"], "normalization_floor": floor})


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
@pytest.mark.parametrize("variant", ["control", "transformed"])
@pytest.mark.parametrize("variable", ["mrro", "snw"])
def test_nonfinite_required_or_optional_outputs_cannot_silently_pass(
    probe, value, variant, variable,
):
    runs = _runs({"mrro": [0.0], "snw": [0.0]}, {"mrro": [0.0], "snw": [0.0]})
    runs[variant].table.loc[0, variable] = value
    with pytest.raises(ValueError, match="invariance needs finite values"):
        invariance(
            runs, probe,
            {"unchanged": ["mrro"], "optional": ["snw"], "normalization_floor": 1.0},
        )


@pytest.mark.parametrize("factor", [np.nan, np.inf, -np.inf])
def test_nonfinite_scaled_factor_cannot_silently_pass(probe, factor):
    runs = _runs({"dis": [0.0]}, {"dis": [0.0]})
    with pytest.raises(ValueError, match="invariance needs finite values"):
        invariance(runs, probe, {"scaled": {"dis": factor}, "normalization_floor": 1.0})


@pytest.mark.parametrize(
    "control, transformed, factor",
    [
        ([1e308, 1e308], [5e307, 5e307], 1.0),  # Mean magnitude overflows.
        ([1e308], [1.0], 2.0),  # Expected scaled value overflows.
        ([-1e308], [1e308], 1.0),  # Pointwise difference overflows.
    ],
)
def test_finite_inputs_with_overflow_cannot_silently_pass(probe, control, transformed, factor):
    runs = _runs({"dis": control}, {"dis": transformed})
    with pytest.raises(ValueError, match="non-finite intermediate or deviation"):
        invariance(runs, probe, {"scaled": {"dis": factor}, "normalization_floor": 1.0})

"""Interception response: the paired static counterfactual."""

import pytest
import pandas as pd

from hydroturing import registry
from hydroturing.criteria import get, is_paired
from hydroturing.harness import run_probe
from hydroturing.protocol import Case, RunResult
from hydroturing.seeds import gate_seeds

PROBE = "mass/interception-response"
POSITIVE = (
    "reference_bucket",
    "flex_lumped",
    "flex_topo",
    "reference_interception_no_canopy",
)
NEGATIVE = (
    "reference_interception_blind",
    "reference_interception_wetter",
    "reference_interception_drier",
    "reference_interception_leaky",
)


def _probe():
    return registry.find_probe(PROBE)


def _seed():
    return gate_seeds(PROBE, 1)[0]


def _params():
    return dict(_probe().criteria[0].params)


def _forcing(pr_last: float = 0.0) -> pd.DataFrame:
    return pd.DataFrame({
        "time": ["2000-01-01", "2000-01-02"],
        "pr": [0.0, pr_last],
        "tas": [10.0, 10.0],
        "pet": [2.0, 2.0],
    })


def _run(cap, *, et=0.0, mrro=0.0, gwex=0.0, storage=0.0,
         forcing=None) -> RunResult:
    """One synthetic run: two rows, the second is the scored window."""
    case = Case(
        probe_id=PROBE,
        seed=1,
        forcing=_forcing() if forcing is None else forcing,
        static={"canopy_capacity_mm": cap},
        spinup_steps=1,
        timestep="PT1D",
    )
    table = pd.DataFrame({
        "evspsbl": [0.0, et],
        "mrro": [0.0, mrro],
        "gwex": [0.0, gwex],
        "mrso": [0.0, storage],
        "canopy": [0.0, 0.0],
    })
    return RunResult(case=case, table=table, meta={}, wall_seconds=0.0)


def _score(low, high):
    return get("interception_response")(
        {"capacity_low": low, "capacity_high": high}, _probe(), _params())


def test_the_criterion_is_paired():
    assert is_paired("interception_response")


@pytest.mark.parametrize("model", POSITIVE)
def test_must_pass_models_respond_and_close(model):
    outcome = run_probe(registry.find_model(model), _probe(), [_seed()])
    assert outcome.verdict == "PASS", outcome.failing


@pytest.mark.parametrize("model", NEGATIVE)
def test_each_negative_control_fails_the_paired_criterion(model):
    outcome = run_probe(registry.find_model(model), _probe(), [_seed()])
    assert outcome.verdict == "FAIL", outcome.failing
    assert "interception_response" in outcome.failing


def test_a_declared_gwex_is_a_source_in_the_paired_budget():
    """gwex is positive into the catchment, so the residual subtracts it."""
    low = _run(2.0)
    high = _run(10.0, et=5.0, mrro=-2.0, gwex=3.0)
    result = _score(low, high)

    assert result.passed, result.message
    assert result.diagnostics["residual_mm"] == pytest.approx(0.0, abs=1e-12)


def test_runoff_rising_is_caught_by_the_direction_check():
    # 0.5 mm rise is 6.25% of the capacity change, over the 5% rule, while the
    # paired residual stays inside its 1 mm allowance (R = 1.0 mm).
    result = _score(_run(2.0), _run(10.0, et=4.0, mrro=0.5, storage=-3.5))
    assert not result.passed
    assert result.diagnostics["failing_check"] == "runoff_rise"


def test_the_residual_floor_is_on_the_paired_residual_only():
    # R = 0.6 mm is 7.5% of the capacity change, over the 5% rule, but inside
    # the 1 mm residual floor. The direction checks still use 5%.
    result = _score(_run(2.0), _run(10.0, et=4.0, mrro=-3.4))
    assert result.passed, result.message
    assert result.diagnostics["residual_mm"] == pytest.approx(0.6)
    assert result.diagnostics["residual_allowance_mm"] == pytest.approx(1.0)


def test_evaporation_falling_is_caught_by_the_direction_check():
    result = _score(_run(2.0), _run(10.0, et=-5.0, mrro=0.0, storage=5.0))
    assert not result.passed
    assert result.diagnostics["failing_check"] == "evap_fall"


def test_a_partial_allocation_is_caught_by_the_residual():
    # 4 mm returned as evaporation, 1 mm less runoff, 1 mm leaked: R = 2 mm
    # against a max(0.2, 1) = 1 mm allowance. The evaporation response itself
    # is present, so criterion 4 passes.
    result = _score(_run(2.0), _run(10.0, et=4.0, mrro=-1.0, storage=-1.0))
    assert not result.passed
    assert result.diagnostics["failing_check"] == "residual"
    assert result.diagnostics["evap_response_share"] >= 0.5


def test_an_ignored_capacity_fails_the_response_floor():
    result = _score(_run(2.0), _run(10.0))
    assert not result.passed
    assert result.diagnostics["failing_check"] == "evaporation_response"


def test_the_controls_pin_which_check_each_trips():
    """One criterion cannot separate them; the diagnostics can, so pin them."""
    diagnostics = {}
    for model in NEGATIVE:
        outcome = run_probe(registry.find_model(model), _probe(), [_seed()])
        diagnostics[model] = outcome.criteria[0].diagnostics

    assert diagnostics["reference_interception_blind"]["failing_check"] == (
        "evaporation_response")
    assert diagnostics["reference_interception_wetter"]["failing_check"] == (
        "runoff_rise")
    assert diagnostics["reference_interception_wetter"]["evap_response_share"] >= 0.5
    assert diagnostics["reference_interception_drier"]["failing_check"] == "evap_fall"
    assert diagnostics["reference_interception_drier"]["evap_response_share"] < 0.5
    assert diagnostics["reference_interception_leaky"]["failing_check"] == "residual"
    assert diagnostics["reference_interception_leaky"]["evap_response_share"] >= 0.5


def test_a_missing_capacity_is_a_configuration_error():
    low = _run(2.0)
    high = _run(10.0)
    high.case.static.pop("canopy_capacity_mm")
    with pytest.raises(ValueError, match="canopy_capacity_mm"):
        _score(low, high)


def test_divergent_weather_is_rejected():
    high_forcing = _forcing(pr_last=1.0)
    with pytest.raises(ValueError, match="differs"):
        _score(_run(2.0), _run(10.0, forcing=high_forcing))


def test_the_capacity_has_to_change():
    with pytest.raises(ValueError, match="must exceed"):
        _score(_run(10.0), _run(10.0))


def test_the_probe_declares_the_pair_and_the_optional_store():
    probe = _probe()
    assert probe.variants == ("capacity_low", "capacity_high")
    assert probe.requires_static == ("canopy_capacity_mm",)
    # canopy is optional: a model without a canopy store is scored rather
    # than rejected as INCOMPLETE.
    assert probe.requires_states == ("mrso",)
    assert [c.name for c in probe.criteria] == ["interception_response"]


def test_a_no_canopy_model_is_scored_not_incomplete():
    outcome = run_probe(
        registry.find_model("reference_interception_no_canopy"),
        _probe(), [_seed()],
    )
    assert outcome.verdict == "PASS", outcome.failing

"""Event closure on the existing hydrological models, through the /io contract."""

from dataclasses import replace

import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.harness import build_case, load_generator
from hydroturing.protocol import Case
from hydroturing.runner import get_runner
from hydroturing.seeds import gate_seeds


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/extreme-event-closure")


@pytest.mark.parametrize("model_name", [
    "reference_bucket", "flex_lumped", "flex_topo", "sacsma_snow17",
])
@pytest.mark.parametrize("overlap", [False, True], ids=["ordinary-weather", "overlapped-weather"])
def test_existing_physical_models_close_each_event(probe, model_name, overlap, tmp_path):
    case = build_case(probe, 20260912)
    if not overlap:
        forcing, static = load_generator(probe).generate_baseline(case.seed)
        case = replace(case, forcing=forcing, static=static)
    model = registry.find_model(model_name)
    run = get_runner(model).run(model, probe, case, tmp_path)
    result = get("event_water_closure")(run, probe, {})
    assert result.passed, result.message
    assert result.diagnostics["event_count"] > 0
    assert result.diagnostics["n_failed_events"] == 0


def test_existing_negative_loses_water_in_a_known_warm_wet_case(probe, tmp_path):
    # Ten 50 mm days wet the bucket without activating its >55 mm defect.
    # A separating dry day is followed by one 120 mm day. This is a unit
    # fixture for the known defect, not a return-period calibration.
    rain = [0.0] + [50.0] * 10 + [0.0, 120.0, 0.0]
    forcing = pd.DataFrame({
        "time": pd.date_range("2001-07-01", periods=len(rain), freq="D"),
        "pr": rain, "tas": 20.0, "pet": 0.0,
    })
    case = Case(
        probe_id=probe.id, seed=1, forcing=forcing,
        static=dict(load_generator(probe).STATIC), spinup_steps=12, timestep="PT1D",
    )
    results = {}
    for name in ("reference_bucket", "reference_in_sample"):
        model = registry.find_model(name)
        run = get_runner(model).run(model, probe, case, tmp_path / name)
        results[name] = get("event_water_closure")(run, probe, {})
    assert results["reference_bucket"].passed
    negative = results["reference_in_sample"]
    assert not negative.passed
    assert negative.diagnostics["event_count"] == 1
    # The saturated bucket has enough surface runoff to lose the full excess.
    assert negative.diagnostics["events"][0]["residual_mm"] == pytest.approx(35.75)
    assert negative.value == pytest.approx(35.75 / 120.0)


@pytest.mark.parametrize("seed", gate_seeds("mass/extreme-event-closure", 5))
def test_existing_negative_passes_aggregate_but_fails_event_closure(probe, seed, tmp_path):
    # The gate requires the declared failure; this regression also requires
    # aggregate closure and the other safeguards to pass on the SAME run.
    # Twenty years or a median-wet selection alone do not guarantee this.
    case = build_case(probe, seed)
    model = registry.find_model("reference_in_sample")
    run = get_runner(model).run(model, probe, case, tmp_path)
    results = {
        spec.name: get(spec.name)(run, probe, dict(spec.params))
        for spec in probe.criteria
    }
    event = results.pop("event_water_closure")
    assert not event.passed, event.message
    assert event.value > event.threshold
    for result in results.values():
        assert result.passed, result.message

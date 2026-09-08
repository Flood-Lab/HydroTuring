"""A clean control must not hide a broken budget in the shifted run."""

from dataclasses import replace

import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import FAIL, PASS, get, is_paired
from hydroturing.harness import evaluate_criteria
from hydroturing.protocol import Case, RunResult
from hydroturing.spec import Criterion


PARAMS = {
    "variants": ["control", "shifted"],
    "denominator": "sum_pr",
    "threshold": 0.05,
    "sinks": ["evspsbl", "mrro"],
}


@pytest.fixture
def probe():
    return replace(
        registry.find_probe("mass/catchment-closure"),
        variants=("control", "shifted"),
        criteria=(Criterion("paired_closure", PARAMS),),
    )


def _run(variant, *, precipitation=10.0, evaporation=3.0, runoff=7.0):
    origin = "2000-01-01" if variant == "control" else "1972-01-01"
    forcing = pd.DataFrame({
        "time": pd.date_range(origin, periods=4, freq="D").strftime("%Y-%m-%d"),
        "pr": precipitation,
    })
    table = forcing.assign(
        evspsbl=evaporation, mrro=runoff, mrso=100.0, snw=0.0, canopy=0.0,
    )
    case = Case(
        probe_id=f"mass/test-paired-closure@{variant}", seed=7,
        forcing=forcing, static={}, spinup_steps=1,
    )
    return RunResult(case=case, table=table, meta={"status": "ok"}, wall_seconds=0.0)


def test_clean_pair_closes_and_reports_each_budget(probe):
    runs = {name: _run(name) for name in probe.variants}
    assert is_paired("paired_closure")
    result, = evaluate_criteria(runs, probe)
    assert result.name == "paired_closure"
    assert result.status == PASS
    assert result.value == 0.0
    assert result.threshold == 0.05
    assert result.diagnostics["failed_variants"] == []
    assert result.diagnostics["suspicious_exact"] is True
    for entry in result.diagnostics["variants"].values():
        assert entry["status"] == PASS
        assert entry["value"] == 0.0
        assert entry["diagnostics"]["denominator_total"] == 30.0


def test_shifted_only_leak_fails_without_changing_control_scoring(probe):
    runs = {"control": _run("control"), "shifted": _run("shifted", runoff=6.0)}
    ordinary = replace(probe, criteria=(Criterion("closure", PARAMS),))
    control_result, = evaluate_criteria(runs, ordinary)
    assert control_result.status == PASS

    result, = evaluate_criteria(runs, probe)
    assert result.status == FAIL
    assert result.value == pytest.approx(0.1)
    assert result.diagnostics["failed_variants"] == ["shifted"]
    assert result.diagnostics["variants"]["control"]["status"] == PASS
    assert result.diagnostics["variants"]["shifted"]["status"] == FAIL
    assert result.diagnostics["variants"]["shifted"]["diagnostics"]["cumulative_residual"] == 3.0
    assert "shifted" in result.message


def test_zero_denominator_shifted_run_is_not_hidden_by_clean_control(probe):
    runs = {
        "control": _run("control"),
        "shifted": _run("shifted", precipitation=0.0, evaporation=0.0, runoff=0.0),
    }
    result, = evaluate_criteria(runs, probe)
    assert result.status == FAIL
    assert result.value == 0.0  # The control is the only defined relative residual.
    assert result.diagnostics["failed_variants"] == ["shifted"]
    assert result.diagnostics["variants"]["shifted"]["value"] is None
    assert "shifted" in result.message and "accumulated to zero" in result.message


def test_default_uses_all_probe_variants(probe):
    runs = {"control": _run("control"), "shifted": _run("shifted", runoff=6.0)}
    params = {key: value for key, value in PARAMS.items() if key != "variants"}
    result = get("paired_closure")(runs, probe, params)
    assert result.status == FAIL
    assert result.diagnostics["failed_variants"] == ["shifted"]


@pytest.mark.parametrize("variants, message", [
    ([], "at least two"),
    (["control"], "at least two"),
    ("control,shifted", "at least two"),
    (None, "at least two"),
    (["control", "control"], "unique"),
    (["control", "missing"], "missing requested variants"),
    (["control", 1], "nonempty strings"),
])
def test_invalid_variant_configuration_is_rejected(probe, variants, message):
    runs = {name: _run(name) for name in probe.variants}
    with pytest.raises(ValueError, match=message):
        get("paired_closure")(runs, probe, {**PARAMS, "variants": variants})


def test_default_rejects_a_probe_without_two_variants(probe):
    with pytest.raises(ValueError, match="at least two"):
        get("paired_closure")({"control": _run("control")}, replace(probe, variants=()), {})

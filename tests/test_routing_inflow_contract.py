"""Contract tests for volumetric routing inflow q_in."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from hydroturing.criteria import get
from hydroturing.protocol import Case, RunResult
from hydroturing.spec import Criterion, ProbeSpec

ROOT = Path(__file__).resolve().parents[1]


def _probe() -> ProbeSpec:
    return ProbeSpec(
        id="momentum/test-routing-inflow",
        title="test",
        law="momentum",
        track="synthetic",
        version=1,
        authors=({"name": "test", "affiliation": "test"},),
        citation="",
        requires_fluxes=("mrro",),
        requires_states=(),
        requires_forcing=("q_in",),
        requires_static=("area_km2",),
        generator="generate.py",
        n_seeds=1,
        timestep="PT1H",
        period_years=0.0,
        period_days=3.0 / 24.0,
        spinup_days=1,
        max_output_mb=1.0,
        max_runtime_s=10.0,
        variants=(),
        criteria=(Criterion("closure", {}),),
        must_pass=("reference_bucket",),
        must_fail={"reference_leaky": "closure"},
        provenance="test",
        path=ROOT / "tests",
    )


def _run(area_km2: float = 1.0, outflow_scale: float = 1.0) -> RunResult:
    time = pd.date_range("2001-01-01", periods=4, freq="h")
    forcing = pd.DataFrame({
        "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "q_in": [1.0] * 4,
    })
    # 1 m3/s over 1 km2 = 86.4 mm/day.  Scale inversely with area.
    mrro = 86.4 / area_km2 * outflow_scale
    table = pd.DataFrame({
        "time": forcing["time"],
        "mrro": [mrro] * 4,
    })
    case = Case(
        probe_id="momentum/test-routing-inflow",
        seed=1,
        forcing=forcing,
        static={"area_km2": area_km2},
        spinup_steps=1,
        timestep="PT1H",
    )
    return RunResult(case, table, {"status": "ok"}, 0.0)


def _params() -> dict:
    return {
        "denominator": "sum_inflow",
        "sinks": ["mrro"],
        "sources": [],
        "states": [],
        "threshold": 1.0e-12,
    }


@pytest.mark.parametrize("area_km2", [1.0, 2.0, 100.0])
def test_sum_inflow_converts_q_in_m3s_to_catchment_depth(area_km2: float):
    result = get("closure")(_run(area_km2=area_km2), _probe(), _params())
    assert result.passed, result.message
    assert result.value < 1.0e-12


def test_sum_inflow_still_detects_a_routing_loss():
    result = get("closure")(_run(outflow_scale=0.9), _probe(), _params())
    assert not result.passed
    assert result.value == pytest.approx(0.1)


def test_sum_inflow_requires_positive_area():
    run = _run()
    run.case.static.pop("area_km2")
    with pytest.raises(ValueError, match="sum_inflow.*area_km2"):
        get("closure")(run, _probe(), _params())

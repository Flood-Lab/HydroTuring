from __future__ import annotations
import numpy as np
import pandas as pd
from hydroturing import registry
from hydroturing.criteria.datum import datum_flux_invariance
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult


def _runs(q_control, q_up, q_down):
    probe = registry.find_probe("mass/groundwater-datum-invariance")
    case = build_case(probe, 7, "control")
    n = len(case.forcing)
    def run(values):
        q=np.asarray(values, dtype=float)
        if len(q)!=n: q=np.resize(q,n)
        table=pd.DataFrame({"gw_sw_exchange":q,"gw_to_sw":np.minimum(q,0),"sw_to_gw":np.maximum(q,0),"gw":np.full(n,2000.)})
        return RunResult(case=case,table=table,meta={},wall_seconds=0.)
    return {"control":run(q_control),"datum_up":run(q_up),"datum_down":run(q_down)}, probe

def test_datum_invariance_passes_identical_active_trajectory():
    n=len(build_case(registry.find_probe("mass/groundwater-datum-invariance"),7,"control").forcing)
    q=np.zeros(n); q[2:6]=1.; q[10:14]=-1.
    runs, probe = _runs(q,q,q)
    params=dict(next(c.params for c in probe.criteria if c.name == "datum_flux_invariance"))
    result=datum_flux_invariance(runs, probe, params)
    assert result.passed

def test_datum_invariance_rejects_zero_activity():
    n=len(build_case(registry.find_probe("mass/groundwater-datum-invariance"),7,"control").forcing)
    q=np.zeros(n)
    runs, probe = _runs(q,q,q)
    params=dict(next(c.params for c in probe.criteria if c.name == "datum_flux_invariance"))
    result=datum_flux_invariance(runs, probe, params)
    assert not result.passed
    assert "activity" in result.message

def test_datum_invariance_uses_absolute_difference_before_cancellation():
    n=len(build_case(registry.find_probe("mass/groundwater-datum-invariance"),7,"control").forcing)
    q=np.zeros(n); q[2:6]=1.; q[10:14]=-1.
    shifted=q.copy(); shifted[2:6]+=0.01; shifted[10:14]-=0.01
    runs, probe = _runs(q,shifted,shifted)
    params=dict(next(c.params for c in probe.criteria if c.name == "datum_flux_invariance"))
    result=datum_flux_invariance(runs, probe, params)
    assert not result.passed

def test_datum_invariance_rejects_nonfinite_transformed_output():
    n=len(build_case(registry.find_probe("mass/groundwater-datum-invariance"),7,"control").forcing)
    q=np.ones(n); shifted=q.copy(); shifted[3]=np.nan
    runs, probe = _runs(q,shifted,shifted)
    params=dict(next(c.params for c in probe.criteria if c.name == "datum_flux_invariance"))
    result=datum_flux_invariance(runs, probe, params)
    assert not result.passed

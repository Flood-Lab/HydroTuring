from __future__ import annotations
import numpy as np, pandas as pd, pytest
from dataclasses import replace
from hydroturing import registry
from hydroturing.criteria.recharge_ordering import aquifer_recharge_ordering
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult

def make(qc, qp, sc, sp, rc=None, rp=None):
 probe=registry.find_probe('mass/aquifer-recharge-ordering'); case=build_case(probe,7,'control'); control_case=case; n=len(case.forcing)
 def run(q,s,r,run_case):
  forcing=run_case.forcing.copy(); forcing['gw_recharge']=np.asarray(r if r is not None else np.zeros(n)); table=pd.DataFrame({'gw_sw_exchange':q,'gw':s,'gw_to_sw':np.minimum(q,0),'sw_to_gw':np.maximum(q,0)}); return RunResult(case=replace(run_case, forcing=forcing),table=table,meta={},wall_seconds=0.)
 control=run(np.asarray(qc),np.asarray(sc),rc if rc is not None else np.zeros(n),control_case)
 # Rebuild the case so the two forcing frames do not alias.
 case=build_case(probe,7,'recharge_added')
 perturbed=run(np.asarray(qp),np.asarray(sp),rp if rp is not None else np.zeros(n),case)
 return {'control':control,'recharge_added':perturbed},probe

def params(probe): return dict(next(c.params for c in probe.criteria if c.name=='aquifer_recharge_ordering'))
def test_recharge_ordering_rejects_exchange_reversal():
 n=len(build_case(registry.find_probe('mass/aquifer-recharge-ordering'),7,'control').forcing); r=np.zeros(n); r[2:4]=5; q=np.zeros(n); qp=q.copy(); qp[2]=-1; s=np.full(n,2000.); sp=s.copy(); sp[2:]=sp[2:]+np.cumsum(r)[2:]
 runs,p=make(q,qp,s,sp,np.zeros(n),r); assert not aquifer_recharge_ordering(runs,p,params(p)).passed
def test_recharge_ordering_rejects_negative_input_increment():
 n=len(build_case(registry.find_probe('mass/aquifer-recharge-ordering'),7,'control').forcing); r=np.zeros(n); r[2]=5; r[3]=-1; q=np.zeros(n); s=np.full(n,2000.); runs,p=make(q,q,s,s,r,np.zeros(n));
 with pytest.raises(ValueError,match='nonnegative'): aquifer_recharge_ordering(runs,p,params(p))


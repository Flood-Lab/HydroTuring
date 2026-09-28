from __future__ import annotations
import numpy as np
from hydroturing.criteria.base import FAIL,PASS,CriterionResult,criterion,make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec
@criterion("aquifer_recharge_ordering",paired=True)
def aquifer_recharge_ordering(runs:dict[str,RunResult],probe:ProbeSpec,params:dict)->CriterionResult:
 c=str(params.get("control","control")); p=str(params.get("perturbed","recharge_added"))
 if c not in runs or p not in runs: raise ValueError("aquifer_recharge_ordering needs control and perturbed variants")
 wc=make_window(runs[c],probe); wp=make_window(runs[p],probe)
 recharge=str(params.get("recharge","gw_recharge")); exchange=str(params.get("exchange","gw_sw_exchange")); storage=str(params.get("storage","gw"))
 for w,n in ((wc,"control"),(wp,"perturbed")):
  if recharge not in w.forcing.columns: raise ValueError(f"{n} forcing needs '{recharge}'")
  for v in (exchange,storage):
   if v not in w.table.columns: raise ValueError(f"{n} result needs '{v}'")
 if len(wc.table)!=len(wp.table): raise ValueError("paired recharge runs have different scored lengths")
 added=wp.volume(wp.forcing[recharge])-wc.volume(wc.forcing[recharge]); I=np.cumsum(added); Istar=float(I[-1])
 if Istar<=0: raise ValueError("aquifer_recharge_ordering needs positive added recharge")
 ds=wp.table[storage].to_numpy(float)-wc.table[storage].to_numpy(float)
 ds-=float(wp.state0[storage])-float(wc.state0[storage])
 d=wc.volume(wc.table[exchange])-wp.volume(wp.table[exchange]); D=np.cumsum(d)
 es=float(np.max(np.maximum(-ds,0.0))); eq=float(np.sum(np.maximum(-d,0.0))); eb=float(np.max(np.abs(ds+D-I)))
 eps=max(float(params.get("rel_tol",0.01))*Istar,float(params.get("abs_tol_mm",1e-3)))
 ratios={"storage":es/eps,"exchange":eq/eps,"budget":eb/eps}; worst_name,worst=max(ratios.items(),key=lambda x:x[1]); ok=worst<=1.0
 return CriterionResult(name="aquifer_recharge_ordering",status=PASS if ok else FAIL,value=worst,threshold=1.0,message=(f"added recharge preserves storage and exchange ordering (I={Istar:.6g} mm)" if ok else f"recharge ordering fails {worst_name}: {max(es,eq,eb):.6g} mm exceeds {eps:.6g} mm"),diagnostics={"added_recharge_mm":Istar,"epsilon_mm":eps,"E_S_mm":es,"E_Q_mm":eq,"E_B_mm":eb,"ratios":ratios})

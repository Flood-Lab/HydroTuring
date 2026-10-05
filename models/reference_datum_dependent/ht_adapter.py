from __future__ import annotations
import argparse, csv, json
import numpy as np
from pathlib import Path
COLUMNS = ["time", "gw_sw_exchange", "gw_to_sw", "sw_to_gw", "gw"]
def simulate(forcing, static, mode="dependent"):
    sy=float(static["aquifer_specific_yield"]); area=float(static["area_km2"])*1e6
    b=sy*1000.0; conductance=float(static["river_conductance_m2_per_day"])/area*1000.0 * np.exp(float(static["aquifer_initial_head_m"])/500.0)
    head=float(static["aquifer_initial_head_m"]); bottom=float(static["aquifer_bottom_m"]); gw=b*(head-bottom)
    rows=[]
    for row in forcing:
        stage=float(row["sw_stage_m"]); recharge=float(row["gw_recharge"])
        # Exact forward linear reservoir over a one-day interval.
        decay=np.exp(-conductance/b)
        equilibrium=stage + recharge/conductance if conductance else head
        new_head=equilibrium+(head-equilibrium)*decay
        q=(new_head-head-recharge/b)*b
        gw=b*(new_head-bottom); head=new_head
        rows.append({"time":row["time"],"gw_sw_exchange":q,"gw_to_sw":min(q,0.0),"sw_to_gw":max(q,0.0),"gw":gw})
    return rows
def main():
    p=argparse.ArgumentParser(); p.add_argument("--request",required=True); a=p.parse_args()
    path=Path(a.request).resolve(); req=json.loads(path.read_text()); d=path.parent
    with open(d/req["input"]["forcing"],newline="") as f: forcing=list(csv.DictReader(f))
    static=json.loads((d/req["input"]["static"]).read_text()); rows=simulate(forcing,static)
    with open(d/req["output"]["table"],"w",newline="") as f: w=csv.DictWriter(f,fieldnames=COLUMNS); w.writeheader(); w.writerows(rows)
    (d/req["output"]["run"]).write_text(json.dumps({"status":"ok"}))
if __name__=="__main__": main()

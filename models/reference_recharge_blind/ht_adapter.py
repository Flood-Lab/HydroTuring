from __future__ import annotations
import argparse,csv,json
from pathlib import Path
MODE = "blind"
COLUMNS=["time","gw_sw_exchange","gw_to_sw","sw_to_gw","gw"]
def simulate(forcing,static,mode="exact"):
 sy=float(static["aquifer_specific_yield"]); area=float(static["area_km2"])*1e6; b=sy*1000.; k=float(static["river_conductance_m2_per_day"])/area*1000.; head=float(static["aquifer_initial_head_m"]); gw=b*head; rows=[]
 for row in forcing:
  stage=float(row["sw_stage_m"]); r=float(row["gw_recharge"]); q=k*(stage-head); new=head+(r+q)/b; q=b*(new-head)-r
  if mode=="blind": r_effect=0.; new=head+(r_effect+q)/b; q=b*(new-head)-r_effect
  elif mode=="overshoot" and r>0: q=-1.2*r; new=head+(r+q)/b
  head=new; gw=b*head
  rows.append({"time":row["time"],"gw_sw_exchange":q,"gw_to_sw":min(q,0.),"sw_to_gw":max(q,0.),"gw":gw})
 return rows
def main():
 p=argparse.ArgumentParser(); p.add_argument("--request",required=True); a=p.parse_args(); path=Path(a.request).resolve(); req=json.loads(path.read_text()); d=path.parent
 with open(d/req["input"]["forcing"],newline="") as f: forcing=list(csv.DictReader(f))
 static=json.loads((d/req["input"]["static"]).read_text()); rows=simulate(forcing,static,MODE)
 with open(d/req["output"]["table"],"w",newline="") as f: w=csv.DictWriter(f,fieldnames=COLUMNS); w.writeheader(); w.writerows(rows)
 (d/req["output"]["run"]).write_text(json.dumps({"status":"ok"}))
if __name__=="__main__": main()

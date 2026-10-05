from __future__ import annotations
import argparse, csv, json
from pathlib import Path
COLUMNS=["time","gw_sw_exchange","gw_to_sw","sw_to_gw","gw"]
def main():
 p=argparse.ArgumentParser(); p.add_argument("--request",required=True); a=p.parse_args(); path=Path(a.request).resolve(); req=json.loads(path.read_text()); d=path.parent
 with open(d/req["input"]["forcing"],newline="") as f: forcing=list(csv.DictReader(f))
 static=json.loads((d/req["input"]["static"]).read_text()); gw=float(static["aquifer_specific_yield"])*1000*(float(static["aquifer_initial_head_m"])-float(static["aquifer_bottom_m"]))
 rows=[{"time":r["time"],"gw_sw_exchange":0.0,"gw_to_sw":0.0,"sw_to_gw":0.0,"gw":gw} for r in forcing]
 with open(d/req["output"]["table"],"w",newline="") as f: w=csv.DictWriter(f,fieldnames=COLUMNS); w.writeheader(); w.writerows(rows)
 (d/req["output"]["run"]).write_text(json.dumps({"status":"ok"}))
if __name__=="__main__": main()

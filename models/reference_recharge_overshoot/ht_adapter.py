"""Recharge-ordering reference: one aquifer cell draining to a river.

The three reference_recharge_* models share this file and differ only in MODE:

exact         the exact solution of S_y dh/dt = r + C/A (stage - h) over each
              day, so the budget closes and the added water drains on
              tau = S_y A / C, 200 days on this probe's case
blind         ignores gw_recharge in its storage and exchange, so the run's
              own budget no longer closes
overshoot     on a recharge day sends 120% of it to the river, so storage
              falls below the control's
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

MODE = "overshoot"
COLUMNS = ["time", "gw_sw_exchange", "gw_to_sw", "sw_to_gw", "gw"]


def simulate(forcing, static, mode=MODE):
    b = float(static["aquifer_specific_yield"]) * 1000.0  # mm of storage per m of head
    k = float(static["river_conductance_m2_per_day"]) / (float(static["area_km2"]) * 1.0e6) * 1000.0  # mm/day per m
    head = float(static["aquifer_initial_head_m"])
    bottom = float(static["aquifer_bottom_m"])
    rows = []
    for row in forcing:
        stage = float(row["sw_stage_m"])
        r = float(row["gw_recharge"])
        if mode == "blind":
            r = 0.0
        if mode == "overshoot" and r > 0.0:
            q = -1.2 * r
        else:
            equilibrium = stage + r / k
            new_head = equilibrium + (head - equilibrium) * math.exp(-k / b)
            q = b * (new_head - head) - r
        head += (r + q) / b
        rows.append({"time": row["time"], "gw_sw_exchange": q, "gw_to_sw": min(q, 0.0), "sw_to_gw": max(q, 0.0), "gw": b * (head - bottom)})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    path = Path(parser.parse_args().request).resolve()
    request = json.loads(path.read_text())
    root = path.parent
    with open(root / request["input"]["forcing"], newline="") as f:
        forcing = list(csv.DictReader(f))
    static = json.loads((root / request["input"]["static"]).read_text())
    rows = simulate(forcing, static)
    with open(root / request["output"]["table"], "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (root / request["output"]["run"]).write_text(json.dumps({"status": "ok"}))


if __name__ == "__main__":
    main()

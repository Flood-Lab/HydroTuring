#!/usr/bin/env python3
"""A confined aquifer with native storage and integrated river exchange.

Positive controls integrate b dh/dt = recharge + c (river_head - h), with
b = 1000 * confined storativity and c = 1000 * conductance / aquifer area.
One exact interval solution and backward Euler use separate flux quadratures.
The two named negative controls change only a head or storage multiplier;
they remain water-conservative and report the same directional conventions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

COLUMNS = ("time", "gw_sw_exchange", "gw_to_sw", "sw_to_gw", "gw")
TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1 / 24, "PT15M": 1 / 96,
                 "PT5M": 1 / 288, "PT1M": 1 / 1440}


def simulate(forcing, static, dt_days=1.0, method="exact", head_gain=1.0, storage_gain=1.0):
    """Integrate the supplied rows without resampling or a fitted equilibrium."""
    if method not in ("exact", "backward_euler"):
        raise ValueError("unknown integration method")
    if not math.isfinite(dt_days) or dt_days <= 0.0:
        raise ValueError("the integration step must be finite and positive")
    area = float(static["area_km2"]) * 1.0e6
    storativity = float(static["aquifer_storage_coefficient"])
    specific_yield = float(static["aquifer_specific_yield"])
    initial_head = float(static["aquifer_initial_head_m"])
    top = float(static["aquifer_top_m"])
    bottom = float(static["aquifer_bottom_m"])
    conductance = float(static["river_conductance_m2_per_day"])
    bed_offset = float(static["river_bottom_offset_m"])
    values = (area, storativity, specific_yield, initial_head, top, bottom,
              conductance, bed_offset, head_gain, storage_gain)
    if not all(math.isfinite(v) for v in values):
        raise ValueError("aquifer parameters must be finite")
    if min(area, storativity, specific_yield, conductance, bed_offset, head_gain, storage_gain) <= 0.0:
        raise ValueError("area, storage, conductance and multipliers must be positive")
    if not bottom < top < initial_head:
        raise ValueError("the reference requires a confined initial state")

    capacity = 1000.0 * storativity * storage_gain
    exchange_coefficient = 1000.0 * conductance / area
    initial_inventory = 1000.0 * specific_yield * (initial_head - bottom)
    head = initial_head
    rows = []
    for step in forcing:
        stage = float(step["sw_stage_m"])
        recharge = float(step["gw_recharge"])
        target = initial_head + head_gain * (stage - initial_head)
        if not all(math.isfinite(v) for v in (stage, recharge, target)):
            raise ValueError("aquifer forcing must be finite")
        if min(head, target) <= top or min(head, target) <= stage - bed_offset:
            raise ValueError("the reference needs confined heads and a connected river")
        if method == "exact":
            rate_step = exchange_coefficient * dt_days / capacity
            response = -math.expm1(-rate_step)
            mean_decay = response / rate_step
            difference = target - head
            # Integrate the constitutive river flux over the same interval;
            # do not infer a missing flux from the reported storage residual.
            q = (exchange_coefficient * difference * mean_decay
                 + recharge * (mean_decay - 1.0))
            head += (difference + recharge / exchange_coefficient) * response
        else:
            increment = dt_days * (recharge + exchange_coefficient * (target - head))
            head += increment / (capacity + exchange_coefficient * dt_days)
            q = exchange_coefficient * (target - head)
        storage = initial_inventory + capacity * (head - initial_head)
        rows.append({
            "time": step["time"], "gw_sw_exchange": q,
            "gw_to_sw": min(q, 0.0), "sw_to_gw": max(q, 0.0), "gw": storage,
        })
    return rows


def main(*, method="exact", head_gain=1.0, storage_gain=1.0):
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    request_path = Path(parser.parse_args().request).resolve()
    root = request_path.parent
    request = json.loads(request_path.read_text())
    with (root / request["input"]["forcing"]).open(newline="") as source:
        forcing = list(csv.DictReader(source))
    static = json.loads((root / request["input"]["static"]).read_text())
    dt_days = TIMESTEP_DAYS[request["timestep"]]
    rows = simulate(forcing, static, dt_days, method, head_gain, storage_gain)
    with (root / request["output"]["table"]).open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (root / request["output"]["run"]).write_text(json.dumps({
        "status": "ok", "n_steps": len(rows), "integration": method,
        "head_multiplier": head_gain, "storage_multiplier": storage_gain,
        "control_volume": "one confined aquifer and a connected prescribed-stage river boundary",
        "time_convention": "interval-mean mm/day exchange; interval-end absolute storage in mm",
        "not_reported": ["soil water", "snow", "canopy", "channel", "catchment runoff", "energy fluxes"],
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()

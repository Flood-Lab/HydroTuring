from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

COLUMNS = ["time", "gw_sw_exchange", "gw_to_sw", "sw_to_gw", "gw"]


def simulate(forcing: list[dict], static: dict) -> list[dict]:
    """A lagged-head linear-reservoir aquifer that closes its own
    groundwater balance exactly and reports its signed exchange components
    honestly, with no injected fault. The physical baseline for
    mass/gw-sw-exchange-consistency: trusted-subprocess, so the gate does
    not need Docker or a MODFLOW 6 build to score this probe.

    The head lags a smoothed river stage rather than equilibrating to
    today's stage, so it genuinely runs both above and below the river
    bottom over the record, giving exchange_directions a real
    reverse-direction magnitude on both sides."""
    conductance = float(static["river_conductance_m2_per_day"])
    area_m2 = float(static["area_km2"]) * 1.0e6
    bottom_offset = float(static["river_bottom_offset_m"])
    initial_head = float(static.get("aquifer_initial_head_m", 10.0))
    specific_yield = float(static["aquifer_specific_yield"])

    lag = 0.985
    smoothed_stage = float(forcing[0]["sw_stage_m"])
    head = initial_head
    gw_mm = specific_yield * head * 1000.0
    rows = []
    for row in forcing:
        recharge_mm = float(row["gw_recharge"])
        stage_m = float(row["sw_stage_m"])
        bottom_m = stage_m - bottom_offset
        smoothed_stage = lag * smoothed_stage + (1.0 - lag) * stage_m
        head = smoothed_stage

        # q_m3_day > 0 is river losing to the aquifer (sw_to_gw), matching
        # MODFLOW's RIV convention and this probe's
        # gw_sw_exchange-positive-into-aquifer rule. Below the river bottom,
        # MODFLOW's RIV package caps the exchange at the bottom elevation
        # instead of the (unreachable) head.
        effective_head = max(head, bottom_m)
        q_m3_day = conductance * (stage_m - effective_head)
        q_mm_day = q_m3_day / area_m2 * 1000.0
        sw_to_gw = max(q_mm_day, 0.0)
        gw_to_sw = min(q_mm_day, 0.0)
        gw_sw_exchange = gw_to_sw + sw_to_gw

        gw_mm += recharge_mm + gw_sw_exchange

        rows.append({
            "time": row["time"],
            "gw_sw_exchange": gw_sw_exchange,
            "gw_to_sw": gw_to_sw,
            "sw_to_gw": sw_to_gw,
            "gw": gw_mm,
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    io_dir = request_path.parent
    request = json.loads(request_path.read_text())
    with open(io_dir / request["input"]["forcing"], newline="") as handle:
        forcing = list(csv.DictReader(handle))
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    rows = simulate(forcing, static)
    with open(io_dir / request["output"]["table"], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(json.dumps({"status": "ok"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

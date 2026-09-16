from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

COLUMNS = ["time", "gw_sw_exchange", "gw_to_sw", "sw_to_gw", "gw"]


def simulate(forcing: list[dict], static: dict) -> list[dict]:
    """A lagged-head bookkeeping reference (not a linear-reservoir aquifer;
    see reference_exchange_exact) that closes its own budget exactly, but
    on every other scored step reports the sign of gw_to_sw flipped: a
    positive aquifer-to-river term (should be <= 0) instead of negating it,
    so gw_to_sw + sw_to_gw != gw_sw_exchange on those steps. The head lags a
    smoothed river stage rather than equilibrating to today's stage, so it
    genuinely runs both above and below the river bottom over the record;
    the untouched steps still carry a real negative gw_to_sw, so
    exchange_directions sees a genuine reverse-direction magnitude and only
    exchange_components, which checks the running sum against
    gw_sw_exchange, catches the error."""
    conductance = float(static["river_conductance_m2_per_day"])
    area_m2 = float(static["area_km2"]) * 1.0e6
    bottom_offset = float(static["river_bottom_offset_m"])
    initial_head = float(static.get("aquifer_initial_head_m", 10.0))
    specific_yield = float(static["aquifer_specific_yield"])

    # The head tracks a heavily smoothed (lagged) river stage: an aquifer
    # that responds far slower than the river, so it is out of phase with
    # today's stage and genuinely sits below the river bottom for part of
    # the cycle and above it for the rest.
    lag = 0.985
    smoothed_stage = float(forcing[0]["sw_stage_m"])
    head = initial_head
    gw_mm = specific_yield * head * 1000.0
    rows = []
    for step, row in enumerate(forcing):
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
        gw_to_sw_true = min(q_mm_day, 0.0)
        gw_sw_exchange = gw_to_sw_true + sw_to_gw

        gw_mm += recharge_mm + gw_sw_exchange

        # Deliberate sign error, on every other step only: report the
        # aquifer-to-river magnitude without negating it, so the reported
        # components no longer sum to the correctly-signed gw_sw_exchange on
        # those steps. The untouched steps keep a real gw_to_sw <= 0, so
        # exchange_directions still sees a genuine reverse-direction
        # magnitude; only exchange_components catches the flipped steps.
        gw_to_sw_reported = abs(gw_to_sw_true) if step % 2 == 0 else gw_to_sw_true

        rows.append({
            "time": row["time"],
            "gw_sw_exchange": gw_sw_exchange,
            "gw_to_sw": gw_to_sw_reported,
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

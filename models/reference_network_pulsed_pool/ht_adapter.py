#!/usr/bin/env python3
"""Conserves routing mass while retaining a downstream pool."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    io_dir = request_path.parent
    request = json.loads(request_path.read_text())

    forcing_path = io_dir / request["input"]["forcing"]
    static_path = io_dir / request["input"]["static"]
    static = json.loads(static_path.read_text())

    reaches = static["routing_network"]["reaches"]
    if set(reaches) != {"A", "B", "C"}:
        raise ValueError("reference_network_pulsed_pool requires reaches A, B, and C")

    with forcing_path.open(newline="") as fh:
        forcing = list(csv.DictReader(fh))

    table_path = io_dir / request["output"]["table"]
    table_path.parent.mkdir(parents=True, exist_ok=True)
    with table_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["time"])
        writer.writeheader()
        for row in forcing:
            writer.writerow({"time": row["time"]})

    routing_path = io_dir / request["output"]["routing"]
    routing_path.parent.mkdir(parents=True, exist_ok=True)
    with routing_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["time", "reach_id", "reach_q_in", "q_out", "channel_storage"],
        )
        writer.writeheader()
        seconds_per_day = 86400.0
        spinup_days = 365
        if len(forcing) != 455:
            raise ValueError("reference_network_pulsed_pool requires 455 daily rows")

        last_inflow = max(
            i for i in range(spinup_days, len(forcing))
            if float(forcing[i]["q_in_A"]) > 0.0
            or float(forcing[i]["q_in_B"]) > 0.0
        )
        pulses = {last_inflow, last_inflow + 31, last_inflow + 62}
        if max(pulses) >= len(forcing):
            raise ValueError("routing case is too short for the pulse control")

        areas = static["routing_network"]["contributing_area_km2"]
        target_a = 0.049 * float(areas["A"]) * 1000.0
        target_b = 0.049 * float(areas["B"]) * 1000.0
        target_c = 0.30 * float(areas["C"]) * 1000.0

        storage_a = 0.0
        storage_b = 0.0
        storage_c = 0.0

        for i, row in enumerate(forcing):
            inflow_a = float(row["q_in_A"])
            inflow_b = float(row["q_in_B"])

            capture_a = 0.0
            capture_b = 0.0
            if spinup_days <= i <= last_inflow:
                capture_a = min(
                    inflow_a * seconds_per_day,
                    max(0.0, target_a - storage_a),
                )
                capture_b = min(
                    inflow_b * seconds_per_day,
                    max(0.0, target_b - storage_b),
                )

            storage_a += capture_a
            storage_b += capture_b

            release_a = 0.0
            release_b = 0.0
            if i in pulses:
                release_a = min(storage_a, target_a / 3.0)
                release_b = min(storage_b, target_b / 3.0)

            storage_a -= release_a
            storage_b -= release_b
            outflow_a = inflow_a + (release_a - capture_a) / seconds_per_day
            outflow_b = inflow_b + (release_b - capture_b) / seconds_per_day

            inflow_c = outflow_a + outflow_b
            capture_c = 0.0
            if spinup_days <= i <= last_inflow:
                capture_c = min(
                    inflow_c * seconds_per_day,
                    max(0.0, target_c - storage_c),
                )

            storage_c += capture_c
            outflow_c = inflow_c - capture_c / seconds_per_day

            inflows = {"A": inflow_a, "B": inflow_b, "C": inflow_c}
            outflows = {"A": outflow_a, "B": outflow_b, "C": outflow_c}
            storages = {"A": storage_a, "B": storage_b, "C": storage_c}

            for reach in reaches:
                writer.writerow(
                    {
                        "time": row["time"],
                        "reach_id": reach,
                        "reach_q_in": inflows[reach],
                        "q_out": outflows[reach],
                        "channel_storage": storages[reach],
                    }
                )

    run_path = io_dir / request["output"]["run"]
    run_path.write_text(
        json.dumps({"status": "ok", "n_steps": len(forcing)}, indent=2)
    )


if __name__ == "__main__":
    main()
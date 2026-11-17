#!/usr/bin/env python3
"""Deliberately retains part of the inflow in channel storage."""

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
        raise ValueError("reference_network_stuck_router requires reaches A, B, and C")

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
        storage_a = 0.0
        for row in forcing:
            inflow_a = float(row["q_in_A"])
            inflow_b = float(row["q_in_B"])

            retained_a = 0.1 * inflow_a
            storage_a += retained_a * 86400.0
            outflow_a = inflow_a - retained_a

            inflows = {
                "A": inflow_a,
                "B": inflow_b,
                "C": outflow_a + inflow_b,
            }
            outflows = {
                "A": outflow_a,
                "B": inflow_b,
                "C": outflow_a + inflow_b,
            }
            storages = {"A": storage_a, "B": 0.0, "C": 0.0}

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
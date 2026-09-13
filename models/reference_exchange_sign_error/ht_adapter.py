from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

COLUMNS = ["time", "gwex", "gw_to_sw", "sw_to_gw", "gw"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    io_dir = request_path.parent
    request = json.loads(request_path.read_text())
    with open(io_dir / request["input"]["forcing"], newline="") as handle:
        forcing = list(csv.DictReader(handle))
    rows = []
    for row in forcing:
        # Deliberate sign error: sw_to_gw must be non-positive but is positive.
        gw_to_sw = 0.01
        sw_to_gw = 0.01
        rows.append({
            "time": row["time"],
            "gwex": gw_to_sw + sw_to_gw,
            "gw_to_sw": gw_to_sw,
            "sw_to_gw": sw_to_gw,
            "gw": 20.0,
        })
    with open(io_dir / request["output"]["table"], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(json.dumps({"status": "ok"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

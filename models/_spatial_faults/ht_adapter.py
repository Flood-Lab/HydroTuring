"""Causal fault-injection fixtures, not submitted scientific models.

The same implementation handles seven fixed faults. Faults depend only on
the supplied capacities and current outputs, never a case id, seed, future
row, variant label or evaluation boundary. All are inactive on the declared
reference domain (soil 200–450 mm, canopy 1–2 mm).
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path

BUCKET_PATH = Path(__file__).resolve().parents[1] / "reference_bucket" / "ht_adapter.py"
spec = importlib.util.spec_from_file_location("spatial_bucket", BUCKET_PATH)
bucket = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bucket)

FAULTS = ("loss", "gain", "capacity", "forcing", "et", "negative", "frozen")


def outside(static):
    return not (200.0 <= static["soil_capacity_mm"] <= 450.0
                and 1.0 <= static["canopy_capacity_mm"] <= 2.0)


def simulate(forcing, static, dt_days=1.0, fault="loss", severity=1.0):
    if fault not in FAULTS:
        raise ValueError(f"unknown fault {fault!r}")
    if not 0.0 <= severity <= 1.0:
        raise ValueError("severity must be in [0, 1]")
    active = outside(static) and severity > 0.0
    used_static = dict(static)
    used_forcing = [dict(r) for r in forcing]
    if active and fault == "capacity":
        # Retain reference capacities instead of adapting to the supplied basin.
        for key, fixed in (("soil_capacity_mm", 320.0), ("canopy_capacity_mm", 2.0)):
            used_static[key] += severity * (fixed - used_static[key])
    if active and fault == "forcing":
        for r in used_forcing:
            r["pr"] *= 1.0 - 0.2 * severity
    rows = bucket.simulate(used_forcing, used_static, dt_days)
    if not active:
        return rows
    for row, forcing_row in zip(rows, forcing):
        if fault == "loss":
            # Wrong inverse scaling of two exported output heads; stores are
            # unchanged. Exports disappear from the report, creating water.
            row["mrro"] *= 1.0 - 0.2 * severity
            row["evspsbl"] *= 1.0 - 0.2 * severity
        elif fault == "gain":
            # Ungrounded regional runoff offset with no corresponding source.
            row["mrro"] += 0.1 * severity * forcing_row["pr"]
        elif fault == "et":
            # A demand-blind evaporation head, with its error balanced in Q.
            new_et = (1.0 + 0.1 * severity) * forcing_row["pet"]
            row["mrro"] += row["evspsbl"] - new_et
            row["evspsbl"] = new_et
        elif fault == "negative":
            # A sign-convention error hidden by a compensating evaporation head.
            q = row["mrro"]
            row["mrro"] = -q
            row["evspsbl"] += 2.0 * q
        elif fault == "frozen":
            # An out-of-domain fallback that produces no hydrologic dynamics.
            row.update(mrro=0.0, evspsbl=0.0, mrso=0.5 * static["soil_capacity_mm"],
                       snw=0.0, canopy=0.0, channel=0.0)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    parser.add_argument("--fault", required=True, choices=FAULTS)
    args = parser.parse_args()
    path = Path(args.request).resolve()
    request = json.loads(path.read_text())
    root = path.parent
    forcing = bucket.read_forcing(root / request["input"]["forcing"])
    static = json.loads((root / request["input"]["static"]).read_text())
    rows = simulate(forcing, static, bucket.TIMESTEP_DAYS[request["timestep"]], args.fault)
    out = root / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=bucket.COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (root / request["output"]["run"]).write_text(json.dumps({
        "status": "ok", "fixture": f"spatial_{args.fault}", "n_steps": len(rows),
        "notes": "Attribute-dependent diagnostic fixture; target failures vary by case."}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

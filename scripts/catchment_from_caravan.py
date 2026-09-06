#!/usr/bin/env python3
"""Write catchments/<id>.json from a local copy of Caravan.

    python scripts/catchment_from_caravan.py /path/to/Caravan/attributes/camels \
        camels_09510200 semi-arid "Semi-arid ..."

Caravan's attributes directory holds, per source dataset, three CSVs
(attributes_hydroatlas_*.csv, attributes_caravan_*.csv,
attributes_other_*.csv) keyed by gauge_id. This copies one basin's rows into
the shape the harness reads, so a new catchment is one command rather than
a hand-edited file.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


def _rows(path: Path, gauge_id: str) -> dict:
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["gauge_id"] == gauge_id:
                return row
    raise SystemExit(f"{gauge_id} not in {path}")


def _num(value: str):
    try:
        f = float(value)
    except ValueError:
        return value
    return int(f) if f.is_integer() and "." not in value else f


def main(argv: list[str]) -> int:
    if len(argv) != 5:
        print(__doc__)
        return 2
    attr_dir, gauge_id, name, description = Path(argv[1]), argv[2], argv[3], argv[4]
    source = attr_dir.name
    hydroatlas = _rows(attr_dir / f"attributes_hydroatlas_{source}.csv", gauge_id)
    caravan = _rows(attr_dir / f"attributes_caravan_{source}.csv", gauge_id)
    other = _rows(attr_dir / f"attributes_other_{source}.csv", gauge_id)
    doc = {
        "id": name,
        "description": description,
        "source": {
            "gauge_id": gauge_id,
            "gauge_name": other["gauge_name"],
            "country": other["country"],
            "dataset": f"Caravan, {source} subset; HydroATLAS attributes aggregated by Caravan",
            "copy": str(attr_dir),
        },
        "area_km2": float(other["area"]),
        "latitude_deg": float(other["gauge_lat"]),
        "longitude_deg": float(other["gauge_lon"]),
        "caravan": {k: _num(v) for k, v in caravan.items() if k != "gauge_id"},
        "hydroatlas": {k: _num(v) for k, v in hydroatlas.items() if k != "gauge_id"},
    }
    out = Path(__file__).resolve().parent.parent / "catchments" / f"{name}.json"
    out.write_text(json.dumps(doc, indent=1))
    print(f"wrote {out} ({other['gauge_name']}, {float(other['area']):.0f} km2)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

#!/usr/bin/env python3
"""Cut the domain the adapter runs on out of upstream's test domain, at build time.

Upstream's `test_domain` is the Mosel basin: 11,636 km2 of real morphological
grids at 500 m. A lumped HydroTuring case is one catchment of a few hundred km2,
so this takes the natural sub-basin draining through one cell of that river
network, follows upstream's own flow directions to find it, and writes every
grid back with everything outside it as nodata. No morphological value is
invented or changed. The one thing that is new is the land surface's phenology
table, made flat over the year (see `flat_lai`), and README.md says why.

Nothing here runs at evaluation time. Its output is /model/domain.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import netCDF4
import numpy as np

UPSTREAM = Path("/model/upstream/test_domain/input")
DOMAIN = Path(sys.argv[1] if len(sys.argv) > 1 else "/model/domain")

# Row and column in test_domain's 500 m grid (row 0 at the north) of the most
# downstream cell of the chosen sub-basin: 286 km2, of the size the probes state
# for their catchment, whose land cover (34% forest, 3% urban, 63% pervious) is
# the nearest of the candidates of that size to the whole basin's. The further
# arguments build the variants README.md reports and nothing else uses.
OUTLET_ROW = int(sys.argv[2]) if len(sys.argv) > 2 else 92
OUTLET_COL = int(sys.argv[3]) if len(sys.argv) > 3 else 65

# Hydrology, routing and forcing all run at this resolution, as upstream's test
# domain has them all equal. Its 24 km would give this sub-basin one or two
# cells, and no river; 4 km gives it a network. README.md has what 8 km moves.
RESOLUTION_M = int(sys.argv[4]) if len(sys.argv) > 4 else 4000
MONTHLY_LAI = len(sys.argv) > 5 and sys.argv[5] == "monthly"
L0_M = 500
FACTOR = RESOLUTION_M // L0_M
GAUGE_ID = 1

# ArcGIS D8 codes, row increasing southwards.
D8 = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}


def read_asc(path: Path) -> tuple[dict, np.ndarray]:
    header = {}
    with open(path) as fh:
        for _ in range(6):
            key, value = fh.readline().split()
            header[key.lower()] = float(value)
    return header, np.loadtxt(path, skiprows=6)


def write_asc(path: Path, header: dict, data: np.ndarray) -> None:
    with open(path, "w") as fh:
        fh.write(f"ncols {data.shape[1]}\nnrows {data.shape[0]}\n")
        fh.write(f"xllcorner {header['xllcorner']:.1f}\nyllcorner {header['yllcorner']:.1f}\n")
        fh.write(f"cellsize {header['cellsize']:.1f}\nNODATA_value -9999\n")
        np.savetxt(fh, data, fmt="%.10g")


def upstream_of(fdir: np.ndarray, mask: np.ndarray, row: int, col: int) -> np.ndarray:
    nrows, ncols = fdir.shape
    into: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for r, c in np.argwhere(mask):
        code = int(fdir[r, c])
        if code in D8:
            rr, cc = r + D8[code][0], c + D8[code][1]
            if 0 <= rr < nrows and 0 <= cc < ncols and mask[rr, cc]:
                into.setdefault((rr, cc), []).append((r, c))
    keep = np.zeros_like(mask)
    stack = [(row, col)]
    while stack:
        cell = stack.pop()
        keep[cell] = True
        stack += into.get(cell, [])
    return keep


def flat_lai(path: Path, out: Path) -> None:
    """Upstream's LAI table with each class at its own annual mean, every month.

    The table is an input, not model code, and upstream's holds a Central European
    seasonal cycle. Nothing in a case says when this catchment leafs out, so
    keeping it would hand the model a calendar the forcing does not have.
    """
    lines = path.read_text().splitlines()
    result = []
    for line in lines:
        fields = line.split()
        if len(fields) == 14 and fields[0].isdigit():
            mean = np.mean([float(v) for v in fields[2:]])
            result.append(" ".join(fields[:2] + [f"{mean:.5f}"] * 12))
        else:
            result.append(line)
    out.write_text("\n".join(result) + "\n")


def main() -> int:
    morph_in = UPSTREAM / "morph"
    morph = DOMAIN / "morph"
    luse = DOMAIN / "luse"
    morph.mkdir(parents=True, exist_ok=True)
    luse.mkdir(parents=True, exist_ok=True)

    header, dem = read_asc(morph_in / "dem.asc")
    _, fdir = read_asc(morph_in / "fdir.asc")
    keep = upstream_of(fdir, dem > -9000, OUTLET_ROW, OUTLET_COL)

    # Crop to a whole number of hydrology cells around it.
    rows, cols = np.nonzero(keep)
    r0 = rows.min() - rows.min() % FACTOR
    c0 = cols.min() - cols.min() % FACTOR
    nrows = int(np.ceil((rows.max() + 1 - r0) / FACTOR)) * FACTOR
    ncols = int(np.ceil((cols.max() + 1 - c0) / FACTOR)) * FACTOR
    r1, c1 = r0 + nrows, c0 + ncols
    full_nrows = dem.shape[0]
    cropped = dict(header)
    cropped["xllcorner"] = header["xllcorner"] + c0 * L0_M
    cropped["yllcorner"] = header["yllcorner"] + (full_nrows - r1) * L0_M
    sub = keep[r0:r1, c0:c1]

    for name in ["dem", "slope", "aspect", "fdir", "facc", "soil_class", "soil_class_horizon_01",
                 "soil_class_horizon_02", "geology_class", "LAI_class"]:
        _, grid = read_asc(morph_in / f"{name}.asc")
        write_asc(morph / f"{name}.asc", cropped, np.where(sub, grid[r0:r1, c0:c1], -9999))
    _, cover = read_asc(UPSTREAM / "luse" / "lc_1991.asc")
    write_asc(luse / "lc_1991.asc", cropped, np.where(sub, cover[r0:r1, c0:c1], -9999))

    gauges = np.full(sub.shape, -9999.0)
    gauges[OUTLET_ROW - r0, OUTLET_COL - c0] = GAUGE_ID
    write_asc(morph / "idgauges.asc", cropped, gauges)

    for name in ["soil_classdefinition.txt", "soil_classdefinition_iFlag_soilDB_1.txt", "geology_classdefinition.txt"]:
        shutil.copy(morph_in / name, morph / name)
    if MONTHLY_LAI:
        shutil.copy(morph_in / "LAI_classdefinition.txt", morph / "LAI_classdefinition.txt")
    else:
        flat_lai(morph_in / "LAI_classdefinition.txt", morph / "LAI_classdefinition.txt")

    # mHM reads latitude and longitude at each level whatever it computes with
    # them. With PET as an input they label the output and set the hemisphere.
    with netCDF4.Dataset(UPSTREAM / "latlon" / "latlon_1.nc") as src:
        lat0 = np.array(src["lat_l0"][:])
        lon0 = np.array(src["lon_l0"][:])
        y0_ascending = bool(src["yc_l0"][0] < src["yc_l0"][-1])
    # Put the upstream arrays north-up to match the rows above.
    if y0_ascending:
        lat0, lon0 = lat0[::-1], lon0[::-1]
    lat_sub, lon_sub = lat0[r0:r1, c0:c1], lon0[r0:r1, c0:c1]
    ny1, nx1 = nrows // FACTOR, ncols // FACTOR
    lat1 = lat_sub.reshape(ny1, FACTOR, nx1, FACTOR).mean(axis=(1, 3))
    lon1 = lon_sub.reshape(ny1, FACTOR, nx1, FACTOR).mean(axis=(1, 3))
    xll, yll = cropped["xllcorner"], cropped["yllcorner"]
    with netCDF4.Dataset(DOMAIN / "latlon.nc", "w") as nc:
        def grid(suffix: str, nx: int, ny: int, cell: float, lat: np.ndarray, lon: np.ndarray) -> None:
            nc.createDimension(f"xc{suffix}", nx)
            nc.createDimension(f"yc{suffix}", ny)
            nc.createVariable(f"xc{suffix}", "f8", (f"xc{suffix}",))[:] = xll + (np.arange(nx) + 0.5) * cell
            nc.createVariable(f"yc{suffix}", "f8", (f"yc{suffix}",))[:] = yll + (ny - 1 - np.arange(ny) + 0.5) * cell
            nc.createVariable(f"lon{suffix}", "f8", (f"yc{suffix}", f"xc{suffix}"))[:] = lon
            nc.createVariable(f"lat{suffix}", "f8", (f"yc{suffix}", f"xc{suffix}"))[:] = lat
        grid("_l0", ncols, nrows, L0_M, lat_sub, lon_sub)
        grid("", nx1, ny1, RESOLUTION_M, lat1, lon1)
        grid("_l11", nx1, ny1, RESOLUTION_M, lat1, lon1)

    # The level-1 area of each cell is its count of level-0 cells: what the
    # adapter weights a cell's depth by. Upstream's outputs are per cell.
    counts = sub.reshape(ny1, FACTOR, nx1, FACTOR).sum(axis=(1, 3))
    meta = {
        "xll": xll,
        "yll": yll,
        "resolution_m": RESOLUTION_M,
        "ncols_l1": int(nx1),
        "nrows_l1": int(ny1),
        "l0_cell_counts_north_up": counts.astype(int).tolist(),
        "area_km2": float(sub.sum()) * (L0_M / 1000.0) ** 2,
        "outlet_x": xll + (OUTLET_COL - c0 + 0.5) * L0_M,
        "outlet_y": yll + (nrows - (OUTLET_ROW - r0) - 0.5) * L0_M,
        "gauge_id": GAUGE_ID,
        "land_cover_shares": {
            name: float((cover[r0:r1, c0:c1][sub] == code).mean())
            for code, name in [(1, "forest"), (2, "impervious"), (3, "pervious")]
        },
    }
    (DOMAIN / "domain.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps({k: v for k, v in meta.items() if k != "l0_cell_counts_north_up"}, indent=2))
    print("level-1 cells", int((counts > 0).sum()), "of", counts.size)
    return 0


if __name__ == "__main__":
    sys.exit(main())

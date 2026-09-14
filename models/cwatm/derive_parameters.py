"""Recompute the parameter table in README.md from CWatM's own 30' input maps.

The adapter takes every quantity CWatM reads from a map, and that the probe
does not give, as the median over the land cells of CWatM's 30-arcminute
input set. This script prints those medians with the mean and the 10th and
90th percentiles, and the area-weighted forest share of forest plus grassland.
It is provenance, not part of the image.

Fetch the maps (Git LFS objects) at the commit the adapter used, then run it
anywhere numpy and netCDF4 are installed, the model image included:

    git clone https://github.com/iiasa/CWatM-Earth-30min
    git -C CWatM-Earth-30min checkout e9dfd99e7feac31c7a03d62cb94de91dd6a28624
    git -C CWatM-Earth-30min lfs pull
    python3 derive_parameters.py CWatM-Earth-30min

A land cell is a cell with a valid local drain direction (1-9). The 10-day
crop-coefficient and interception-capacity maps are averaged over the year
per cell before the median is taken.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def grid(path: Path, var: str | None = None) -> np.ndarray:
    with Dataset(path) as nc:
        if var is None:
            names = [v for v in nc.variables if v not in ("lat", "lon", "x", "y", "time", "crs")]
            var = names[-1]
        return np.ma.filled(nc.variables[var][:].astype(float), np.nan)


def stats(a: np.ndarray, mask: np.ndarray) -> dict:
    if a.ndim == 3:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # ocean cells are empty in every period
            a = np.nanmean(a, axis=0)
    v = a[mask & np.isfinite(a)]
    return {"median": float(np.median(v)), "mean": float(np.mean(v)),
            "p10": float(np.percentile(v, 10)), "p90": float(np.percentile(v, 90)), "n": int(v.size)}


def main(root: Path) -> dict:
    ldd = grid(root / "routing/ldd.nc")
    mask = np.isfinite(ldd) & (ldd >= 1) & (ldd <= 9)
    out: dict = {"land_cells": int(mask.sum())}

    soil = {}
    for layer in (1, 2, 3):
        for name in ("ksat", "alpha", "lambda", "thetas", "thetar"):
            soil[f"{name}{layer}"] = stats(grid(root / f"soil/{name}{layer}.nc"), mask)
            if layer < 3:
                soil[f"forest_{name}{layer}"] = stats(grid(root / f"soil/forest_{name}{layer}.nc"), mask)
    for name in ("percolationImp", "storageDepth1", "storageDepth2", "cropgrp"):
        soil[name] = stats(grid(root / f"soil/{name}.nc"), mask)
    out["soil"] = soil

    cover = {}
    for lc, cap in (("forest", "Forest"), ("grassland", "Grassland")):
        cover[lc] = {
            "cropCoefficient": stats(grid(root / f"landcover/{lc}/cropCoefficient{cap}_10days.nc"), mask),
            "interceptCap": stats(grid(root / f"landcover/{lc}/interceptCap{cap}_10days.nc"), mask),
            "maxRootDepth": stats(grid(root / f"landcover/{lc}/maxRootDepth.nc"), mask),
            "rootFraction1": stats(grid(root / f"landcover/{lc}/rootFraction1.nc"), mask),
        }
    for lc in ("irrPaddy", "irrNonPaddy"):
        cover[lc] = {"maxRootDepth": stats(grid(root / f"landcover/{lc}/maxRootDepth.nc"), mask),
                     "rootFraction1": stats(grid(root / f"landcover/{lc}/rootFraction1.nc"), mask)}
    out["landcover"] = cover

    out["groundwater"] = {n: stats(grid(root / f"groundwater/{n}.nc"), mask) for n in ("recessionCoeff", "specificYield")}
    out["topo"] = {"elvstd": stats(grid(root / "landsurface/topo/elvstd.nc"), mask),
                   "tanslope": stats(grid(root / "landsurface/topo/tanslope.nc"), mask)}
    dz = root / "landsurface/topo/dzRel_hydro1k.nc"
    with Dataset(dz) as nc:
        names = [v for v in nc.variables if v.startswith("dzRel")]
    out["topo"]["dzRel"] = {n: stats(grid(dz, n), mask) for n in names}

    cellarea = grid(root / "routing/cellarea.nc")
    with Dataset(root / "landsurface/fractionLandcover.nc") as nc:
        f = np.ma.filled(nc.variables["fracforest"][-1].astype(float), np.nan)
        g = np.ma.filled(nc.variables["fracgrassland"][-1].astype(float), np.nan)
    ok = mask & np.isfinite(f) & np.isfinite(g) & ((f + g) > 0) & np.isfinite(cellarea)
    out["forest_share_of_forest_plus_grassland"] = float(
        np.sum(f[ok] * cellarea[ok]) / np.sum((f[ok] + g[ok]) * cellarea[ok]))
    return out


if __name__ == "__main__":
    print(json.dumps(main(Path(sys.argv[1] if len(sys.argv) > 1 else "CWatM-Earth-30min")), indent=1))

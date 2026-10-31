#!/usr/bin/env python3
"""HydroTuring adapter for the Entity-Aware LSTM of Kratzert et al. (2019).

The network, its weights and its statistics are upstream's and are used
unmodified: `main.Model` and `papercode.ealstm.EALSTM` from
kratzert/ealstm_regional_modeling, the published `run_1906_1004_seed111`
checkpoint from HydroShare (EA-LSTM trained with the NSE loss), the input and
output statistics hard-coded in `papercode.datautils.SCALER`, and the attribute
statistics that `main.evaluate` takes over the 531 basins in the run's own
`attributes.db`. What this file adds is the mapping from a probe's forcing to
the inputs the network was trained on.

The model predicts streamflow and nothing else, so it is declared as `mrro`
plus `dis`, with `states: []`, and is N/A (INCOMPLETE) on every budget probe.

How the forcing is mapped
-------------------------
The network reads five Maurer forcings a probe does not all generate. Only
precipitation is given directly. Everything else is built from the row's own
pr, tas and pet, with no calendar dependence, and is listed under
`mock_inputs` in run.json:

* Tmax and Tmin are tas plus and minus half a fixed diurnal range, the
  difference between the Tmax and Tmin training means in SCALER.
* Vp follows FAO-56 eq. 48, the saturation vapour pressure at Tmin, in Pa.
* SRAD inverts Priestley-Taylor from pet, adds FAO-56 net longwave, removes
  albedo, and doubles the 24-hour mean to a daylight mean under a fixed 12-hour
  day. The doubling is because the training statistics (372 W/m2 on average)
  are a daylight mean, not a 24-hour one; the true day length follows the
  calendar and so cannot be used.

The 27 static attributes: the eight climate indices with an unambiguous CAMELS
definition (Addor et al. 2017) are derived from the first year of the forcing
the model is given, over days at any step, since that is how they are defined. The
other 19, p_seasonality included, are held at the
training mean, which is zero after standardisation.

How the record is simulated
---------------------------
As `main.evaluate` does it: many-to-one, a fresh state per row, and a window of
`seq_length` (270) rows ending on that row. Rows with less history behind them
use the history that exists. Predictions are de-standardised with upstream's
`rescale_features`, and negative values are set to zero, as `evaluate_basin`
does. Rows are fed as given at whatever step the case runs: there is no
resampling, so at PT1H the 270-row window covers 11 days and each hourly rate
is read as a daily one. That is what the resolution probe exists to measure.

Sensitivity switches, for the table in README.md and never set in the image:
EALSTM_HOLD lists attribute groups to hold at the training mean, or "all";
EALSTM_MOCK=mean holds SRAD, Tmax, Tmin and Vp at their training means;
EALSTM_AREA=static feeds the case's area to area_gages2.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch

UPSTREAM = Path(os.environ.get("EALSTM_UPSTREAM", "/model/upstream"))
RUN_DIR = Path(os.environ.get("EALSTM_RUN_DIR", "/model/run"))
sys.path.insert(0, str(UPSTREAM))

from main import Model  # noqa: E402  upstream's network wrapper
from papercode.datautils import SCALER, load_attributes, normalize_features, rescale_features  # noqa: E402
from papercode.utils import get_basin_list  # noqa: E402

MODEL = {"name": "ealstm", "version": "0.1.2-d118158-seed111"}
COLUMNS = ["time", "mrro", "dis"]
WEIGHTS = "model_epoch30.pt"
THREADS = int(os.environ.get("EALSTM_THREADS", "4"))
BATCH = 256

TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0, "PT15M": 1.0 / 96.0, "PT5M": 1.0 / 288.0, "PT1M": 1.0 / 1440.0}
CLIMATOLOGY_DAYS = 365.0
DAYS_PER_YEAR = 365.25

# Which attributes the first year of forcing can give, by group, for the
# sensitivity table. Every other attribute stays at the training mean.
ATTRIBUTE_GROUPS = {
    "precipitation": ["p_mean", "high_prec_freq", "high_prec_dur", "low_prec_freq", "low_prec_dur"],
    "demand": ["pet_mean", "aridity"],
    "snow": ["frac_snow"],
}

# Mock forcing. Order as in `SCALER`: PRCP, SRAD, Tmax, Tmin, Vp.
ALBEDO = 0.23
PRIESTLEY_TAYLOR_ALPHA = 1.26
LATENT_HEAT = 2.45  # MJ kg-1
STEFAN_BOLTZMANN = 4.903e-9  # MJ K-4 m-2 day-1
MJ_PER_DAY_TO_W = 1.0e6 / 86400.0
CLEARNESS_DRY, CLEARNESS_WET = 0.70 / 0.75, 0.45 / 0.75
WET_DAY_MM = 1.0
DAYLIGHT_FRACTION = 0.5
MOCK_INPUT_NOTES = {
    "Tmax, Tmin": "tas plus and minus half the difference of the Tmax and Tmin training means; degC",
    "Vp": "FAO-56 eq. 48, saturation vapour pressure at Tmin; Pa",
    "SRAD": (
        f"Priestley-Taylor (alpha {PRIESTLEY_TAYLOR_ALPHA}) inverted from pet, plus FAO-56 net "
        f"longwave, over 1 - albedo ({ALBEDO}), doubled to a daylight mean under a fixed 12-hour "
        "day because the training statistics are a daylight mean; W/m2"
    ),
}


def mock_inputs(pr: np.ndarray, tas: np.ndarray, pet: np.ndarray, elevation_m: float) -> np.ndarray:
    """The five Maurer inputs per row, from that row's forcing and nothing else."""
    half_range = 0.5 * float(SCALER["input_means"][2] - SCALER["input_means"][3])
    tmax, tmin = tas + half_range, tas - half_range
    vp_kpa = 0.6108 * np.exp(17.27 * tmin / (tmin + 237.3))

    pressure = 101.3 * ((293.0 - 0.0065 * elevation_m) / 293.0) ** 5.26
    saturation_kpa = 0.6108 * np.exp(17.27 * tas / (tas + 237.3))
    slope = 4098.0 * saturation_kpa / (tas + 237.3) ** 2
    gamma = 0.000665 * pressure
    net_radiation = pet * LATENT_HEAT * (slope + gamma) / (PRIESTLEY_TAYLOR_ALPHA * slope)
    clearness = np.where(pr >= WET_DAY_MM, CLEARNESS_WET, CLEARNESS_DRY)
    net_longwave = (
        STEFAN_BOLTZMANN * (tas + 273.16) ** 4
        * (0.34 - 0.14 * np.sqrt(vp_kpa))
        * (1.35 * clearness - 0.35)
    )
    incoming = np.maximum(net_radiation + net_longwave, 0.0) / (1.0 - ALBEDO)
    srad = incoming * MJ_PER_DAY_TO_W / DAYLIGHT_FRACTION

    if os.environ.get("EALSTM_MOCK") == "mean":
        means = SCALER["input_means"]
        return np.column_stack([pr, np.full_like(pr, means[1]), np.full_like(pr, means[2]),
                                np.full_like(pr, means[3]), np.full_like(pr, means[4])])
    return np.column_stack([pr, srad, tmax, tmin, vp_kpa * 1000.0])


def _mean_run_length(mask: np.ndarray) -> float:
    """Average length of the runs of True, in rows."""
    runs, length = [], 0
    for flag in mask:
        if flag:
            length += 1
        elif length:
            runs.append(length)
            length = 0
    if length:
        runs.append(length)
    return float(np.mean(runs)) if runs else 0.0


def climate_attributes(pr: np.ndarray, tas: np.ndarray, pet: np.ndarray, dt_days: float) -> dict[str, float]:
    """The CAMELS climate indices with one definition (Addor et al. 2017).

    Frequencies are days per year and durations are days, the units of the
    training table, whose means are 20.8 and 252 days per year. Those indices
    and the snow fraction are defined over days, so a sub-daily record is
    grouped into whole days here. This is the one place rows are grouped, and
    nothing the network reads per row comes from it; applied to the rows
    instead, the same weather would describe a different catchment at each step.
    """
    rows_per_day = max(1, int(round(1.0 / dt_days)))
    n_days = len(pr) // rows_per_day
    if n_days == 0:  # less than a day of record: the rows are all there is
        rows_per_day, n_days = 1, len(pr)
    keep = n_days * rows_per_day
    daily_pr = pr[:keep].reshape(n_days, rows_per_day).mean(axis=1)
    daily_tas = tas[:keep].reshape(n_days, rows_per_day).mean(axis=1)
    p_mean, pet_mean = float(pr.mean()), float(pet.mean())
    high = daily_pr >= 5.0 * p_mean
    low = daily_pr < 1.0
    total = float(daily_pr.sum())
    return {
        "p_mean": p_mean,
        "pet_mean": pet_mean,
        "aridity": pet_mean / max(p_mean, 1e-9),
        "frac_snow": float(daily_pr[daily_tas < 0.0].sum() / total) if total > 0 else 0.0,
        "high_prec_freq": float(high.mean()) * DAYS_PER_YEAR,
        "high_prec_dur": _mean_run_length(high),
        "low_prec_freq": float(low.mean()) * DAYS_PER_YEAR,
        "low_prec_dur": _mean_run_length(low),
    }


def load_model():
    """As `main.evaluate` builds it, from the run's own configuration."""
    cfg = json.loads((RUN_DIR / "cfg.json").read_text())
    model = Model(
        input_size_dyn=5,
        input_size_stat=27,
        hidden_size=cfg["hidden_size"],
        dropout=cfg["dropout"],
        concat_static=cfg["concat_static"],
        no_static=cfg["no_static"],
    )
    model.load_state_dict(torch.load(RUN_DIR / WEIGHTS, map_location="cpu", weights_only=True))
    model.eval()
    return cfg, model


@torch.no_grad()
def predict(model: Model, x_d: np.ndarray, x_s: torch.Tensor, seq_length: int) -> np.ndarray:
    """One prediction per row, each from a fresh state, as `evaluate_basin` does."""
    n = len(x_d)
    out = np.full(n, np.nan)
    by_length: dict[int, list[int]] = {}
    for row in range(n):
        by_length.setdefault(min(row + 1, seq_length), []).append(row)
    tensor = torch.from_numpy(x_d.astype(np.float32))
    for length, rows in sorted(by_length.items()):
        for start in range(0, len(rows), BATCH):
            ends = torch.as_tensor(rows[start : start + BATCH])
            index = ends[:, None] - (length - 1) + torch.arange(length)[None, :]
            pred = model(tensor[index], x_s.expand(len(ends), -1))[0]
            out[ends.numpy()] = pred.numpy().flatten()
    return out


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet"):
            row[key] = float(row[key])
    return rows


def simulate(forcing: list[dict], static: dict, seed: int, timestep: str) -> tuple[list[dict], dict]:
    torch.set_num_threads(THREADS)
    torch.manual_seed(seed)  # the network draws nothing in eval mode
    cfg, model = load_model()

    pr = np.array([r["pr"] for r in forcing])
    tas = np.array([r["tas"] for r in forcing])
    pet = np.array([r["pet"] for r in forcing])
    dt_days = TIMESTEP_DAYS[timestep]

    # main.evaluate's own statistics, in the column order of its own table.
    attributes = load_attributes(str(RUN_DIR / "attributes.db"), get_basin_list(), drop_lat_lon=True)
    means, stds = attributes.mean(), attributes.std()

    n_clim = min(len(pr), int(round(CLIMATOLOGY_DAYS / dt_days)))
    derived = climate_attributes(pr[:n_clim], tas[:n_clim], pet[:n_clim], dt_days)
    hold = {g for g in os.environ.get("EALSTM_HOLD", "").split(",") if g}
    if "all" in hold:
        derived = {}
    for group in hold & set(ATTRIBUTE_GROUPS):
        for name in ATTRIBUTE_GROUPS[group]:
            derived.pop(name, None)
    if os.environ.get("EALSTM_AREA") == "static":
        derived["area_gages2"] = float(static["area_km2"])
    values = np.array([derived.get(name, means[name]) for name in attributes.columns])
    x_s = torch.from_numpy(((values - means.values) / stds.values).astype(np.float32))[None, :]

    x_d = normalize_features(mock_inputs(pr, tas, pet, float(means["elev_mean"])), variable="inputs")
    scaled = predict(model, x_d, x_s, int(cfg["seq_length"]))
    depth = rescale_features(scaled[:, None], variable="output").flatten()
    depth[depth < 0] = 0  # as evaluate_basin does

    discharge = depth * 1.0e-3 * float(static["area_km2"]) * 1.0e6 / 86400.0
    rows = [
        {"time": step["time"], "mrro": float(q), "dis": float(d)}
        for step, q, d in zip(forcing, depth, discharge)
    ]
    notes = {
        "checkpoint": Path(cfg["run_dir"]).name,
        "weights": WEIGHTS,
        "prediction": "many-to-one from a fresh state per row, negative values set to zero, as evaluate_basin",
        "seq_length_rows": int(cfg["seq_length"]),
        "climatology_rows": int(n_clim),
        "mock_inputs": MOCK_INPUT_NOTES,
        "derived_static_attributes": sorted(derived),
        "static_attributes_at_training_mean": [a for a in attributes.columns if a not in derived],
    }
    return rows, notes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    started = time.monotonic()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent

    forcing = read_forcing(io_dir / request["input"]["forcing"])
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    rows, notes = simulate(forcing, static, int(request.get("seed", 0)), str(request.get("timestep", "PT1D")))

    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    (io_dir / request["output"]["run"]).write_text(
        json.dumps(
            {
                "status": "ok",
                "model": MODEL,
                "n_steps": len(rows),
                "wall_seconds": round(time.monotonic() - started, 2),
                "not_reported": {"states": "the network has no stores to report", "evspsbl": "the network predicts no evaporation"},
                "notes": notes,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

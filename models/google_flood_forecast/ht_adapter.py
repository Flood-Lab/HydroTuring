#!/usr/bin/env python3
"""HydroTuring adapter for Google's operational flood forecasting model.

The model is the mean-embedding forecast LSTM behind Google Flood Hub
(Gauch et al. 2025, HESS; the earlier state-handoff LSTM is Nearing et al.
2024, Nature), run through the open-source OpenHydroNet package with the
published `google-floodhub-settings-55-epochs` weights. It predicts
streamflow and nothing else, so it is declared as `mrro` plus `dis` and is
scored INCOMPLETE on any budget probe. This adapter exists so that the
discharge it produces is available to probes that can use it, and so that
the model is run as it is operated rather than approximated.

How the probe's forcing is mapped onto the model's inputs
---------------------------------------------------------
The model was trained on four weather products (HRES, GraphCast, IMERG, CPC)
and 84 Caravan/HydroATLAS catchment attributes. The probe hands it daily
precipitation, air temperature and potential ET for one lumped catchment.

* Precipitation and temperature are given to every product that carries
  them, each standardised with that product's own training statistics. All
  products are in mm/day and degC, as the forcing is.
* HRES also needs net solar and thermal radiation and surface pressure,
  which the probe cannot supply. The whole HRES product is marked missing
  (NaN), which is the model's own documented path for an unavailable input
  product: product embeddings are combined with a NaN-aware mean. Nothing
  is invented to fill the gap.
* The climate attributes that Caravan derives from forcing are derived here
  from the forcing the model is given, using Caravan's definitions. Every
  other attribute (land cover, terrain, soils, human footprint, ...) has no
  counterpart in a synthetic lumped catchment and is set to its training
  mean, i.e. zero after standardisation.

How the record is simulated
---------------------------
The operational model issues a forecast every day from a 365-day hindcast
window. The row reported for each day is the day-0 member of that day's
forecast: the prediction for the issue day itself, given all forcing up to
and including it. Days with fewer than 365 days of history behind them use
the history that exists. Forecast lead days are not needed for a nowcast
and are not computed. The point prediction is the median of the model's own
CMAL mixture samples, as in the package's tester, seeded from the request so
that a case reproduces exactly.

The catchment area from the static attributes converts the model's
streamflow depth (mm/day, reported as `mrro`) into discharge (`dis`, m3/s).
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

MODEL = {"name": "google_flood_forecast", "version": "0.1.0-828dfc5"}
COLUMNS = ["time", "mrro", "dis"]

RUN_DIR = Path(os.environ.get("GFF_RUN_DIR", "/model/run"))
WEIGHTS = "model_epoch055.pt"
THREADS = int(os.environ.get("GFF_THREADS", "8"))
BATCH = 64

# Model input feature -> forcing column. None means the probe cannot supply
# it, and the product it belongs to is then treated as missing in full.
FEATURE_SOURCE = {
    "hres_total_precipitation": "pr",
    "hres_temperature_2m": "tas",
    "hres_surface_net_solar_radiation": None,
    "hres_surface_net_thermal_radiation": None,
    "hres_surface_pressure": None,
    "graphcast_total_precipitation": "pr",
    "graphcast_temperature_2m": "tas",
    "imerg_precipitation": "pr",
    "cpc_precipitation": "pr",
}

TARGET = "streamflow"  # Caravan streamflow, mm/day
DAYS_PER_YEAR = 365.25


# --- catchment attributes ---------------------------------------------------


def _mean_run_length(mask: np.ndarray) -> float:
    """Average length of the runs of True in a boolean series."""
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


def climate_attributes(pr: np.ndarray, tas: np.ndarray, pet: np.ndarray,
                       months: np.ndarray) -> dict[str, float]:
    """The attributes Caravan computes from a catchment's own forcing.

    Definitions follow Caravan (Kratzert et al. 2023), which follows CAMELS
    (Addor et al. 2017) for the precipitation indices and Knoben et al.
    (2018) for the moisture index and its seasonality. Frequencies are
    fractions of days, durations are in days. The four HydroATLAS climate
    fields that are plain annual totals of the same quantities are filled
    from the forcing as well.
    """
    eps = 1e-9
    p_mean = float(pr.mean())
    pet_mean = float(pet.mean())
    total_pr = float(pr.sum())

    monthly_p = np.array([pr[months == m].mean() for m in range(1, 13) if (months == m).any()])
    monthly_e = np.array([pet[months == m].mean() for m in range(1, 13) if (months == m).any()])
    moisture = np.where(
        monthly_p > monthly_e,
        1.0 - monthly_e / np.maximum(monthly_p, eps),
        np.where(monthly_p < monthly_e, monthly_p / np.maximum(monthly_e, eps) - 1.0, 0.0),
    )

    high = pr >= 5.0 * p_mean
    low = pr < 1.0
    return {
        "p_mean": p_mean,
        "pet_mean_ERA5_LAND": pet_mean,
        "aridity_ERA5_LAND": pet_mean / max(p_mean, eps),
        "frac_snow": float(pr[tas < 0.0].sum() / total_pr) if total_pr > 0 else 0.0,
        "moisture_index_ERA5_LAND": float(moisture.mean()),
        "seasonality_ERA5_LAND": float(moisture.max() - moisture.min()),
        "high_prec_freq": float(high.mean()),
        "high_prec_dur": _mean_run_length(high),
        "low_prec_freq": float(low.mean()),
        "low_prec_dur": _mean_run_length(low),
        # HydroATLAS: annual precipitation and PET in mm, aridity index as
        # 100 * P / PET, annual mean temperature in tenths of a degree.
        "pre_mm_syr": p_mean * DAYS_PER_YEAR,
        "pet_mm_syr": pet_mean * DAYS_PER_YEAR,
        "ari_ix_sav": 100.0 * p_mean / max(pet_mean, eps),
        "tmp_dc_syr": 10.0 * float(tas.mean()),
    }


# --- the model --------------------------------------------------------------


def load_model():
    """Build the network from its own run configuration and load the weights."""
    from ruamel.yaml import YAML
    from googlehydrology.modelzoo.mean_embedding_forecast_lstm import MeanEmbeddingForecastLSTM
    from googlehydrology.utils.config import Config

    raw = YAML(typ="safe").load((RUN_DIR / "config.yml").read_text())
    # The training paths in the shipped configuration point at the machine
    # the model was trained on. Inference reads nothing but the scaler, which
    # lives beside the weights.
    raw["run_dir"] = str(RUN_DIR)
    raw["device"] = "cpu"
    cfg = Config(raw)

    model = MeanEmbeddingForecastLSTM(cfg)
    state = torch.load(RUN_DIR / WEIGHTS, map_location="cpu", weights_only=True)
    state = {key.removeprefix("_orig_mod."): value for key, value in state.items()}
    model.load_state_dict(state)
    model.eval()
    return cfg, model


def scaler_stats(model, names: list[str]) -> tuple[dict[str, float], dict[str, float]]:
    """Centre and scale of each feature, from the scaler the model shipped with."""
    table = model._scaler.scaler
    center = {n: float(table[n].sel(parameter="center")) for n in names}
    scale = {n: float(table[n].sel(parameter="scale")) for n in names}
    return center, scale


def _flatten(groups: dict[str, list[str]] | list) -> list[str]:
    if isinstance(groups, dict):
        return [f for features in groups.values() for f in features]
    return list(groups)


def dynamic_inputs(cfg, forcing: dict[str, np.ndarray], center, scale) -> tuple[dict[str, np.ndarray], list[str]]:
    """Standardised input series per model feature; NaN for a missing product."""
    products: dict[str, list[str]] = {}
    for group in (cfg.hindcast_inputs, cfg.forecast_inputs):
        for product, features in group.items():
            products.setdefault(product, [])
            products[product] += [f for f in features if f not in products[product]]

    n = len(next(iter(forcing.values())))
    series: dict[str, np.ndarray] = {}
    missing: list[str] = []
    for product, features in products.items():
        available = all(FEATURE_SOURCE.get(f) is not None for f in features)
        if not available:
            missing.append(product)
        for f in features:
            if available:
                values = (forcing[FEATURE_SOURCE[f]] - center[f]) / scale[f]
            else:
                values = np.full(n, np.nan)
            series[f] = values.astype(np.float32)
    return series, missing


@torch.no_grad()
def nowcast(model, cfg, x_s: torch.Tensor, series: dict[str, np.ndarray],
            seed: int, seq_length: int) -> np.ndarray:
    """Day-0 prediction for every day, in the model's standardised units.

    Each day is its own forecast issue: a window of up to `seq_length` days
    of history ending on that day, processed from a fresh state exactly as
    the operational model does. Days are batched by history length.
    """
    hindcast = _flatten(cfg.hindcast_inputs)
    forecast = _flatten(cfg.forecast_inputs)
    n = len(next(iter(series.values())))

    # Nowcast only: the sequence ends on the issue day and the head's last
    # output is the day-0 prediction, so no forecast lead days are built.
    model.lead_time = 0
    torch.manual_seed(seed)

    by_length: dict[int, list[int]] = {}
    for day in range(n):
        by_length.setdefault(min(day + 1, seq_length), []).append(day)

    out = np.full(n, np.nan)
    for length, days in sorted(by_length.items()):
        model.seq_length = length
        for start in range(0, len(days), BATCH):
            issue = np.asarray(days[start : start + BATCH])
            rows = issue[:, None] - (length - 1) + np.arange(length)[None, :]

            def window(feature: str) -> torch.Tensor:
                return torch.from_numpy(series[feature][rows][..., None])

            data = {
                "x_s": x_s.expand(len(issue), -1),
                "x_d_hindcast": {f: window(f) for f in hindcast},
                "x_d_forecast": {f: window(f) for f in forecast},
                # The sampler reads the batch size off the target tensor.
                "y": torch.zeros(len(issue), 1, 1),
            }
            samples = model.sample(data, cfg.n_samples)["y_hat"]  # [batch, steps, target, sample]
            out[issue] = samples[:, -1, 0, :].median(dim=-1).values.numpy()
    return out


# --- contract plumbing ------------------------------------------------------


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet"):
            row[key] = float(row[key])
    return rows


def simulate(forcing: list[dict], static: dict, seed: int) -> tuple[list[dict], dict]:
    torch.set_num_threads(THREADS)
    cfg, model = load_model()

    pr = np.array([r["pr"] for r in forcing])
    tas = np.array([r["tas"] for r in forcing])
    pet = np.array([r["pet"] for r in forcing])
    months = np.array([int(str(r["time"])[5:7]) for r in forcing])

    names = list(cfg.static_attributes) + list(FEATURE_SOURCE) + [TARGET]
    center, scale = scaler_stats(model, names)

    derived = climate_attributes(pr, tas, pet, months)
    x_s = torch.tensor(
        [[(derived.get(a, center[a]) - center[a]) / scale[a] for a in cfg.static_attributes]],
        dtype=torch.float32,
    )
    series, missing_products = dynamic_inputs(cfg, {"pr": pr, "tas": tas}, center, scale)

    scaled = nowcast(model, cfg, x_s, series, seed, int(cfg.seq_length))
    depth = np.maximum(scaled * scale[TARGET] + center[TARGET], 0.0)  # mm/day
    area_m2 = float(static["area_km2"]) * 1.0e6
    discharge = depth * 1.0e-3 * area_m2 / 86400.0  # m3/s

    rows = [
        {"time": step["time"], "mrro": float(q), "dis": float(d)}
        for step, q, d in zip(forcing, depth, discharge)
    ]
    notes = {
        "weights": WEIGHTS,
        "prediction": "day-0 member of each day's forecast, median of CMAL samples",
        "n_samples": int(cfg.n_samples),
        "hindcast_days": int(cfg.seq_length),
        "missing_products": missing_products,
        "derived_static_attributes": sorted(a for a in derived if a in cfg.static_attributes),
        "static_attributes_at_training_mean": [
            a for a in cfg.static_attributes if a not in derived
        ],
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

    rows, notes = simulate(forcing, static, int(request.get("seed", 0)))

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
                "notes": notes,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

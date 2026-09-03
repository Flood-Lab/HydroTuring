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
  which the probe does not generate. For the synthetic mass-balance test
  these are mocked from the forcing with textbook (FAO-56) relations:
  extraterrestrial radiation at the catchment's latitude, attenuated on wet
  days, an albedo of 0.23, net longwave from air temperature at 70 percent
  humidity, and surface pressure from the elevation the static attributes
  assume. They are labelled as mock inputs in run.json. Set
  GFF_MOCK_HRES=0 to mark the HRES product missing instead, which is the
  model's own path for an unavailable product (a NaN-aware mean over
  product embeddings).
* The climate attributes that Caravan derives from forcing are derived here
  from the forcing the model is given, using Caravan's definitions. Every
  other attribute (land cover, terrain, soils, human footprint, ...) has no
  counterpart in a synthetic lumped catchment and is set to its training
  mean, i.e. zero after standardisation: the mock catchment is an average
  Caravan basin at the probe's latitude.

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
import datetime as dt
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
MOCK_HRES = os.environ.get("GFF_MOCK_HRES", "1") != "0"
BATCH = 64

# Model input feature -> forcing column, or a mock series built from the
# forcing. A product with any input that is neither is treated as missing
# in full.
FEATURE_SOURCE = {
    "hres_total_precipitation": "pr",
    "hres_temperature_2m": "tas",
    "hres_surface_net_solar_radiation": "mock_net_solar_radiation",
    "hres_surface_net_thermal_radiation": "mock_net_thermal_radiation",
    "hres_surface_pressure": "mock_surface_pressure",
    "graphcast_total_precipitation": "pr",
    "graphcast_temperature_2m": "tas",
    "imerg_precipitation": "pr",
    "cpc_precipitation": "pr",
}

TARGET = "streamflow"  # Caravan streamflow, mm/day
DAYS_PER_YEAR = 365.25

# Mock HRES fields (units as in the Caravan-MultiMet training data: W/m2
# for radiation, kPa for pressure). Textbook FAO-56 daily relations.
ALBEDO = 0.23
CLOUD_FACTOR_DRY, CLOUD_FACTOR_WET = 0.70, 0.45  # fraction of extraterrestrial radiation
WET_DAY_MM = 1.0
RELATIVE_HUMIDITY = 0.70
SOLAR_CONSTANT = 0.0820  # MJ m-2 min-1
STEFAN_BOLTZMANN = 4.903e-9  # MJ K-4 m-2 day-1
MJ_PER_DAY_TO_W = 1.0e6 / 86400.0
MOCK_INPUT_NOTES = {
    "mock_net_solar_radiation": (
        "FAO-56 extraterrestrial radiation at the catchment latitude, times "
        f"{CLOUD_FACTOR_DRY} on dry and {CLOUD_FACTOR_WET} on wet days, "
        f"times (1 - albedo {ALBEDO}); W/m2"
    ),
    "mock_net_thermal_radiation": (
        "FAO-56 net longwave from air temperature at "
        f"{int(RELATIVE_HUMIDITY * 100)} percent relative humidity and the "
        "same cloudiness, reported as a negative (outgoing) flux; W/m2"
    ),
    "mock_surface_pressure": (
        "FAO-56 standard atmosphere at the elevation the static attributes "
        "assume (the training-mean ele_mt_sav); kPa, constant"
    ),
}


def mock_hres_fields(pr: np.ndarray, tas: np.ndarray, doy: np.ndarray,
                     latitude_deg: float, elevation_m: float) -> dict[str, np.ndarray]:
    """Plausible daily radiation and pressure for a catchment the probe only
    describes by its weather, latitude and an assumed elevation."""
    phi = np.radians(latitude_deg)
    j = doy.astype(float)
    dr = 1.0 + 0.033 * np.cos(2.0 * np.pi * j / 365.0)
    delta = 0.409 * np.sin(2.0 * np.pi * j / 365.0 - 1.39)
    omega = np.arccos(np.clip(-np.tan(phi) * np.tan(delta), -1.0, 1.0))
    ra = (24.0 * 60.0 / np.pi) * SOLAR_CONSTANT * dr * (
        omega * np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.sin(omega)
    )  # MJ m-2 day-1

    cloud = np.where(pr >= WET_DAY_MM, CLOUD_FACTOR_WET, CLOUD_FACTOR_DRY)
    rs = cloud * ra
    rso = (0.75 + 2.0e-5 * elevation_m) * ra
    clearness = np.clip(rs / np.maximum(rso, 1e-6), 0.3, 1.0)

    net_solar = (1.0 - ALBEDO) * rs
    saturation_kpa = 0.6108 * np.exp(17.27 * tas / (tas + 237.3))
    ea = RELATIVE_HUMIDITY * saturation_kpa
    net_longwave = (
        STEFAN_BOLTZMANN * (tas + 273.16) ** 4
        * (0.34 - 0.14 * np.sqrt(ea))
        * (1.35 * clearness - 0.35)
    )
    pressure = 101.3 * ((293.0 - 0.0065 * elevation_m) / 293.0) ** 5.26

    return {
        "mock_net_solar_radiation": net_solar * MJ_PER_DAY_TO_W,
        "mock_net_thermal_radiation": -net_longwave * MJ_PER_DAY_TO_W,
        "mock_surface_pressure": np.full(len(pr), pressure),
    }


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
        available = all(FEATURE_SOURCE.get(f) in forcing for f in features)
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
    dates = [dt.date.fromisoformat(str(r["time"])[:10]) for r in forcing]
    months = np.array([d.month for d in dates])
    doy = np.array([d.timetuple().tm_yday for d in dates])

    names = list(cfg.static_attributes) + list(FEATURE_SOURCE) + [TARGET]
    center, scale = scaler_stats(model, names)

    derived = climate_attributes(pr, tas, pet, months)
    x_s = torch.tensor(
        [[(derived.get(a, center[a]) - center[a]) / scale[a] for a in cfg.static_attributes]],
        dtype=torch.float32,
    )

    inputs = {"pr": pr, "tas": tas}
    mocked: dict[str, str] = {}
    if MOCK_HRES:
        inputs.update(mock_hres_fields(
            pr, tas, doy,
            latitude_deg=float(static.get("latitude_deg", 40.0)),
            elevation_m=float(center.get("ele_mt_sav", 0.0)),
        ))
        mocked = {
            feature: MOCK_INPUT_NOTES[source]
            for feature, source in FEATURE_SOURCE.items()
            if source in MOCK_INPUT_NOTES
        }
    series, missing_products = dynamic_inputs(cfg, inputs, center, scale)

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
        "mock_inputs": mocked,
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

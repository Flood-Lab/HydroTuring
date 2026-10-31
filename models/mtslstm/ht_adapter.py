#!/usr/bin/env python3
"""HydroTuring adapter for the Multi-Timescale LSTM of Gauch et al. (2021).

The network and its weights are upstream's and unmodified: `MTSLSTM` from
neuralhydrology at commit 9626578, built from the run's own config.yml, with
the published `mtslstm_seed333` epoch-30 checkpoint from Zenodo. The MTS-LSTM
archive holds no scaler, so the statistics are not this run's own: they are the
training scaler of a sister sMTS-LSTM run on the same data and period, exported to
JSON at image build time because the pickled xarray objects no longer open on a
current stack. With it this driver reproduces this checkpoint's published simulation
to single precision (see README.md).

The model predicts streamflow and nothing else, so it is declared as `mrro`
plus `dis`, with `states: []`, and is N/A (INCOMPLETE) on every budget probe.

Which branch answers which step
-------------------------------
The network has a daily branch and an hourly branch, and the step in the
request says which one is reported. Fluxes leave in mm/day at either step;
both heads predict mm/hour, the unit of the training target.

* PT1D. The daily LSTM and its head, on the daily rows. In `MTSLSTM.forward`
  the daily branch never sees the hourly inputs, and splitting it at the
  state-transfer step changes nothing, so this is the model's own daily
  prediction. There is no resampling anywhere on this path.
* PT1H. The hourly LSTM and its head, on the hourly rows, started from the
  state the daily branch hands over at the transfer step. Each day of rows is
  one sample and its last 24 hours are kept, as upstream's tester does.

The daily branch has to read something, and upstream's own data pipeline gives
it the daily mean of the hourly record (`BaseDataset`, `resample(freq).mean()`).
The adapter does the same for that branch alone, and only for days wholly before
the hourly window. This is the one place a model input is aggregated, and it is
flagged because AGENTS.md says not to: the rule is about reaching a step the model
prefers, and here no row is reported from anything but the hourly head fed the
hourly rows. MTSLSTM_DAILY_BRANCH=none skips the aggregation and starts the
hourly branch from the transfer of a zero state instead; README.md gives what
that moves.

Any other step has no branch and is refused rather than fed to the wrong one.

How the forcing is mapped
-------------------------
Eleven NLDAS inputs, from the row's own pr, tas and pet with no calendar
dependence, standardised with the run's scaler, listed in run.json:

* total_precipitation and potential_evaporation are pr and pet, in mm/hour.
* temperature is tas.
* shortwave_radiation inverts Priestley-Taylor from pet, adds FAO-56 net
  longwave and removes albedo, in W/m2 as a mean over the row.
* longwave_radiation is Brutsaert's clear-sky emissivity at the row's vapour
  pressure, in W/m2.
* specific_humidity is from air temperature at 70 percent relative humidity.
* pressure, convective_fraction, potential_energy, wind_u and wind_v have no
  counterpart in the row and are held at the training mean, zero after
  standardisation.

The 26 static attributes: the eight climate indices with an unambiguous CAMELS
definition are derived from the first year of the forcing the model is given,
over days at any step, since that is how they are defined; the other 18 are held
at the training mean.

Sensitivity switches, for the table in README.md and never set in the image:
MTSLSTM_HOLD lists attribute groups to hold at the training mean, or "all";
MTSLSTM_MOCK=mean holds every input but precipitation at the training mean;
MTSLSTM_AREA=static feeds the case's area to area_gages2.
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

RUN_DIR = Path(os.environ.get("MTSLSTM_RUN_DIR", "/model/run"))
THREADS = int(os.environ.get("MTSLSTM_THREADS", "4"))
DAILY_BRANCH = os.environ.get("MTSLSTM_DAILY_BRANCH", "aggregate")

MODEL = {"name": "mtslstm", "version": "0.1.3-9626578-seed333"}
COLUMNS = ["time", "mrro", "dis"]
WEIGHTS = "model_epoch030.pt"
BATCH = 512

TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0}
CLIMATOLOGY_DAYS = 365.0
DAYS_PER_YEAR = 365.25
HOURS_PER_DAY = 24

ATTRIBUTE_GROUPS = {
    "precipitation": ["p_mean", "high_prec_freq", "high_prec_dur", "low_prec_freq", "low_prec_dur"],
    "demand": ["pet_mean", "aridity"],
    "snow": ["frac_snow"],
}

ALBEDO = 0.23
PRIESTLEY_TAYLOR_ALPHA = 1.26
LATENT_HEAT = 2.45  # MJ kg-1
STEFAN_BOLTZMANN = 4.903e-9  # MJ K-4 m-2 day-1
MJ_PER_DAY_TO_W = 1.0e6 / 86400.0
CLEARNESS_DRY, CLEARNESS_WET = 0.70 / 0.75, 0.45 / 0.75
WET_RATE_MM_DAY = 1.0
RELATIVE_HUMIDITY = 0.70
MOCK_INPUT_NOTES = {
    "total_precipitation": "pr / 24, mm/hour",
    "potential_evaporation": "pet / 24, mm/hour",
    "temperature": "tas",
    "shortwave_radiation": (
        f"Priestley-Taylor (alpha {PRIESTLEY_TAYLOR_ALPHA}) inverted from pet, plus FAO-56 net "
        f"longwave, over 1 - albedo ({ALBEDO}); W/m2"
    ),
    "longwave_radiation": f"Brutsaert clear-sky emissivity at {int(RELATIVE_HUMIDITY * 100)} percent humidity; W/m2",
    "specific_humidity": f"air temperature at {int(RELATIVE_HUMIDITY * 100)} percent humidity; kg/kg",
    "pressure, convective_fraction, potential_energy, wind_u, wind_v": "training mean",
}


def physical_inputs(pr, tas, pet, pressure_pa: float) -> dict[str, np.ndarray]:
    """The inputs a row can speak to; anything not returned stays at the mean."""
    saturation = 0.6108 * np.exp(17.27 * tas / (tas + 237.3))
    vapour = RELATIVE_HUMIDITY * saturation
    slope = 4098.0 * saturation / (tas + 237.3) ** 2
    gamma = 0.000665 * pressure_pa / 1000.0
    net_radiation = pet * LATENT_HEAT * (slope + gamma) / (PRIESTLEY_TAYLOR_ALPHA * slope)
    clearness = np.where(pr >= WET_RATE_MM_DAY, CLEARNESS_WET, CLEARNESS_DRY)
    kelvin = tas + 273.16
    net_longwave = STEFAN_BOLTZMANN * kelvin**4 * (0.34 - 0.14 * np.sqrt(vapour)) * (1.35 * clearness - 0.35)
    shortwave = np.maximum(net_radiation + net_longwave, 0.0) / (1.0 - ALBEDO) * MJ_PER_DAY_TO_W
    emissivity = 1.24 * (10.0 * vapour / kelvin) ** (1.0 / 7.0)
    longwave = emissivity * STEFAN_BOLTZMANN * kelvin**4 * MJ_PER_DAY_TO_W
    humidity = 0.622 * vapour / (pressure_pa / 1000.0 - 0.378 * vapour)
    return {
        "total_precipitation": pr / HOURS_PER_DAY,
        "potential_evaporation": pet / HOURS_PER_DAY,
        "temperature": tas,
        "shortwave_radiation": shortwave,
        "longwave_radiation": longwave,
        "specific_humidity": humidity,
    }


def _mean_run_length(mask: np.ndarray) -> float:
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


def climate_attributes(pr, tas, pet, dt_days: float) -> dict[str, float]:
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
    """Built from the run's own configuration, as neuralhydrology's tester does."""
    from neuralhydrology.modelzoo.mtslstm import MTSLSTM
    from neuralhydrology.utils.config import Config

    cfg = Config(RUN_DIR / "config.yml")
    model = MTSLSTM(cfg=cfg)
    model.load_state_dict(torch.load(RUN_DIR / WEIGHTS, map_location="cpu", weights_only=True))
    model.eval()
    return cfg, model


def _windows(n_rows: int, seq_length: int, ends: list[int]):
    """Row indices of windows ending on `ends`, grouped by length."""
    by_length: dict[int, list[int]] = {}
    for end in ends:
        by_length.setdefault(min(end + 1, seq_length), []).append(end)
    for length, group in sorted(by_length.items()):
        for start in range(0, len(group), BATCH):
            chunk = torch.as_tensor(group[start : start + BATCH])
            yield chunk, chunk[:, None] - (length - 1) + torch.arange(length)[None, :]


def _with_static(x_d: torch.Tensor, x_s: torch.Tensor) -> torch.Tensor:
    """[seq, batch, features + attributes], as `MTSLSTM._prepare_inputs` builds it."""
    return torch.cat([x_d, x_s[None, None, :].repeat(x_d.shape[0], x_d.shape[1], 1)], dim=-1)


@torch.no_grad()
def predict_daily(model, x: torch.Tensor, x_s: torch.Tensor, seq_length: int) -> np.ndarray:
    """The daily branch alone: a fresh state and up to a year of rows per row."""
    n = x.shape[0]
    out = np.full(n, np.nan)
    for ends, index in _windows(n, seq_length, list(range(n))):
        inputs = _with_static(x[index].transpose(0, 1), x_s)
        states = (inputs.new_zeros((1, inputs.shape[1], model._hidden_size["1D"])),) * 2
        lstm_output, _ = model.lstms["1D"](inputs, states)
        head = model.heads["1D"](model.dropout(lstm_output.transpose(0, 1)))["y_hat"]
        out[ends.numpy()] = head[:, -1, 0].numpy()
    return out


@torch.no_grad()
def predict_hourly(model, x: torch.Tensor, x_s: torch.Tensor, seq_daily: int, seq_hourly: int) -> np.ndarray:
    """The hourly branch per day of rows, started from the daily branch's handover.

    The window of hourly rows ends on the day's last row. The daily branch
    reads the daily means of whole days before that window, and its final state
    goes through the transfer layers, as in `MTSLSTM.forward`.
    """
    n = x.shape[0]
    out = np.full(n, np.nan)
    hidden = model._hidden_size["1D"]
    for day_end in range(HOURS_PER_DAY - 1, n + HOURS_PER_DAY - 1, HOURS_PER_DAY):
        end = min(day_end, n - 1)
        # From the nominal day's end, or a record ending mid-day would leave rows that
        # neither branch reads, and extending the record would change earlier history.
        first = max(0, day_end - seq_hourly + 1)
        hourly = _with_static(x[first : end + 1][:, None, :], x_s)

        h = c = hourly.new_zeros((1, 1, hidden))
        n_days = first // HOURS_PER_DAY
        if DAILY_BRANCH == "aggregate" and n_days > 0:
            days = x[: n_days * HOURS_PER_DAY].reshape(n_days, HOURS_PER_DAY, -1).mean(dim=1)
            # The upstream window is seq_daily days ending on the sample's day, of
            # which the last seq_hourly / 24 are past the handover.
            keep = seq_daily - seq_hourly // HOURS_PER_DAY
            daily = _with_static(days[-keep:][:, None, :], x_s)
            _, (h, c) = model.lstms["1D"](daily, (h, c))
        h = model.transfer_fcs["h_1D"](h)
        c = model.transfer_fcs["c_1D"](c)

        lstm_output, _ = model.lstms["1H"](hourly, (h, c))
        head = model.heads["1H"](model.dropout(lstm_output.transpose(0, 1)))["y_hat"]
        start = day_end - HOURS_PER_DAY + 1
        out[start : end + 1] = head[0, -(end + 1 - start) :, 0].numpy()
    return out


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet"):
            row[key] = float(row[key])
    return rows


def simulate(forcing: list[dict], static: dict, seed: int, timestep: str) -> tuple[list[dict], dict]:
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"MTS-LSTM has a daily and an hourly branch, and none for {timestep}")
    torch.set_num_threads(THREADS)
    torch.manual_seed(seed)  # the network draws nothing in eval mode
    cfg, model = load_model()
    scaler = json.loads((RUN_DIR / "scaler.json").read_text())
    means, stds = scaler["xarray_means"], scaler["xarray_stds"]
    attr_means, attr_stds = scaler["attribute_means"], scaler["attribute_stds"]

    pr = np.array([r["pr"] for r in forcing])
    tas = np.array([r["tas"] for r in forcing])
    pet = np.array([r["pet"] for r in forcing])
    dt_days = TIMESTEP_DAYS[timestep]

    physical = physical_inputs(pr, tas, pet, means["pressure"])
    if os.environ.get("MTSLSTM_MOCK") == "mean":
        physical = {"total_precipitation": physical["total_precipitation"]}
    columns = [
        (physical[name] - means[name]) / stds[name] if name in physical else np.zeros(len(pr))
        for name in cfg.dynamic_inputs
    ]
    x = torch.from_numpy(np.column_stack(columns).astype(np.float32))

    # A record shorter than a year has no climatology to give, and all of it may be
    # scored, so the attributes stay at the training means rather than look ahead.
    n_clim = int(round(CLIMATOLOGY_DAYS / dt_days))
    derived = climate_attributes(pr[:n_clim], tas[:n_clim], pet[:n_clim], dt_days) if len(pr) >= n_clim else {}
    hold = {g for g in os.environ.get("MTSLSTM_HOLD", "").split(",") if g}
    if "all" in hold:
        derived = {}
    for group in hold & set(ATTRIBUTE_GROUPS):
        for name in ATTRIBUTE_GROUPS[group]:
            derived.pop(name, None)
    if os.environ.get("MTSLSTM_AREA") == "static":
        derived["area_gages2"] = float(static["area_km2"])
    # In the order the scaler holds them: the order the network was trained on.
    names = list(attr_means)
    x_s = torch.tensor(
        [(derived.get(a, attr_means[a]) - attr_means[a]) / attr_stds[a] for a in names],
        dtype=torch.float32,
    )

    if timestep == "PT1D":
        scaled = predict_daily(model, x, x_s, cfg.seq_length["1D"])
        branch = "1D"
    else:
        scaled = predict_hourly(model, x, x_s, cfg.seq_length["1D"], cfg.seq_length["1H"])
        branch = "1H"
    target = cfg.target_variables[0]
    per_hour = scaled * stds[target] + means[target]
    depth = np.maximum(per_hour, 0.0) * HOURS_PER_DAY  # clip_targets_to_zero, mm/hour to mm/day
    discharge = depth * 1.0e-3 * float(static["area_km2"]) * 1.0e6 / 86400.0

    rows = [
        {"time": step["time"], "mrro": float(q), "dis": float(d)}
        for step, q, d in zip(forcing, depth, discharge)
    ]
    notes = {
        "checkpoint": f"{cfg.experiment_name} {WEIGHTS}",
        "branch_reported": branch,
        "daily_branch_input_at_hourly_step": (
            ("daily means of hourly rows before the hourly window, as upstream's pipeline; "
             "no row is reported from them") if DAILY_BRANCH == "aggregate" else "none"
        ) if branch == "1H" else "the daily rows",
        "climatology_rows": int(n_clim),
        "mock_inputs": MOCK_INPUT_NOTES,
        "derived_static_attributes": sorted(derived),
        "static_attributes_at_training_mean": [a for a in names if a not in derived],
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

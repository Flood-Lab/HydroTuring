#!/usr/bin/env python3
"""HydroTuring adapter for δHBV 2.0 (mhpi/dhbv2).

δHBV 2.0 (Song et al. 2025, WRR) is a differentiable HBV: an LSTM writes
three HBV parameters day by day from the forcing and the catchment
attributes, an MLP writes the rest from the attributes alone, and the HBV
bucket model runs with them. The package this is built from is the NextGen
module mhpi/dhbv2, whose BMI exposes one output, routed runoff. This adapter
does not go through the BMI. It builds the same δMG model the BMI builds,
from the same configuration and weights, and reads the whole output
dictionary, because the point of submitting this model is that it carries
explicit stores and an explicit evaporation and can be scored on closure.

What is reported
----------------
Fluxes, as rates in mm per day: `pr` echoed as given, `evspsbl` the HBV
actual evaporation, `snm` the rain passing through the snow module plus
liquid water released from it (`RAIN + tosoil`), `mrro` the routed
streamflow, `dis` that streamflow over the catchment area in m3/s.
States, absolute, in mm, each averaged
over the HBV components the model runs in parallel:

* `snw`     snowpack plus the liquid water held in it (SNOWPACK + MELTWATER)
* `mrso`    the soil moisture box (SM)
* `gw`      the two groundwater boxes (SUZ + SLZ)
* `channel` runoff that has been generated but not yet released by the unit
            hydrograph: cumulative unrouted minus cumulative routed flow
* `canopy`  zero. HBV has no interception store, so its canopy storage is
            identically zero. This is a statement about the model's
            structure, not a value invented to fill a column.
* `gwex`    the regional groundwater exchange HBV 2.0 adds to its lower
            store every day, declared as the flux it is (positive into the
            catchment, mm per day, mean over components). It is computed
            from the parameters the MLP wrote, parRT and parAC, and the
            upstream area, exactly as the model applies it; where the term
            would drain the store below zero the model clamps it and the
            declared flux is clamped the same way. Declared, the budget can
            close; the residual the closure probe found before this was
            this term, hidden.

How the probe's forcing is mapped onto the model's inputs
---------------------------------------------------------
The model takes precipitation, air temperature and potential evaporation,
which is exactly what the probe generates, so nothing has to be mocked.
The probe's `pet` is handed to the model directly; the module's own
Hargreaves routine is not used. The adapter does not resample: the rows go
in at the step they came at, which is what a resolution probe measures.

Units at a step other than a day
--------------------------------
The networks were trained to read forcing in mm per day and to write HBV
parameters in daily units: recession coefficients per day, a percolation
cap and a degree-day factor and a groundwater exchange in mm per day, a
unit hydrograph in days. Run at an hourly step, the physics would drain
its stores twenty-four times too fast unless those parameters are put in
the units of the step. So the networks are fed rates, as they were
trained, the HBV core is fed depths per row (rate times step length), and
every parameter with a time in its units is rescaled to the step before
the core uses it: a per-day fraction k becomes 1 - (1 - k)^dt, a per-day
amount becomes amount * dt, the unit hydrograph keeps its shape in days by
scaling its time constant and its length. At a daily step the scaling is
the identity. What remains step-dependent after that is the LSTM's own
recurrence, whose memory was learned with one step meaning one day, and
that is the part of the model the resolution probe then measures.

Of the 28 catchment attributes the parameterisation network reads, the
ones with an unambiguous definition are derived from the first year of the
forcing the model is given (or all of it when the record is shorter), so
that they are a climatology the model has already seen and nothing later
in the record can reach back through them: mean annual precipitation and
potential evaporation, their ratio, mean temperature, and the share of
precipitation falling below the snow threshold. Upstream area is the
catchment area. Glaciers and permafrost are zero. Every other attribute is
set to its training mean, i.e. zero after standardisation: soil texture,
porosity, permeability, slope, elevation, NDVI and free water because a
synthetic lumped catchment has no counterpart for them, and the two
seasonality indices and the snow-cover fraction because the definitions
the training set used could not be established, and a guessed definition
that lands outside the training range moves the parameters the network
writes more than any input the probe controls (README.md has the numbers).

Structure read from the checkpoint
----------------------------------
The configuration shipped with the weights declares four parallel HBV
components and no routing. The weights disagree: their output layers are
sized for three components and a two-parameter unit hydrograph, and the
model cannot be loaded as configured. The adapter reads the number of
components and the presence of routing from the checkpoint's own shapes
and records what it found in run.json. Nothing else is overridden.

The model is deterministic: dropout is off in evaluation and the
parameter-dropout rate in the shipped configuration is zero. The request
seed is recorded and otherwise unused.
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
import yaml

MODEL = {"name": "dhbv2", "version": "0.5.4-hbv2ep100.4"}
COLUMNS = ["time", "pr", "evspsbl", "snm", "mrro", "dis", "gwex", "mrso", "snw", "canopy", "gw", "channel"]

MODEL_DIR = Path(os.environ.get("DHBV_MODEL_DIR", "/model/dhbv_2"))
THREADS = int(os.environ.get("DHBV_THREADS", "2"))
PHY_MODEL = "Hbv_2"
CLIMATOLOGY_DAYS = 365.0  # the attributes come from the first year the model is shown
DAYS_PER_YEAR = 365.25
TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0, "PT15M": 1.0 / 96.0, "PT5M": 1.0 / 288.0, "PT1M": 1.0 / 1440.0}
EPS = 1e-6  # the module's own normalisation epsilon


# --- the catchment ------------------------------------------------------------


def climate_attributes(pr: np.ndarray, tas: np.ndarray, pet: np.ndarray,
                       snow_threshold: float) -> dict[str, float]:
    """The attributes a catchment's own forcing determines. Rates are mm/day."""
    total_pr = float(pr.sum())
    cold = tas < snow_threshold
    return {
        "meanP": float(pr.mean()) * DAYS_PER_YEAR,
        "ETPOT_Hargr": float(pet.mean()) * DAYS_PER_YEAR,
        "aridity": float(pet.mean()) / max(float(pr.mean()), 1e-9),
        "meanTa": float(tas.mean()),
        "snowfall_fraction": float(pr[cold].sum() / total_pr) if total_pr > 0 else 0.0,
    }


# --- the model ----------------------------------------------------------------


class StepScaledHbv2:
    """Mixin swapped onto the loaded HBV core: puts the daily parameters the
    networks wrote into the units of the case's step before the physics
    runs. See the module docstring. Installed by `load_model` with
    `__class__` assignment so the trained weights and configuration are
    untouched; `dt_days` is set per case."""

    dt_days = 1.0
    ROUTING_LENGTH_DAYS = 15  # hydrodl2's lenF at the daily step

    PER_DAY_FRACTION = ("parK0", "parK1", "parK2", "parC")
    PER_DAY_AMOUNT = ("parPERC", "parCFMAX", "parRT")

    @classmethod
    def _to_step(cls, name: str, value, dt: float):
        if name in cls.PER_DAY_FRACTION:
            return 1.0 - (1.0 - value) ** dt
        if name in cls.PER_DAY_AMOUNT:
            return value * dt
        return value

    def _PBM(self, forcing, Ac, Elevation, states, phy_dy_params_dict, phy_static_params_dict):
        dt = float(self.dt_days)
        if dt != 1.0:
            phy_dy_params_dict = {k: self._to_step(k, v, dt) for k, v in phy_dy_params_dict.items()}
            phy_static_params_dict = {k: self._to_step(k, v, dt) for k, v in phy_static_params_dict.items()}
            if self.routing:
                # uh_gamma uses theta = relu(route_b) + 0.5 in steps over lenF
                # steps; keep the hydrograph's shape in days.
                theta_days = torch.relu(self.routing_param_dict["route_b"]) + 0.5
                self.routing_param_dict["route_b"] = theta_days / dt - 0.5
                self.lenF = int(round(self.ROUTING_LENGTH_DAYS / dt))
        return super()._PBM(forcing, Ac, Elevation, states, phy_dy_params_dict, phy_static_params_dict)


def checkpoint_structure(state: dict, n_phy: int, n_dynamic: int) -> tuple[int, bool]:
    """Components and routing as the trained weights actually have them.

    The LSTM's output layer holds one value per dynamic parameter per
    component; the MLP's holds the static parameters per component plus the
    unit hydrograph's two parameters when the model was trained with routing.
    """
    n_lstm_out = int(state["nn_model.lstminv.linear_out.bias"].shape[0])
    n_mlp_out = int(state["nn_model.ann.h2o.bias"].shape[0])
    if n_lstm_out % n_dynamic:
        raise RuntimeError(
            f"checkpoint LSTM output {n_lstm_out} is not a multiple of "
            f"{n_dynamic} dynamic parameters"
        )
    nmul = n_lstm_out // n_dynamic
    n_route = n_mlp_out - (n_phy - n_dynamic) * nmul
    if n_route not in (0, 2):
        raise RuntimeError(
            f"checkpoint MLP output {n_mlp_out} fits neither {nmul} components "
            "without routing nor with a two-parameter unit hydrograph"
        )
    return nmul, n_route == 2


def load_model() -> tuple[dict, torch.nn.Module, dict, dict]:
    """Build the δMG model as the module's BMI does, from its own files."""
    from dmg import ModelHandler
    from dmg.core import Dates

    config = yaml.safe_load((MODEL_DIR / "config.yaml").read_text())
    config["device"] = "cpu"
    config["model_dir"] = str(MODEL_DIR)
    sim = Dates(config["sim"], config["model"]["rho"])
    config["sim_time"] = [sim.start_time, sim.end_time]
    for key in ("plot_dir", "sim_dir", "log_dir"):
        config[key] = ""
    config["model"]["phy"]["nearzero"] = float(config["model"]["phy"]["nearzero"])
    # One batched pass over the whole record; nothing is carried between calls.
    config["cache_states"] = False
    config["model"]["phy"]["cache_states"] = False
    config["model"]["nn"]["cache_states"] = False

    epoch = config["test"]["test_epoch"]
    weights = MODEL_DIR / f"{PHY_MODEL.lower()}_ep{epoch}.pt"
    state = torch.load(weights, map_location="cpu", weights_only=True)
    from hydrodl2.models.hbv.hbv_2 import Hbv_2
    n_phy = len(Hbv_2().parameter_bounds)
    n_dynamic = len(config["model"]["phy"]["dynamic_params"][PHY_MODEL])
    nmul, routing = checkpoint_structure(state, n_phy, n_dynamic)
    declared = {"nmul": config["model"]["phy"]["nmul"], "routing": config["model"]["phy"]["routing"]}
    config["model"]["phy"]["nmul"] = nmul
    config["model"]["phy"]["routing"] = routing

    handler = ModelHandler(config, device="cpu")
    handler.eval()
    model = handler.model_dict[PHY_MODEL]
    core = model.phy_model
    core.__class__ = type("StepScaledHbv2", (StepScaledHbv2, core.__class__), {})
    norm = json.loads((MODEL_DIR / "normalization_statistics.json").read_text())
    structure = {
        "weights": weights.name,
        "components_nmul": nmul,
        "unit_hydrograph_routing": routing,
        "as_declared_in_config": declared,
    }
    return config, model, norm, structure


def regional_exchange(core, static_params, slz: np.ndarray, percolation: np.ndarray,
                      upstream_area_km2: float, dt: float) -> np.ndarray:
    """The lateral groundwater term HBV 2.0 adds to its lower store each step.

    Reconstructed from the parameters the network wrote, as the core applies
    them: parRT (mm/day) times a factor in [-1, 1] that depends on the
    upstream area against parAC, per component; scaled to the step as the
    core's own parameters are. The core clamps the store at zero, so a
    draining term can remove at most what the store holds after that step's
    percolation; the same clamp is applied here using the reported store of
    the previous step plus the mean percolation, which is exact when the
    term is a source and approximate when it is a sink. Returned as a depth
    per step, mean over components.
    """
    names = [p for p in core.phy_param_names if p not in core.dynamic_params]
    nmul = core.nmul
    raw = static_params[0, : len(names) * nmul].view(len(names), nmul).double().numpy()
    b = core.parameter_bounds
    def descale(name):
        lo, hi = b[name]
        return lo + (hi - lo) * raw[names.index(name)]
    rt = descale("parRT") * dt
    ac = descale("parAC")
    a = float(upstream_area_km2)
    if a < 2500.0:
        factor = np.clip((a - ac) / 1000.0, -1.0, 1.0)
    else:
        factor = np.exp(np.clip(-(a - 2500.0) / 50.0, -10.0, 0.0))
    lf = rt * factor  # per component, per step
    n = slz.shape[0]
    added = np.zeros(n)
    for t in range(n):
        before = (slz[t - 1] if t > 0 else np.full(nmul, 0.001)) + percolation[t]
        added[t] = float(np.maximum(lf, -before).mean())
    return added

def snow_module_outflow(core, static_params, forcing: torch.Tensor,
                        elevation_m: float, tosoil: torch.Tensor) -> np.ndarray:
    """Liquid water leaving HBV's snow module each step.

    HBV partitions precipitation into RAIN and SNOW using parTT, with a
    fixed 4 degC threshold above 2000 m. RAIN bypasses snow storage while
    tosoil is liquid water released from MELTWATER. Their sum is therefore
    the liquid flux delivered from the snow module to the soil system.
    """
    names = [p for p in core.phy_param_names if p not in core.dynamic_params]
    nmul = core.nmul
    raw = static_params[0, : len(names) * nmul].view(len(names), nmul)

    lo, hi = core.parameter_bounds["parTT"]
    par_tt = lo + (hi - lo) * raw[names.index("parTT")]

    elevation = torch.full_like(par_tt, float(elevation_m))
    threshold = torch.where(
        elevation >= 2000.0,
        torch.full_like(par_tt, 4.0),
        par_tt,
    )

    pr = forcing[:, 0, 0].unsqueeze(-1)
    tas = forcing[:, 0, 1].unsqueeze(-1)
    rain = pr * (tas >= threshold.unsqueeze(0)).to(pr.dtype)

    return (rain.mean(-1) + tosoil[:, 0, 0]).double().numpy()

def standardise(values: np.ndarray, names: list[str], norm: dict) -> np.ndarray:
    mean = np.array([norm[v][2] for v in names], dtype=np.float64)
    std = np.array([norm[v][3] for v in names], dtype=np.float64)
    return (values - mean) / (std + EPS)


# --- the record ---------------------------------------------------------------


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet"):
            row[key] = float(row[key])
    return rows


def simulate(forcing: list[dict], static: dict, timestep: str) -> tuple[list[dict], dict]:
    torch.set_num_threads(THREADS)
    config, model, norm, structure = load_model()
    dt = TIMESTEP_DAYS.get(timestep, 1.0)

    pr_rate = np.array([r["pr"] for r in forcing])
    tas = np.array([r["tas"] for r in forcing])
    pet_rate = np.array([r["pet"] for r in forcing])
    n = len(pr_rate)

    # The physics wants depths per row; the networks read rates per day, as
    # they were trained to, and write daily parameters that the core rescales.
    pr = pr_rate * dt
    pet = pet_rate * dt
    model.phy_model.dt_days = dt

    # Attributes from the first year only: a climatology the model has seen,
    # through which nothing later in the record can reach back.
    n_clim = min(n, int(round(CLIMATOLOGY_DAYS / dt)))
    derived = climate_attributes(
        pr_rate[:n_clim], tas[:n_clim], pet_rate[:n_clim],
        float(static.get("snow_threshold_degC", 0.0)),
    )
    derived["uparea"] = float(static["area_km2"])
    derived["glaciers"] = 0.0
    derived["permafrost"] = 0.0
    if "elevation_m" in static:
        derived["meanelevation"] = float(static["elevation_m"])

    attr_names = list(config["model"]["nn"]["attributes"])
    attrs = np.array([[derived.get(a, norm[a][2]) for a in attr_names]])
    c_norm = standardise(attrs, attr_names, norm)

    forcing_names = list(config["model"]["phy"]["forcings"])  # P, T, PET
    x = np.stack([pr, tas, pet], axis=-1)[:, None, :]  # [time, 1, 3], depths per row
    x_rates = np.stack([pr_rate, tas, pet_rate], axis=-1)[:, None, :]  # mm/day, for the LSTM
    x_norm = standardise(x_rates, forcing_names, norm)

    data = {
        "xc_nn_norm": torch.tensor(
            np.concatenate([x_norm, np.repeat(c_norm[None], n, axis=0)], axis=-1),
            dtype=torch.float32,
        ),
        "c_nn_norm": torch.tensor(c_norm, dtype=torch.float32),
        "x_phy": torch.tensor(x, dtype=torch.float32),
        "ac_all": torch.tensor([derived["uparea"]], dtype=torch.float32),
        "elev_all": torch.tensor([derived.get("meanelevation", norm["meanelevation"][2])], dtype=torch.float32),
        "areas": torch.tensor([derived["uparea"]], dtype=torch.float32),
    }
    with torch.no_grad():
        out = model(data)
        _, static_params = model.nn_model(data["xc_nn_norm"], data["c_nn_norm"])
    states = model.phy_model.get_states()
    snowpack, meltwater, sm, suz, slz = (s[:, 0, :].mean(-1).double().numpy() for s in states)
    snm = snow_module_outflow(
        model.phy_model,
        static_params,
        data["x_phy"],
        float(derived.get("meanelevation", norm["meanelevation"][2])),
        out["tosoil"],
    )

    snow_storage = snowpack + meltwater
    pr_depth = np.asarray(
        [float(step["pr"]) * dt for step in forcing],
        dtype=np.float64,
    )
    # The model exposes end-of-step snow states but not the state immediately
    # before row 0, so validate every transition that can be checked directly.
    snow_residual = np.diff(snow_storage) - (pr_depth[1:] - snm[1:])
    snow_abs_residual = np.abs(snow_residual)

    # HBV exposes its native snow states and fluxes as float32. At deep
    # snowpacks their representable spacing can exceed the fixed absolute
    # tolerance used for ordinary states. Retain the 1e-4 mm floor, but
    # allow several native float32 ULPs for roundoff.
    native_scale = np.maximum.reduce([
        np.abs(snowpack[:-1]),
        np.abs(snowpack[1:]),
        np.abs(meltwater[:-1]),
        np.abs(meltwater[1:]),
        np.abs(snm[1:]),
        np.ones_like(snow_residual),
    ])
    native_ulp = np.spacing(
        native_scale.astype(np.float32)
    ).astype(np.float64)
    snow_tolerance = np.maximum(
        1.0e-4,
        4.0 * native_ulp,
    )

    snow_module_max_residual = (
        float(np.max(snow_abs_residual))
        if snow_abs_residual.size
        else 0.0
    )
    snow_module_max_tolerance = (
        float(np.max(snow_tolerance))
        if snow_tolerance.size
        else 1.0e-4
    )

    bad = snow_abs_residual > snow_tolerance
    if np.any(bad):
        worst = int(np.argmax(snow_abs_residual / snow_tolerance))
        raise RuntimeError(
            "reconstructed HBV rain partition is inconsistent with "
            "SNOWPACK + MELTWATER: step residual "
            f"{snow_abs_residual[worst]:.6g} mm exceeds "
            f"precision-aware tolerance {snow_tolerance[worst]:.6g} mm"
        )

    gwex = regional_exchange(model.phy_model, static_params, states[4][:, 0, :].double().numpy(),
                             out["percolation"][:, 0, 0].double().numpy(), derived["uparea"], dt)
    routed = out["streamflow"][:, 0, 0].double().numpy()
    unrouted = out["streamflow_no_rout"][:, 0, 0].double().numpy()
    aet = out["AET_hydro"][:, 0, 0].double().numpy()
    channel = np.cumsum(unrouted) - np.cumsum(routed)

    area_m2 = float(static["area_km2"]) * 1.0e6
    mrro = routed / dt  # depth per row back to a rate per day
    rows = [
        {
            "time": step["time"],
            "pr": step["pr"],
            "evspsbl": float(aet[i] / dt),
            "snm": float(snm[i] / dt),
            "mrro": float(mrro[i]),
            "dis": float(mrro[i] * 1.0e-3 * area_m2 / 86400.0),
            "gwex": float(gwex[i] / dt),
            "mrso": float(sm[i]),
            "snw": float(snowpack[i] + meltwater[i]),
            "canopy": 0.0,
            "gw": float(suz[i] + slz[i]),
            "channel": float(channel[i]),
        }
        for i, step in enumerate(forcing)
    ]
    notes = {
        **structure,
        "timestep": timestep,
        "parameter_units": (
            "daily, as written by the networks" if dt == 1.0 else
            f"rescaled to a {timestep} step: per-day fractions as 1-(1-k)^dt, "
            "per-day amounts as amount*dt, unit hydrograph time constant and length in days"
        ),
        "climatology_rows": int(n_clim),
        "derived_static_attributes": sorted(a for a in derived if a in attr_names),
        "static_attributes_at_training_mean": [a for a in attr_names if a not in derived],
        "gwex": "regional groundwater exchange parRT * clamp((Ac - parAC)/1000, -1, 1), "
                "reconstructed from the network's parameters; declared as a source",
        "snm": (
            "RAIN + tosoil: liquid precipitation bypassing snow storage plus "
            "model-reported liquid water released from MELTWATER, mean over components"
        ),
        "snow_module_max_step_residual_mm": snow_module_max_residual,
        "snow_module_max_step_tolerance_mm": snow_module_max_tolerance,
        "states": {
            "snw": "SNOWPACK + MELTWATER, mean over components",
            "mrso": "SM, mean over components",
            "gw": "SUZ + SLZ, mean over components",
            "channel": "cumulative unrouted minus routed flow: water inside the unit hydrograph",
            "canopy": "identically zero; HBV has no interception store",
        },
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

    rows, notes = simulate(forcing, static, str(request.get("timestep", "PT1D")))
    notes["seed"] = int(request.get("seed", 0))

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

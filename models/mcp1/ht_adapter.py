#!/usr/bin/env python3
"""HydroTuring adapter for MCPModel1, the single-node Mass-Conserving Perceptron.

MC{O_sigma L_sigma^con} of Wang & Gupta (2024a, WRR, e2023WR036461), from
YuanHWang/5MCPmodelnet_HydroTuring at the pinned commit. One store C and three
gates. At row t, with the store as it stands before that row's rain:

    Q_t = oo_t * C_t            oo: gate on the store
    E_t = olc_t * C_t           olc = min(ol, PET/C), ol: gate on PET
    C_{t+1} = C_t - Q_t - E_t + P_t

The daily path below runs the repository's own class,
`MCPBRNN_lib_tools.MCNZoo.MCPBRNN_Generic_PETconstraint_Scaling`, with the
repository's checkpoint and the fixed constants from its evaluation script.

What is reported
----------------
`pr` echoed, `evspsbl` the constrained loss E, `mrro` the output Q, both in mm
per day, and `mrso` the store C at the END of the row. The upstream class saves
the store before the row's update (its `c_n[t]` is C_t), so C_{t+1} is what the
harness must difference against the row's fluxes; it is read from `c_n[t+1]` by
giving the model one extra row of no weather at the end of the record, whose
own outputs are thrown away. No budget term is derived from a residual.

The model has no snow, no canopy, no groundwater and no routing store, and
reads neither temperature nor anything about the catchment. Snow and canopy are
reported as identically zero under `zero_stores` in run.json; the rest are named
under `not_reported`.

Spin-up
-------
The model starts from an empty store. Upstream's published recipe repeats the
first water year three times in front of the record, and it is followed here:
the first 365 days of the rows handed over (all of them, for a shorter record)
run three times through the model before row 0. That is the only use
of later rows, and it reaches no further than the first year, as AGENTS.md
allows a climatology to.

Units at a step other than a day
--------------------------------
The upstream class knows one step, one row meaning one day. A resolution probe
runs it at the next finer step anyway, and that path (`forward_step_scaled`) is
an adapter-side copy of the upstream loop rather than the repository's own
code, with three changes and no others: rain is a depth per row, the gates that
are per-day fractions become 1 - (1 - k)^dt, and the PET cap becomes a depth
per row. The gates still read rates in mm per day, as they were trained. At
dt = 1 the copy reproduces the upstream forward to single-precision round-off
(checked in prototyping: at most 1.2e-4 mm on the store over 3000 rows), and the
daily path does not touch it.

The model is deterministic. The request seed is recorded and otherwise unused.
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
from torch import nn

UPSTREAM = Path(os.environ.get("MCP_UPSTREAM", "/model/upstream"))
sys.path.insert(0, str(UPSTREAM))

from MCPBRNN_lib_tools.MCNZoo import MCPBRNN_Generic_PETconstraint_Scaling  # noqa: E402

MODEL = {"name": "mcp1", "version": "5944b5f-epoch27.3"}
COLUMNS = ["time", "pr", "evspsbl", "mrro", "mrso", "snw", "canopy"]
CHECKPOINT = UPSTREAM / "MCPmodel1" / "model_epoch27.pt"

# Defaults of evaluate_model1_single_node.py at the pinned commit.
C_MEAN = 412.9139085
C_STD = 77.49943867
# Inside the upstream forward; only the step-scaled copy needs them here.
PET_MEAN = 2.9086
PET_STD = 1.8980

SPINUP_REPEATS = int(os.environ.get("MCP_SPINUP_REPEATS", "3"))  # upstream's recipe
SPINUP_DAYS = 365
TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0, "PT15M": 1.0 / 96.0, "PT5M": 1.0 / 288.0, "PT1M": 1.0 / 1440.0}

# Stores this model does not have. They hold no water, so they are reported as
# identically zero rather than left out, and the budget probes can ask the question.
ZERO_STORES = {
    "snw": "no snow store: rain and snow are one input, temperature is never read",
    "canopy": "no interception store",
}

NOT_REPORTED = {
    "gw": "no groundwater store: one undifferentiated store, reported as mrso",
    "channel": "no routing store: the output leaves the store directly",
    "dis": "no river: the model knows no catchment area",
    "gwex": "no exchange with the outside: the store update conserves mass by construction",
    "energy": "hfls, hfss, hfg, rlus and the other energy outputs: the model has no energy balance",
}


class Node(MCPBRNN_Generic_PETconstraint_Scaling):
    """The upstream class, with the step-scaled loop beside its own forward."""

    @torch.no_grad()
    def forward_step_scaled(self, x: torch.Tensor, cmean: float, cstd: float, dt: float) -> dict:
        """The upstream loop with its per-day quantities put in the units of the step."""
        n = x.shape[0]
        c_state = x.new_zeros(1, 1)
        h_n, c_n, lc_n = x.new_zeros(n, 1), x.new_zeros(n, 1), x.new_zeros(n, 1)

        exp_yom = torch.exp(self.weight_r_yom)
        exp_ylm = torch.exp(self.weight_r_ylm)
        exp_yfm = torch.exp(self.weight_r_yfm)
        gate_sum = exp_yom + exp_ylm + exp_yfm
        oo_scale = exp_yom / gate_sum
        ol_scale = exp_ylm / gate_sum
        bias_b0_yom = self.bias_b0_yom.unsqueeze(0)
        bias_b0_ylm = self.bias_b0_ylm.unsqueeze(0)

        for b in range(n):
            precipitation = x[b, 0, 0].reshape(1, 1) * dt
            pet_rate = x[b, 0, 1].reshape(1, 1)
            oo = oo_scale * torch.sigmoid(
                torch.addmm(bias_b0_yom, (c_state - cmean) / cstd, self.weight_b1_yom)
            )
            ol = ol_scale * torch.sigmoid(
                torch.addmm(bias_b0_ylm, (pet_rate - PET_MEAN) / PET_STD, self.weight_b2_ylm)
            )
            oo = 1.0 - (1.0 - oo) ** dt
            ol = 1.0 - (1.0 - ol) ** dt
            if c_state.item() > 0:
                ol_constraint = ol - self.relu_l(ol - pet_rate * dt / c_state)
            else:
                ol_constraint = ol
            f = 1.0 - oo - ol_constraint
            h_n[b, :] = (oo * c_state)[0]
            c_n[b, :] = c_state[0]
            lc_n[b, :] = (ol_constraint * c_state)[0]
            c_state = f * c_state + precipitation
        return {"q": h_n, "c": c_n, "et": lc_n}


class Model(nn.Module):
    """As evaluate_model1_single_node.py builds it, so the checkpoint's keys fit.

    spinLen and traintimeLen only slice the observations for a diagnostic the
    adapter does not read; there are no observations here.
    """

    def __init__(self, n_rows: int):
        super().__init__()
        self.MCPBRNNNode = Node(
            input_size=1, hidden_size=1, gate_dim=1,
            spinLen=0, traintimeLen=n_rows, initial_forget_bias=0,
        )


def load_model(n_rows: int) -> Model:
    model = Model(n_rows)
    model.load_state_dict(torch.load(CHECKPOINT, map_location="cpu"), strict=True)
    return model.eval()


def run_model(model: Model, x: torch.Tensor, dt: float) -> dict:
    """Depths per row: `q` and `et` per row, `c` the store before the row."""
    node = model.MCPBRNNNode
    if dt != 1.0:
        return node.forward_step_scaled(x, C_MEAN, C_STD, dt)
    y_obs = x.new_zeros(x.shape[0], 1)
    with torch.no_grad():
        h_n, c_n, _l, lc_n, *_rest = node(x, 0, 0, y_obs, C_MEAN, C_STD)
    return {"q": h_n, "c": c_n, "et": lc_n}


def spinup_rows(pr: np.ndarray, pet: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """The first year of the record, or all of it when shorter, repeated."""
    first = min(int(round(SPINUP_DAYS / dt)), len(pr))
    return np.tile(pr[:first], SPINUP_REPEATS), np.tile(pet[:first], SPINUP_REPEATS)


def simulate(forcing: list[dict], timestep: str) -> tuple[list[dict], dict]:
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    dt = TIMESTEP_DAYS[timestep]
    pr = np.array([r["pr"] for r in forcing], dtype=np.float64)
    pet = np.array([r["pet"] for r in forcing], dtype=np.float64)
    n = len(pr)

    spin_pr, spin_pet = spinup_rows(pr, pet, dt)
    n_spin = len(spin_pr)
    # One extra row of no weather so that the store after the last real row
    # comes back in c_n; its outputs are discarded.
    x = np.stack(
        [np.concatenate([spin_pr, pr, [0.0]]), np.concatenate([spin_pet, pet, [0.0]])], axis=-1
    )[:, None, :]
    x_t = torch.tensor(x, dtype=torch.float32)

    torch.set_num_threads(1)
    model = load_model(x_t.shape[0])
    out = run_model(model, x_t, dt)

    q = out["q"][n_spin:n_spin + n, 0].double().numpy()
    et = out["et"][n_spin:n_spin + n, 0].double().numpy()
    store_end = out["c"][n_spin + 1:n_spin + n + 1, 0].double().numpy()

    rows = [
        {
            "time": step["time"],
            "pr": step["pr"],
            "evspsbl": float(et[i] / dt),
            "mrro": float(q[i] / dt),
            "mrso": float(store_end[i]),
            "snw": 0.0,
            "canopy": 0.0,
        }
        for i, step in enumerate(forcing)
    ]
    notes = {
        "checkpoint": str(CHECKPOINT.relative_to(UPSTREAM)),
        "timestep": timestep,
        "code_path": (
            "upstream MCPBRNN_Generic_PETconstraint_Scaling.forward" if dt == 1.0 else
            "adapter step-scaled copy of the upstream loop: per-day gates as 1-(1-k)^dt, "
            "rain and the PET cap as depths per row, gates read rates in mm/day"
        ),
        "constants": {"c_mean": C_MEAN, "c_std": C_STD},
        "spinup": f"first {SPINUP_DAYS} days of the rows, {SPINUP_REPEATS} times, from an empty store",
        "states": {"mrso": "the single store C at the end of the row (upstream c_n[t+1])"},
    }
    return rows, notes


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "pet"):
            row[key] = float(row[key])
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    started = time.monotonic()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent

    forcing = read_forcing(io_dir / request["input"]["forcing"])
    rows, notes = simulate(forcing, str(request.get("timestep", "PT1D")))
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
                "zero_stores": ZERO_STORES,
                "not_reported": NOT_REPORTED,
                "notes": notes,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

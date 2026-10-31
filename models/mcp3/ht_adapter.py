#!/usr/bin/env python3
"""HydroTuring adapter for MCPModel3, the MA_5 three-node MCP architecture.

MA_5 of Wang & Gupta (2024b, WRR, e2024WR037224), from
YuanHWang/5MCPmodelnet_HydroTuring at the pinned commit. A primary store C
receives the rain and has three outflows; a surface node S and a groundwater
node G route two of them, and the discharge is what those two release:

    Qm_t = oo_t * C_t        R_t = oog_t * C_t      E_t = olc_t * C_t
    C_{t+1} = C_t - Qm_t - R_t - E_t + P_t          olc = min(ol, PET/C)
    Qs_t = os_t * S_t        S_{t+1} = S_t - Qs_t + Qm_t
    Qg_t = og_t * G_t        G_{t+1} = G_t - Qg_t + R_t
    Q_t = Qs_t + Qg_t

The daily path below runs the repository's own three classes from
`MCPBRNN_lib_tools.MCNZoo`, wired as its evaluation script wires them, with the
repository's checkpoint and the fixed constants of that script (normalisation
of each store, and G's initial value of 1.6871063 mm).

What is reported
----------------
`pr` echoed, `evspsbl` the constrained loss E and `mrro` the total Q in mm per
day, and three stores at the END of the row: `mrso` the primary store C, `gw`
the groundwater node G, `channel` the surface node S, which holds water already
generated as runoff that its routing has not released. The upstream classes
save each store before the row's update (`c_n[t]`), so the store after the row
is read from `c_n[t+1]` by giving the model one extra row of no weather at the
end of the record, whose own outputs are thrown away. No budget term is derived
from a residual.

The model has no snow and no canopy store and reads neither temperature nor
anything about the catchment. Those columns are reported as identically zero
under `zero_stores` in run.json.

Spin-up
-------
The primary and surface stores start empty, G at its published value. Upstream's
published recipe repeats the first water year three times in front of the
record, and it is followed here: the first 365 days of the rows handed over
(all of them, for a shorter record) run three times through the model
before row 0. That is the only use of later rows, and it reaches no further than
the first year, as AGENTS.md allows a climatology to.

Units at a step other than a day
--------------------------------
The upstream classes know one step, one row meaning one day. A resolution probe
runs them at the next finer step anyway, and that path (`forward_step_scaled`
on each node) is an adapter-side copy of the upstream loops rather than the
repository's own code, with three changes and no others: rain is a depth per
row, the gates that are per-day fractions become 1 - (1 - k)^dt, and the PET cap
becomes a depth per row. The gates still read rates in mm per day, as they were
trained. At dt = 1 the copy reproduces the upstream forward to single-precision
round-off (checked in prototyping: at most 1.1e-4 mm on the stores over 3000
rows), and the daily path does not touch it.

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

from MCPBRNN_lib_tools.MCNZoo import (  # noqa: E402
    MCPBRNN_Generic_PETconstraint_Two_VariantOutputGate,
    MCPBRNN_GWVariant_Routing,
    MCPBRNN_SW_Variant_Routing,
)

MODEL = {"name": "mcp3", "version": "5944b5f-epoch619.3"}
COLUMNS = ["time", "pr", "evspsbl", "mrro", "mrso", "snw", "canopy", "gw", "channel"]
CHECKPOINT = UPSTREAM / "MCPmodel3" / "model_epoch619.pt"

# Defaults of evaluate_model3_mca5.py at the pinned commit.
C_MEAN = 412.9139085
C_STD = 77.49943867
SW_MEAN = 3.051248249
SW_STD = 3.365261828
GW_MEAN = 1.87176748
GW_STD = 2.085571043
INI_C = 1.6871063
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
    "dis": "no river: the model knows no catchment area",
    "gwex": "no exchange with the outside: every node's update conserves mass by construction",
    "energy": "hfls, hfss, hfg, rlus and the other energy outputs: the model has no energy balance",
}


def to_step(gate: torch.Tensor, dt: float) -> torch.Tensor:
    """A per-day fraction as the fraction of one row."""
    return 1.0 - (1.0 - gate) ** dt


class PrimaryNode(MCPBRNN_Generic_PETconstraint_Two_VariantOutputGate):
    """The upstream class, with the step-scaled loop beside its own forward."""

    @torch.no_grad()
    def forward_step_scaled(self, x: torch.Tensor, p_mean: float, p_std: float, dt: float) -> dict:
        n = x.shape[0]
        c_state = x.new_zeros(1, 1)
        h_n, gw_n, c_n, lc_n = (x.new_zeros(n, 1) for _ in range(4))
        bias_b0_yom = self.bias_b0_yom.unsqueeze(0)
        bias_b0_yom_gw = self.bias_b0_yom_gw.unsqueeze(0)
        bias_b0_ylm = self.bias_b0_ylm.unsqueeze(0)
        exp_yom = torch.exp(self.weight_r_yom)
        exp_gw = torch.exp(self.weight_r_yom_gw)
        exp_ylm = torch.exp(self.weight_r_ylm)
        exp_yfm = torch.exp(self.weight_r_yfm)
        gate_sum = exp_yom + exp_gw + exp_ylm + exp_yfm

        for b in range(n):
            precipitation = x[b, 0, 0].reshape(1, 1) * dt
            pet_rate = x[b, 0, 1].reshape(1, 1)
            gate_surface = (exp_yom / gate_sum) * torch.sigmoid(
                torch.addmm(bias_b0_yom, (c_state - p_mean) / p_std, self.weight_b1_yom)
            )
            gate_groundwater = (exp_gw / gate_sum) * torch.sigmoid(
                torch.addmm(bias_b0_yom_gw, (c_state - p_mean) / p_std, self.weight_b1_yom_gw)
            )
            gate_loss = (exp_ylm / gate_sum) * torch.sigmoid(
                torch.addmm(bias_b0_ylm, (pet_rate - PET_MEAN) / PET_STD, self.weight_b2_ylm)
            )
            gate_surface = to_step(gate_surface, dt)
            gate_groundwater = to_step(gate_groundwater, dt)
            gate_loss = to_step(gate_loss, dt)
            if c_state.item() > 0:
                constrained_loss = gate_loss - self.relu_l(gate_loss - pet_rate * dt / c_state)
            else:
                constrained_loss = gate_loss
            gate_retention = 1.0 - gate_surface - gate_groundwater - constrained_loss

            h_n[b, :] = (gate_surface * c_state)[0]
            gw_n[b, :] = (gate_groundwater * c_state)[0]
            c_n[b, :] = c_state[0]
            lc_n[b, :] = (constrained_loss * c_state)[0]
            c_state = gate_retention * c_state + precipitation
        return {"q": h_n, "recharge": gw_n, "c": c_n, "et": lc_n}


def routing_step_scaled(node, inflow: torch.Tensor, c_mean: float, c_std: float,
                        initial: float, dt: float) -> tuple[torch.Tensor, torch.Tensor]:
    """The loop the upstream routing nodes share. `inflow` is a depth per row already."""
    n = inflow.shape[0]
    c_state = inflow.new_ones(1, 1) * initial
    h_n, c_n = inflow.new_zeros(n, 1), inflow.new_zeros(n, 1)
    bias_b0_yom = node.bias_b0_yom.unsqueeze(0)
    exp_yom = torch.exp(node.weight_r_yom)
    exp_yfm = torch.exp(node.weight_r_yfm)
    for b in range(n):
        gate_output = (exp_yom / (exp_yom + exp_yfm)) * torch.sigmoid(
            torch.addmm(bias_b0_yom, (c_state - c_mean) / c_std, node.weight_b1_yom)
        )
        gate_output = to_step(gate_output, dt)
        h_n[b, :] = (gate_output * c_state)[0]
        c_n[b, :] = c_state[0]
        c_state = (1.0 - gate_output) * c_state + inflow[b, 0].reshape(1, 1)
    return h_n, c_n


class Model(nn.Module):
    """As evaluate_model3_mca5.py builds it, so the checkpoint's keys fit.

    spinLen and traintimeLen only slice the observations for a diagnostic the
    adapter does not read; there are no observations here.
    """

    def __init__(self, n_rows: int):
        super().__init__()
        common = dict(input_size=1, hidden_size=1, gate_dim=1, spinLen=0, traintimeLen=n_rows,
                      initial_forget_bias=0)
        self.MCPBRNNNode = PrimaryNode(**common)
        self.MCPBRNNNode_SW = MCPBRNN_SW_Variant_Routing(**common)
        self.MCPBRNNNode_GW = MCPBRNN_GWVariant_Routing(**common)


def load_model(n_rows: int) -> Model:
    model = Model(n_rows)
    model.load_state_dict(torch.load(CHECKPOINT, map_location="cpu"), strict=True)
    return model.eval()


def run_model(model: Model, x: torch.Tensor, dt: float) -> dict:
    """Depths per row. `q` is the total; stores are as they stand before the row."""
    if dt != 1.0:
        with torch.no_grad():
            main = model.MCPBRNNNode.forward_step_scaled(x, C_MEAN, C_STD, dt)
            sw_q, sw_c = routing_step_scaled(model.MCPBRNNNode_SW, main["q"], SW_MEAN, SW_STD, 0.0, dt)
            gw_q, gw_c = routing_step_scaled(model.MCPBRNNNode_GW, main["recharge"], GW_MEAN, GW_STD, INI_C, dt)
        return {"q": sw_q + gw_q, "c": main["c"], "sw": sw_c, "gw": gw_c, "et": main["et"]}

    y_obs = x.new_zeros(x.shape[0], 1)
    with torch.no_grad():
        main_q, main_c, _l, lc_n, _bp, recharge, *_rest = model.MCPBRNNNode(x, 0, 0, y_obs, C_MEAN, C_STD)
        sw_q, sw_c, *_ = model.MCPBRNNNode_SW(main_q, 0, 0, y_obs, SW_MEAN, SW_STD)
        gw_q, gw_c, *_ = model.MCPBRNNNode_GW(recharge, 0, 0, y_obs, GW_MEAN, GW_STD, INI_C)
    return {"q": sw_q + gw_q, "c": main_c, "sw": sw_c, "gw": gw_c, "et": lc_n}


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
    # One extra row of no weather so that the stores after the last real row
    # come back in c_n; its outputs are discarded.
    x = np.stack(
        [np.concatenate([spin_pr, pr, [0.0]]), np.concatenate([spin_pet, pet, [0.0]])], axis=-1
    )[:, None, :]
    x_t = torch.tensor(x, dtype=torch.float32)

    torch.set_num_threads(1)
    model = load_model(x_t.shape[0])
    out = run_model(model, x_t, dt)

    q = out["q"][n_spin:n_spin + n, 0].double().numpy()
    et = out["et"][n_spin:n_spin + n, 0].double().numpy()
    main_end = out["c"][n_spin + 1:n_spin + n + 1, 0].double().numpy()
    sw_end = out["sw"][n_spin + 1:n_spin + n + 1, 0].double().numpy()
    gw_end = out["gw"][n_spin + 1:n_spin + n + 1, 0].double().numpy()

    rows = [
        {
            "time": step["time"],
            "pr": step["pr"],
            "evspsbl": float(et[i] / dt),
            "mrro": float(q[i] / dt),
            "mrso": float(main_end[i]),
            "snw": 0.0,
            "canopy": 0.0,
            "gw": float(gw_end[i]),
            "channel": float(sw_end[i]),
        }
        for i, step in enumerate(forcing)
    ]
    notes = {
        "checkpoint": str(CHECKPOINT.relative_to(UPSTREAM)),
        "timestep": timestep,
        "code_path": (
            "upstream MCNZoo classes, wired as evaluate_model3_mca5.py wires them" if dt == 1.0 else
            "adapter step-scaled copy of the upstream loops: per-day gates as 1-(1-k)^dt, "
            "rain and the PET cap as depths per row, gates read rates in mm/day"
        ),
        "constants": {
            "c_mean": C_MEAN, "c_std": C_STD, "sw_mean": SW_MEAN, "sw_std": SW_STD,
            "gw_mean": GW_MEAN, "gw_std": GW_STD, "initial_gw_store": INI_C,
        },
        "spinup": f"first {SPINUP_DAYS} days of the rows, {SPINUP_REPEATS} times; "
                  "primary and surface stores from empty, groundwater from its published value",
        "states": {
            "mrso": "primary store C at the end of the row",
            "gw": "groundwater routing node G at the end of the row",
            "channel": "surface routing node S at the end of the row",
        },
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

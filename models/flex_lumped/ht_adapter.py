#!/usr/bin/env python3
"""HydroTuring adapter for the lumped FLEX/HBV model of chrimerss/HydrologicModels.

The equations are those of `lumped_model/HBVMod.py` at commit cc0aa6f
(the semi-distributed directory's copy of the same file, which clips every
outflow at its store), rewritten so that the step is a parameter rather
than the constant `dt = 1` of the original, and with two deviations that
are stated here so they can be argued with:

* Transpiration in the original is `Ep * Su / (Sumax * Ce)` with no upper
  limit, so for `Ce < 1` a wet soil transpires above the potential rate.
  FLEX's published form is `Ep * min(1, Su / (Ce * Sumax))`; that limit is
  applied here. A physical reference model must not evaporate more than
  the atmosphere asks for, and without the limit the et_plausible
  criterion catches it (README.md has the number).
* The catchment's soil and canopy capacities come from static.json rather
  than from the calibrated parameter set, because a physical model is told
  its catchment; the remaining parameters keep their Wark values.

Stores reported: canopy (Si), soil (Su), groundwater (the slow reservoir
Ss), and the fast reservoir plus the water inside the triangular lag as
`channel`, since both hold runoff that has been generated and not yet
released. No snow module: snow is identically zero and precipitation
below freezing is treated as rain, which conserves water and is wrong
about timing, a limitation of the model rather than the adapter.

Rates with a time in their units are rescaled to the step exactly as any
submitted model's must be: per-day fractions as 1 - (1 - k)^dt, per-day
amounts as amount * dt, the lag in days.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLUMNS = ["time", "pr", "evspsbl", "mrro", "gwex", "mrso", "snw", "canopy", "gw", "channel"]
MODEL = {"name": "flex_lumped", "version": "1.0.0"}

TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0, "PT15M": 1.0 / 96.0, "PT5M": 1.0 / 288.0, "PT1M": 1.0 / 1440.0}

# lumped_model/A_MC_HBV.py's feasible ranges, and a set inside them near the
# repository's best Wark fit (C_run_model_lumped.py): Imax Ce Sumax beta Pmax Tlag Kf Ks.
PARAMS = {
    "Imax": 2.0,      # mm, replaced by canopy_capacity_mm when the catchment gives one
    "Ce": 0.68,       # soil moisture fraction at which transpiration reaches its potential
    "Sumax": 90.0,    # mm, replaced by soil_capacity_mm when the catchment gives one
    "beta": 1.85,     # partition curvature
    "Pmax": 0.09,     # mm/day, percolation at a full unsaturated store
    "Tlag": 1.1,      # days, triangular lag
    "Kf": 0.1,        # 1/day, fast reservoir
    "Ks": 0.008,      # 1/day, slow reservoir
}


def lag_weights(tlag_steps: float) -> list[float]:
    """Weigfun.py: a triangle of base Tlag, discretised to steps, summing to one."""
    nmax = int(-(-tlag_steps // 1))  # ceil
    if nmax <= 1:
        return [1.0]
    w = [0.0] * nmax
    th = tlag_steps / 2.0
    nh = int(th // 1)
    for i in range(nh):
        w[i] = ((i + 1) - 0.5) / th
    i = nh
    w[i] = (1 + ((i + 1) - 1) / th) * (th - (th // 1)) / 2 + (1 + (tlag_steps - (i + 1)) / th) * ((th // 1) + 1 - th) / 2
    for i in range(nh + 1, int(tlag_steps // 1)):
        w[i] = (tlag_steps - (i + 1) + 0.5) / th
    if tlag_steps > tlag_steps // 1:
        w[int(tlag_steps // 1)] = (tlag_steps - (tlag_steps // 1)) ** 2 / (2 * th)
    total = sum(w)
    return [x / total for x in w]


def per_step(fraction_per_day: float, dt: float) -> float:
    return 1.0 - (1.0 - fraction_per_day) ** dt


def simulate(forcing: list[dict], static: dict, dt: float) -> list[dict]:
    p = dict(PARAMS)
    if "soil_capacity_mm" in static:
        p["Sumax"] = float(static["soil_capacity_mm"])
    if "canopy_capacity_mm" in static:
        p["Imax"] = float(static["canopy_capacity_mm"])
    kf, ks = per_step(p["Kf"], dt), per_step(p["Ks"], dt)
    pmax = p["Pmax"] * dt
    weights = lag_weights(p["Tlag"] / dt)

    si = 0.0
    su = 0.5 * p["Sumax"]
    sf = 0.0
    ss = 0.0
    generated: list[float] = []  # unrouted runoff per step, for the lag
    rows = []
    for step in forcing:
        pr_rate, pet_rate = step["pr"], step["pet"]
        P = pr_rate * dt
        Ep = pet_rate * dt

        # The human term, first call on the store: a prescribed withdrawal
        # (`abstr`, mm/day net of return flow) is taken from the unsaturated
        # store, and whatever the store cannot supply later from the day's
        # generated runoff. Absent the column nothing changes.
        want = max(step.get("abstr", 0.0), 0.0) * dt
        removed = min(su, want)
        su -= removed

        # Interception store; evaporation from it only on rainless steps.
        if P > 0.0:
            si += P
            pe = max(0.0, si - p["Imax"])
            si -= pe
            ei = 0.0
        else:
            pe = 0.0
            ei = min(Ep, si)
            si -= ei

        # Unsaturated store: a share rho of effective rain goes to fast runoff.
        if pe > 0.0:
            rho = (su / p["Sumax"]) ** p["beta"]
            su += (1.0 - rho) * pe
            quf = rho * pe
        else:
            quf = 0.0

        # Transpiration, limited by soil moisture and by demand.
        ep_left = max(0.0, Ep - ei)
        ea = ep_left * min(1.0, su / (p["Sumax"] * p["Ce"]))
        ea = min(ea, su)
        su -= ea

        # Percolation to the slow reservoir.
        qus = min(pmax * su / p["Sumax"], su)
        su -= qus

        # Fast and slow linear reservoirs.
        sf += quf
        qf = min(kf * sf, sf)
        sf -= qf
        ss += qus
        qs = min(ks * ss, ss)
        ss -= qs

        # The human term, second call: the day's generated runoff.
        rest = want - removed
        take = min(qf, rest)
        qf -= take
        removed += take
        rest -= take
        take = min(qs, rest)
        qs -= take
        removed += take

        generated.append(qf + qs)
        n = len(generated)
        routed = sum(weights[k] * generated[n - 1 - k] for k in range(min(len(weights), n)))

        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (ei + ea) / dt,
            "mrro": routed / dt,
            "gwex": -removed / dt,
            "mrso": su,
            "snw": 0.0,
            "canopy": si,
            "gw": ss,
            "channel": sf,  # completed below with the lag's contents
        })
    # Water inside the lag: cumulative generated minus cumulative routed.
    cum_gen = 0.0
    cum_out = 0.0
    for row, g in zip(rows, generated):
        cum_gen += g
        cum_out += row["mrro"] * dt
        row["channel"] += cum_gen - cum_out
    return rows


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet", "abstr"):
            if key in row:
                row[key] = float(row[key])
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent
    forcing = read_forcing(io_dir / request["input"]["forcing"])
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    timestep = request.get("timestep", "PT1D")
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    rows = simulate(forcing, static, TIMESTEP_DAYS[timestep])

    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows),
                    "notes": {"snw": "identically zero; the model has no snow module",
                              "channel": "fast reservoir plus water inside the triangular lag"}}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

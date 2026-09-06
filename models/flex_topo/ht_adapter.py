#!/usr/bin/env python3
"""HydroTuring adapter for FLEX-Topo from chrimerss/HydrologicModels.

The equations are those of `semi-distributed_model/{FLEXtopo,plateau,
hillslope,wetland}.py` at commit cc0aa6f, rewritten with the step as a
parameter and area-weighted so that every reported quantity is a depth
over the whole catchment. Three landscape units run side by side on the
same weather:

* plateau:   interception, unsaturated store with beta partition, fast
             reservoir, percolation `Pmax * Su/Sumax` to the shared slow
             reservoir;
* hillslope: the same, but percolation is a share `D` of the fast
             partition (preferential flow) rather than a soil-moisture
             function;
* wetland:   no percolation; capillary rise `Cmax * (1 - Su/Sumax)` from
             the shared slow reservoir into its unsaturated store.

The slow reservoir drains linearly (`Ks`); total runoff is the slow
outflow plus each unit's fast outflow times its fraction, passed through a
triangular lag.

Deviations from the repository code, stated so they can be argued with:

* Transpiration is `Ep * min(1, Su / (Ce * Sumax))`, FLEX's published form;
  the original omits the `min(1, .)` and can transpire above demand.
* On the hillslope the preferential share `D` of the fast partition goes
  to the slow reservoir and the remainder to the fast reservoir. The
  original adds the share to the slow reservoir, drops it from the fast
  partition, and also removes it from the unsaturated store, so that share
  of every storm is lost twice over; here it is moved once.
* The wetland's capillary rise is limited by, and removed from, the slow
  reservoir using the wetland fraction. The original passes the plateau
  fraction into the wetland routine (`landscapes[2]`), which creates or
  destroys water whenever the two fractions differ. With the Wark
  fractions used here they differ by a factor of eight.
* Soil and canopy capacities are told to the model by static.json: the
  three units' `Sumax` are scaled so that their area-weighted mean equals
  the catchment's, keeping their ratios, and every `Imax` is the canopy
  capacity.

Stores reported, all area-weighted: canopy (Si), soil (Su), groundwater
(the shared Ss), and the fast reservoirs plus the water inside the lag as
`channel`. No snow module: snow is identically zero.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLUMNS = ["time", "pr", "evspsbl", "mrro", "mrso", "snw", "canopy", "gw", "channel"]
MODEL = {"name": "flex_topo", "version": "1.0.0"}
TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0, "PT15M": 1.0 / 96.0, "PT5M": 1.0 / 288.0, "PT1M": 1.0 / 1440.0}

# Wark catchment landscape fractions from the repository's HAND, slope and
# basin grids with A_landscapes.py's rule (hillslope: slope > 11; plateau:
# HAND > 5 and slope < 11; wetland: 0 < HAND <= 5), normalised over the
# classified cells (5.3 percent of the basin, the channel cells at HAND 0,
# fall in no class).
_CLASSIFIED_CELLS = {"plateau": 6494, "hillslope": 6449, "wetland": 844}  # of 14550 basin cells
FRACTIONS = {k: v / sum(_CLASSIFIED_CELLS.values()) for k, v in _CLASSIFIED_CELLS.items()}

# B_run_model.py's parameter sets for the Wark catchment.
UNITS = {
    "plateau":   {"Imax": 3.2,  "Ce": 0.50, "Sumax": 17.40,  "beta": 0.95, "Pmax": 1.76, "Kf": 0.91},
    "hillslope": {"Imax": 3.25, "Ce": 0.50, "Sumax": 321.99, "beta": 0.99, "D": 0.4,     "Kf": 0.97},
    "wetland":   {"Imax": 9.94, "Ce": 0.50, "Sumax": 53.25,  "beta": 0.70, "Cmax": 0.65, "Kf": 0.45},
}
CATCHMENT = {"Ks": 0.0281, "Tlag": 2.21}


def lag_weights(tlag_steps: float) -> list[float]:
    nmax = int(-(-tlag_steps // 1))
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


def scaled_parameters(static: dict, dt: float) -> tuple[dict, dict]:
    units = {name: dict(p) for name, p in UNITS.items()}
    if "soil_capacity_mm" in static:
        mean = sum(FRACTIONS[n] * units[n]["Sumax"] for n in units)
        factor = float(static["soil_capacity_mm"]) / mean
        for p in units.values():
            p["Sumax"] *= factor
    if "canopy_capacity_mm" in static:
        for p in units.values():
            p["Imax"] = float(static["canopy_capacity_mm"])
    for p in units.values():
        p["Kf"] = per_step(p["Kf"], dt)
        if "Pmax" in p:
            p["Pmax"] *= dt
        if "Cmax" in p:
            p["Cmax"] *= dt
    catchment = {"Ks": per_step(CATCHMENT["Ks"], dt), "Tlag": CATCHMENT["Tlag"] / dt}
    return units, catchment


def simulate(forcing: list[dict], static: dict, dt: float) -> list[dict]:
    units, catchment = scaled_parameters(static, dt)
    weights = lag_weights(catchment["Tlag"])
    frac = FRACTIONS

    si = {n: 0.0 for n in units}
    su = {n: 0.5 * units[n]["Sumax"] for n in units}
    sf = {n: 0.0 for n in units}
    ss = 0.0
    generated: list[float] = []
    rows = []
    for step in forcing:
        pr_rate, pet_rate = step["pr"], step["pet"]
        P = pr_rate * dt
        Ep = pet_rate * dt

        evap = 0.0
        fast_out = 0.0
        to_slow = 0.0
        for name, p in units.items():
            # Interception; evaporation from it on rainless steps only.
            if P > 0.0:
                si[name] += P
                pe = max(0.0, si[name] - p["Imax"])
                si[name] -= pe
                ei = 0.0
            else:
                pe = 0.0
                ei = min(Ep, si[name])
                si[name] -= ei
            # Unsaturated store and beta partition.
            if pe > 0.0:
                rho = (su[name] / p["Sumax"]) ** p["beta"]
                su[name] += (1.0 - rho) * pe
                quf = rho * pe
            else:
                quf = 0.0
            # Transpiration.
            ep_left = max(0.0, Ep - ei)
            ea = min(ep_left * min(1.0, su[name] / (p["Sumax"] * p["Ce"])), su[name])
            su[name] -= ea
            # Percolation or preferential flow to the slow reservoir. On the
            # hillslope a share D of the fast partition recharges the slow
            # reservoir directly; the rest goes to the fast reservoir.
            if name == "plateau":
                qus = min(p["Pmax"] * su[name] / p["Sumax"], su[name])
                su[name] -= qus
            elif name == "hillslope":
                qus = p["D"] * quf
                quf -= qus
            else:
                qus = 0.0
            to_slow += frac[name] * qus
            # Wetland capillary rise from the shared slow reservoir.
            if name == "wetland":
                qr = p["Cmax"] * (1.0 - su[name] / p["Sumax"])
                qr = min(qr, ss / frac[name])
                qr = min(qr, p["Sumax"] - su[name])
                qr = max(qr, 0.0)
                su[name] += qr
                ss -= qr * frac[name]
            # Fast reservoir.
            sf[name] += quf
            qf = min(p["Kf"] * sf[name], sf[name])
            sf[name] -= qf
            fast_out += frac[name] * qf
            evap += frac[name] * (ei + ea)

        ss += to_slow
        qs = min(catchment["Ks"] * ss, ss)
        ss -= qs

        generated.append(qs + fast_out)
        n = len(generated)
        routed = sum(weights[k] * generated[n - 1 - k] for k in range(min(len(weights), n)))
        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": evap / dt,
            "mrro": routed / dt,
            "mrso": sum(frac[n_] * su[n_] for n_ in units),
            "snw": 0.0,
            "canopy": sum(frac[n_] * si[n_] for n_ in units),
            "gw": ss,
            "channel": sum(frac[n_] * sf[n_] for n_ in units),
        })
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
        for key in ("pr", "tas", "pet"):
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
                    "notes": {"landscape_fractions": FRACTIONS,
                              "snw": "identically zero; the model has no snow module",
                              "channel": "fast reservoirs plus water inside the triangular lag, area-weighted"}}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

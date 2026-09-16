#!/usr/bin/env python3
"""HydroTuring adapter for the lumped FLEX/HBV model of chrimerss/HydrologicModels.

The equations are those of `lumped_model/HBVMod.py` at commit cc0aa6f
(the semi-distributed directory's copy of the same file, which clips every
outflow at its store), rewritten so that the step is a parameter rather
than the constant `dt = 1` of the original, and with three deviations that
are stated here so they can be argued with:

* Transpiration in the original is `Ep * Su / (Sumax * Ce)` with no upper
  limit, so for `Ce < 1` a wet soil transpires above the potential rate.
  FLEX's published form is `Ep * min(1, Su / (Ce * Sumax))`; that limit is
  applied here. A physical reference model must not evaporate more than
  the atmosphere asks for, and without the limit the et_plausible
  criterion catches it (README.md has the number).
* The catchment's soil and canopy capacities come from static.json rather
  than from the calibrated parameter set, because a physical model is told
  its catchment; the remaining parameters keep their Wark values unless the
  following geometry mapping applies.
* When complete channel geometry is supplied, the model's existing triangular
  lag is configured independently of the probe's Snyder criterion. Travel time
  is the centroid-to-outlet channel distance divided by a fixed, literature-
  anchored flood-wave celerity. FLEX's `Tlag` is a triangle base. The travel
  time is assigned to that routing kernel using a row-centred discrete-time
  convention. Without both lengths, it retains the calibrated Wark value.

Stores reported: canopy (Si), soil (Su), groundwater (the slow reservoir
Ss), and the fast reservoir plus the water inside the triangular lag as
`channel`, since both hold runoff that has been generated and not yet
released. No snow module: snow is identically zero and precipitation
below freezing is treated as rain, which conserves water and is wrong
about timing, a limitation of the model rather than the adapter.

Rates with a time in their units are rescaled to the step exactly as any
submitted model's must be: per-day fractions as 1 - (1 - k)^dt, per-day
amounts as amount * dt, the lag in days.

A `stage` is also reported, as a diagnostic rather than a store: the depth
Manning's normal-depth relation gives the reach's own discharge, so that
`momentum/stage-discharge-monotonic` has a gauge to read. A `dis` in m3/s is
reported alongside it — the same flow over the catchment area — so the rating
can be drawn against discharge rather than against a store. Neither is
differenced into any budget.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "gwex", "mrso", "snw", "canopy", "gw", "channel", "stage"]
MODEL = {"name": "flex_lumped", "version": "1.1.0"}

# Default reach geometry, used when the catchment does not hand one over.
# These are the same defaults the reference rating adapters carry, so every
# model judged on this probe is read in the same channel.
DEFAULT_WIDTH_M = 18.0
DEFAULT_SLOPE = 0.0015
DEFAULT_MANNING_N = 0.035
DEFAULT_REACH_LENGTH_M = 4500.0
SECONDS_PER_DAY = 86400.0

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

GEOMETRY_KEYS = (
    "area_km2",
    "main_channel_length_km",
    "centroid_channel_length_km",
)
# Fixed before gate evaluation. Beven (2020, Appendix equations A18--A20 and
# Figure A3) separates water velocity from kinematic-wave celerity and uses
# 1 m/s for an upland-channel example; see also Beven (1979), WRR 15(5).
# This is a synthetic control, not a universal constant or a probe-fitted value.
FLOOD_WAVE_CELERITY_M_S = 1.0
KM_PER_DAY_PER_M_S = 86.4


def routing_parameters(static: dict, dt_days: float) -> tuple[float, dict]:
    """Configure FLEX's triangular lag from geometry when it is available.

    ``Tlag`` in the original FLEX ``Weigfun`` is the full base of a symmetric
    triangle. The routed source depths are interval totals, while the legacy
    kernel bins are indexed from each source row's start. The independent
    channel travel time is the centroid-to-outlet distance divided by a fixed
    flood-wave celerity. To represent that travel time from the source interval
    centre, the continuous mode from the row start is the travel time plus half
    a model step, and ``Tlag`` is twice that mode. Cases without both channel
    lengths retain the calibrated Wark ``Tlag``; an area supplied on its own is
    still checked but cannot define a travel path.
    """
    if not math.isfinite(dt_days) or dt_days <= 0.0:
        raise ValueError("model timestep must be finite and positive")

    area = None
    if "area_km2" in static:
        try:
            area = float(static["area_km2"])
        except (TypeError, ValueError):
            raise ValueError("area_km2 must be numeric") from None
        if not math.isfinite(area) or area <= 0.0:
            raise ValueError("area_km2 must be finite and positive")

    has_length = "main_channel_length_km" in static
    has_centroid = "centroid_channel_length_km" in static
    if not has_length and not has_centroid:
        return float(PARAMS["Tlag"]), {
            "parameter_source": "calibrated Wark fallback",
            "triangle_base_days": float(PARAMS["Tlag"]),
            "triangle_mode_from_source_row_start_days": (
                0.5 * float(PARAMS["Tlag"])
            ),
        }

    missing = [name for name in GEOMETRY_KEYS if name not in static]
    if missing:
        raise ValueError(
            "routing geometry must supply area and both channel lengths together; "
            f"missing {missing}"
        )
    try:
        length = float(static["main_channel_length_km"])
        centroid_length = float(static["centroid_channel_length_km"])
    except (TypeError, ValueError):
        raise ValueError("channel lengths must be numeric") from None
    if not all(math.isfinite(value) for value in (length, centroid_length)):
        raise ValueError("channel lengths must be finite")
    if length <= 0.0 or centroid_length <= 0.0:
        raise ValueError("channel lengths must be positive")
    if centroid_length > length:
        raise ValueError("centroid channel length cannot exceed main-channel length")

    channel_travel_days = centroid_length / (
        FLOOD_WAVE_CELERITY_M_S * KM_PER_DAY_PER_M_S
    )
    if not math.isfinite(channel_travel_days) or channel_travel_days <= 0.0:
        raise ValueError("channel geometry produced an invalid travel time")
    source_centroid_offset_days = 0.5 * dt_days
    triangle_mode_days = channel_travel_days + source_centroid_offset_days
    triangle_base_days = 2.0 * triangle_mode_days
    return triangle_base_days, {
        "parameter_source": "constant-celerity travel time from supplied channel geometry",
        "area_km2": area,
        "main_channel_length_km": length,
        "centroid_channel_length_km": centroid_length,
        "travel_distance_km": centroid_length,
        "distance_definition": (
            "centroid-to-outlet channel path (centroid_channel_length_km)"
        ),
        "flood_wave_celerity_m_s": FLOOD_WAVE_CELERITY_M_S,
        "celerity_reference": "https://doi.org/10.5194/hess-24-2655-2020",
        "celerity_assumption": (
            "fixed first-order upland-channel benchmark; not site-specific"
        ),
        "channel_travel_time_days": channel_travel_days,
        "travel_time_origin": "generated-runoff interval centroid",
        "source_interval_centroid_offset_days": source_centroid_offset_days,
        "triangle_base_days": triangle_base_days,
        "triangle_mode_from_source_row_start_days": triangle_mode_days,
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


def discharge_m3s(runoff_mm_per_day: float, static: dict) -> float:
    """The reach's discharge, in m3/s, from the runoff over the catchment.

    `mrro` is the model's own outflow — the fast reservoir's release plus the
    slow reservoir's, both already passed through the lag — so it is the whole
    of what the reach is carrying and nothing else needs adding to it.
    """
    area_km2 = float(static.get("area_km2", 0.0))
    return max(runoff_mm_per_day, 0.0) * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY


def manning_depth(q_m3s: float, static: dict) -> float:
    """The depth a steady flow makes in the reach's cross-section, in metres.

    A stage is a *length* read off a staff gauge in a cross-section, and the
    length is set by the flow passing through it. Manning's normal depth is
    the honest bridge between the two:

        Q   = w * h * (1/n) * h^(2/3) * S^(1/2)
        h   = ( Q * n / (w * sqrt(S)) )^(3/5)

    The flow is the reach's own discharge and nothing else. An earlier version
    added the groundwater store divided by the step, which was wrong twice
    over: that store's *release* is already inside `mrro` (the slow reservoir
    drains into the runoff this model reports), so the term counted the same
    water twice, and dividing the store rather than its release made the gauge
    read the whole reservoir instead of the water leaving it — six times the
    real flow at a daily step, and twenty-four times more again at an hourly
    one, because the divisor shrinks with the step.

    The depth is a function of the discharge at the same step, so this gauge
    reports a single-valued rating: the model has one state carrying the storm
    and the gauge reads it directly.
    """
    width_m = float(static.get("width_m", DEFAULT_WIDTH_M))
    slope = float(static.get("slope", DEFAULT_SLOPE))
    manning_n = float(static.get("manning_n", DEFAULT_MANNING_N))
    if q_m3s <= 0.0 or width_m <= 0.0 or slope <= 0.0:
        return 0.0
    return (q_m3s * manning_n / (width_m * slope ** 0.5)) ** 0.6


def simulate(forcing: list[dict], static: dict, dt: float) -> list[dict]:
    p = dict(PARAMS)
    if "soil_capacity_mm" in static:
        p["Sumax"] = float(static["soil_capacity_mm"])
    if "canopy_capacity_mm" in static:
        p["Imax"] = float(static["canopy_capacity_mm"])
    kf, ks = per_step(p["Kf"], dt), per_step(p["Ks"], dt)
    pmax = p["Pmax"] * dt
    p["Tlag"], _routing = routing_parameters(static, dt)
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
        # The discharge the gauge reads is the reach's own outflow — the whole
        # of what it is carrying, since the slow reservoir's drainage is
        # already inside `mrro` — reported in m3/s so
        # `momentum/stage-discharge-monotonic` can score the rating against
        # discharge rather than against the store.
        row["dis"] = discharge_m3s(row["mrro"], static)
        row["stage"] = manning_depth(row["dis"], static)
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
    _tlag_days, routing = routing_parameters(static, TIMESTEP_DAYS[timestep])

    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows),
                    "notes": {"snw": "identically zero; the model has no snow module",
                              "channel": "fast reservoir plus water inside the triangular lag",
                              "routing": routing}}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

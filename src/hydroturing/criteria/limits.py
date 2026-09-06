"""Limits a budget must respect even when the budget is not reported.

Three statements that hold for any catchment and need no store to be
reported, or only one. Water that did not run off either evaporated or is
still there, and evaporation cannot exceed demand: so runoff is bounded
below by rain minus demand minus what the catchment can hold, and above by
rain plus what it could have released. Evaporation must answer demand when
the soil is wet and be limited by water when it is dry. And a routing store
holds runoff for at most as long as its hydrograph, so it cannot grow past
the flow that feeds it.

The first is the probe that makes a streamflow-only model falsifiable on
mass. The second is an energy statement made with water variables. The
third is what `channel` was added to the contract for.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


def _capacity(run: RunResult, keys: list[str]) -> float:
    static = run.case.static
    missing = [k for k in keys if k not in static]
    if missing:
        raise ValueError(f"the catchment attributes lack {missing}, needed for the storage bound")
    return float(sum(float(static[k]) for k in keys))


@criterion("runoff_bounds")
def runoff_bounds(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Integrated runoff must lie between rain minus demand and rain, give or
    take what the catchment can store.

    Over a long record, P = ET + Q + dS, and 0 <= ET <= PET, so

        P - PET - S_max <= Q <= P + S_max

    where S_max is the most the catchment could have held at the start: its
    bounded capacities, plus whatever it reports in stores that have no
    capacity to name. A model that reports no stores is held to the
    capacities alone. Both bounds are loosened by a share of the rain, the
    engineering rule again. Neither bound needs evaporation or storage to be
    reported, which is the point: the runoff-only models that the closure
    probe can only call INCOMPLETE are asked a mass question here and have
    to answer it.
    """
    keys = list(params.get("capacity", ["soil_capacity_mm", "canopy_capacity_mm"]))
    slack = float(params.get("tolerance", 0.02))
    w = make_window(run, probe)
    for col in ("pr", "pet"):
        if col not in w.forcing.columns:
            raise ValueError(f"runoff_bounds needs '{col}' in the forcing")
    if "mrro" not in w.table.columns:
        raise ValueError("runoff_bounds needs 'mrro' in the model result")

    rain = float(w.volume(w.forcing["pr"]).sum())
    demand = float(w.volume(w.forcing["pet"]).sum())
    runoff = float(w.volume(w.table["mrro"]).sum())
    # A declared exchange with the outside (`gwex`) is supply the model
    # admits to; both bounds move with it, so a declared source is judged as
    # a source and a hidden one as runoff from nowhere.
    declared = float(w.volume(w.table["gwex"]).sum()) if "gwex" in w.table.columns else 0.0
    stored = _capacity(run, keys)
    for var in ("snw", "gw", "channel"):
        if var in w.state0.index:
            stored += max(0.0, float(w.state0[var]))
    if rain <= 0:
        return CriterionResult(name="runoff_bounds", status=FAIL, message="no rain fell; the case is degenerate")

    lower = rain + declared - demand - stored - slack * rain
    upper = rain + declared + stored + slack * rain
    ratio = runoff / rain
    failures = []
    if runoff < lower:
        failures.append(
            f"only {runoff:.0f} mm ran off of {rain:.0f} mm of rain; even evaporating at the "
            f"full {demand:.0f} mm of demand and filling {stored:.0f} mm of storage leaves "
            f"{lower:.0f} mm that had to run off"
        )
    if runoff > upper:
        failures.append(
            f"{runoff:.0f} mm ran off of {rain:.0f} mm of rain; at most {stored:.0f} mm could "
            f"have been released from storage, so no more than {upper:.0f} mm could have run off"
        )
    ok = not failures
    return CriterionResult(
        name="runoff_bounds",
        status=PASS if ok else FAIL,
        value=ratio,
        threshold=None,
        message=(
            f"runoff is {ratio:.2f} of the rain, inside [{max(lower, 0.0) / rain:.2f}, {upper / rain:.2f}] "
            f"(demand {demand / rain:.2f} of rain, storage {stored:.0f} mm"
            + (f", declared exchange {declared / rain:+.2f} of rain" if declared else "") + ")"
            if ok else "; ".join(failures)
        ),
        diagnostics={"rain_mm": rain, "demand_mm": demand, "runoff_mm": runoff, "declared_mm": declared,
                     "storage_bound_mm": stored, "lower_mm": lower, "upper_mm": upper},
    )


@criterion("demand_consistency")
def demand_consistency(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Evaporation answers demand when water is plentiful and is limited by
    water when it is not.

    The regimes are the model's own: the steps on which its reported soil
    store is in the top share of its range are its wet regime, the bottom
    share its dry one. In the wet regime evaporation must reach at least a
    share of potential and may not exceed it; in the dry regime the ratio
    to potential must be lower than in the wet one by a margin, because
    evaporation efficiency does not fall as the soil gets wetter; and over
    the whole window evaporation may not exceed potential. This is the
    surface energy balance stated with the water variables every model
    that reports evaporation has: the available energy sets the ceiling,
    the available water decides how much of it is used.
    """
    wet_q = float(params.get("wet_quantile", 0.8))
    dry_q = float(params.get("dry_quantile", 0.2))
    min_wet = float(params.get("min_wet_ratio", 0.7))
    min_contrast = float(params.get("min_contrast", 0.05))
    ceiling = float(params.get("max_ratio", 1.0))
    slack = float(params.get("tolerance", 0.01))

    w = make_window(run, probe)
    for col in ("evspsbl", "mrso"):
        if col not in w.table.columns:
            raise ValueError(f"demand_consistency needs '{col}' in the model result")
    if "pet" not in w.forcing.columns:
        raise ValueError("demand_consistency needs 'pet' in the forcing")

    et = w.volume(w.table["evspsbl"])
    pet = w.volume(w.forcing["pet"])
    soil = w.table["mrso"].to_numpy(dtype=float)
    demand_on = pet > 1e-9
    if demand_on.sum() < 10:
        raise ValueError("demand_consistency needs steps with positive potential evaporation")
    hi, lo = np.quantile(soil[demand_on], [wet_q, dry_q])
    wet = (soil >= hi) & demand_on
    dry = (soil <= lo) & demand_on
    wet_ratio = float(et[wet].sum() / pet[wet].sum())
    dry_ratio = float(et[dry].sum() / pet[dry].sum())
    total = float(et.sum() / max(pet.sum(), 1e-12))

    failures = []
    if wet_ratio < min_wet:
        failures.append(
            f"with the soil in the wettest {1 - wet_q:.0%} of its range ({hi:.0f} mm and above), "
            f"evaporation was only {wet_ratio:.2f} of demand (at least {min_wet:g} expected)"
        )
    if wet_ratio > ceiling + slack:
        failures.append(
            f"with the soil in the wettest {1 - wet_q:.0%} of its range, evaporation was "
            f"{wet_ratio:.2f} of demand: more than the atmosphere asked for"
        )
    if dry_ratio > wet_ratio - min_contrast:
        failures.append(
            f"evaporation reached {dry_ratio:.2f} of demand with the soil in the driest "
            f"{dry_q:.0%} of its range ({lo:.0f} mm and below) against {wet_ratio:.2f} in the "
            f"wettest; it should be water-limited when dry"
        )
    if total > ceiling + slack:
        failures.append(f"cumulative evaporation is {total:.3f} of potential (limit {ceiling:g})")

    ok = not failures
    return CriterionResult(
        name="demand_consistency",
        status=PASS if ok else FAIL,
        value=wet_ratio,
        threshold=min_wet,
        message=(
            f"evaporation follows demand when the soil is wet ({wet_ratio:.2f} of potential) and "
            f"is water-limited when it is dry ({dry_ratio:.2f}); {total:.2f} of potential over the record"
            if ok else "; ".join(failures)
        ),
        diagnostics={"wet_ratio": wet_ratio, "dry_ratio": dry_ratio, "total_ratio": total,
                     "wet_threshold_mm": float(hi), "dry_threshold_mm": float(lo),
                     "wet_steps": int(wet.sum()), "dry_steps": int(dry.sum())},
    )


@criterion("routing_conservation")
def routing_conservation(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """A routing store holds runoff for as long as its hydrograph and no longer.

    `channel` is what a model has generated as runoff and not yet released.
    It can never be negative, and it can never hold more than the flow that
    feeds it could have put there over the length of the hydrograph: at any
    step, channel <= max_lag_days times the largest runoff rate seen over
    the preceding window. A kernel that does not sum to one, or a store that
    leaks or accumulates, breaks one of the two.
    """
    max_lag = float(params.get("max_lag_days", 15.0))
    lookback = float(params.get("lookback_days", 2.0 * max_lag))
    slack = float(params.get("tolerance", 0.05))

    w = make_window(run, probe)
    for col in ("mrro", "channel"):
        if col not in w.table.columns:
            raise ValueError(f"routing_conservation needs '{col}' in the model result")
    channel = w.table["channel"].to_numpy(dtype=float)
    runoff = w.table["mrro"].to_numpy(dtype=float)
    n_back = max(1, int(round(lookback / w.dt_days)))
    peak = np.array([runoff[max(0, i - n_back): i + 1].max() for i in range(len(runoff))])
    allowed = max_lag * peak * (1.0 + slack) + 1e-6
    excess = channel - allowed

    failures = []
    if channel.min() < -1e-6:
        failures.append(f"the channel store goes negative ({channel.min():.3g} mm)")
    if excess.max() > 0:
        i = int(excess.argmax())
        failures.append(
            f"the channel holds {channel[i]:.1f} mm on {w.forcing['time'].iloc[i]} while the "
            f"largest runoff of the preceding {lookback:g} days is {peak[i]:.2f} mm/day; a "
            f"{max_lag:g}-day hydrograph cannot hold more than {allowed[i]:.1f} mm of that"
        )
    ok = not failures
    ratio = float((channel / np.maximum(allowed, 1e-9)).max())
    return CriterionResult(
        name="routing_conservation",
        status=PASS if ok else FAIL,
        value=ratio,
        threshold=1.0,
        message=(
            f"the channel store stays non-negative and within {ratio:.2f} of what a "
            f"{max_lag:g}-day hydrograph can hold"
            if ok else "; ".join(failures)
        ),
        diagnostics={"max_channel_mm": float(channel.max()), "max_ratio": ratio},
    )

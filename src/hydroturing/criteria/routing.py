"""Reach-indexed routing criteria."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from hydroturing.criteria.base import (
    FAIL,
    PASS,
    CriterionResult,
    criterion,
    make_window,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


SECONDS_PER_DAY = 86_400.0
TIME_COL = "time"
REACH_COL = "reach_id"
INFLOW_COL = "reach_q_in"
OUTFLOW_COL = "q_out"
STORAGE_COL = "channel_storage"


@dataclass(frozen=True)
class RoutingWindow:
    """Scored routing rows and the preceding storage for every reach."""

    table: pd.DataFrame
    state0: dict[str, float]
    reaches: tuple[str, ...]
    network: dict[str, Any]
    dt_seconds: float


def _reject_unknown(name: str, params: dict, allowed: set[str]) -> None:
    unknown = set(params) - allowed
    if unknown:
        raise ValueError(f"{name}: unknown parameters {sorted(unknown)}")


def _finite_float(
    name: str,
    value: Any,
    *,
    positive: bool = False,
    nonnegative: bool = False,
) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite; got {value!r}")
    if positive and result <= 0.0:
        raise ValueError(f"{name} must be positive; got {result}")
    if nonnegative and result < 0.0:
        raise ValueError(f"{name} must be non-negative; got {result}")
    return result


def _routing_window(run: RunResult) -> RoutingWindow:
    """Return routing rows after spin-up, preserving each reach's prior state."""

    if run.routing is None:
        raise ValueError("routing criterion needs routing output")

    network = run.case.static.get("routing_network")
    if not isinstance(network, dict):
        raise ValueError("routing criterion needs static.json routing_network")

    raw_reaches = network.get("reaches")
    if not isinstance(raw_reaches, list) or not raw_reaches:
        raise ValueError(
            "static.json routing_network.reaches must be a non-empty list"
        )
    reaches = tuple(str(reach) for reach in raw_reaches)

    start = run.case.spinup_steps
    scored_parts: list[pd.DataFrame] = []
    state0: dict[str, float] = {}

    for reach in reaches:
        reach_table = (
            run.routing.loc[run.routing[REACH_COL] == reach]
            .reset_index(drop=True)
        )
        if start >= len(reach_table):
            raise ValueError(
                f"spinup consumes routing reach {reach!r}; nothing left to score"
            )

        prior = max(start - 1, 0)
        state0[reach] = float(reach_table.loc[prior, STORAGE_COL])
        scored_parts.append(reach_table.iloc[start:].reset_index(drop=True))

    return RoutingWindow(
        table=pd.concat(scored_parts, ignore_index=True),
        state0=state0,
        reaches=reaches,
        network=network,
        dt_seconds=run.case.dt_days * SECONDS_PER_DAY,
    )


def _contributing_areas(window: RoutingWindow) -> dict[str, float]:
    raw = window.network.get("contributing_area_km2")
    if not isinstance(raw, dict):
        raise ValueError(
            "routing criterion needs routing_network.contributing_area_km2"
        )

    areas: dict[str, float] = {}
    for reach in window.reaches:
        if reach not in raw:
            raise ValueError(
                f"routing_network.contributing_area_km2 is missing {reach!r}"
            )
        areas[reach] = _finite_float(
            f"contributing area for reach {reach!r}",
            raw[reach],
            positive=True,
        )
    return areas


@criterion("routing_input_fidelity")
def routing_input_fidelity(
    run: RunResult,
    probe: ProbeSpec,
    params: dict,
) -> CriterionResult:
    """Reported headwater inflows must match their prescribed forcing."""

    _reject_unknown(
        "routing_input_fidelity",
        params,
        {"headwater_forcing", "rtol", "atol_m3_s"},
    )

    mapping = params.get("headwater_forcing")
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError(
            "routing_input_fidelity needs a non-empty headwater_forcing map"
        )

    rtol = _finite_float(
        "routing_input_fidelity rtol",
        params.get("rtol", 1.0e-6),
        nonnegative=True,
    )
    atol = _finite_float(
        "routing_input_fidelity atol_m3_s",
        params.get("atol_m3_s", 1.0e-9),
        nonnegative=True,
    )

    ordinary = make_window(run, probe)
    routing = _routing_window(run)

    worst_ratio = 0.0
    worst_error = 0.0
    worst_reach: str | None = None
    worst_forcing: str | None = None
    worst_time: str | None = None

    for raw_reach, raw_forcing in mapping.items():
        reach = str(raw_reach)
        forcing = str(raw_forcing)

        if reach not in routing.reaches:
            raise ValueError(
                f"routing_input_fidelity: unknown headwater reach {reach!r}"
            )
        if forcing not in ordinary.forcing.columns:
            raise ValueError(
                f"routing_input_fidelity: forcing column {forcing!r} is missing"
            )

        reach_table = (
            routing.table.loc[routing.table[REACH_COL] == reach]
            .reset_index(drop=True)
        )
        reported = reach_table[INFLOW_COL].to_numpy(dtype=float)
        prescribed = ordinary.forcing[forcing].to_numpy(dtype=float)

        if len(reported) != len(prescribed):
            raise ValueError(
                f"routing_input_fidelity: reach {reach!r} has "
                f"{len(reported)} scored rows, expected {len(prescribed)}"
            )

        error = np.abs(reported - prescribed)
        allowance = atol + rtol * np.abs(prescribed)
        ratio = error / np.maximum(allowance, np.finfo(float).tiny)
        index = int(np.argmax(ratio))

        if float(ratio[index]) > worst_ratio:
            worst_ratio = float(ratio[index])
            worst_error = float(error[index])
            worst_reach = reach
            worst_forcing = forcing
            worst_time = str(reach_table.loc[index, TIME_COL])

    ok = worst_ratio <= 1.0
    return CriterionResult(
        name="routing_input_fidelity",
        status=PASS if ok else FAIL,
        value=worst_ratio,
        threshold=1.0,
        message=(
            "reported headwater inflows match the prescribed forcing"
            if ok
            else (
                f"reach {worst_reach!r} differs from {worst_forcing!r} "
                f"at {worst_time}: absolute error {worst_error:.3e} m3/s "
                f"is {worst_ratio:.3g} times its tolerance"
            )
        ),
        diagnostics={
            "worst_reach": worst_reach,
            "worst_forcing": worst_forcing,
            "worst_time": worst_time,
            "worst_absolute_error_m3_s": worst_error,
            "rtol": rtol,
            "atol_m3_s": atol,
        },
    )


@criterion("routing_reach_closure")
def routing_reach_closure(
    run: RunResult,
    probe: ProbeSpec,
    params: dict,
) -> CriterionResult:
    """Each reach must conserve inflow, outflow, and channel storage."""

    _reject_unknown(
        "routing_reach_closure",
        params,
        {
            "threshold",
            "aggregate",
            "minimum_inflow_volume_mm",
        },
    )

    threshold = _finite_float(
        "routing_reach_closure threshold",
        params.get("threshold", 0.05),
        nonnegative=True,
    )
    minimum_mm = _finite_float(
        "routing_reach_closure minimum_inflow_volume_mm",
        params.get("minimum_inflow_volume_mm", 1.0),
        nonnegative=True,
    )
    aggregate = params.get(
        "aggregate",
        "sum_absolute_interval_residual",
    )
    if aggregate != "sum_absolute_interval_residual":
        raise ValueError(
            "routing_reach_closure aggregate must be "
            "'sum_absolute_interval_residual'"
        )

    routing = _routing_window(run)
    areas = _contributing_areas(routing)

    worst_ratio = 0.0
    worst_reach: str | None = None
    diagnostics: dict[str, dict[str, float]] = {}

    for reach in routing.reaches:
        reach_table = (
            routing.table.loc[routing.table[REACH_COL] == reach]
            .reset_index(drop=True)
        )
        inflow = reach_table[INFLOW_COL].to_numpy(dtype=float)
        outflow = reach_table[OUTFLOW_COL].to_numpy(dtype=float)
        storage = reach_table[STORAGE_COL].to_numpy(dtype=float)
        previous_storage = np.concatenate(
            ([routing.state0[reach]], storage[:-1])
        )

        inflow_volume = float(np.sum(inflow) * routing.dt_seconds)
        minimum_volume = minimum_mm * areas[reach] * 1_000.0

        if inflow_volume < minimum_volume:
            return CriterionResult(
                name="routing_reach_closure",
                status=FAIL,
                value=0.0,
                threshold=threshold,
                message=(
                    f"reach {reach!r} is degenerate: incoming volume "
                    f"{inflow_volume:.3g} m3 is below the required "
                    f"{minimum_volume:.3g} m3"
                ),
                diagnostics={
                    "degenerate_reach": reach,
                    "inflow_volume_m3": inflow_volume,
                    "minimum_inflow_volume_m3": minimum_volume,
                },
            )

        residual = (
            (inflow - outflow) * routing.dt_seconds
            - (storage - previous_storage)
        )
        absolute_residual = float(np.sum(np.abs(residual)))
        ratio = absolute_residual / inflow_volume

        diagnostics[reach] = {
            "inflow_volume_m3": inflow_volume,
            "sum_absolute_residual_m3": absolute_residual,
            "ratio": ratio,
        }

        if ratio > worst_ratio:
            worst_ratio = ratio
            worst_reach = reach

    ok = worst_ratio <= threshold
    return CriterionResult(
        name="routing_reach_closure",
        status=PASS if ok else FAIL,
        value=worst_ratio,
        threshold=threshold,
        message=(
            "every reach closes its interval water budget"
            if ok
            else (
                f"reach {worst_reach!r} has cumulative absolute residual "
                f"{worst_ratio:.2%} of incoming volume "
                f"(limit {threshold:.2%})"
            )
        ),
        diagnostics={
            "worst_reach": worst_reach,
            "by_reach": diagnostics,
        },
    )


@criterion("routing_junction_closure")
def routing_junction_closure(
    run: RunResult,
    probe: ProbeSpec,
    params: dict,
) -> CriterionResult:
    """Every declared junction must conserve its reported routed flow."""

    _reject_unknown(
        "routing_junction_closure",
        params,
        {
            "threshold",
            "aggregate",
            "minimum_throughflow_volume_mm",
        },
    )

    threshold = _finite_float(
        "routing_junction_closure threshold",
        params.get("threshold", 0.05),
        nonnegative=True,
    )
    minimum_mm = _finite_float(
        "routing_junction_closure minimum_throughflow_volume_mm",
        params.get("minimum_throughflow_volume_mm", 1.0),
        nonnegative=True,
    )
    aggregate = params.get(
        "aggregate",
        "sum_absolute_interval_residual",
    )
    if aggregate != "sum_absolute_interval_residual":
        raise ValueError(
            "routing_junction_closure aggregate must be "
            "'sum_absolute_interval_residual'"
        )

    routing = _routing_window(run)
    areas = _contributing_areas(routing)

    raw_junctions = routing.network.get("junctions")
    if not isinstance(raw_junctions, list) or not raw_junctions:
        raise ValueError(
            "routing_junction_closure needs a non-empty "
            "routing_network.junctions list"
        )

    by_reach = {
        reach: (
            routing.table.loc[routing.table[REACH_COL] == reach]
            .reset_index(drop=True)
        )
        for reach in routing.reaches
    }

    worst_ratio = 0.0
    worst_junction: str | None = None
    diagnostics: dict[str, dict[str, float]] = {}

    for index, raw_junction in enumerate(raw_junctions):
        if not isinstance(raw_junction, dict):
            raise ValueError(
                "every routing_network.junctions entry must be an object"
            )

        junction_id = raw_junction.get("id", f"junction_{index}")
        upstream = raw_junction.get("upstream")
        downstream = raw_junction.get("downstream")

        if (
            not isinstance(junction_id, str)
            or not junction_id
            or not isinstance(upstream, list)
            or not upstream
            or any(
                not isinstance(reach, str) or not reach
                for reach in upstream
            )
            or not isinstance(downstream, str)
            or not downstream
        ):
            raise ValueError(
                "each routing junction needs a non-empty string id, "
                "a non-empty upstream list, and a downstream reach"
            )

        unknown_reaches = sorted(
            (set(upstream) | {downstream}) - set(routing.reaches)
        )
        if unknown_reaches:
            raise ValueError(
                f"junction {junction_id!r} refers to unknown reaches "
                f"{unknown_reaches}"
            )

        upstream_outflow = np.sum(
            [
                by_reach[reach][OUTFLOW_COL].to_numpy(dtype=float)
                for reach in upstream
            ],
            axis=0,
        )
        downstream_inflow = by_reach[downstream][INFLOW_COL].to_numpy(
            dtype=float
        )

        throughflow_volume = float(
            np.sum(upstream_outflow) * routing.dt_seconds
        )
        minimum_volume = minimum_mm * areas[downstream] * 1_000.0

        if throughflow_volume < minimum_volume:
            return CriterionResult(
                name="routing_junction_closure",
                status=FAIL,
                value=threshold + 1.0,
                threshold=threshold,
                message=(
                    f"junction {junction_id!r} is degenerate: throughflow "
                    f"{throughflow_volume:.3g} m3 is below the required "
                    f"{minimum_volume:.3g} m3"
                ),
                diagnostics={
                    "degenerate_junction": junction_id,
                    "throughflow_volume_m3": throughflow_volume,
                    "minimum_throughflow_volume_m3": minimum_volume,
                },
            )

        residual = (
            upstream_outflow - downstream_inflow
        ) * routing.dt_seconds
        absolute_residual = float(np.sum(np.abs(residual)))
        ratio = absolute_residual / throughflow_volume

        diagnostics[junction_id] = {
            "throughflow_volume_m3": throughflow_volume,
            "sum_absolute_residual_m3": absolute_residual,
            "ratio": ratio,
        }

        if ratio > worst_ratio:
            worst_ratio = ratio
            worst_junction = junction_id

    ok = worst_ratio <= threshold
    return CriterionResult(
        name="routing_junction_closure",
        status=PASS if ok else FAIL,
        value=worst_ratio,
        threshold=threshold,
        message=(
            "every routing junction conserves its reported transfer"
            if ok
            else (
                f"junction {worst_junction!r} has cumulative absolute "
                f"residual {worst_ratio:.2%} of throughflow "
                f"(limit {threshold:.2%})"
            )
        ),
        diagnostics={
            "worst_junction": worst_junction,
            "by_junction": diagnostics,
        },
    )


@criterion("routing_recession_bound")
def routing_recession_bound(
    run: RunResult,
    probe: ProbeSpec,
    params: dict,
) -> CriterionResult:
    """A reach may retain only a bounded amount after its inflow recedes."""

    _reject_unknown(
        "routing_recession_bound",
        params,
        {
            "apply_to",
            "max_lag_days",
            "lookback_days",
            "tolerance",
            "min_allowance_mm",
        },
    )

    apply_to = params.get("apply_to", "each_reach")
    if apply_to != "each_reach":
        raise ValueError(
            "routing_recession_bound apply_to must be 'each_reach'"
        )

    max_lag_days = _finite_float(
        "routing_recession_bound max_lag_days",
        params.get("max_lag_days", 15.0),
        positive=True,
    )
    lookback_days = _finite_float(
        "routing_recession_bound lookback_days",
        params.get("lookback_days", 2.0 * max_lag_days),
        positive=True,
    )
    tolerance = _finite_float(
        "routing_recession_bound tolerance",
        params.get("tolerance", 0.05),
        nonnegative=True,
    )
    min_allowance_mm = _finite_float(
        "routing_recession_bound min_allowance_mm",
        params.get("min_allowance_mm", 1e-6),
        nonnegative=True,
    )

    routing = _routing_window(run)
    areas = _contributing_areas(routing)
    n_back = max(
        1,
        int(round(lookback_days / run.case.dt_days)),
    )

    worst_ratio = 0.0
    worst_reach: str | None = None
    worst_time: str | None = None
    failures: list[str] = []
    diagnostics: dict[str, dict[str, float | str]] = {}

    for reach in routing.reaches:
        reach_table = (
            routing.table.loc[routing.table[REACH_COL] == reach]
            .reset_index(drop=True)
        )
        inflow = reach_table[INFLOW_COL].to_numpy(dtype=float)
        storage = reach_table[STORAGE_COL].to_numpy(dtype=float)

        peak = np.array(
            [
                inflow[max(0, index - n_back) : index + 1].max()
                for index in range(len(inflow))
            ],
            dtype=float,
        )
        allowed = (
            max_lag_days
            * peak
            * SECONDS_PER_DAY
            * (1.0 + tolerance)
            + min_allowance_mm * areas[reach] * 1_000.0
        )
        excess = storage - allowed
        ratio = float(
            np.max(storage / np.maximum(allowed, 1e-12))
        )
        index = int(np.argmax(excess))
        maximum_excess = float(excess[index])

        diagnostics[reach] = {
            "max_storage_to_allowance_ratio": ratio,
            "max_excess_m3": maximum_excess,
            "time_of_max_excess": str(reach_table.loc[index, TIME_COL]),
        }

        if ratio > worst_ratio:
            worst_ratio = ratio
            worst_reach = reach
            worst_time = str(reach_table.loc[index, TIME_COL])

        if maximum_excess > 0.0:
            failures.append(
                f"reach {reach!r} holds {storage[index]:.3g} m3 on "
                f"{reach_table.loc[index, TIME_COL]!r}, above its "
                f"{allowed[index]:.3g} m3 allowance"
            )

    ok = not failures
    return CriterionResult(
        name="routing_recession_bound",
        status=PASS if ok else FAIL,
        value=worst_ratio,
        threshold=1.0,
        message=(
            "every reach stays within its recession-storage allowance"
            if ok
            else "; ".join(failures)
        ),
        diagnostics={
            "worst_reach": worst_reach,
            "worst_time": worst_time,
            "max_ratio": worst_ratio,
            "by_reach": diagnostics,
        },
    )


@criterion("routing_state_bounds")
def routing_state_bounds(
    run: RunResult,
    probe: ProbeSpec,
    params: dict,
) -> CriterionResult:
    """Reported routing inflow, outflow, and storage must remain physical."""

    _reject_unknown(
        "routing_state_bounds",
        params,
        {
            "storage_lower_m3",
            "discharge_lower_m3_s",
            "storage_atol_m3",
            "discharge_atol_m3_s",
        },
    )

    storage_lower = _finite_float(
        "routing_state_bounds storage_lower_m3",
        params.get("storage_lower_m3", 0.0),
    )
    discharge_lower = _finite_float(
        "routing_state_bounds discharge_lower_m3_s",
        params.get("discharge_lower_m3_s", 0.0),
    )
    storage_atol = _finite_float(
        "routing_state_bounds storage_atol_m3",
        params.get("storage_atol_m3", 1e-6),
        nonnegative=True,
    )
    discharge_atol = _finite_float(
        "routing_state_bounds discharge_atol_m3_s",
        params.get("discharge_atol_m3_s", 1e-9),
        nonnegative=True,
    )

    routing = _routing_window(run)

    minimum_storage = float(routing.table[STORAGE_COL].min())
    minimum_inflow = float(routing.table[INFLOW_COL].min())
    minimum_outflow = float(routing.table[OUTFLOW_COL].min())

    storage_shortfall = max(
        0.0,
        storage_lower - storage_atol - minimum_storage,
    )
    inflow_shortfall = max(
        0.0,
        discharge_lower - discharge_atol - minimum_inflow,
    )
    outflow_shortfall = max(
        0.0,
        discharge_lower - discharge_atol - minimum_outflow,
    )
    maximum_shortfall = max(
        storage_shortfall,
        inflow_shortfall,
        outflow_shortfall,
    )

    failures: list[str] = []
    if storage_shortfall > 0.0:
        failures.append(
            f"channel storage falls to {minimum_storage:.3g} m3"
        )
    if inflow_shortfall > 0.0:
        failures.append(
            f"reach inflow falls to {minimum_inflow:.3g} m3/s"
        )
    if outflow_shortfall > 0.0:
        failures.append(
            f"reach outflow falls to {minimum_outflow:.3g} m3/s"
        )

    ok = not failures
    return CriterionResult(
        name="routing_state_bounds",
        status=PASS if ok else FAIL,
        value=maximum_shortfall,
        threshold=0.0,
        message=(
            "routing discharge and channel storage remain non-negative"
            if ok
            else "; ".join(failures)
        ),
        diagnostics={
            "minimum_channel_storage_m3": minimum_storage,
            "minimum_reach_q_in_m3_s": minimum_inflow,
            "minimum_q_out_m3_s": minimum_outflow,
            "maximum_shortfall": maximum_shortfall,
        },
    )
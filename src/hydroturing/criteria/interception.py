"""Interception response: a larger canopy must give the water back as evaporation.

The usual paired probe changes the weather. This one changes a declared static
parameter, `canopy_capacity_mm`, and leaves the weather byte-identical. The
question is not whether the budget closes — a model can close it while ignoring
the parameter — but whether the two runs differ, and whether the difference is
allocated as physics requires:

  * runoff may not rise with a larger canopy;
  * evaporation may not fall;
  * the paired budget must account for the difference, using the fluxes and
    every store the model reports;
  * the response must exist and must be an evaporation response, so a model
    that ignores the parameter (two identical runs) cannot pass by doing
    nothing, and a model that parks the water permanently in storage fails
    the same test.

The normalized denominator is the capacity change the case itself prescribes,
read from the two runs' `case.static` rather than from a number in the probe
file, so the criterion cannot drift from the generator.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL,
    PASS,
    CriterionResult,
    criterion,
    make_window,
    reported_states,
)
from hydroturing.criteria.response import pick
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


def _capacity(run: RunResult, parameter: str) -> float:
    static = run.case.static
    if parameter not in static:
        raise ValueError(
            f"interception_response needs '{parameter}' in the case static; "
            "the generator has to declare the capacity it varies"
        )
    try:
        value = float(static[parameter])
    except (TypeError, ValueError):
        raise ValueError(f"case static '{parameter}' must be numeric") from None
    if not np.isfinite(value):
        raise ValueError(f"case static '{parameter}' must be finite")
    return value


def _same_weather(low, high) -> None:
    """Only the declared static may differ; otherwise the residual cannot be attributed."""
    if len(low.forcing) != len(high.forcing):
        raise ValueError(
            "the variants produced forcing records of different lengths "
            f"({len(low.forcing)} against {len(high.forcing)})"
        )
    for column in low.forcing.columns:
        if column not in high.forcing.columns:
            raise ValueError(
                f"forcing column '{column}' is missing from the high variant"
            )
        a = low.forcing[column].to_numpy()
        b = high.forcing[column].to_numpy()
        if a.dtype.kind in "fc" or b.dtype.kind in "fc":
            if not np.allclose(
                np.asarray(a, dtype=float), np.asarray(b, dtype=float),
                rtol=0, atol=1e-9, equal_nan=True,
            ):
                raise ValueError(
                    f"forcing column '{column}' differs between the variants; "
                    "only the declared static parameter may change"
                )
        elif not np.array_equal(a.astype(str), b.astype(str)):
            raise ValueError(
                f"forcing column '{column}' differs between the variants; "
                "only the declared static parameter may change"
            )


@criterion("interception_response", paired=True)
def interception_response(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    parameter = str(params.get("parameter", "canopy_capacity_mm"))
    tolerance = float(params.get("tolerance", 0.05))
    response_floor = float(params.get("response_floor", 0.5))
    floor_mm = float(params.get("floor_mm", 1.0))

    low_run = pick(runs, params, "low", "capacity_low")
    high_run = pick(runs, params, "high", "capacity_high")
    low = make_window(low_run, probe)
    high = make_window(high_run, probe)
    _same_weather(low, high)

    cap_low = _capacity(low_run, parameter)
    cap_high = _capacity(high_run, parameter)
    dcap = cap_high - cap_low
    if dcap <= 0:
        raise ValueError(
            f"the high variant's {parameter} ({cap_high:g}) must exceed the "
            f"low variant's ({cap_low:g}); the generator has to make the two "
            "variants actually differ"
        )

    for var in ("evspsbl", "mrro"):
        if var not in low.table.columns or var not in high.table.columns:
            raise ValueError(f"interception_response needs '{var}' in the model result")

    d_evspsbl = float(
        high.volume(high.table["evspsbl"]).sum()
        - low.volume(low.table["evspsbl"]).sum()
    )
    d_mrro = float(
        high.volume(high.table["mrro"]).sum()
        - low.volume(low.table["mrro"]).sum()
    )
    d_gwex = 0.0
    if "gwex" in low.table.columns and "gwex" in high.table.columns:
        d_gwex = float(
            high.volume(high.table["gwex"]).sum()
            - low.volume(low.table["gwex"]).sum()
        )

    states = reported_states(low, probe)
    if states:
        d_storage = float(
            high.storage(states)[-1] - high.storage_initial(states)
        ) - float(
            low.storage(states)[-1] - low.storage_initial(states)
        )
    else:
        d_storage = 0.0

    residual = d_evspsbl + d_mrro - d_gwex + d_storage
    residual_limit = max(tolerance, floor_mm / dcap)

    runoff_rise = max(0.0, d_mrro) / dcap
    evap_fall = max(0.0, -d_evspsbl) / dcap
    residual_share = abs(residual) / dcap
    evap_response = d_evspsbl / dcap

    # Each failure carries its own normalized exceedance, so the worst seed is
    # the worst check and the reported value is the quantity that failed.
    ranked: list[tuple[str, float, str]] = []
    if runoff_rise > tolerance:
        ranked.append((
            "runoff_rise", runoff_rise / tolerance - 1,
            f"runoff rose with the larger canopy ({d_mrro:+.2f} mm, "
            f"{runoff_rise:.1%} of the capacity change, limit {tolerance:.1%})",
        ))
    if evap_fall > tolerance:
        ranked.append((
            "evap_fall", evap_fall / tolerance - 1,
            f"evaporation fell with the larger canopy ({d_evspsbl:+.2f} mm, "
            f"{evap_fall:.1%}, limit {tolerance:.1%})",
        ))
    if residual_share > residual_limit:
        ranked.append((
            "residual", residual_share / residual_limit - 1,
            f"the paired budget does not account for the difference "
            f"({residual:+.2f} mm, {residual_share:.1%}, limit {residual_limit:.1%})",
        ))
    if evap_response < response_floor:
        ranked.append((
            "evaporation_response", 1 - evap_response / response_floor,
            f"the evaporation response is too small or absent "
            f"({d_evspsbl:+.2f} mm, {evap_response:.1%} of the capacity change, "
            f"floor {response_floor:.0%})",
        ))

    ok = not ranked
    if ok:
        value, threshold = residual_share, residual_limit
        failing_check = None
    else:
        failing_check, value, _ = max(ranked, key=lambda item: item[1])
        threshold = 0.0
    message = (
        f"capacity {cap_low:g} -> {cap_high:g} mm: evaporation {d_evspsbl:+.2f} mm, "
        f"runoff {d_mrro:+.2f} mm, storage {d_storage:+.2f} mm, "
        f"residual {residual:.2f} mm of {dcap:g} mm ({residual_share:.1%})"
        if ok
        else f"capacity {cap_low:g} -> {cap_high:g} mm: " + "; ".join(
            message for _, _, message in ranked
        )
    )
    return CriterionResult(
        name="interception_response",
        status=PASS if ok else FAIL,
        value=value,
        threshold=threshold,
        message=message,
        diagnostics={
            "cap_low_mm": cap_low,
            "cap_high_mm": cap_high,
            "delta_cap_mm": dcap,
            "d_evspsbl_mm": d_evspsbl,
            "d_mrro_mm": d_mrro,
            "d_gwex_mm": d_gwex,
            "d_storage_mm": d_storage,
            "residual_mm": residual,
            "runoff_rise_share": runoff_rise,
            "evap_fall_share": evap_fall,
            "residual_share": residual_share,
            "evap_response_share": evap_response,
            "tolerance": tolerance,
            "response_floor": response_floor,
            "floor_mm": floor_mm,
            "residual_limit": residual_limit,
            "failing_check": failing_check,
        },
    )

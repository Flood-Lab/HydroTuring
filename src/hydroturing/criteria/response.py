"""Counterfactual response: the budget has to move when the forcing moves.

This is the structural answer to closure by construction. A model that solves
for one budget term as the residual closes perfectly and stays inside its
storage bounds if it is careful, and no single-run criterion can distinguish
it from physics with certainty. Run the same seed twice, once with more rain,
and the question changes from "does the budget close" to "where did the extra
water go" — which a model has to actually partition to answer.

The three failure modes this separates:

  nothing responds        the model ignores the perturbation
  one term absorbs it     evaporation swallows every extra millimetre, which
                          is the degenerate partition dressed up
  the shares do not sum   the model invents or destroys water in the
                          difference between the two runs, even though each
                          run closes on its own
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL,
    PASS,
    CriterionResult,
    criterion,
    make_window,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


def pick(runs: dict[str, RunResult], params: dict, key: str, default: str) -> RunResult:
    """Resolve a variant named in the probe to the run the model produced for it."""
    name = str(params.get(key, default))
    if name not in runs:
        raise ValueError(
            f"variant '{name}' is not one of this probe's variants "
            f"({sorted(runs)}); check `case.variants` and the criterion's "
            f"'{key}' parameter"
        )
    return runs[name]


@criterion("counterfactual_response", paired=True)
def counterfactual_response(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    """Added water must be partitioned, and no single term may take it all."""
    driver = str(params.get("driver", "pr"))
    terms = list(params.get("terms", ["evspsbl", "mrro"]))
    min_share = float(params.get("min_share", 0.05))
    max_share = float(params.get("max_share", 0.90))
    sum_tolerance = float(params.get("sum_tolerance", 0.10))

    control = make_window(pick(runs, params, "control", "control"), probe)
    perturbed = make_window(pick(runs, params, "perturbed", "perturbed"), probe)

    if driver not in control.forcing.columns:
        raise ValueError(f"counterfactual_response needs '{driver}' in the forcing")

    added = float(
        perturbed.volume(perturbed.forcing[driver]).sum()
        - control.volume(control.forcing[driver]).sum()
    )
    if added <= 0:
        raise ValueError(
            f"the perturbed variant adds no {driver} ({added:.4g} mm); the "
            "generator has to make the two variants actually differ"
        )

    shares: dict[str, float] = {}
    for var in terms:
        if var not in control.table.columns:
            raise ValueError(f"counterfactual_response needs '{var}' in the model result")
        change = float(
            perturbed.volume(perturbed.table[var]).sum()
            - control.volume(control.table[var]).sum()
        )
        shares[var] = change / added

    states = probe.requires_states
    if states:
        def storage_change(w):
            return float(w.storage(states)[-1]) - w.storage_initial(states)

        shares["storage"] = (storage_change(perturbed) - storage_change(control)) / added

    failures = []
    for var, share in shares.items():
        if var == "storage":
            continue
        if abs(share) < min_share:
            failures.append(
                f"{var} barely responds ({share:+.3f} of the added {driver}, "
                f"minimum {min_share:g})"
            )

    hog, hog_share = max(shares.items(), key=lambda kv: abs(kv[1]))
    if abs(hog_share) > max_share:
        failures.append(
            f"{hog} absorbs {abs(hog_share):.3f} of the added {driver} "
            f"(limit {max_share:g})"
        )

    accounted = float(sum(shares.values()))
    if abs(accounted - 1.0) > sum_tolerance:
        failures.append(
            f"the responses account for {accounted:.3f} of the added {driver}, "
            f"not 1.000 (tolerance {sum_tolerance:g})"
        )

    ok = not failures
    detail = ", ".join(f"{k} {v:+.3f}" for k, v in shares.items())
    return CriterionResult(
        name="counterfactual_response",
        status=PASS if ok else FAIL,
        value=accounted,
        threshold=1.0 + sum_tolerance,
        message=(
            f"added {driver} is partitioned ({detail})" if ok else "; ".join(failures)
        ),
        diagnostics={
            "added_mm": added,
            "shares": {k: float(v) for k, v in shares.items()},
            "largest_share": hog,
            "accounted": accounted,
        },
    )


@criterion("response_sign", paired=True)
def response_sign(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    """Perturb one driver, hold the rest, and require the response to point
    the way the physics says.

    Warm the air over the same rain and more water evaporates, so less runs
    off. That is not a closure statement; it is a statement about the sign
    of an internal relationship, and a model can close its budget perfectly
    while getting it backwards, or while ignoring the driver altogether. Each
    expected variable must move in its declared direction by at least
    `min_share` of the change in the driver's integrated demand, and by no
    more than `max_share` of it: a model that loses more runoff than the
    warming asked for has invented a loss.

    Only the driver may differ between the variants. The precipitation is
    checked to be identical, so that the response cannot be to a change in
    the rain that slipped in with the temperature.
    """
    driver = str(params.get("driver", "pet"))
    fixed = list(params.get("fixed", ["pr"]))
    expect = dict(params.get("expect", {"mrro": "decrease"}))
    min_share = float(params.get("min_share", 0.1))
    max_share = float(params.get("max_share", 1.1))

    control = make_window(pick(runs, params, "control", "control"), probe)
    perturbed = make_window(pick(runs, params, "perturbed", "perturbed"), probe)
    if len(control.table) != len(perturbed.table):
        raise ValueError(
            f"the two variants produced windows of different length "
            f"({len(control.table)} and {len(perturbed.table)}); a driver "
            "perturbation must preserve the number of steps"
        )
    for column in (driver, *fixed):
        if column not in control.forcing.columns or column not in perturbed.forcing.columns:
            raise ValueError(f"response_sign needs '{column}' in the forcing")
    for column in fixed:
        a = np.asarray(control.forcing[column], dtype=float)
        b = np.asarray(perturbed.forcing[column], dtype=float)
        if not np.allclose(a, b, rtol=0, atol=1e-9):
            raise ValueError(
                f"'{column}' differs between the variants; only '{driver}' may "
                "change, or the response cannot be attributed"
            )

    added = float(
        perturbed.volume(perturbed.forcing[driver]).sum()
        - control.volume(control.forcing[driver]).sum()
    )
    if added <= 0:
        raise ValueError(
            f"the perturbed variant adds no {driver} ({added:.4g} mm); the "
            "generator has to raise the driver, not lower it"
        )

    shares: dict[str, float] = {}
    failures: list[str] = []
    for var, direction in expect.items():
        if var not in control.table.columns or var not in perturbed.table.columns:
            continue
        change = float(
            perturbed.volume(perturbed.table[var]).sum()
            - control.volume(control.table[var]).sum()
        )
        share = change / added
        shares[var] = share
        signed = -share if direction == "decrease" else share
        if direction not in ("decrease", "increase"):
            raise ValueError(f"response_sign expects 'increase' or 'decrease' for {var}, not {direction!r}")
        if signed < min_share:
            failures.append(
                f"{var} should {direction} by at least {min_share:.0%} of the added "
                f"{driver} demand; it moved {share:+.3f} of it"
            )
        elif signed > max_share:
            failures.append(
                f"{var} {direction}s by {abs(share):.3f} of the added {driver} demand, "
                f"more than the {max_share:.0%} the perturbation can account for"
            )
    if not shares:
        raise ValueError(
            f"response_sign found none of {list(expect)} in the model result"
        )

    ok = not failures
    detail = ", ".join(f"{k} {v:+.3f}" for k, v in shares.items())
    lead = next(iter(shares))
    return CriterionResult(
        name="response_sign",
        status=PASS if ok else FAIL,
        value=shares[lead],
        threshold=min_share,
        message=(
            f"the response to more {driver} points the right way ({detail}, "
            f"as shares of the added demand)"
            if ok
            else "; ".join(failures)
        ),
        diagnostics={
            "driver": driver,
            "added_demand_mm": added,
            "shares": {k: float(v) for k, v in shares.items()},
        },
    )

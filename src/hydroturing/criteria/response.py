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
    reported_states,
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

    states = reported_states(control, probe)
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
    the way the physics says, in every direction the driver is pushed.

    Warm the air over the same rain and more water evaporates, so less runs
    off; cool it and the reverse must happen. The expectation is written as
    a derivative, the change in each variable per unit change in the
    driver's integrated demand: `decrease` means that derivative is
    negative, `increase` that it is positive, whichever way the driver
    moved. Each perturbed variant is scored against the control on its own,
    because a model can get one direction right and the other wrong: a
    learned rule that only fires when it is hotter, a response that
    saturates, an output clipped at zero.

    Each expected variable must move by at least `min_share` of the change
    in demand, and by no more than `max_share` of it: a model that loses
    more runoff than the warming asked for has invented a loss. The value
    reported is the weakest response found, the one closest to failing.

    Only the driver may differ between the variants. The precipitation is
    checked to be identical, so that the response cannot be to a change in
    the rain that slipped in with the temperature.
    """
    driver = str(params.get("driver", "pet"))
    fixed = list(params.get("fixed", ["pr"]))
    expect = dict(params.get("expect", {"mrro": "decrease"}))
    min_share = float(params.get("min_share", 0.1))
    max_share = float(params.get("max_share", 1.1))
    perturbed_names = params.get("perturbed", "perturbed")
    if isinstance(perturbed_names, str):
        perturbed_names = [perturbed_names]
    for direction in expect.values():
        if direction not in ("decrease", "increase"):
            raise ValueError(
                f"response_sign expects 'increase' or 'decrease', not {direction!r}"
            )

    control = make_window(pick(runs, params, "control", "control"), probe)
    for column in (driver, *fixed):
        if column not in control.forcing.columns:
            raise ValueError(f"response_sign needs '{column}' in the forcing")

    per_variant: dict[str, dict] = {}
    failures: list[str] = []
    weakest: tuple[float, float] | None = None  # (signed margin, share) of the lead variable
    lead = next(iter(expect))

    for name in perturbed_names:
        if name not in runs:
            raise ValueError(
                f"variant '{name}' is not one of this probe's variants ({sorted(runs)}); "
                "check `case.variants` and the criterion's 'perturbed' parameter"
            )
        perturbed = make_window(runs[name], probe)
        if len(control.table) != len(perturbed.table):
            raise ValueError(
                f"variant '{name}' produced a window of a different length "
                f"({len(perturbed.table)} against {len(control.table)}); a driver "
                "perturbation must preserve the number of steps"
            )
        for column in fixed:
            a = np.asarray(control.forcing[column], dtype=float)
            b = np.asarray(perturbed.forcing[column], dtype=float)
            if not np.allclose(a, b, rtol=0, atol=1e-9):
                raise ValueError(
                    f"'{column}' differs between control and '{name}'; only '{driver}' "
                    "may change, or the response cannot be attributed"
                )

        delta_driver = float(
            perturbed.volume(perturbed.forcing[driver]).sum()
            - control.volume(control.forcing[driver]).sum()
        )
        if abs(delta_driver) < 1e-9:
            raise ValueError(
                f"variant '{name}' does not change {driver}; the generator has to "
                "move the driver for the response to be attributable"
            )
        pushed = "more" if delta_driver > 0 else "less"

        shares: dict[str, float] = {}
        for var, direction in expect.items():
            if var not in control.table.columns or var not in perturbed.table.columns:
                continue
            change = float(
                perturbed.volume(perturbed.table[var]).sum()
                - control.volume(control.table[var]).sum()
            )
            share = change / delta_driver
            shares[var] = share
            signed = -share if direction == "decrease" else share
            moved = "fell" if change < 0 else ("rose" if change > 0 else "did not move")
            if signed < min_share:
                failures.append(
                    f"{name} ({pushed} {driver}): {var} {moved} by {abs(change):.1f} mm, "
                    f"{share:+.3f} per unit of demand; a {direction} of at least "
                    f"{min_share:g} per unit is expected"
                )
            elif signed > max_share:
                failures.append(
                    f"{name} ({pushed} {driver}): {var} moved {abs(share):.3f} per unit "
                    f"of demand, more than the {max_share:g} the perturbation can "
                    "account for"
                )
            if var == lead and (weakest is None or signed < weakest[0]):
                weakest = (signed, share)
        if not shares:
            raise ValueError(
                f"response_sign found none of {list(expect)} in the model result"
            )
        per_variant[name] = {"driver_change_mm": delta_driver, "shares": shares}

    ok = not failures
    detail = "; ".join(
        f"{name} ({'+' if v['driver_change_mm'] > 0 else ''}{v['driver_change_mm']:.0f} mm "
        f"of {driver}): "
        + ", ".join(f"{var} {share:+.3f}" for var, share in v["shares"].items())
        for name, v in per_variant.items()
    )
    both_ways = " in both directions" if len(per_variant) > 1 else ""
    return CriterionResult(
        name="response_sign",
        status=PASS if ok else FAIL,
        value=None if weakest is None else weakest[1],
        threshold=min_share,
        message=(
            f"the response to {driver} points the right way{both_ways} "
            f"(per unit of demand: {detail})"
            if ok
            else "; ".join(failures)
        ),
        diagnostics={
            "driver": driver,
            "variants": per_variant,
        },
    )

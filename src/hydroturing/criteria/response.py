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

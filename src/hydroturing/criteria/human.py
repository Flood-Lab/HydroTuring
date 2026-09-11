"""Human abstraction: water the model was told to remove must leave the budget.

Irrigation withdrawal is the oldest unaccounted sink in hydrology. A model
that never reads the prescribed abstraction runs the same with it or without
it — and every single-run criterion passes, because nothing leaked: the water
was never taken. The only way to catch the omission is to run the same
weather twice, once with the human term and once without, and ask where the
abstracted water went.

The driver is prescribed in the forcing (a visible `abstr` column, mm/day),
not added to the rain, which is what separates this from
`counterfactual_response`: the perturbed variant does not receive more water,
it is told to lose some.

The three failure modes this separates:

  abstraction ignored     the model never reads the driver, so both runs are
                          identical and the residual is the whole abstracted
                          volume, not a fraction of it
  unreported sink         the model removes the water from nothing it
                          reports, so each run may close while the difference
                          does not move by the prescribed amount
  inconsistent removal    the model removes the right total on one seed and
                          the wrong total on the next; every seed must pass
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


@criterion("human_abstraction", paired=True)
def human_abstraction(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    """The budget difference between the paired runs must equal the
    prescribed abstraction, neither more nor less."""
    driver = str(params.get("driver", "abstr"))
    terms = list(params.get("terms", ["evspsbl", "mrro"]))
    tolerance = float(params.get("tolerance", 0.05))
    floor_mm = float(params.get("floor_mm", 5.0))
    weather = list(params.get("weather", ["pr", "tas", "pet"]))

    control = make_window(pick(runs, params, "control", "natural"), probe)
    perturbed = make_window(pick(runs, params, "perturbed", "irrigated"), probe)

    if driver not in perturbed.forcing.columns:
        raise ValueError(
            f"human_abstraction needs a '{driver}' column in the forcing of the "
            "perturbed variant; the generator has to prescribe the abstraction"
        )
    if len(control.table) != len(perturbed.table):
        raise ValueError(
            f"the variants produced windows of different lengths "
            f"({len(control.table)} against {len(perturbed.table)}); prescribing "
            "an abstraction must not change the number of steps"
        )

    # Only the human term may differ between the variants. If the weather
    # moved too, the budget difference can no longer be attributed to the
    # abstraction.
    for column in weather:
        if column not in control.forcing.columns:
            continue
        a = np.asarray(control.forcing[column], dtype=float)
        b = np.asarray(perturbed.forcing[column], dtype=float)
        if not np.allclose(a, b, rtol=0, atol=1e-9):
            raise ValueError(
                f"'{column}' differs between control and perturbed; only "
                f"'{driver}' may change, or the residual cannot be attributed"
            )

    abstracted = float(perturbed.volume(perturbed.forcing[driver]).sum())
    if driver in control.forcing.columns:
        abstracted -= float(control.volume(control.forcing[driver]).sum())
    if abstracted <= 0:
        raise ValueError(
            f"the perturbed variant prescribes no {driver} ({abstracted:.4g} mm); "
            "the generator has to make the two variants actually differ"
        )

    deltas: dict[str, float] = {}
    for var in terms:
        if var not in control.table.columns or var not in perturbed.table.columns:
            raise ValueError(f"human_abstraction needs '{var}' in the model result")
        deltas[var] = float(
            perturbed.volume(perturbed.table[var]).sum()
            - control.volume(control.table[var]).sum()
        )

    states = reported_states(control, probe)
    if states:
        d_perturbed = float(perturbed.storage(states)[-1]) - perturbed.storage_initial(states)
        d_control = float(control.storage(states)[-1]) - control.storage_initial(states)
        deltas["storage"] = d_perturbed - d_control

    # The budget terms, taken as losses from the catchment, must account for
    # exactly the prescribed removal: outflows fall and/or storage ends lower
    # by sum(A), so their signed difference plus sum(A) is zero.
    residual = float(sum(deltas.values())) + abstracted
    relative = abs(residual) / abstracted

    # One effective limit: the relative tolerance, or the absolute floor
    # expressed as a share of the abstracted volume when that floor is the more
    # permissive of the two. Reporting it as `threshold` keeps value and
    # threshold on one scale, so a pass never reads as value above threshold.
    effective = max(tolerance, floor_mm / abstracted)
    ok = relative <= effective
    detail = ", ".join(f"{k} {v:+.1f} mm" for k, v in deltas.items())
    return CriterionResult(
        name="human_abstraction",
        status=PASS if ok else FAIL,
        value=relative,
        threshold=effective,
        message=(
            f"the prescribed {abstracted:.0f} mm of {driver} left the budget "
            f"({detail}); residual {relative:.2%} of it (limit {effective:.2%})"
            if ok
            else f"of the prescribed {abstracted:.0f} mm of {driver}, the budget "
            f"difference accounts for all but {residual:+.1f} mm "
            f"({relative:.2%}, limit {effective:.2%}): {detail}"
        ),
        diagnostics={
            "abstracted_mm": abstracted,
            "residual_mm": residual,
            "floor_mm": floor_mm,
            "effective_limit": effective,
            "deltas_mm": {k: float(v) for k, v in deltas.items()},
        },
    )

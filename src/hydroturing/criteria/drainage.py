"""Recession drainage: a reach that is not being fed cannot fill.

A channel store is water generated as runoff and not yet released. When it
rains, the store fills; when the rain stops, the only thing it can do is
drain. That gives a statement about the store that no water budget makes,
because the budget is differenced over *every* store at once and a reach
that quietly manufactures water inside its routing is invisible to it.

The statement is deliberately one-sided and weak, so that every honest
routing passes it:

    on a step where no new water enters the reach, the store must not rise.

"No new water enters" is read off the forcing, not off the model. A step is
in recession when it and the few steps before it carried no rain
(`settle_days`). Rain that fell earlier may still be draining through the
reach, and that drainage is *exactly* what the criterion allows: an honest
store falls while it releases what it holds. What it may not do is rise.

The measure is a **frequency**, not a worst case. An honest model can show
a single small rise: a discrete unit hydrograph renormalised over a partial
history, or the first step of a record where the hydrograph is still
filling from a zero state, moves the store by a fraction of a percent of
its own scale. Taking the worst step would fail such a model on one
startup sample. Taking the *share of recession steps on which the store
rose* separates the two populations cleanly instead: an honest router sits
near zero, a router that loses or invents water rises on nearly every
recession step, because the error is applied every step and nothing is
entering to mask it.

A rise is only counted when it is larger than a relative floor, so that the
floating-point dust of a per-step renormalisation is not scored as a rise.
The floor scales with the store's own maximum, so it reads as a fraction of
the reach rather than as an absolute depth.

The store is `channel`. A model that does not route holds nothing in transit
and reports a channel of zero at every step; a constant store never rises
and passes trivially, which is correct: the criterion asserts only that a
store which exists must not fill from nothing.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

# The default store the criterion watches. `channel` is the reach's water in
# transit, and it is the one the roadmap entry this probe claims is about.
DEFAULT_STORE = "channel"
# A step is in recession when this many consecutive steps, counting the step
# itself, carried no rain. Three days lets the fastest surface response of a
# typical event clear before the step is scored, so a small rise caused by the
# tail of a storm that has not quite finished draining is not counted.
DEFAULT_SETTLE_DAYS = 3
# Rain below this rate, in mm/day, counts as no rain. The generators produce
# exactly zero on dry days, so this only forgives rounding.
DRY_PR_MM_PER_DAY = 0.05
# A rise has to clear this share of the store's own maximum before it counts.
# An honest router renormalising a discrete kernel over a partial history
# moves the store by ~1e-5 of its scale; a real leak moves it by ~1e-3. The
# floor sits between them with an order of magnitude on each side.
DEFAULT_RISE_FRACTION = 1e-4
# The share of recession steps that may rise before the store is called
# self-filling. Honest routers sit under 0.04; a router that loses a fixed
# share of every step is near 0.9, since the loss is applied on every step and
# nothing enters to hide it. A quarter sits with an order of magnitude of room
# on the honest side and a factor of three on the failing one.
DEFAULT_MAX_RISING_FRACTION = 0.25
# The first steps of the scored window are not scored: a unit hydrograph
# starting from a zero store fills for the first few steps whatever the
# weather, and that is arithmetic, not a violation. Excluding them keeps the
# frequency honest without weakening the test anywhere the model is settled.
DEFAULT_SETTLE_STEPS = 30


@criterion("recession_drainage")
def recession_drainage(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """On recession steps the reach store must not rise.

    Reports the share of recession steps on which the store rose more than a
    relative floor, and fails when that share passes `max_rising_fraction`.
    """
    store_var = str(params.get("store", DEFAULT_STORE))
    settle_days = int(params.get("settle_days", DEFAULT_SETTLE_DAYS))
    rise_fraction = float(params.get("rise_fraction", DEFAULT_RISE_FRACTION))
    max_fraction = float(params.get("max_rising_fraction", DEFAULT_MAX_RISING_FRACTION))
    settle_steps = int(params.get("settle_steps", DEFAULT_SETTLE_STEPS))

    w = make_window(run, probe)
    # A missing column is a failure, not a crash and not a pass: the criterion
    # cannot make its statement, and one that cannot be scored must not be
    # reported as satisfied. The harness normally short-circuits a model that
    # does not report the state before any criterion runs, so this is a second
    # line of defence rather than the first.
    if store_var not in w.table.columns:
        return CriterionResult(
            name="recession_drainage",
            status=FAIL,
            message=(
                f"the run reports no '{store_var}' column, so there is no store "
                "to watch on the recession steps"
            ),
        )
    if "pr" not in w.forcing.columns:
        return CriterionResult(
            name="recession_drainage",
            status=FAIL,
            message=(
                "the forcing carries no 'pr' column, so a recession step cannot "
                "be located"
            ),
        )

    store = np.asarray(w.table[store_var], dtype=float)
    pr = np.asarray(w.forcing["pr"], dtype=float)
    n = len(store)

    # A step is in recession when it and the `settle_days` steps before it were
    # all rainless. The comparison is against the forcing, so a model cannot
    # declare itself in recession to dodge the test.
    dry = pr <= DRY_PR_MM_PER_DAY
    recession = np.copy(dry)
    for k in range(1, settle_days + 1):
        recession[k:] &= dry[:-k]
    # The startup transient of the scored window is excluded, not scored.
    recession[: max(0, settle_steps)] = False

    if not recession.any():
        return CriterionResult(
            name="recession_drainage",
            status=FAIL,
            message=(
                "the record contains no recession steps to score "
                f"(no rainless run of {settle_days + 1} steps outside the first "
                f"{settle_steps}); the case is degenerate for this criterion"
            ),
            diagnostics={"recession_steps": 0},
        )

    change = np.diff(store, prepend=store[0])
    scale = float(np.max(np.abs(store))) if n else 0.0
    floor = rise_fraction * max(scale, 1e-12)
    rising = change[recession] > floor
    fraction = float(rising.mean())
    worst = float(change[recession].max())
    n_rising = int(rising.sum())
    n_scored = int(recession.sum())

    ok = fraction <= max_fraction
    return CriterionResult(
        name="recession_drainage",
        status=PASS if ok else FAIL,
        value=fraction,
        threshold=max_fraction,
        message=(
            f"the {store_var} store drains on recession steps "
            f"({n_rising} of {n_scored} rise above {floor:.2e} mm, "
            f"worst {worst:+.4f} mm)"
            if ok
            else (
                f"the {store_var} store fills from nothing on "
                f"{fraction:.1%} of recession steps ({n_rising} of {n_scored}, "
                f"limit {max_fraction:.0%}, worst rise {worst:+.4f} mm); a reach "
                f"with no water entering it cannot rise"
            )
        ),
        diagnostics={
            "store": store_var,
            "recession_steps": n_scored,
            "rising_steps": n_rising,
            "rising_fraction": fraction,
            "worst_rise_mm": worst,
            "rise_floor_mm": floor,
            "store_max_mm": scale,
        },
    )

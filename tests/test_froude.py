"""Counterexamples for the Froude criterion, on cross-sections built by hand.

`froude_subcritical` reads a `stage` diagnostic against a discharge and asks
whether the pair is one the declared cross-section could carry. The gate
exercises it on four honest gauges and one whose rating is drawn for a wider
reach, but a gate only ever drives the paths its baselines happen to take, and
three things about this criterion are easy to get wrong while staying green:

* the discharge can come from `dis` *or* from `mrro` over the catchment area,
  and the fallback is the only path a model that reports no `dis` can take —
  a fallback that forgets the 86400 would still fail a wrong gauge, so the
  gate cannot tell it apart from a correct one;
* the width of the section is a static attribute, not a constant, so a
  criterion that hard-codes the width passes every gate while refusing to
  judge a model in the channel the case declared;
* a step is excused at one end of the depth range and refused at the other, and
  both doors have to be read from the flow rather than from the depth: a step
  the reach carries no water through is not scored, while a *shallow* step
  carrying water is the fastest in the record and is scored at the floor —
  excusing it is how a model buys a pass by reporting its flood peaks as
  sub-centimetre, and the direction of the error is the tell, because making
  those peaks *more* wrong moves the verdict from a failure to a pass;
* `stage` is an elevation on the case's fixed vertical datum, so the depth is
  `stage - bed_elevation_m` and the criterion may not infer a datum: a case
  that declares none is refused, a reading taken from another zero reaches the
  harness as N/A (INCOMPATIBLE) rather than as a violation, a non-finite
  reading is refused rather than quietly dropped out of the mask, and shifting
  the datum and the level together must leave the verdict exactly where it was.

These tests build the frames directly so each of those is pinned: a
Manning-consistent pair, a rating drawn for a wider reach, the `mrro`
fallback, a width that changes the verdict, dry steps, the depth ceiling and
the two tolerance knobs, and the columns the criterion cannot do without —
missing ones must fail it rather than crash it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing import harness, registry
from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS, CriterionIncompatibleError
from hydroturing.harness import build_case, compatibility_issues, run_probe
from hydroturing.protocol import Case, RunResult
from hydroturing.scoring import INCOMPATIBLE, NOT_SCORED
from hydroturing.seeds import gate_seeds

GRAVITY = 9.81
WIDTH_M = 18.0
SLOPE = 0.0015
MANNING_N = 0.035
AREA_KM2 = 250.0
SECONDS_PER_DAY = 86400.0
# The fixed vertical datum the case declares. `stage` is a water-surface
# elevation on it, so a depth the tests reason about is reported above it.
BED_M = 120.0
FROUDE = get("froude_subcritical")


def build(q=None, stage=None, mrro=None, static=None, n=None, discharge=True, spinup=0):
    """A RunResult carrying the columns the Froude criterion reads.

    `q` is `dis` in m3/s and `stage` is the flow **depth** in metres, because
    that is what these tests reason about. The frame is built with the case
    declaring `BED_M` and the depth reported above it, because the contract
    makes `stage` a water-surface elevation on the case's fixed datum and the
    criterion forms the depth as `stage - bed_elevation_m`. Passing `mrro`
    instead of `q` and `discharge=False` builds the frame a model that reports
    only a runoff depth rate produces, which is the frame the fallback has to
    score.
    """
    if q is None:
        q = np.linspace(1.0, 50.0, n or 400)
    q = np.asarray(q, dtype=float)
    n_steps = len(q)
    if stage is None:
        stage = manning_depth(q)
    depth = np.broadcast_to(np.asarray(stage, dtype=float), (n_steps,))
    table = {"time": pd.date_range("2001-01-01", periods=n_steps, freq="D"),
             "stage": BED_M + depth}
    if discharge:
        table["dis"] = q
    if mrro is not None:
        table["mrro"] = np.broadcast_to(np.asarray(mrro, dtype=float), (n_steps,))
    case = Case(
        probe_id="momentum/froude-regime", seed=1,
        forcing=pd.DataFrame({"time": table["time"]}),
        static=dict(static or {}, width_m=WIDTH_M, area_km2=AREA_KM2,
                    slope=SLOPE, manning_n=MANNING_N, bed_elevation_m=BED_M),
        spinup_steps=spinup, timestep="PT1D",
    )
    return RunResult(case, pd.DataFrame(table), {}, 0.0)


def manning_depth(q, width_m=WIDTH_M):
    """The depth a reach of `width_m` at `SLOPE` needs to carry `q`.

    This is the honest gauge when the width is the declared one: the depth and
    the discharge are two readings of one cross-section, so given one the
    other follows, and the Froude number of the pair is a property of the
    section rather than of the flow. Passing a wider section is the same law
    applied to a reach that is not the one the case declares.
    """
    q = np.maximum(np.asarray(q, dtype=float), 1e-9)
    return (q * MANNING_N / (width_m * SLOPE**0.5)) ** 0.6


def test_a_manning_consistent_pair_is_subcritical_at_every_flow():
    """A gauge that reads the flow it is carrying never outruns its own wave."""
    result = FROUDE(build(), None, {})
    assert result.status == PASS
    assert result.diagnostics["max_froude"] < 1.0
    assert result.diagnostics["exceeding_steps"] == 0
    # Well under the limit, not scraping it: on a mild slope the honest pair
    # sits near 0.4 at every flow the case produces.
    assert result.diagnostics["max_froude"] < 0.6
    assert "stays subcritical" in result.message


def test_a_gauge_pinned_at_one_depth_goes_supercritical():
    """Ten centimetres whatever the flow: the velocity absorbs all of it."""
    q = np.linspace(1.0, 50.0, 400)
    result = FROUDE(build(q=q, stage=0.10), None, {})
    assert result.status == FAIL
    assert result.diagnostics["exceed_fraction"] > 0.5
    assert result.diagnostics["max_froude"] > 1.05
    assert "supercritical" in result.message


def test_the_pinned_gauge_fails_only_where_the_flow_is_large():
    """The violation is on the rising limb, not at baseflow.

    A pinned depth is plausible at low flow — a real reach is ten centimetres
    deep much of the time — so a criterion that failed every step would be
    reading the depth rather than the pair.
    """
    q = np.full(400, 0.05)  # baseflow: 0.05 m3/s through 18 m is very slow
    result = FROUDE(build(q=q, stage=0.10), None, {})
    assert result.status == PASS
    assert result.diagnostics["max_froude"] < 1.0


def test_a_rating_drawn_for_a_wider_reach_goes_supercritical():
    """The right law on the wrong section: depth is what velocity divides by.

    This is the shape the must-fail model has, and the shape the criterion has
    to catch alone. The pair is derived from the flow rather than pinned, so
    nothing about the frame is degenerate — it is the rating of a reach five
    times wider than the one declared, which reports a depth `5**0.6` times
    too shallow and so a Froude number `5**0.9` (~4.3) times too large. The
    honest pair over the same flow is subcritical at every step.
    """
    q = np.linspace(1.0, 50.0, 400)
    honest = FROUDE(build(q=q), None, {})
    wider = FROUDE(build(q=q, stage=manning_depth(q, width_m=5.0 * WIDTH_M)), None, {})
    assert honest.status == PASS
    assert wider.status == FAIL
    assert wider.diagnostics["exceed_fraction"] > 0.8
    assert wider.diagnostics["max_froude"] > 1.05


def test_the_wider_rating_is_not_degenerate_and_not_non_monotone():
    """Why no other criterion in the suite can see the wider-reach rating.

    A gauge that does not move is caught by `non_degenerate` on
    `momentum/stage-discharge-monotonic`, so the must-fail could not be that
    and demonstrate anything about this criterion. This one varies exactly as
    much as the honest gauge — scaling a depth leaves its coefficient of
    variation alone — and rises with the flow, so every variability and
    monotonicity test passes while the pair is supercritical.
    """
    q = np.linspace(1.0, 50.0, 400)
    honest = manning_depth(q)
    wider = manning_depth(q, width_m=5.0 * WIDTH_M)
    assert np.all(np.diff(wider) >= 0.0)
    assert wider.std() / wider.mean() == pytest.approx(honest.std() / honest.mean())
    assert wider.std() / wider.mean() > 0.1


def test_a_model_that_reports_only_runoff_is_scored_through_the_area():
    """The `mrro` fallback has to convert, not merely rescale.

    A depth rate over the catchment is a volume rate, and the conversion
    carries the 86400 seconds of a day. Dropping it makes every discharge
    5 orders of magnitude too large, which would fail this frame — so the
    honest gauge is what separates a correct fallback from a wrong one.
    """
    mrro = np.linspace(0.5, 12.0, 400)
    result = FROUDE(build(mrro=mrro, discharge=False), None, {})
    assert result.status == PASS
    assert result.diagnostics["discharge_source"] == "mrro"
    assert result.diagnostics["max_froude"] < 0.6


def test_the_fallback_and_the_declared_discharge_agree():
    """The two sources are interchangeable, not merely comparable."""
    q = np.linspace(1.0, 50.0, 400)
    mrro = q * SECONDS_PER_DAY / (1e-3 * AREA_KM2 * 1e6)
    from_dis = FROUDE(build(q=q), None, {})
    from_mrro = FROUDE(build(mrro=mrro, discharge=False), None, {})
    assert from_dis.diagnostics["discharge_source"] == "dis"
    assert from_mrro.diagnostics["discharge_source"] == "mrro"
    assert from_dis.diagnostics["max_froude"] == pytest.approx(
        from_mrro.diagnostics["max_froude"], rel=1e-9)


def test_the_width_of_the_section_is_read_from_the_case_not_assumed():
    """A wider section carries the same flow more slowly, and that decides it.

    A criterion that hard-codes the width would fail both of these, because
    the pair is identical and only the channel differs.
    """
    q = np.linspace(1.0, 50.0, 400)
    pinned = dict(q=q, stage=0.10)
    narrow = FROUDE(build(**pinned), None, {"width_m": WIDTH_M})
    wide = FROUDE(build(**pinned), None, {"width_m": 1000.0})
    assert narrow.status == FAIL
    assert wide.status == PASS
    assert wide.diagnostics["width_m"] == 1000.0


def test_a_dry_reach_has_no_regime_to_read():
    """Steps that only clear the depth floor are not scored.

    A receding flow reaches the floor by arithmetic, not by hydraulics, and
    the velocity computed on it is an artefact of the floor.
    """
    q = np.full(400, 0.0)
    result = FROUDE(build(q=q, stage=0.0), None, {})
    assert result.status == FAIL
    assert result.diagnostics["scored_steps"] == 0
    assert "degenerate" in result.message


def test_a_step_with_no_flow_is_still_excused():
    """Dryness is read from the flow, so a step with no water in it leaves the
    sample.

    The floor is a floor on the division rather than a dryness test. A step the
    reach carries nothing through has no velocity to read in it, and scoring it
    would fail an honest receding flow on arithmetic instead of on its
    hydraulics — so the one thing that still excuses a step is the absence of
    flow, which is what keeps this from failing every record at baseflow.
    """
    q = np.linspace(1.0, 50.0, 400)
    stage = manning_depth(q)
    stage[:200] = 0.0  # the first half of the record is dry
    q = q.copy()
    q[:200] = 0.0      # ... and carries nothing, which is what makes it dry
    result = FROUDE(build(q=q, stage=stage), None, {})
    assert result.status == PASS
    assert result.diagnostics["scored_steps"] == 200


def _runoff_for(q):
    """The `mrro` in mm/day that is the same outflow as `q` in m3/s."""
    return np.asarray(q, dtype=float) * SECONDS_PER_DAY / (AREA_KM2 * 1e6) * 1e3


def test_a_zero_in_dis_with_runoff_behind_it_is_scored_on_the_runoff():
    """A step is still only when every flow the model reports says so.

    `dis` is preferred, so with the scored set decided on it alone a model that
    reports both columns picks which steps are scored: the wider-reach rating,
    with `dis` zeroed on exactly the steps it goes supercritical on and `mrro`
    still carrying the water, passed — on the probe's own gate seeds too, where
    it left 144 of 1460 steps scored. Under the contract the two columns are two
    readings of one outflow, so the silent steps are scored on the runoff.
    """
    # Baseflow low enough that even the wider rating stays subcritical, as it
    # does on the gate record, then the floods it fails on.
    q = np.concatenate([np.full(60, 0.2), np.linspace(1.0, 50.0, 340)])
    stage = manning_depth(q, width_m=5 * WIDTH_M)
    fr = q / (WIDTH_M * stage ** 1.5 * GRAVITY ** 0.5)
    lying = q.copy()
    lying[fr > 1.05] = 0.0
    n_zeroed = int((lying == 0.0).sum())
    assert n_zeroed == 340  # every flood step, which is the point

    result = FROUDE(build(q=lying, stage=stage, mrro=_runoff_for(q)), None, {})
    assert result.status == FAIL
    assert result.diagnostics["scored_steps"] == 400
    assert result.diagnostics["second_reading"] == "mrro"
    assert result.diagnostics["steps_read_from_second_reading"] == n_zeroed
    assert f"{n_zeroed} of the scored steps report no 'dis'" in result.message

    # With only `dis` to read, the same record is bounded by the floor alone:
    # the criterion's own default lets it through, and the floor this probe
    # declares refuses it as degenerate rather than passing it.
    dis_only = build(q=lying, stage=stage)
    assert FROUDE(dis_only, None, {}).status == PASS
    declared = FROUDE(dis_only, None, {"min_scored_fraction": 0.9})
    assert declared.status == FAIL
    assert "dry or flat" in declared.message


def test_a_step_both_readings_call_still_is_excused():
    """The rule reads both columns and changes nothing where they agree: a
    step with no flow in either is still dry, and an honest record that reports
    both is scored exactly as it was, with no clause added to its message."""
    q = np.linspace(1.0, 50.0, 400)
    stage = manning_depth(q)
    q[:200] = 0.0
    runoff = _runoff_for(q)
    result = FROUDE(build(q=q, stage=stage, mrro=runoff), None, {})
    assert result.status == PASS
    assert result.diagnostics["scored_steps"] == 200
    assert result.diagnostics["steps_read_from_second_reading"] == 0
    assert "report no 'dis'" not in result.message


def test_a_bad_value_in_the_reading_a_silent_step_falls_back_to_is_refused():
    """The non-finite refusal covers the flow as it is divided by. A `NaN` in
    `mrro` behind a zero in `dis` would otherwise read as stillness and leave the
    scored set by the door the refusal exists to close."""
    q = np.linspace(1.0, 50.0, 400)
    stage = manning_depth(q)
    runoff = _runoff_for(q)
    q[150] = 0.0
    runoff[150] = np.nan
    result = FROUDE(build(q=q, stage=stage, mrro=runoff), None, {})
    assert result.status == FAIL
    assert result.diagnostics["non_finite_steps"] == 1


def test_the_probe_declares_a_scored_floor_the_honest_models_clear():
    """The floor is what bounds a model that writes a zero into every flow
    column it reports. Ninety per cent: the must-pass baselines score every step
    on the gate seeds and at least 96% over seeds 0-99."""
    from hydroturing.registry import find_probe

    probe = find_probe("momentum/froude-regime")
    (froude,) = [c for c in probe.criteria if c.name == "froude_subcritical"]
    assert froude.params["min_scored_fraction"] == 0.9


def test_a_shallow_step_carrying_water_is_scored_not_excused():
    """The door the `NaN` fix did not close: depth decides dryness, so a model
    that reports its worst steps as sub-centimetre takes them out of the sample.

    A step with discharge behind it and a depth under the floor is not a dry
    reach — at a given discharge it is the *fastest* step the record can hold,
    because `Fr` grows without bound as the depth falls. The verdict was
    therefore monotone in the wrong direction: peaks at a fifth of their honest
    depth failed, and the same model made more wrong still — peaks at five
    millimetres — passed, because the steps it would have failed on left the
    scored set while the remaining 70% of the record stayed comfortably above
    `min_scored_fraction` and nothing reported the gap.

    The discharge is untouched throughout, so what this pins is that the depth
    alone may not excuse a step.
    """
    q = np.linspace(1.0, 50.0, 400)
    stage = manning_depth(q)
    peaks = q >= np.quantile(q, 0.70)   # the 120 largest discharges
    stage[peaks] = 0.005                # sub-centimetre, discharge unchanged

    # The same record with its honest depths, so what fails is the depth.
    assert FROUDE(build(q=q, stage=manning_depth(q)), None, {}).status == PASS

    result = FROUDE(build(q=q, stage=stage), None, {})
    assert result.status == FAIL
    assert result.diagnostics["scored_steps"] == 400
    assert result.diagnostics["exceeding_steps"] == int(peaks.sum())
    assert result.diagnostics["max_froude"] > 100.0
    assert "supercritical" in result.message


def test_the_clamped_floor_is_what_a_shallow_flow_is_scored_at():
    """A flowing step has a number to divide by, and the number is the floor.

    `Fr = Q / (w d**1.5 sqrt(g))` with `d` at zero is not defined, so the floor
    is doing its job when a millimetre of reported depth is read as the one
    centimetre the division needs — and the result is enormous, which is the
    point: that pair is not a reach, it is a discharge with no water under it.
    """
    q = np.linspace(1.0, 50.0, 400)
    result = FROUDE(build(q=q, stage=np.full(400, 0.005)), None, {})
    assert result.status == FAIL
    assert result.diagnostics["scored_steps"] == 400
    floor_fr = q / (WIDTH_M * 0.01**1.5 * GRAVITY**0.5)
    assert result.diagnostics["max_froude"] == pytest.approx(float(floor_fr.max()))


def test_a_negative_reading_is_water_moving_rather_than_an_excuse():
    """A discharge is a magnitude here, so neither the mask nor the quotient may
    read the sign.

    Flipping the sign of the steps a model would fail on is the same move as
    reporting them shallow, one step further out; and a numerator that keeps the
    sign computes a negative quotient, which is under the limit at every step —
    a reading that cannot fail rather than one that passes.
    """
    q = np.linspace(1.0, 50.0, 400)
    honest_depth = manning_depth(q)
    result = FROUDE(build(q=-q, stage=honest_depth), None, {})
    assert result.diagnostics["scored_steps"] == 400

    # The same pair with the sign taken off is the frame every other test here
    # uses, so the magnitude is what the criterion reads.
    flipped = FROUDE(build(q=q, stage=honest_depth), None, {})
    assert result.diagnostics["max_froude"] == pytest.approx(
        flipped.diagnostics["max_froude"])
    assert result.status == flipped.status == PASS

    # And a pinning that only a positive quotient would notice: the depth is a
    # tenth of a metre at every flow, so the pair is supercritical on most of the
    # record whether the discharge is reported up or down.
    pinned = np.full(400, 0.10)
    up = FROUDE(build(q=q, stage=pinned), None, {})
    down = FROUDE(build(q=-q, stage=pinned), None, {})
    assert up.status == down.status == FAIL
    assert down.diagnostics["max_froude"] == pytest.approx(
        up.diagnostics["max_froude"])
    assert down.diagnostics["max_froude"] > 1.05


def test_a_stage_on_the_wrong_datum_is_a_convention_mismatch():
    """The depth is `stage - bed_elevation_m`, so a model that reports the depth
    itself — the convention this probe used to carry — is a convention mismatch
    rather than a deep reach.

    `depth_series` raises `CriterionIncompatibleError` for it and the criterion
    must let that exception through: the classification is the harness's to
    make, and a criterion that caught it would report a conservation violation
    for a model that merely reported a different quantity. The verdict that
    classification produces, end to end, is pinned by the test below.
    """
    q = np.linspace(1.0, 50.0, 400)
    run = build(q=q, stage=manning_depth(q))
    # Report the depth itself, as a model that never read the declared datum would.
    run.table["stage"] = run.table["stage"] - BED_M
    with pytest.raises(CriterionIncompatibleError, match="not depth above the bed"):
        FROUDE(run, None, {})


def test_the_wrong_datum_control_is_incompatible_rather_than_a_violation(monkeypatch):
    """The same reading, taken through the harness, has to arrive as N/A.

    A criterion that swallowed the exception and returned a failed result would
    put a conservation violation in the archive for a control that reported a
    level on another zero, and would take those seeds out of the harness's own
    account of how many of them were scored — which is the machinery that lets a
    genuinely incompatible model be recorded as incompatible and a genuinely
    failing one as a failure.
    """
    probe = registry.find_probe("momentum/froude-regime")
    model = registry.find_model("reference_bucket")
    real = harness.get_runner(model)

    class DepthAsStage:
        """`reference_bucket`'s own run, with the declared bed taken back out."""

        def run(self, model, probe, case, io_dir):
            run = real.run(model, probe, case, io_dir)
            table = run.table.copy()
            table["stage"] = table["stage"] - float(case.static["bed_elevation_m"])
            return RunResult(
                case=run.case, table=table, meta=run.meta, wall_seconds=0.0
            )

    monkeypatch.setattr(harness, "get_runner", lambda _model: DepthAsStage())
    outcome = run_probe(model, probe, gate_seeds(probe.id, 3))

    assert (outcome.verdict, outcome.reason) == (NOT_SCORED, INCOMPATIBLE)
    assert outcome.criteria == []
    assert len(outcome.incompatible) == 3
    assert "not depth above the bed" in outcome.incompatible[0]


def test_a_non_finite_reading_cannot_leave_the_scored_set():
    """`NaN` compares false against both bounds, so the masks would drop the step
    rather than score it: one step fewer, a share measured on what is left, and
    nothing said about the value that went missing.

    A criterion that lets that happen reports a pass on a record it could not
    read, so the pair is refused and the count of bad steps is reported instead.
    """
    q = np.linspace(1.0, 50.0, 400)
    stage = manning_depth(q)
    peak = int(np.argmax(q))  # the step a mask is most likely to lose unnoticed

    broken_q = q.copy()
    broken_q[peak] = np.nan
    result = FROUDE(build(q=broken_q, stage=stage), None, {})
    assert result.status == FAIL
    assert result.diagnostics["non_finite_steps"] == 1
    assert "non-finite" in result.message

    # The depth enters the other half of the same mask, so it is the same hole.
    broken_stage = stage.copy()
    broken_stage[peak] = np.nan
    result = FROUDE(build(q=q, stage=broken_stage), None, {})
    assert result.status == FAIL
    assert result.diagnostics["non_finite_steps"] == 1

    # The same pair without the bad value passes, so what failed is the value
    # and not the record.
    assert FROUDE(build(q=q, stage=stage), None, {}).status == PASS


def test_a_case_that_declares_no_datum_is_refused():
    """The criterion may not assume a datum, so a case with no bed is refused
    rather than scored against zero — defaulting the bed to zero would let a
    model that never read the key score exactly as one that did."""
    run = build()
    del run.case.static["bed_elevation_m"]
    result = FROUDE(run, None, {})
    assert result.status == FAIL
    assert "bed elevation" in result.message


def test_shifting_the_datum_moves_the_level_and_not_the_verdict():
    """A fixed datum cancels from a difference: raising the bed and the level
    together leaves the depth, and the verdict, exactly where they were."""
    q = np.linspace(1.0, 50.0, 400)
    base = FROUDE(build(q=q), None, {})
    run = build(q=q)
    run.case.static["bed_elevation_m"] = BED_M + 250.0
    run.table["stage"] = run.table["stage"] + 250.0
    shifted = FROUDE(run, None, {})
    assert shifted.status == base.status == PASS
    assert shifted.diagnostics["max_froude"] == pytest.approx(
        base.diagnostics["max_froude"]
    )
    assert shifted.diagnostics["scored_steps"] == base.diagnostics["scored_steps"]


def test_the_ceiling_leaves_a_depth_the_section_can_hold_alone():
    """It refuses what is not a depth, not what is merely deep."""
    q = np.linspace(1.0, 50.0, 400)
    result = FROUDE(build(q=q), None, {"max_depth_m": 10.0})
    assert result.status == PASS
    assert result.diagnostics["too_deep_steps"] == 0
    assert result.diagnostics["max_depth_m"] == 10.0
    assert result.diagnostics["scored_steps"] == 400


def test_the_ceiling_still_refuses_where_the_floor_no_longer_does():
    """The floor left the mask; the ceiling stays in it.

    Both bounds used to be applied the same way — to the depth — and they are
    not the same kind of statement. A shallow step is a real step and the floor
    only gives the division a number; a step deeper than the section can hold is
    not a reading of that section at all. Clamping it the way the floor clamps
    would put the offset back inside the range the bound exists to exclude, and
    an offset only ever adds to the depth, so a clamped offset is a comfortable
    pass rather than a refusal.
    """
    q = np.linspace(1.0, 50.0, 400)
    offset = manning_depth(q) + 25.0
    refused = FROUDE(build(q=q, stage=offset), None, {"max_depth_m": 10.0})
    assert refused.status == FAIL
    assert refused.diagnostics["scored_steps"] == 0
    assert refused.diagnostics["too_deep_steps"] == 400
    assert "ceiling" in refused.message

    # The same depth with the ceiling lifted is the pass a clamp would have
    # bought: an offset only adds to the depth, and a deeper depth is a smaller
    # Froude number.
    quiet = FROUDE(build(q=q, stage=offset), None, {})
    assert quiet.status == PASS
    assert quiet.diagnostics["max_froude"] < 0.1


def test_there_is_no_ceiling_unless_the_probe_sets_one():
    """How deep is too deep is a statement about the case's reach, not about
    `stage`, so the shared criterion does not guess at one."""
    q = np.linspace(1.0, 50.0, 400)
    result = FROUDE(build(q=q, stage=manning_depth(q) + 100.0), None, {})
    assert result.diagnostics["max_depth_m"] is None
    assert result.diagnostics["too_deep_steps"] == 0


def test_the_probe_declares_the_section_and_a_model_that_never_read_it_is_not_judged():
    """`requires.static` is what keeps a model carrying its own river geometry
    out of the archive, and nothing else in this file would notice if the
    declaration were dropped."""
    probe = registry.find_probe("momentum/froude-regime")
    assert probe.requires_static == ("width_m", "slope", "manning_n", "bed_elevation_m")
    assert probe.requires_diagnostics == ("stage",)

    case = build_case(probe, 0)
    flat = registry.find_model("reference_flat_stage")
    assert compatibility_issues(flat, probe, case) != []
    outcome = run_probe(flat, probe, gate_seeds(probe.id, 1))
    assert (outcome.verdict, outcome.reason) == (NOT_SCORED, INCOMPATIBLE)
    assert outcome.incompatible == [
        "model does not declare that it consumes static width_m, slope, manning_n"
    ]

    # The five baselines do declare all three, so the gate still separates in
    # the section the case declares.
    for name in ("reference_bucket", "flex_lumped", "flex_topo", "sacsma_snow17",
                 "reference_shallow_rating"):
        issues = compatibility_issues(registry.find_model(name), probe, case)
        assert issues == [], f"{name}: {issues}"


def test_the_tolerance_is_read_and_can_forgive_a_violation():
    """The margin covers the rectangular approximation, not a wrong velocity,
    so it has to be a parameter rather than a constant baked into the test."""
    q = np.linspace(1.0, 50.0, 400)
    strict = FROUDE(build(q=q, stage=0.10), None, {"tolerance": 0.05})
    lax = FROUDE(build(q=q, stage=0.10), None, {"tolerance": 1000.0})
    assert strict.status == FAIL
    assert lax.status == PASS


def test_the_exceed_fraction_limit_is_read():
    """Half the record supercritical is a verdict the probe can set either
    way, which is why the measure is a frequency and not a worst case."""
    q = np.linspace(1.0, 50.0, 400)
    stage = manning_depth(q)
    stage[200:] = 0.10  # the second half is pinned
    strict = FROUDE(build(q=q, stage=stage), None, {"max_exceed_fraction": 0.0})
    lenient = FROUDE(build(q=q, stage=stage), None, {"max_exceed_fraction": 0.6})
    assert strict.status == FAIL
    assert lenient.status == PASS
    assert 0.4 < strict.diagnostics["exceed_fraction"] < 0.6


def test_a_missing_stage_is_not_a_pass():
    """A missing column is a failure, not a crash and not a pass.

    The criterion cannot make its statement without a water surface, and one
    that cannot be scored must not be reported as satisfied.
    """
    run = build()
    del run.table["stage"]
    result = FROUDE(run, None, {})
    assert result.status == FAIL
    assert "stage" in result.message


def test_a_runoff_with_no_area_to_convert_it_is_not_a_pass():
    """The fallback needs the area the case declares; without it there is no
    discharge to compute a Froude number from, and silently reading the depth
    rate as a discharge would be a different error wearing the same name."""
    run = build(mrro=np.linspace(0.5, 12.0, 400), discharge=False)
    del run.case.static["area_km2"]
    result = FROUDE(run, None, {})
    assert result.status == FAIL
    assert "area_km2" in result.message


def test_a_model_that_reports_neither_discharge_nor_runoff_is_not_a_pass():
    result = FROUDE(build(discharge=False), None, {})
    assert result.status == FAIL
    assert "no flow" in result.message


def test_a_reach_with_no_width_is_not_a_pass():
    """The section the case declares has to exist for the question to be asked."""
    result = FROUDE(build(), None, {"width_m": 0.0})
    assert result.status == FAIL
    assert "width" in result.message


def test_a_case_that_declares_no_width_is_not_scored_in_an_invented_one():
    """The section is the case's to declare; with none there is nothing to read.

    Falling back to a width would score every model in a channel no case ever
    declared, and the message would name a number that appears nowhere in the
    run — which the archive cannot tell apart from a real violation.
    """
    run = build()
    del run.case.static["width_m"]
    result = FROUDE(run, None, {})
    assert result.status == FAIL
    assert "width_m" in result.message
    assert "invent" not in result.message  # it names the missing key, not a width


def test_turning_the_scored_floor_off_does_not_divide_by_zero():
    """`min_scored_fraction: 0.0` is how a probe disables the floor.

    `0 / n < 0.0` is False, so a guard that only tests the fraction falls
    through to `n_exceed / n_scored` with no scored steps and raises — an
    ERROR row for a broken adapter, when the cause is a dry record.
    """
    q = np.full(400, 0.0)
    result = FROUDE(build(q=q, stage=0.0), None, {"min_scored_fraction": 0.0})
    assert result.status == FAIL
    assert result.diagnostics["scored_steps"] == 0
    assert "degenerate" in result.message


def test_a_lenient_limit_does_not_report_that_no_step_was_supercritical():
    """The sentence is archived, so it may not state the opposite of the record.

    `max_exceed_fraction` is a knob a future case can raise, and the message
    has to carry the share in the passing branch too — otherwise a model that
    went supercritical on 5% of steps is recorded as having gone supercritical
    on none.
    """
    q = np.linspace(1.0, 50.0, 400)
    stage = manning_depth(q)
    stage[200:] = 0.10  # the second half is pinned
    result = FROUDE(build(q=q, stage=stage), None, {"max_exceed_fraction": 0.6})
    assert result.status == PASS
    assert result.diagnostics["exceeding_steps"] > 0
    assert "every scored step" not in result.message
    assert "supercritical" in result.message

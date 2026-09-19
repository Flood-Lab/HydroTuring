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
* steps where the reach is dry are excluded rather than scored, and a
  criterion that scores them fails an honest receding flow on the depth
  floor instead of on its hydraulics;
* `stage` is a depth above the reach bed by contract, and nothing in the
  criterion can convert a reading taken from another zero — an offset only
  adds to the depth, and the depth enters at the three-halves power, so a
  level a hundred metres above the bed passes comfortably unless the probe's
  ceiling refuses it.

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

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS
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
FROUDE = get("froude_subcritical")


def build(q=None, stage=None, mrro=None, static=None, n=None, discharge=True, spinup=0):
    """A RunResult carrying the columns the Froude criterion reads.

    `q` is `dis` in m3/s and `stage` is the reported depth in metres above the
    reach bed, which is what the contract says it is. Passing `mrro` instead of
    `q` and `discharge=False` builds the frame a model that reports only a
    runoff depth rate produces, which is the frame the fallback has to score.
    """
    if q is None:
        q = np.linspace(1.0, 50.0, n or 400)
    q = np.asarray(q, dtype=float)
    n_steps = len(q)
    if stage is None:
        stage = manning_depth(q)
    stage = np.broadcast_to(np.asarray(stage, dtype=float), (n_steps,))
    table = {"time": pd.date_range("2001-01-01", periods=n_steps, freq="D"),
             "stage": stage}
    if discharge:
        table["dis"] = q
    if mrro is not None:
        table["mrro"] = np.broadcast_to(np.asarray(mrro, dtype=float), (n_steps,))
    case = Case(
        probe_id="momentum/froude-regime", seed=1,
        forcing=pd.DataFrame({"time": table["time"]}),
        static=dict(static or {}, width_m=WIDTH_M, area_km2=AREA_KM2,
                    slope=SLOPE, manning_n=MANNING_N),
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


def test_nearly_dry_steps_are_excluded_from_the_scored_set():
    """A stage at the floor is a dry step wearing a number."""
    q = np.linspace(1.0, 50.0, 400)
    stage = manning_depth(q)
    stage[:200] = 0.0  # the first half of the record is dry
    result = FROUDE(build(q=q, stage=stage), None, {})
    assert result.status == PASS
    assert result.diagnostics["scored_steps"] == 200


def test_an_absolute_stage_is_refused_when_the_probe_sets_a_ceiling():
    """`stage` is a depth above the bed by contract, so a level on another zero
    is not convertible — and only a ceiling can notice it.

    An offset only adds to the depth, and the depth enters at the three-halves
    power, so a model reporting an absolute level looks *more* subcritical the
    further its zero sits from the bed: the offset is the cheapest way to pass.
    With the probe's ceiling the reading is refused and the record fails naming
    it; without one the criterion has no way to tell it from a depth.
    """
    q = np.linspace(1.0, 50.0, 400)
    level = manning_depth(q) + 100.0
    refused = FROUDE(build(q=q, stage=level), None, {"max_depth_m": 10.0})
    assert refused.status == FAIL
    assert refused.diagnostics["scored_steps"] == 0
    assert refused.diagnostics["too_deep_steps"] == 400
    assert "ceiling" in refused.message
    # The same reading with no ceiling is a hundred-metre-deep reach that passes.
    unscored = FROUDE(build(q=q, stage=level), None, {})
    assert unscored.status == PASS
    assert unscored.diagnostics["max_froude"] < 0.01


def test_the_ceiling_leaves_a_depth_the_section_can_hold_alone():
    """It refuses what is not a depth, not what is merely deep."""
    q = np.linspace(1.0, 50.0, 400)
    result = FROUDE(build(q=q), None, {"max_depth_m": 10.0})
    assert result.status == PASS
    assert result.diagnostics["too_deep_steps"] == 0
    assert result.diagnostics["max_depth_m"] == 10.0
    assert result.diagnostics["scored_steps"] == 400


def test_there_is_no_ceiling_unless_the_probe_sets_one():
    """How deep is too deep is a statement about the case's reach, not about
    `stage`, so the shared criterion does not guess at one."""
    q = np.linspace(1.0, 50.0, 400)
    result = FROUDE(build(q=q, stage=manning_depth(q) + 100.0), None, {})
    assert result.diagnostics["max_depth_m"] is None
    assert result.diagnostics["too_deep_steps"] == 0


def test_a_stage_below_the_section_is_still_read_as_it_stands():
    """With the bed gone there is nothing to subtract: the depth is the reading.

    The criterion used to support a case-declared datum, which is exactly the
    licence a model reporting an absolute level needs. The contract now fixes
    the datum at the bed, so a reading is a depth — shallow or deep — and no
    number the case publishes can move it.
    """
    q = np.linspace(1.0, 50.0, 400)
    run = build(q=q, stage=manning_depth(q) + 100.0)
    run.case.static["bed_elevation_m"] = 100.0
    result = FROUDE(run, None, {})
    assert result.status == PASS
    assert result.diagnostics["max_froude"] < 0.01


def test_the_probe_declares_the_section_and_a_model_that_never_read_it_is_not_judged():
    """`requires.static` is what keeps a model carrying its own river geometry
    out of the archive, and nothing else in this file would notice if the
    declaration were dropped."""
    probe = registry.find_probe("momentum/froude-regime")
    assert probe.requires_static == ("width_m", "slope", "manning_n")
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

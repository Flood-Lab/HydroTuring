"""What a report says about a probe that passed.

A failing probe's `detail` names the criteria that failed. A passing one has
to name what the probe exists to score. On mass/human-abstraction that is the
paired `human_abstraction`: `closure` there checks only that the natural run
closes before the two runs are compared, and its residual, standing alone in
the archive, reads as the probe's result when it is not.
"""

from __future__ import annotations

import re

import pytest

from hydroturing import registry
from hydroturing.harness import run_model
from hydroturing.report import to_csv_rows, to_markdown
from hydroturing.scoring import (
    FAIL,
    OK,
    PASS,
    VIOLATION,
    CriterionOutcome,
    ModelReport,
    ProbeOutcome,
)
from hydroturing.seeds import gate_seeds

# What each merged probe reports a pass by. Pinned in full because a change to
# a probe's criteria or to its must_fail baselines moves this, and that change
# should show up in review as a change to what the archive will say.
HEADLINES = {
    "energy/evaporative-partition": ("partition_shift",),
    "energy/latent-heat-et-consistency": ("energy_closure", "flux_identity"),
    "energy/pet-consistency": ("demand_consistency",),
    "energy/surface-energy-closure": ("energy_closure_by_phase",),
    "mass/antecedent-monotonicity": ("antecedent_monotonicity",),
    "mass/area-invariance": ("invariance",),
    "mass/catchment-closure": ("closure", "state_bounds", "non_degenerate"),
    "mass/causality": ("causality",),
    "mass/dry-down": ("dry_down",),
    "mass/extreme-rain": ("monotone_response",),
    "mass/human-abstraction": ("human_abstraction",),
    "mass/phase-counterfactual": ("phase_invariance",),
    "mass/precipitation-counterfactual": ("counterfactual_response", "monotone_response"),
    "mass/resolution-invariance": ("resolution_invariance",),
    "mass/response-nonnegativity": ("response_nonnegativity",),
    "mass/runoff-bounds": ("runoff_bounds",),
    "mass/steady-state": ("steady_state",),
    "mass/time-origin-invariance": ("invariance",),
    "mass/warming-response": ("response_sign",),
    "momentum/routing-conservation": ("routing_conservation",),
}


def _run(model_name: str, probe_id: str):
    probe = registry.find_probe(probe_id)
    return run_model(registry.find_model(model_name), [probe], gate_seeds(probe.id, 1))


def _criterion(name: str, status: str = "pass", message: str = "", value: float | None = None):
    return CriterionOutcome(name=name, status=status, message=message, value=value)


def _report(*criteria: CriterionOutcome, headline: tuple[str, ...] = ()) -> ModelReport:
    """One probe's outcome built by hand, for the cases no merged probe produces."""
    failing = [c.name for c in criteria if not c.passed]
    outcome = ProbeOutcome(
        probe_id="mass/example",
        law="mass",
        verdict=FAIL if failing else PASS,
        reason=VIOLATION if failing else OK,
        criteria=list(criteria),
        headline=list(headline),
    )
    return ModelReport(model_name="example", model_version="1", suite_version="0", probes=[outcome])


@pytest.fixture(scope="module")
def abstraction():
    return _run("reference_bucket", "mass/human-abstraction")


@pytest.fixture(scope="module")
def latent_heat():
    return _run("reference_coupled", "energy/latent-heat-et-consistency")


@pytest.mark.parametrize("probe_id, headline", sorted(HEADLINES.items()))
def test_headline_is_what_the_probe_exists_to_score(probe_id, headline):
    assert registry.find_probe(probe_id).headline == headline


def test_every_probe_has_a_headline():
    """A probe with none could say nothing about a pass but that it passed."""
    for spec in registry.all_probes():
        assert spec.headline, f"{spec.id} names no criterion to report a pass by"


def test_a_passing_paired_probe_reports_its_paired_criterion(abstraction):
    outcome = abstraction.probes[0]
    assert outcome.verdict == PASS, outcome.failing
    paired = next(c for c in outcome.criteria if c.name == "human_abstraction")

    (row,) = to_csv_rows(abstraction)
    assert row["detail"] == f"human_abstraction: {paired.message}"
    # The control run's closure passed as well, and is not what the row reports.
    assert "closure" not in row["detail"]
    assert f"**human_abstraction**: {paired.message}" in to_markdown(abstraction)


def test_a_passing_single_run_probe_reports_what_its_gate_catches(latent_heat):
    """Closure is declared first on this probe, as the precondition for the
    identity, so neither its position nor its name can stand in for the point."""
    outcome = latent_heat.probes[0]
    assert outcome.verdict == PASS, outcome.failing
    by_name = {c.name: c for c in outcome.criteria}

    (row,) = to_csv_rows(latent_heat)
    assert row["detail"] == (
        f"energy_closure: {by_name['energy_closure'].message}; "
        f"flux_identity: {by_name['flux_identity'].message}"
    )


def test_a_pipe_in_a_message_does_not_split_the_markdown_cell(latent_heat):
    # energy_closure measures its residual against accumulated |rn|.
    (row,) = [
        line for line in to_markdown(latent_heat).splitlines()
        if line.startswith("| `energy/latent-heat-et-consistency` |")
    ]
    assert "\\|rn\\|" in row
    assert len(re.findall(r"(?<!\\)\|", row)) == 6  # the five columns' borders


def test_a_failing_probe_still_reports_only_what_failed():
    report = _run("reference_abstraction_blind", "mass/human-abstraction")
    outcome = report.probes[0]
    assert outcome.verdict == FAIL
    failing = [c for c in outcome.criteria if not c.passed]

    (row,) = to_csv_rows(report)
    assert row["detail"] == "; ".join(f"{c.name}: {c.message}" for c in failing)
    assert row["detail"].startswith("human_abstraction: of the prescribed")


def test_several_failures_are_all_named_and_the_headline_is_not():
    report = _report(
        _criterion("closure", "fail", "cumulative residual 15% of sum_pr"),
        _criterion("human_abstraction", "pass", "the prescribed 380 mm left the budget"),
        _criterion("state_bounds", "fail", "mrso leaves [0, 300] on 4 steps"),
        headline=("human_abstraction",),
    )
    (row,) = to_csv_rows(report)
    assert row["detail"] == (
        "closure: cumulative residual 15% of sum_pr; state_bounds: mrso leaves [0, 300] on 4 steps"
    )


def test_a_long_passing_detail_is_cut_where_a_failing_one_is():
    report = _report(_criterion("human_abstraction", message="x" * 500), headline=("human_abstraction",))
    (row,) = to_csv_rows(report)
    assert len(row["detail"]) == 400
    assert row["detail"].startswith("human_abstraction: xxx")


def test_an_outcome_without_a_headline_does_not_fall_back_to_closure():
    """An outcome built outside run_probe carries no headline, and a bare
    closure residual must not come back as its result."""
    report = _report(_criterion("closure", message="cumulative residual 0.0210% of sum_pr", value=2.1e-4))
    (row,) = to_csv_rows(report)
    assert row["detail"] == "all criteria pass"

"""The second round of probes: runoff bounds, area invariance, step-by-step
response, antecedent memory, phase, demand consistency, routing.

Each is checked the way the gate checks it, on one gate seed: the physical
models pass and the named broken model trips the named criterion. The full
gate does this over every seed; these keep a regression close to the code."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.harness import build_case
from hydroturing.runner import get_runner
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.seeds import gate_seeds

ROUND_2 = {
    "mass/runoff-bounds": ("runoff_bounds", "reference_overflowing"),
    "mass/area-invariance": ("invariance", "reference_area_leak"),
    "mass/response-nonnegativity": ("response_nonnegativity", "reference_overshooting"),
    "mass/antecedent-monotonicity": ("antecedent_monotonicity", "reference_cheater"),
    "mass/phase-counterfactual": ("phase_invariance", "reference_sublimating"),
    "energy/pet-consistency": ("demand_consistency", "reference_thirsty"),
    "momentum/routing-conservation": ("routing_conservation", "reference_stuck_router"),
}


def _score(probe, model_name, seed, criterion):
    model = registry.find_model(model_name)
    runs = {}
    for variant in probe.variants or (None,):
        case = build_case(probe, seed, variant)
        runs[variant or "_"] = get_runner(model).run(model, probe, case, Path(tempfile.mkdtemp()))
    params = dict(next(c.params for c in probe.criteria if c.name == criterion))
    fn = get(criterion)
    subject = runs if probe.variants else runs["_"]
    return fn(subject, probe, params)


@pytest.mark.parametrize("probe_id", sorted(ROUND_2))
def test_physical_models_pass_and_the_broken_one_fails(probe_id):
    probe = registry.find_probe(probe_id)
    criterion, broken = ROUND_2[probe_id]
    seeds = gate_seeds(probe.id, probe.n_seeds)
    for name in probe.must_pass:
        result = _score(probe, name, seeds[0], criterion)
        assert result.status == PASS, f"{name} on {probe_id}: {result.message}"
    # As the gate does: the broken model must trip on at least one gate seed.
    # A fault that only shows on some weather is still a fault.
    tripped = [seed for seed in seeds if _score(probe, broken, seed, criterion).status == FAIL]
    assert tripped, f"{broken} should trip {criterion} on {probe_id} for some gate seed"


def test_declared_exchange_closes_the_budget():
    """A model that declares its exchange with the outside as `gwex` closes;
    the same model with the column hidden is a residual of exactly that size."""
    probe = registry.find_probe("mass/catchment-closure")
    model = registry.find_model("reference_bucket")
    case = build_case(probe, gate_seeds(probe.id, 1)[0])
    run = get_runner(model).run(model, probe, case, Path(tempfile.mkdtemp()))
    params = dict(next(c.params for c in probe.criteria if c.name == "closure"))

    from hydroturing.protocol import RunResult
    leaky = run.table.copy()
    leaky["mrro"] = leaky["mrro"] + 1.0          # a millimetre a day from nowhere
    hidden = RunResult(case=run.case, table=leaky, meta=run.meta, wall_seconds=0.0)
    assert get("closure")(hidden, probe, params).status == FAIL

    declared = leaky.copy()
    declared["gwex"] = 1.0                        # ... now declared as a source
    told = RunResult(case=run.case, table=declared, meta=run.meta, wall_seconds=0.0)
    assert get("closure")(told, probe, params).status == PASS

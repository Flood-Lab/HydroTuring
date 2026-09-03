"""Orchestration: generate a case, run the model, apply the criteria."""

from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd

from hydroturing import SUITE_VERSION, criteria as criteria_mod
from hydroturing.criteria.base import CriterionResult
from hydroturing.protocol import Case, ProtocolError, RunResult
from hydroturing.runner import RunnerError, get_runner
from hydroturing.scoring import (
    FAIL,
    PASS,
    CriterionOutcome,
    ModelReport,
    ProbeOutcome,
    reason_for,
)
from hydroturing.seeds import eval_seeds, gate_seeds
from hydroturing.spec import ModelManifest, ProbeSpec


def load_generator(probe: ProbeSpec):
    """Import the probe's generator module by path.

    Generators are code, not data. That is the point: because forcing is made
    fresh at run time, a reviewer reviews a short deterministic script instead
    of a binary blob, and the repository stays small.
    """
    path = probe.generator_path
    spec = importlib.util.spec_from_file_location(f"hydroturing_gen_{probe.slug}", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import generator {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "generate"):
        raise AttributeError(f"{path} must define generate(seed) -> (DataFrame, dict)")
    return module


def build_case(probe: ProbeSpec, seed: int) -> Case:
    module = load_generator(probe)
    forcing, static = module.generate(seed)

    if not isinstance(forcing, pd.DataFrame):
        raise TypeError(f"{probe.generator} returned {type(forcing)}, expected a DataFrame")
    if "time" not in forcing.columns:
        raise ValueError(f"{probe.generator} must produce a 'time' column")

    expected = int(round(probe.period_years * 365)) + probe.spinup_days
    if len(forcing) != expected:
        raise ValueError(
            f"{probe.generator} produced {len(forcing)} steps, expected {expected} "
            f"({probe.period_years:g} years plus {probe.spinup_days} spinup days)"
        )

    return Case(
        probe_id=probe.id,
        seed=seed,
        forcing=forcing,
        static=dict(static),
        spinup_days=probe.spinup_days,
    )


def evaluate_criteria(run: RunResult, probe: ProbeSpec) -> list[CriterionResult]:
    results = []
    for criterion in probe.criteria:
        fn = criteria_mod.get(criterion.name)
        results.append(fn(run, probe, dict(criterion.params)))
    return results


def run_probe(
    model: ModelManifest,
    probe: ProbeSpec,
    seeds: list[int] | None = None,
    workdir: Path | None = None,
) -> ProbeOutcome:
    """Run one probe across its seeds. Every seed must pass."""
    missing = model.missing_for(probe)
    if missing:
        return ProbeOutcome(
            probe_id=probe.id,
            law=probe.law,
            verdict=FAIL,
            reason=reason_for([], missing, None),
            missing=missing,
        )

    seeds = seeds if seeds is not None else eval_seeds(probe.n_seeds)
    runner = get_runner(model)
    per_criterion: dict[str, list[tuple[int, CriterionResult]]] = {
        c.name: [] for c in probe.criteria
    }
    flags: list[str] = []

    tmp_root = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="hydroturing-"))
    tmp_root.mkdir(parents=True, exist_ok=True)

    for seed in seeds:
        io_dir = tmp_root / f"{model.name}__{probe.slug}__{seed}"
        try:
            run = runner.run(model, probe, build_case(probe, seed), io_dir)
        except (RunnerError, ProtocolError) as exc:
            return ProbeOutcome(
                probe_id=probe.id, law=probe.law, verdict=FAIL,
                reason=reason_for([], [], str(exc)),
                seeds=seeds, error=str(exc),
            )
        for result in evaluate_criteria(run, probe):
            per_criterion[result.name].append((seed, result))
            if result.diagnostics.get("suspicious_exact"):
                flags.append(f"suspicious_exact:{probe.id}")

    outcomes = []
    for name, pairs in per_criterion.items():
        # The worst seed decides. A model that passes four seeds and fails the
        # fifth has not shown conservation, it has shown luck.
        worst_seed, worst = min(pairs, key=lambda pair: (pair[1].passed, -abs(pair[1].value or 0.0)))
        outcomes.append(
            CriterionOutcome(
                name=name,
                status=worst.status,
                message=worst.message,
                value=worst.value,
                threshold=worst.threshold,
                worst_seed=worst_seed,
                per_seed=[
                    {"seed": s, "status": r.status, "value": r.value} for s, r in pairs
                ],
                diagnostics=worst.diagnostics,
            )
        )

    failing = [o.name for o in outcomes if not o.passed]
    return ProbeOutcome(
        probe_id=probe.id,
        law=probe.law,
        verdict=PASS if not failing else FAIL,
        reason=reason_for(failing, [], None),
        seeds=seeds,
        criteria=outcomes,
        flags=sorted(set(flags)),
    )


def run_model(
    model: ModelManifest,
    probes: list[ProbeSpec],
    seeds: list[int] | None = None,
    gate: bool = False,
    workdir: Path | None = None,
) -> ModelReport:
    report = ModelReport(
        model_name=model.name,
        model_version=model.version,
        suite_version=SUITE_VERSION,
    )
    for probe in probes:
        probe_seeds = seeds
        if probe_seeds is None and gate:
            probe_seeds = gate_seeds(probe.id, probe.n_seeds)
        outcome = run_probe(model, probe, probe_seeds, workdir=workdir)
        report.probes.append(outcome)
        report.flags.extend(outcome.flags)
    report.flags = sorted(set(report.flags))
    return report

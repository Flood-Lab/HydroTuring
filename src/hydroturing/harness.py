"""Orchestration: generate a case, run the model, apply the criteria."""

from __future__ import annotations

import importlib.util
import inspect
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


def build_case(probe: ProbeSpec, seed: int, variant: str | None = None) -> Case:
    module = load_generator(probe)

    if variant is None:
        forcing, static = module.generate(seed)
    else:
        if "variant" not in inspect.signature(module.generate).parameters:
            raise TypeError(
                f"{probe.generator} must define generate(seed, variant) for a "
                f"probe that declares case.variants {list(probe.variants)}"
            )
        forcing, static = module.generate(seed, variant=variant)

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
        probe_id=probe.id if variant is None else f"{probe.id}@{variant}",
        seed=seed,
        forcing=forcing,
        static=dict(static),
        spinup_days=probe.spinup_days,
    )


def evaluate_criteria(runs: dict[str, RunResult], probe: ProbeSpec) -> list[CriterionResult]:
    """Score one seed.

    `runs` is keyed by variant, with a single entry under the control name for
    an ordinary probe. Paired criteria are handed the whole mapping because
    what they assert is a relationship between the runs; every other criterion
    sees the control run alone, so that adding a variant to a probe never
    silently changes what its existing criteria measure.
    """
    control = runs[probe.control] if probe.control else next(iter(runs.values()))
    results = []
    for criterion in probe.criteria:
        fn = criteria_mod.get(criterion.name)
        subject = runs if criteria_mod.is_paired(criterion.name) else control
        results.append(fn(subject, probe, dict(criterion.params)))
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
            authors=list(probe.authors),
        )

    seeds = seeds if seeds is not None else eval_seeds(probe.n_seeds)
    runner = get_runner(model)
    per_criterion: dict[str, list[tuple[int, CriterionResult]]] = {
        c.name: [] for c in probe.criteria
    }
    flags: list[str] = []

    tmp_root = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="hydroturing-"))
    tmp_root.mkdir(parents=True, exist_ok=True)

    # An ordinary probe has one unnamed case per seed. A paired probe runs the
    # model once per variant on the same seed, which is what makes a
    # counterfactual askable at all.
    variants: tuple[str | None, ...] = probe.variants or (None,)

    for seed in seeds:
        runs: dict[str, RunResult] = {}
        try:
            for variant in variants:
                suffix = f"__{variant}" if variant else ""
                io_dir = tmp_root / f"{model.name}__{probe.slug}__{seed}{suffix}"
                case = build_case(probe, seed, variant)
                runs[variant or "_"] = runner.run(model, probe, case, io_dir)
        except (RunnerError, ProtocolError) as exc:
            return ProbeOutcome(
                probe_id=probe.id, law=probe.law, verdict=FAIL,
                reason=reason_for([], [], str(exc)),
                seeds=seeds, error=str(exc), authors=list(probe.authors),
            )
        for result in evaluate_criteria(runs, probe):
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
        authors=list(probe.authors),
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

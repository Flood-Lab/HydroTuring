"""Orchestration: generate a case, run the model, apply the criteria."""

from __future__ import annotations

import importlib.util
import inspect
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd

from hydroturing import SUITE_VERSION, criteria as criteria_mod
from hydroturing.criteria.base import CriterionResult
from hydroturing.protocol import Case, RunResult
from hydroturing.runner import get_runner
from hydroturing.scoring import (
    FAIL,
    NOT_SCORED,
    PASS,
    CriterionOutcome,
    ModelReport,
    ProbeOutcome,
    reason_for,
)
from hydroturing.seeds import eval_seeds, gate_seeds
from hydroturing.spec import (
    DEFAULT_WINDOW_DAYS,
    FULL_WINDOW,
    TIMESTEP_DAYS,
    TRUSTED_SUBPROCESS_MODELS,
    ModelManifest,
    ProbeSpec,
)

# What locates a flood event when no reference model can say: the stretch
# with the most precipitation. See `event_signal` for why that is the
# fallback rather than the rule.
WINDOW_DRIVER = "pr"


class WindowError(ValueError):
    """A window cannot be cut from this case without breaking the probe."""


class IncompatibleError(Exception):
    """The model cannot consume this probe, which is N/A (INCOMPATIBLE) for it.

    Raised before the adapter runs, so that a check that never put the probe
    to the model is not mistaken for an adapter that broke the contract.
    """

    def __init__(self, issues: list[str]):
        super().__init__("; ".join(issues))
        self.issues = list(issues)


def resolve_window_days(
    model: ModelManifest, probe: ProbeSpec, override: int | str | None = None
) -> int | None:
    """How many days of the scored record this model is evaluated on.

    None means the whole record. The command line wins over the manifest,
    the manifest over the default, and the default depends on what kind of
    model it is: a submitted model gets a flood-event window unless it asks
    for otherwise, a reference model always gets the full record because
    the acceptance gate is defined on it.
    """
    choice = override if override is not None else model.window_days
    if choice is None:
        if model.name in TRUSTED_SUBPROCESS_MODELS:
            return None
        # The default follows the step the model is actually run at, which
        # on a multi-step probe is its own step, not the probe's finest.
        control = select_variants(model, probe)[0]
        choice = DEFAULT_WINDOW_DAYS[probe.timestep_for(control)]
    if choice == FULL_WINDOW:
        return None
    days = int(choice)
    if days < 1:
        raise ValueError(f"an evaluation window must be at least one day, not {days}")
    # A probe whose expectation only holds over a long enough stretch widens
    # the window to it. The model's own preference is a floor, not a cap.
    return max(days, probe.min_window_days)


def event_signal(case: Case, probe: ProbeSpec) -> np.ndarray | None:
    """A per-step series whose largest accumulation marks the flood event.

    The probe's own physical baseline decides. Every probe must name a
    reference model it passes, and that model's runoff is the probe's
    definition of what this catchment does with its weather: where the
    exact physics puts the flood is where the flood is. Precipitation alone
    is a poor guide in a catchment with a snowpack, where the wettest month
    is winter accumulation and the flood is the melt three months later.

    Falls back to precipitation when no trusted reference model reports
    runoff, and to None, meaning take the first stretch, when there is no
    precipitation column either.
    """
    from hydroturing import registry  # noqa: PLC0415 - registry imports spec, not this module

    for name in probe.must_pass:
        try:
            reference = registry.find_model(name)
        except KeyError:
            continue
        if reference.runner != "subprocess" or "mrro" not in reference.emitted:
            continue
        io_dir = Path(tempfile.mkdtemp(prefix=f"hydroturing-window-{reference.name}-"))
        result = get_runner(reference).run(reference, probe, case, io_dir)
        return result.table["mrro"].to_numpy(dtype=float)

    if WINDOW_DRIVER in case.forcing.columns:
        return case.forcing[WINDOW_DRIVER].to_numpy(dtype=float)
    return None


@dataclass(frozen=True)
class WindowBounds:
    """Where the evaluation window sits, in time rather than rows.

    Time, so that the same stretch of weather can be cut from variants that
    run at different steps: an offset into the scored record and a length,
    both in days.
    """

    offset_days: float
    days: int

    def rows(self, dt_days: float) -> tuple[int, int]:
        """Start offset and length in rows at a given step."""
        return int(round(self.offset_days / dt_days)), max(1, int(round(self.days / dt_days)))


def select_window(case: Case, probe: ProbeSpec, days: int) -> WindowBounds:
    """The `days`-long flood event after spinup.

    The event is the stretch over which `event_signal` accumulates the most,
    so a window lands on the largest multi-day flood of the record rather
    than on a single wet day. A window at least as long as the scored record
    starts at the record's start and covers all of it.
    """
    _, rows = WindowBounds(0.0, days).rows(case.dt_days)
    scored_rows = case.n_steps - case.spinup_steps
    if rows >= scored_rows:
        return WindowBounds(0.0, days)

    offset = 0
    signal = event_signal(case, probe)
    if signal is not None:
        totals = np.convolve(signal[case.spinup_steps :], np.ones(rows), mode="valid")
        offset = int(np.argmax(totals))

    # Snap the start to the coarsest step the probe runs at, so that every
    # variant is cut at the same instant. A window that began in the middle
    # of an hour would hand the hourly variant a different half hour of rain
    # than the minute variant, and the comparison would be off by a storm.
    grain = max(TIMESTEP_DAYS[step] for step in probe.timesteps)
    offset_days = np.floor(offset * case.dt_days / grain + 1e-9) * grain
    return WindowBounds(float(offset_days), days)


def window_case(case: Case, bounds: WindowBounds) -> Case:
    """Cut a case down to a window, keeping the spinup that leads into it.

    The scored stretch starts `bounds.offset_days` into the scored record
    and runs for `bounds.days`, converted to rows at this case's own step.
    The model receives the spinup rows before it as well, so its own spinup
    is the same length as it would be on the full record.

    A generator that labels the record (any `_`-prefixed column) has told the
    criteria that particular stretches matter. If the window drops a label
    that the full scored record carries, the probe cannot be scored on it and
    that is reported as incompatible rather than quietly evaluated on the
    wrong stretch.
    """
    offset, rows = bounds.rows(case.dt_days)
    start = case.spinup_steps + offset
    stop = min(start + rows, case.n_steps)
    if start >= case.n_steps:
        raise WindowError(
            f"a window {bounds.offset_days:g} days into the scored record falls "
            "outside a record that short"
        )
    forcing = case.forcing.iloc[start - case.spinup_steps : stop].reset_index(drop=True)

    for column in case.forcing.columns:
        if not column.startswith("_"):
            continue
        full = set(case.forcing[column].iloc[case.spinup_steps :].astype(str))
        kept = set(forcing[column].iloc[case.spinup_steps :].astype(str))
        lost = sorted(full - kept)
        if lost:
            raise WindowError(
                f"a {bounds.days}-day window drops the {lost} stretch of '{column}' "
                "that this probe scores; evaluate the full record instead "
                "(window_days: full)"
            )

    scored = forcing.iloc[case.spinup_steps :]
    window = {
        "days": bounds.days,
        "rows": stop - start,
        "start": str(scored["time"].iloc[0]),
        "end": str(scored["time"].iloc[-1]),
    }
    return Case(
        probe_id=case.probe_id,
        seed=case.seed,
        forcing=forcing,
        static=dict(case.static),
        spinup_steps=case.spinup_steps,
        timestep=case.timestep,
        window=window,
    )


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

    timestep = probe.timestep_for(variant)
    expected = probe.n_steps_for(variant)
    if len(forcing) != expected:
        which = "" if variant is None else f" for variant '{variant}'"
        raise ValueError(
            f"{probe.generator} produced {len(forcing)} steps{which}, expected "
            f"{expected} ({probe.period_days:g} days plus {probe.spinup_days} "
            f"spinup days at {timestep})"
        )

    return Case(
        probe_id=probe.id if variant is None else f"{probe.id}@{variant}",
        seed=seed,
        forcing=forcing,
        static=dict(static),
        spinup_steps=probe.spinup_steps_for(variant),
        timestep=timestep,
    )


def evaluate_criteria(
    runs: dict[str, RunResult], probe: ProbeSpec, control: str | None = None
) -> list[CriterionResult]:
    """Score one seed.

    `runs` is keyed by variant, with a single entry under the control name for
    an ordinary probe. Paired criteria are handed the whole mapping because
    what they assert is a relationship between the runs; every other criterion
    sees the control run alone, so that adding a variant to a probe never
    silently changes what its existing criteria measure. `control` names the
    control run when it is not the probe's first variant, as on a probe that
    selects variants per model.
    """
    if control is not None:
        control_run = runs[control]
    elif probe.control in runs:
        control_run = runs[probe.control]
    else:
        control_run = next(iter(runs.values()))
    control = control_run
    results = []
    for criterion in probe.criteria:
        fn = criteria_mod.get(criterion.name)
        subject = runs if criteria_mod.is_paired(criterion.name) else control
        results.append(fn(subject, probe, dict(criterion.params)))
    return results


def compatibility_issues(
    model: ModelManifest,
    probe: ProbeSpec,
    case: Case | None = None,
    *,
    check_perturbation: bool = True,
) -> list[str]:
    """Explain why a model cannot be meaningfully run on a probe.

    On a single-step probe a step the model does not declare is an
    incompatibility: the probe measures something else, closure say, and
    cannot measure it on a model it cannot feed. A probe that runs the same
    weather at several steps is asking how the model copes with the step,
    so the manifest's declaration is not consulted: the model is run at the
    steps the probe selects and what it does there is measured.
    """
    issues: list[str] = []
    if len(probe.timesteps) == 1:
        needed = [case.timestep] if case is not None else [probe.timestep]
        unsupported = [t for t in needed if not model.supports_timestep(t)]
        if unsupported:
            issues.append(
                f"model timestep {'/'.join(model.timesteps)} does not cover the "
                f"probe's {'/'.join(unsupported)}"
            )
    if check_perturbation and probe.variants and not model.supports_perturbation:
        issues.append("model does not declare support for paired perturbation cases")
    if case is not None:
        visible = {c for c in case.forcing.columns if not c.startswith("_")}
        missing = [v for v in model.needs_forcing if v not in visible]
        if missing:
            issues.append("forcing does not provide " + ", ".join(missing))
    return issues


def select_variants(model: ModelManifest, probe: ProbeSpec) -> list[str | None]:
    """The variants a model is run on for one seed, the control first.

    An ordinary probe runs every variant, control first. A probe declaring
    `variant_selection: native_and_finer` serves the same weather at several
    steps and runs each model at two of them: the one at the model's own
    step, or the nearest coarser step when the probe has no variant at the
    native step, and the next finer one. The finest step, having nothing
    finer, pairs with the next coarser. The variant at the model's step is
    the control, so single-run criteria judge the model where it lives and
    the paired criterion measures what the step does to it.

    The manifest's list of steps is not consulted here. How the model copes
    with a step it did not declare is precisely what such a probe measures,
    and a declaration is not a measurement.
    """
    if not probe.variants:
        return [None]
    if probe.variant_selection != "native_and_finer":
        return list(probe.variants)

    by_step = sorted(probe.variants, key=lambda v: TIMESTEP_DAYS[probe.timestep_for(v)])
    native_dt = TIMESTEP_DAYS[model.timestep]
    index = next(
        (i for i, v in enumerate(by_step) if TIMESTEP_DAYS[probe.timestep_for(v)] >= native_dt),
        len(by_step) - 1,
    )
    if index == 0:
        return [by_step[0], by_step[1]]
    return [by_step[index], by_step[index - 1]]


def verify_adapter_contract(
    model: ModelManifest,
    probe: ProbeSpec,
    seed: int,
    workdir: Path | None = None,
    window: int | str | None = None,
) -> RunResult:
    """Invoke an adapter once and validate only the outputs it declares.

    Contract verification must not short-circuit merely because a scientific
    probe needs variables the model does not produce. That limitation belongs
    to the probe's later N/A (INCOMPLETE), not to this smoke test.

    A probe the model cannot consume is another matter: a step it does not
    declare, a forcing it needs and the probe does not generate, a window
    that drops a stretch the probe scores. There is nothing to run the
    adapter on, so IncompatibleError is raised before it is invoked, on the
    same grounds that make `run_probe` call the probe N/A (INCOMPATIBLE).

    The case is cut to the model's evaluation window, as it will be in the
    real run, so a model that only fits its time budget on the window is
    smoke-tested under the same conditions it is scored under.
    """
    # The control variant is the one at the model's own step, so a daily
    # model is smoke-tested on daily rows rather than on a month of minutes.
    case = build_case(probe, seed, select_variants(model, probe)[0])
    issues = compatibility_issues(model, probe, case, check_perturbation=False)
    if issues:
        raise IncompatibleError(issues)

    days = resolve_window_days(model, probe, window)
    if days is not None:
        bounds = select_window(case, probe, days)
        try:
            case = window_case(case, bounds)
        except WindowError as exc:
            raise IncompatibleError([str(exc)]) from exc

    smoke_probe = replace(
        probe,
        requires_fluxes=model.emits_fluxes,
        requires_states=model.emits_states,
        requires_diagnostics=model.emits_diagnostics,
        variants=(),
        criteria=(),
    )
    if workdir is None:
        io_dir = Path(tempfile.mkdtemp(prefix=f"hydroturing-verify-{model.name}-"))
    else:
        io_dir = Path(workdir) / f"{model.name}__adapter_verification"
    return get_runner(model).run(model, smoke_probe, case, io_dir)


def run_probe(
    model: ModelManifest,
    probe: ProbeSpec,
    seeds: list[int] | None = None,
    workdir: Path | None = None,
    window: int | str | None = None,
) -> ProbeOutcome:
    """Run one probe across its seeds. Every seed must pass.

    `window` overrides the model's evaluation window: a number of days, or
    "full" for the whole record. Left as None, the manifest decides.
    """
    missing = model.missing_for(probe)
    incompatible = compatibility_issues(model, probe)
    if missing or incompatible:
        return ProbeOutcome(
            probe_id=probe.id,
            law=probe.law,
            verdict=NOT_SCORED,
            reason=reason_for([], missing, None, incompatible),
            missing=missing,
            incompatible=incompatible,
            authors=list(probe.authors),
        )

    seeds = seeds if seeds is not None else eval_seeds(probe.n_seeds)
    days = resolve_window_days(model, probe, window)
    runner = get_runner(model)
    per_criterion: dict[str, list[tuple[int, CriterionResult]]] = {
        c.name: [] for c in probe.criteria
    }
    flags: list[str] = []
    windows: list[dict] = []

    tmp_root = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="hydroturing-"))
    tmp_root.mkdir(parents=True, exist_ok=True)

    # An ordinary probe has one unnamed case per seed. A paired probe runs the
    # model once per variant on the same seed, which is what makes a
    # counterfactual askable at all. The control comes first.
    variants = select_variants(model, probe)

    def incompatible_outcome(issues: list[str]) -> ProbeOutcome:
        return ProbeOutcome(
            probe_id=probe.id,
            law=probe.law,
            verdict=NOT_SCORED,
            reason=reason_for([], [], None, issues),
            incompatible=issues,
            seeds=seeds,
            authors=list(probe.authors),
            window_days=days,
        )

    def error_outcome(exc: Exception) -> ProbeOutcome:
        # A bare assert or an exhausted next() has no message. It is still the
        # machinery failing, and an empty message must not read as reason OK.
        message = str(exc) or type(exc).__name__
        return ProbeOutcome(
            probe_id=probe.id, law=probe.law, verdict=FAIL,
            reason=reason_for([], [], message),
            seeds=seeds, error=message, authors=list(probe.authors),
            window_days=days, windows=windows,
        )

    # Every seed's cases are built, checked and cut to their windows before
    # the model runs on any of them. A case the model cannot consume makes
    # the probe N/A, and N/A has to mean the probe asked the model nothing:
    # found on a later seed, after earlier seeds had been scored, it would
    # discard failures that were already measured.
    prepared: list[tuple[int, list[tuple[str | None, Case]]]] = []
    try:
        for seed in seeds:
            # One window per seed, chosen on the control variant and applied
            # to every variant in time, so a paired probe still compares the
            # same stretch of the same weather under its two treatments, even
            # when the treatments run at different steps.
            bounds: WindowBounds | None = None
            cases: list[tuple[str | None, Case]] = []
            for variant in variants:
                case = build_case(probe, seed, variant)
                issues = compatibility_issues(model, probe, case)
                if issues:
                    return incompatible_outcome(issues)
                if days is not None:
                    if bounds is None:
                        bounds = select_window(case, probe, days)
                    try:
                        case = window_case(case, bounds)
                    except WindowError as exc:
                        return incompatible_outcome([str(exc)])
                    if variant in (None, variants[0]):
                        windows.append({"seed": seed, **case.window})
                cases.append((variant, case))
            prepared.append((seed, cases))
    except Exception as exc:  # noqa: BLE001 - a generator or window-locating failure is an ERROR, not a traceback
        return error_outcome(exc)

    for seed, cases in prepared:
        runs: dict[str, RunResult] = {}
        try:
            for variant, case in cases:
                suffix = f"__{variant}" if variant else ""
                io_dir = tmp_root / f"{model.name}__{probe.slug}__{seed}{suffix}"
                runs[variant or "_"] = runner.run(model, probe, case, io_dir)
            for result in evaluate_criteria(runs, probe, control=variants[0] or "_"):
                per_criterion[result.name].append((seed, result))
                if result.diagnostics.get("suspicious_exact"):
                    flags.append(f"suspicious_exact:{probe.id}")
        except Exception as exc:  # noqa: BLE001 - a runner, protocol or criterion failure is an ERROR
            return error_outcome(exc)

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
        headline=list(probe.headline),
        flags=sorted(set(flags)),
        authors=list(probe.authors),
        window_days=days,
        windows=windows,
    )


def run_model(
    model: ModelManifest,
    probes: list[ProbeSpec],
    seeds: list[int] | None = None,
    gate: bool = False,
    workdir: Path | None = None,
    window: int | str | None = None,
) -> ModelReport:
    report = ModelReport(
        model_name=model.name,
        model_version=model.version,
        suite_version=SUITE_VERSION,
        runner=model.runner,
    )
    for probe in probes:
        probe_seeds = seeds
        if probe_seeds is None and gate:
            probe_seeds = gate_seeds(probe.id, probe.n_seeds)
        outcome = run_probe(model, probe, probe_seeds, workdir=workdir, window=window)
        report.probes.append(outcome)
        report.flags.extend(outcome.flags)
    report.flags = sorted(set(report.flags))
    return report

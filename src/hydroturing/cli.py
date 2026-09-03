"""Command line interface.

    ht init-probe [--template KIND]      start a new probe from a template
    ht init-model --name NAME            start a new model from the template
    ht list                              what probes and models exist
    ht validate                          schema-check every probe and model
    ht verify-adapter --model NAME       does the adapter honour the contract
    ht run --model NAME [--probe ID]     evaluate a model
    ht gate [--probe ID]                 the probe acceptance gate used by CI
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace
import tempfile
from pathlib import Path

from hydroturing import SUITE_VERSION, __version__
from hydroturing import registry
from hydroturing.harness import (
    resolve_window_days,
    run_model,
    run_probe,
    verify_adapter_contract,
)
from hydroturing.report import (
    append_csv,
    append_csv_rows,
    contract_row,
    mark,
    prefix,
    to_markdown,
    to_text,
    write_json,
)
from hydroturing.scaffold import (
    available_templates,
    scaffold_model,
    scaffold_probe,
    write_draft,
)
from hydroturing.scoring import ERROR, FAIL, PASS
from hydroturing.seeds import gate_seeds
from hydroturing.spec import (
    FULL_WINDOW,
    TRUSTED_SUBPROCESS_MODELS,
    SpecError,
    load_model,
    load_probe,
)


def _window_arg(text: str) -> int | str:
    """`--window 30` or `--window full`."""
    if text.strip().lower() == FULL_WINDOW:
        return FULL_WINDOW
    try:
        days = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"expected a number of days or '{FULL_WINDOW}', got {text!r}"
        ) from None
    if days < 1:
        raise argparse.ArgumentTypeError("the window must be at least one day")
    return days


# One line each, shown by `ht init-probe --list-templates`. A template with no
# entry still works; it just goes undescribed.
TEMPLATE_BLURB = {
    "default": "a conservation budget over one generated case",
    "extrapolation-space": "catchments outside the hull models are fitted to",
    "extrapolation-time": "conditions outside anything earlier in the record",
    "counterfactual": "same seed twice, perturbed: where did the extra water go",
    "invariance": "a transform the physics does not depend on must change nothing",
}


def _probe_roots(args) -> list[Path] | None:
    extra = getattr(args, "probe_root", None) or []
    if not extra:
        return None
    paths = [Path(p).resolve() for p in extra]
    missing = [str(path) for path in paths if not path.is_dir()]
    if missing:
        raise SpecError("probe roots do not exist: " + ", ".join(missing))
    return [registry.PROBES_DIR, *paths]


def _probe_paths(args) -> list[Path]:
    roots = _probe_roots(args)
    return registry.probe_paths() if roots is None else [
        path for root in roots for path in registry.probe_paths(root)
    ]


def _all_probes(args):
    return registry.all_probes(_probe_roots(args))


def _probes(args):
    roots = _probe_roots(args)
    return [registry.find_probe(args.probe, roots)] if args.probe else registry.all_probes(roots)


def cmd_init_probe(args) -> int:
    """Two steps, so the contributor never faces a blank page.

    Without --from, hand them a template to fill in. With it, turn the filled
    template into a directory that already validates, leaving only the part
    they alone can write.
    """
    if args.list_templates:
        print("templates (ht init-probe --template <name>):\n")
        for name in sorted(available_templates()):
            print(f"    {name:<22} {TEMPLATE_BLURB.get(name, '')}")
        return 0

    if not args.from_file:
        dest = write_draft(args.output, kind=args.template)
        print(f"wrote {dest}  (template: {args.template})\n")
        print("Fill it in, then run:")
        print(f"    ht init-probe --from {dest}\n")
        print("Read docs/writing-a-probe.md first. The question to answer before")
        print("you write anything: how would a model pass your probe while")
        print("understanding no physics?")
        return 0

    target = scaffold_probe(args.from_file, force=args.force)
    rel = target.relative_to(Path.cwd()) if target.is_relative_to(Path.cwd()) else target
    probe = load_probe(target)
    print(f"created {rel}/\n")
    for name in sorted(p.name for p in target.iterdir()):
        print(f"    {name}")
    print(f"\nNext:")
    print(f"  1. write the forcing in {rel}/{probe.generator}")
    print(f"  2. complete {rel}/README.md")
    print(f"  3. ht gate --probe {probe.id}")
    print("\nThe gate is what proves your probe discriminates. It is not a")
    print("formality: a probe that cannot separate the reference models is")
    print("measuring nothing.")
    return 0


def cmd_init_model(args) -> int:
    target = scaffold_model(args.name, force=args.force)
    rel = target.relative_to(Path.cwd()) if target.is_relative_to(Path.cwd()) else target
    print(f"created {rel}/\n")
    print("Next:")
    print(f"  1. declare what the model honestly emits in {rel}/model.yaml")
    print(f"  2. implement simulate() in {rel}/ht_adapter.py")
    print(f"  3. install the model in {rel}/Dockerfile")
    print(f"  4. ht verify-adapter --model {args.name}")
    print("\nSee AGENTS.md for the contract, including the three rules that are")
    print("easy to get wrong.")
    return 0


def cmd_list(args) -> int:
    probes = _all_probes(args)
    models = registry.all_models()
    print(f"probes ({len(probes)}):")
    for p in probes:
        req = ", ".join(p.required_vars)
        steps = "/".join(p.timesteps)
        print(f"  {p.id:<34} {p.law:<9} {p.track:<10} {steps:<10} seeds={p.n_seeds}  needs: {req}")
    print(f"\nmodels ({len(models)}):")
    for m in models:
        matching = next((p for p in probes if p.timestep == m.timestep), None)
        days = resolve_window_days(m, matching) if matching else m.window_days
        window = "full" if days is None else f"{days}d"
        steps = "/".join(m.timesteps)
        print(
            f"  {m.name:<24} v{m.version:<8} runner={m.runner:<11} {steps:<22} "
            f"window={window:<5} emits: {', '.join(m.emitted)}"
        )
    return 0


def cmd_validate(args) -> int:
    problems = []
    probe_paths = _probe_paths(args)
    for path in probe_paths:
        try:
            load_probe(path)
        except Exception as exc:  # noqa: BLE001 - collect every problem, do not stop at the first
            problems.append(f"probe {path.name}: {exc}")
    for path in registry.model_paths():
        try:
            load_model(path)
        except Exception as exc:  # noqa: BLE001
            problems.append(f"model {path.name}: {exc}")

    if not problems:
        try:
            _all_probes(args)  # also checks duplicate ids across public/private roots
        except Exception as exc:  # noqa: BLE001
            problems.append(str(exc))

    if problems:
        print(f"{prefix(False)}validation failed:")
        for problem in problems:
            print(f"  - {problem}")
        return 1
    n_p, n_m = len(probe_paths), len(registry.model_paths())
    print(f"{prefix(True)}validation passed: {n_p} probe(s), {n_m} model(s), "
          f"suite {SUITE_VERSION}")
    return 0


def cmd_verify_adapter(args) -> int:
    """Smoke test before any physics: does the adapter obey the contract?"""
    model = registry.find_model(args.model)
    probes = _all_probes(args)
    if not probes:
        print("no probes available to smoke test against")
        return 1
    roots = _probe_roots(args)
    if args.probe:
        probe = registry.find_probe(args.probe, roots)
    else:
        probe = next(
            (p for p in probes if all(model.supports_timestep(t) for t in p.timesteps)),
            probes[0],
        )

    seed = gate_seeds(probe.id, 1)[0]
    try:
        result = verify_adapter_contract(model, probe, seed, window=args.window)
    except Exception as exc:  # report a clean smoke-test failure, never a traceback
        print(f"{prefix(False)}adapter contract FAILED for {model.name}:\n  {exc}")
        if args.csv:
            _archive_contract(args.csv, model, probe, seed, error=str(exc))
        return 1
    if args.csv:
        path = _archive_contract(args.csv, model, probe, seed, result=result)
        print(f"appended 1 row to {path}")
    window = result.case.window
    stretch = (
        f" ({window['days']}-day flood event, {window['start']} to {window['end']}, "
        f"{result.case.n_steps} rows with spinup)"
        if window
        else f" (full record, {result.case.n_steps} rows)"
    )
    print(
        f"{prefix(True)}adapter contract OK for {model.name} on {probe.id}{stretch}; "
        f"{result.wall_seconds:.1f}s"
    )
    return 0


def _archive_contract(path, model, probe, seed, *, result=None, error=None) -> Path:
    row = contract_row(
        model.name, model.version, SUITE_VERSION, model.runner, probe.id, seed,
        result=result, error=error,
    )
    return append_csv_rows([row], path)


def cmd_run(args) -> int:
    model = registry.find_model(args.model)
    if args.runner:
        if args.runner == "subprocess" and model.name not in TRUSTED_SUBPROCESS_MODELS:
            raise SpecError(
                "runner 'subprocess' is reserved for trusted reference models; "
                "submitted models must run in Docker"
            )
        model = replace(model, runner=args.runner)
    probes = _probes(args)
    seeds = [args.seed] if args.seed is not None else None
    workdir = Path(args.workdir) if args.workdir else None

    report = run_model(
        model, probes, seeds=seeds, gate=args.gate_seeds, workdir=workdir,
        window=args.window,
    )
    print(to_text(report) if not args.markdown else to_markdown(report))
    if args.json:
        path = write_json(report, args.json)
        print(f"\nwrote {path}")
    if args.csv:
        path = append_csv(report, args.csv)
        print(f"appended {len(report.probes)} row(s) to {path}")
    if report.reason == ERROR:
        return 2
    return 0 if report.verdict == PASS else 1


def cmd_gate(args) -> int:
    """The acceptance gate: a probe is only merged if it discriminates.

    Every probe must pass its reference physical model and must fail each
    deliberately broken one, on the specific criterion that is supposed to do
    the catching. A probe that cannot tell them apart is measuring nothing,
    and a tolerance tightened past what an exact model can meet breaks here
    rather than silently mislabelling honest models later.
    """
    failures = []
    for probe in _probes(args):
        seeds = gate_seeds(probe.id, probe.n_seeds)
        print(f"\n{probe.id}  (seeds {seeds})")

        for name in probe.must_pass:
            outcome = run_probe(registry.find_model(name), probe, seeds)
            ok = outcome.verdict == PASS
            note = "" if ok else "BROKEN: "
            print(f"  {mark(ok)}  must_pass  {name:<24} {note}{_why(outcome)}")
            if not ok:
                failures.append(f"{probe.id}: {name} was expected to PASS")

        for name, expected in probe.must_fail.items():
            outcome = run_probe(registry.find_model(name), probe, seeds)
            tripped = outcome.failing
            ok = outcome.verdict == FAIL and expected in tripped
            note = "" if ok else "BROKEN: "
            print(
                f"  {mark(ok)}  must_fail  {name:<24} {note}"
                f"expected `{expected}`, tripped {tripped or 'nothing'}"
            )
            if not ok:
                failures.append(
                    f"{probe.id}: {name} was expected to FAIL on '{expected}', "
                    f"tripped {tripped or 'nothing'}"
                )

    print()
    if failures:
        print(f"{prefix(False)}GATE FAILED ({len(failures)} problem(s)):")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print(f"{prefix(True)}GATE PASSED: every probe separates the reference models "
          "as declared.")
    return 0


def _why(outcome) -> str:
    if outcome.error:
        return outcome.error.splitlines()[0][:90]
    if outcome.missing:
        return "missing " + ", ".join(outcome.missing)
    return "; ".join(f"{c.name}={c.status}" for c in outcome.criteria)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ht", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version",
                        version=f"hydroturing {__version__} (suite {SUITE_VERSION})")
    sub = parser.add_subparsers(dest="command", required=True)

    init_probe = sub.add_parser("init-probe", help="start a new probe from the template")
    init_probe.add_argument("--from", dest="from_file",
                            help="a filled-in draft; creates probes/<law>/<slug>/")
    init_probe.add_argument("-o", "--output", help="where to write the blank draft")
    init_probe.add_argument("--template", default="default",
                            help="which template to start from (see --list-templates)")
    init_probe.add_argument("--list-templates", action="store_true",
                            help="show the available probe templates and exit")
    init_probe.add_argument("--force", action="store_true", help="overwrite an existing probe")
    init_probe.set_defaults(fn=cmd_init_probe)

    init_model = sub.add_parser("init-model", help="start a new model from the template")
    init_model.add_argument("--name", required=True)
    init_model.add_argument("--force", action="store_true")
    init_model.set_defaults(fn=cmd_init_model)

    def add_probe_roots(command) -> None:
        command.add_argument(
            "--probe-root",
            action="append",
            metavar="PATH",
            help=(
                "add a probe tree outside the repository (repeatable; useful for "
                "private evaluation suites)"
            ),
        )

    list_cmd = sub.add_parser("list", help="show probes and models")
    add_probe_roots(list_cmd)
    list_cmd.set_defaults(fn=cmd_list)

    validate = sub.add_parser("validate", help="schema-check probes and models")
    add_probe_roots(validate)
    validate.set_defaults(fn=cmd_validate)

    def add_window(command) -> None:
        command.add_argument(
            "--window",
            type=_window_arg,
            metavar="DAYS|full",
            help=(
                "days of the scored record to evaluate, taken at the largest flood "
                "event, or 'full' for the whole record; overrides window_days in model.yaml"
            ),
        )

    verify = sub.add_parser("verify-adapter", help="check a model honours the /io contract")
    verify.add_argument("--model", required=True)
    verify.add_argument("--probe")
    verify.add_argument("--csv", metavar="PATH",
                        help="append the contract check as a row to this CSV archive")
    add_window(verify)
    add_probe_roots(verify)
    verify.set_defaults(fn=cmd_verify_adapter)

    run = sub.add_parser("run", help="evaluate a model against the suite")
    run.add_argument("--model", required=True)
    run.add_argument("--probe", help="restrict to one probe")
    run.add_argument("--seed", type=int, help="run a single specific seed (for reproducing a failure)")
    run.add_argument("--gate-seeds", action="store_true", help="use the deterministic gate seeds")
    run.add_argument("--json", help="write the machine-readable report here")
    run.add_argument("--csv", metavar="PATH",
                     help="append one row per probe to this CSV archive (e.g. models/result.csv)")
    run.add_argument("--markdown", action="store_true", help="print the PR-comment table")
    run.add_argument("--workdir", help="keep the /io directories here instead of a temp dir")
    run.add_argument("--runner", choices=["subprocess", "docker"],
                     help="override the runner declared in model.yaml (used by CI to exercise the container path)")
    add_window(run)
    add_probe_roots(run)
    run.set_defaults(fn=cmd_run)

    gate = sub.add_parser("gate", help="run the probe acceptance gate")
    gate.add_argument("--probe", help="restrict to one probe")
    add_probe_roots(gate)
    gate.set_defaults(fn=cmd_gate)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.fn(args)
    except (KeyError, SpecError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

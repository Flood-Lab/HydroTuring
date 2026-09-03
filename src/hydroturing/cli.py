"""Command line interface.

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
from hydroturing.harness import run_model, run_probe
from hydroturing.report import to_markdown, to_text, write_json
from hydroturing.scoring import FAIL, PASS
from hydroturing.seeds import gate_seeds
from hydroturing.spec import SpecError, load_model, load_probe


def _probes(args):
    return [registry.find_probe(args.probe)] if args.probe else registry.all_probes()


def cmd_list(args) -> int:
    probes = registry.all_probes()
    models = registry.all_models()
    print(f"probes ({len(probes)}):")
    for p in probes:
        req = ", ".join(p.required_vars)
        print(f"  {p.id:<34} {p.law:<9} {p.track:<10} seeds={p.n_seeds}  needs: {req}")
    print(f"\nmodels ({len(models)}):")
    for m in models:
        print(f"  {m.name:<24} v{m.version:<8} runner={m.runner:<11} emits: {', '.join(m.emitted)}")
    return 0


def cmd_validate(args) -> int:
    problems = []
    for path in registry.probe_paths():
        try:
            load_probe(path)
        except Exception as exc:  # noqa: BLE001 - collect every problem, do not stop at the first
            problems.append(f"probe {path.name}: {exc}")
    for path in registry.model_paths():
        try:
            load_model(path)
        except Exception as exc:  # noqa: BLE001
            problems.append(f"model {path.name}: {exc}")

    if problems:
        print("validation failed:")
        for problem in problems:
            print(f"  - {problem}")
        return 1
    n_p, n_m = len(registry.probe_paths()), len(registry.model_paths())
    print(f"validation passed: {n_p} probe(s), {n_m} model(s), suite {SUITE_VERSION}")
    return 0


def cmd_verify_adapter(args) -> int:
    """Smoke test before any physics: does the adapter obey the contract?"""
    model = registry.find_model(args.model)
    probes = registry.all_probes()
    if not probes:
        print("no probes available to smoke test against")
        return 1
    probe = registry.find_probe(args.probe) if args.probe else probes[0]

    missing = model.missing_for(probe)
    if missing:
        print(f"{model.name} does not emit {', '.join(missing)}; using it anyway to test the contract")

    outcome = run_probe(model, probe, seeds=[gate_seeds(probe.id, 1)[0]])
    if outcome.error:
        print(f"adapter contract FAILED for {model.name}:\n  {outcome.error}")
        return 1
    print(f"adapter contract OK for {model.name} on {probe.id}")
    return 0


def cmd_run(args) -> int:
    model = registry.find_model(args.model)
    if args.runner:
        model = replace(model, runner=args.runner)
    probes = _probes(args)
    seeds = [args.seed] if args.seed is not None else None
    workdir = Path(args.workdir) if args.workdir else None

    report = run_model(model, probes, seeds=seeds, gate=args.gate_seeds, workdir=workdir)
    print(to_text(report) if not args.markdown else to_markdown(report))
    if args.json:
        path = write_json(report, args.json)
        print(f"\nwrote {path}")
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
            print(f"  must_pass  {name:<24} {'ok' if ok else 'BROKEN'}  {_why(outcome)}")
            if not ok:
                failures.append(f"{probe.id}: {name} was expected to PASS")

        for name, expected in probe.must_fail.items():
            outcome = run_probe(registry.find_model(name), probe, seeds)
            tripped = outcome.failing
            ok = outcome.verdict == FAIL and expected in tripped
            print(
                f"  must_fail  {name:<24} {'ok' if ok else 'BROKEN'}  "
                f"expected `{expected}`, tripped {tripped or 'nothing'}"
            )
            if not ok:
                failures.append(
                    f"{probe.id}: {name} was expected to FAIL on '{expected}', "
                    f"tripped {tripped or 'nothing'}"
                )

    print()
    if failures:
        print(f"GATE FAILED ({len(failures)} problem(s)):")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print("GATE PASSED: every probe separates the reference models as declared.")
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

    sub.add_parser("list", help="show probes and models").set_defaults(fn=cmd_list)
    sub.add_parser("validate", help="schema-check probes and models").set_defaults(fn=cmd_validate)

    verify = sub.add_parser("verify-adapter", help="check a model honours the /io contract")
    verify.add_argument("--model", required=True)
    verify.add_argument("--probe")
    verify.set_defaults(fn=cmd_verify_adapter)

    run = sub.add_parser("run", help="evaluate a model against the suite")
    run.add_argument("--model", required=True)
    run.add_argument("--probe", help="restrict to one probe")
    run.add_argument("--seed", type=int, help="run a single specific seed (for reproducing a failure)")
    run.add_argument("--gate-seeds", action="store_true", help="use the deterministic gate seeds")
    run.add_argument("--json", help="write the machine-readable report here")
    run.add_argument("--markdown", action="store_true", help="print the PR-comment table")
    run.add_argument("--workdir", help="keep the /io directories here instead of a temp dir")
    run.add_argument("--runner", choices=["subprocess", "docker"],
                     help="override the runner declared in model.yaml (used by CI to exercise the container path)")
    run.set_defaults(fn=cmd_run)

    gate = sub.add_parser("gate", help="run the probe acceptance gate")
    gate.add_argument("--probe", help="restrict to one probe")
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

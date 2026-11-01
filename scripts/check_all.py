#!/usr/bin/env python3
"""Run the local test suite and build every packaged model, in one command.

This is the entrypoint of the root Dockerfile, and runs the same way on a host
that has the harness installed. It does what CI's gate job does (validate,
pytest, the acceptance gate), then builds each model's own Dockerfile and runs
its adapter once through `ht verify-adapter`, so a binary that compiles but
cannot run is caught too.

    python3 scripts/check_all.py                    # everything
    python3 scripts/check_all.py --skip-tests       # models only
    python3 scripts/check_all.py --skip-models      # tests only
    python3 scripts/check_all.py --models vic summa --jobs 2

Every step runs even when an earlier one fails, so one run reports everything
that is wrong. Exit 0 if all of it passed, 1 if not. `verify-adapter` exiting 1
means no probe can consume the model, so its adapter was never run. That is
N/A rather than a failure, as it is in the model workflow.

Nothing is written to models/result.csv.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from hydroturing import registry
from hydroturing.runner.base import RunnerError
from hydroturing.runner.docker_runner import BUILD_TIMEOUT_S, build, image_tag, require_docker
from hydroturing.spec import REPO_ROOT, ModelManifest

TEST_STEPS = {
    "validate": ["ht", "validate"],
    "pytest": ["pytest", "-q"],
    "gate": ["ht", "gate"],
}

# Built and run first because it takes seconds. If the daemon cannot mount the
# temp directory, every adapter fails the same way, and this finds that out
# before hours of compiling rather than after.
CANARY = "reference_bucket"


def run_tests(results: dict[str, dict]) -> None:
    # TMPDIR is for the model containers' /io mounts, which live on the host.
    # Nothing in these steps starts a container, so pytest's temp trees stay
    # inside this container instead of accumulating on the host.
    env = {k: v for k, v in os.environ.items() if k != "TMPDIR"}
    for name, argv in TEST_STEPS.items():
        print(f"\n=== {name}: {' '.join(argv)}", flush=True)
        start = time.monotonic()
        code = subprocess.run(argv, cwd=REPO_ROOT, env=env).returncode
        results[name] = {
            "status": "PASS" if code == 0 else f"FAIL (exit {code})",
            "seconds": round(time.monotonic() - start),
        }


def build_one(model: ModelManifest, timeout: float) -> dict:
    print(f"build  start {model.name}", flush=True)
    start = time.monotonic()
    entry = {"image": image_tag(model)}
    try:
        build(model, timeout=timeout)
        entry["build"] = "PASS"
    except RunnerError as exc:
        entry["build"] = "FAIL"
        entry["detail"] = str(exc)
    entry["build_seconds"] = round(time.monotonic() - start)
    print(f"build  {entry['build']:<5}{model.name} ({entry['build_seconds']} s)", flush=True)
    return entry


def verify_one(model: ModelManifest, runner: str | None) -> dict:
    argv = ["ht", "verify-adapter", "--model", model.name]
    if runner:
        argv += ["--runner", runner]
    start = time.monotonic()
    # verify-adapter keeps its /io directory for inspection. Here nobody will
    # look, and it would pile up on the host, so each run gets one to discard.
    with tempfile.TemporaryDirectory() as tmp:
        proc = subprocess.run(
            argv, cwd=REPO_ROOT, capture_output=True, text=True,
            env={**os.environ, "TMPDIR": tmp},
        )
    verdict = {0: "PASS", 1: "N/A"}.get(proc.returncode, f"FAIL (exit {proc.returncode})")
    entry = {"verify": verdict, "verify_seconds": round(time.monotonic() - start)}
    if verdict.startswith("FAIL"):
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-20:]
        entry["detail"] = "\n".join(tail)
    print(f"verify {verdict:<5} {model.name} ({entry['verify_seconds']} s)", flush=True)
    return entry


def check_models(names: list[str] | None, jobs: int, timeout: float) -> dict[str, dict]:
    models = registry.all_models()
    if names:
        unknown = set(names) - {m.name for m in models}
        if unknown:
            raise SystemExit(f"unknown models: {', '.join(sorted(unknown))}")
        models = [m for m in models if m.name in names]

    containers = [m for m in models if (m.path / "Dockerfile").exists()]
    # A model with no Dockerfile runs as a subprocess and has nothing to build.
    # The broken reference controls are the gate's to judge, and some of them
    # break the contract on purpose.
    plain = [
        m for m in models
        if m not in containers and not m.name.startswith("reference_")
    ]
    results: dict[str, dict] = {}

    if containers:
        require_docker()
        canary = registry.find_model(CANARY)
        problem = build_one(canary, timeout)
        if problem["build"] == "PASS":
            problem = verify_one(canary, "docker")
        if problem.get("detail"):
            raise SystemExit(
                f"{CANARY} failed, so no model would pass:\n{problem['detail']}\n\n"
                "Inside a container, TMPDIR must be a directory mounted at the "
                "same path on the host: -v \"$DIR:$DIR\" -e TMPDIR=\"$DIR\"."
            )

    print(f"\n=== building {len(containers)} model images, {jobs} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        built = pool.map(lambda m: build_one(m, timeout), containers)
        for model, entry in zip(containers, built):
            results[model.name] = entry

    # One at a time: a case's time budget assumes the model has the machine,
    # and a build competing for the CPUs must not be what runs it over.
    print(f"\n=== running {len(containers) + len(plain)} adapters", flush=True)
    for model in containers:
        if results[model.name]["build"] == "PASS":
            results[model.name].update(verify_one(model, "docker"))
    for model in plain:
        results[model.name] = {"build": "none", **verify_one(model, None)}
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--skip-tests", action="store_true")
    parser.add_argument("--skip-models", action="store_true")
    parser.add_argument("--models", nargs="+", metavar="NAME", help="only these models")
    parser.add_argument("--jobs", type=int, default=4, help="image builds at once (default 4)")
    parser.add_argument("--build-timeout", type=float, default=BUILD_TIMEOUT_S,
                        help=f"seconds per image (default {BUILD_TIMEOUT_S})")
    parser.add_argument("--json", type=Path, help="also write the summary here")
    args = parser.parse_args()

    tests: dict[str, dict] = {}
    models: dict[str, dict] = {}
    if not args.skip_tests:
        run_tests(tests)
    if not args.skip_models:
        models = check_models(args.models, args.jobs, args.build_timeout)

    print("\n=== summary")
    for name, entry in tests.items():
        print(f"{name:<24} {entry['status']}  ({entry['seconds']} s)")
    for name, entry in models.items():
        print(f"{name:<24} build {entry['build']:<5} verify {entry.get('verify', 'skipped')}")
    failures = [n for n, e in tests.items() if e["status"] != "PASS"]
    failures += [
        n for n, e in models.items()
        if e["build"] == "FAIL" or e.get("verify", "").startswith("FAIL")
    ]
    for name in failures:
        if detail := models.get(name, {}).get("detail"):
            print(f"\n--- {name}\n{detail}")
    if args.json:
        args.json.write_text(json.dumps({"tests": tests, "models": models}, indent=2) + "\n")
    print(f"\n{'all passed' if not failures else 'FAILED: ' + ', '.join(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

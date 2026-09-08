"""Local runner for the in-repo reference models."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from hydroturing.runner.base import Runner, RunnerError
from hydroturing.spec import ModelManifest, ProbeSpec


class SubprocessRunner(Runner):
    name = "subprocess"

    @staticmethod
    def resolve_entrypoint(entrypoint: tuple[str, ...] | list[str]) -> list[str]:
        """An entrypoint that names a Python interpreter runs on the one the
        harness itself runs on.

        Manifests say `python3` because that is what a container has. The
        host is not a container: Windows installs Python as `python` or
        `py` and has no `python3`, and a virtual environment's interpreter
        is the one with the harness's dependencies. A trusted in-repo model
        is run by the interpreter already running, whatever it is called.
        """
        argv = list(entrypoint)
        if argv and argv[0] in ("python", "python3", "py"):
            argv[0] = sys.executable
        return argv

    def invoke(self, model: ModelManifest, probe: ProbeSpec, io_dir: Path, request_path: Path) -> None:
        argv = [*self.resolve_entrypoint(model.entrypoint), "--request", str(request_path)]
        env = dict(os.environ)
        env["PYTHONPATH"] = str(model.path)
        env.pop("PYTHONDONTWRITEBYTECODE", None)

        try:
            proc = subprocess.run(
                argv,
                cwd=model.path,
                env=env,
                capture_output=True,
                text=True,
                timeout=probe.max_runtime_s,
            )
        except FileNotFoundError:
            raise RunnerError(
                f"{model.name}: entrypoint {argv[0]!r} not found in {model.path}"
            ) from None
        except subprocess.TimeoutExpired:
            raise RunnerError(
                f"{model.name}: exceeded the {probe.max_runtime_s:g}s budget for {probe.id}"
            ) from None

        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-12:]
            raise RunnerError(
                f"{model.name}: adapter exited {proc.returncode}\n"
                + "\n".join("    " + line for line in tail)
            )

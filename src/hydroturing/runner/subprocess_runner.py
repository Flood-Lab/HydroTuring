"""Local runner for the in-repo reference models."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from hydroturing.runner.base import Runner, RunnerError
from hydroturing.spec import ModelManifest, ProbeSpec


class SubprocessRunner(Runner):
    name = "subprocess"

    def invoke(self, model: ModelManifest, probe: ProbeSpec, io_dir: Path, request_path: Path) -> None:
        argv = [*model.entrypoint, "--request", str(request_path)]
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

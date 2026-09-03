"""Container runner for submitted models.

Isolation is part of the benchmark, not an operational detail. The model runs
with no network, a memory and CPU cap, and a mount that exposes only the case.
It cannot fetch the probe definition, read the tolerance, or phone home.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from hydroturing.runner.base import Runner, RunnerError
from hydroturing.spec import ModelManifest, ProbeSpec

IMAGE_PREFIX = "hydroturing"


def image_tag(model: ModelManifest) -> str:
    return f"{IMAGE_PREFIX}/{model.name}:{model.version}"


def require_docker() -> str:
    """Fail with the actual reason, not a generic build error.

    A missing binary and an unreachable daemon are different problems with
    different fixes, and conflating them wastes a contributor's afternoon.
    """
    docker = shutil.which("docker")
    if docker is None:
        raise RunnerError(
            "docker is not installed. Reference models declare "
            "runner: subprocess and do not need it; submitted models do."
        )
    probe = subprocess.run(
        [docker, "info", "--format", "{{.ServerVersion}}"],
        capture_output=True, text=True,
    )
    if probe.returncode != 0:
        raise RunnerError(
            "the docker CLI is installed but the daemon is not reachable. "
            "Start Docker, or use --runner subprocess for an in-repo model."
        )
    return docker


def build(model: ModelManifest, quiet: bool = True) -> str:
    docker = require_docker()
    dockerfile = model.path / "Dockerfile"
    if not dockerfile.exists():
        raise RunnerError(f"{model.name}: no Dockerfile at {dockerfile}")

    tag = image_tag(model)
    argv = [docker, "build", "-t", tag, "-f", str(dockerfile), str(model.path)]
    if quiet:
        argv.insert(2, "--quiet")
    proc = subprocess.run(argv, capture_output=True, text=True)
    if proc.returncode != 0:
        tail = (proc.stderr or "").strip().splitlines()[-20:]
        raise RunnerError(
            f"{model.name}: image build failed\n"
            + "\n".join("    " + line for line in tail)
        )
    return tag


class DockerRunner(Runner):
    name = "docker"

    @staticmethod
    def command(docker: str, tag: str, model: ModelManifest, io_dir: Path) -> list[str]:
        """Build the run command.

        Isolation is part of the benchmark, not an operational detail, so this
        is asserted by the test suite rather than trusted. No network means the
        model cannot phone home; the single bind mount means it sees the case
        and nothing else, so it cannot read the probe definition or the
        tolerance it is judged against.
        """
        resources = model.resources or {}
        return [
            docker, "run", "--rm",
            "--network", "none",
            "--mount", f"type=bind,source={io_dir.resolve()},target=/io",
            "--memory", f"{resources.get('memory_gb', 4)}g",
            "--cpus", str(resources.get("cpu", 2)),
            "--workdir", "/io",
            tag,
            *model.entrypoint, "--request", "/io/request.json",
        ]

    def invoke(self, model: ModelManifest, probe: ProbeSpec, io_dir: Path, request_path: Path) -> None:
        docker = require_docker()
        tag = build(model)
        argv = self.command(docker, tag, model, io_dir)
        try:
            proc = subprocess.run(
                argv, capture_output=True, text=True, timeout=probe.max_runtime_s + 60
            )
        except subprocess.TimeoutExpired:
            raise RunnerError(
                f"{model.name}: container exceeded the time budget for {probe.id}"
            ) from None

        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-12:]
            raise RunnerError(
                f"{model.name}: container exited {proc.returncode}\n"
                + "\n".join("    " + line for line in tail)
            )

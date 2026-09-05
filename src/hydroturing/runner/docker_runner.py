"""Container runner for submitted models.

Isolation is part of the benchmark, not an operational detail. The model runs
with no network, a memory and CPU cap, and a mount that exposes only the case.
It cannot fetch the probe definition, read the tolerance, or phone home.
"""

from __future__ import annotations

import os
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


def build(model: ModelManifest, quiet: bool = True, docker: str | None = None) -> str:
    docker = docker or require_docker()
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

    def __init__(self) -> None:
        self._tag: str | None = None

    @staticmethod
    def user_flags() -> list[str]:
        """`--user uid:gid` of the invoking user, where the platform has one."""
        if not hasattr(os, "getuid"):
            return []
        return ["--user", f"{os.getuid()}:{os.getgid()}"]

    @staticmethod
    def command(docker: str, tag: str, model: ModelManifest, io_dir: Path) -> list[str]:
        """Build the run command.

        Isolation is part of the benchmark, not an operational detail, so this
        is asserted by the test suite rather than trusted. Inputs and the
        request are mounted read-only; only the output directory is writable.
        The image's own working directory is preserved so a relative entrypoint
        resolves exactly as it did when the image was built.
        """
        resources = model.resources or {}
        request_path = (io_dir / "request.json").resolve()
        input_dir = (io_dir / "input").resolve()
        output_dir = (io_dir / "output").resolve()
        # A manifest asks for the CPUs the model would like; the host has what
        # it has. Docker refuses a quota above the host's count outright, so
        # the request is capped rather than passed through, and a model that
        # asked for eight runs on a four-core runner instead of not at all.
        cpus = min(int(resources.get("cpu", 2)), os.cpu_count() or 1)
        # The container runs as the user invoking the harness, not as root.
        # With every capability dropped, root inside the container has no
        # DAC override, so it cannot write into an output directory the host
        # user owns; on a Linux host (the CI runner) that is a PermissionError
        # on result.csv. As the host user it can, and it can touch nothing
        # else on the host either. HOME points at the tmpfs so a library that
        # wants a cache directory has one that vanishes with the container.
        return [
            docker, "run", "--rm",
            *DockerRunner.user_flags(),
            "--env", "HOME=/tmp",
            "--network", "none",
            "--read-only",
            "--cap-drop", "ALL",
            "--security-opt", "no-new-privileges",
            "--pids-limit", "256",
            "--tmpfs", "/tmp:rw,nosuid,nodev,size=64m",
            "--mount", f"type=bind,source={request_path},target=/io/request.json,readonly",
            "--mount", f"type=bind,source={input_dir},target=/io/input,readonly",
            "--mount", f"type=bind,source={output_dir},target=/io/output",
            "--memory", f"{resources.get('memory_gb', 4)}g",
            "--cpus", str(cpus),
            tag,
            *model.entrypoint, "--request", "/io/request.json",
        ]

    def invoke(self, model: ModelManifest, probe: ProbeSpec, io_dir: Path, request_path: Path) -> None:
        docker = require_docker()
        if self._tag is None:
            self._tag = build(model, docker=docker)
        tag = self._tag
        argv = self.command(docker, tag, model, io_dir)
        try:
            proc = subprocess.run(
                argv, capture_output=True, text=True, timeout=probe.max_runtime_s
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

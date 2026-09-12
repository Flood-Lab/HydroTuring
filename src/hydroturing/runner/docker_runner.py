"""Container runner for submitted models.

Isolation is part of the benchmark, not an operational detail. The model runs
with no network, a memory and CPU cap, and a mount that exposes only the case.
It cannot fetch the probe definition, read the tolerance, or phone home.
"""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path

from hydroturing.runner.base import Runner, RunnerError
from hydroturing.spec import ModelManifest, ProbeSpec

IMAGE_PREFIX = "hydroturing"

# Twice the fifteen minutes a first build of a heavy image has taken, and
# inside the model workflow's 45-minute job, so CI reports a stuck build
# rather than cancelling it. A cached rebuild takes seconds; this is for a
# build that never exits, not a budget for a slow one.
BUILD_TIMEOUT_S = 1800

# How long `docker kill`, and then `docker rm -f`, may take to end a container
# that ran past its budget. Either normally returns within seconds; the limit
# is for a daemon that has stopped answering, so the budget error still
# reaches the caller instead of a second hang.
KILL_TIMEOUT_S = 30


def image_tag(model: ModelManifest) -> str:
    return f"{IMAGE_PREFIX}/{model.name}:{model.version}"


def container_name(model: ModelManifest) -> str:
    """A name unique to one run, so two sessions running the same model on
    one daemon never collide."""
    return f"{IMAGE_PREFIX}-{model.name}-{uuid.uuid4().hex}"


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


def build(
    model: ModelManifest,
    quiet: bool = True,
    docker: str | None = None,
    timeout: float = BUILD_TIMEOUT_S,
) -> str:
    docker = docker or require_docker()
    dockerfile = model.path / "Dockerfile"
    if not dockerfile.exists():
        raise RunnerError(f"{model.name}: no Dockerfile at {dockerfile}")

    tag = image_tag(model)
    argv = [docker, "build", "-t", tag, "-f", str(dockerfile), str(model.path)]
    if quiet:
        argv.insert(2, "--quiet")
    # The log goes to a file, not a pipe. On Docker Desktop the build starts
    # `docker-credential-desktop get`, and that helper can be orphaned still
    # holding the stderr it inherited. A pipe reaches end of file only when
    # every holder has closed it, so reading one waited on the orphan and a
    # run hung after its image was built. With a file the wait is on docker
    # alone. The build keeps the caller's session: in a new one, Ctrl+C at
    # the terminal would no longer reach docker.
    with tempfile.TemporaryFile() as log:
        try:
            proc = subprocess.run(
                argv, stdout=subprocess.DEVNULL, stderr=log, timeout=timeout
            )
        except subprocess.TimeoutExpired:
            by_hand = shlex.join(["docker", "build", "-t", tag, str(model.path)])
            raise RunnerError(
                f"{model.name}: building {tag} did not finish in {timeout:g} s. "
                "A credential helper that never answers stalls a build this way; "
                "on Docker Desktop, look for a `docker-credential-desktop` process "
                "and end it. A first build of a large image can also take this "
                f"long: run `{by_hand}` once by hand, "
                "and the cached rebuild here takes seconds."
            ) from None
        if proc.returncode != 0:
            log.seek(0)
            tail = log.read().decode(errors="replace").strip().splitlines()[-20:]
            raise RunnerError(
                f"{model.name}: image build failed\n"
                + "\n".join("    " + line for line in tail)
            )
    return tag


def kill_container(docker: str, name: str) -> None:
    """End a container whose `docker run` was cut off.

    `docker kill` stops a running container, and --rm then removes it. A
    container that is not running, such as one created but never started,
    fails the kill, and `docker rm -f` removes it instead. A start the client
    sent before it died holds the container's lock in the daemon, the lock
    both commands take, so that start either finishes and is killed or is
    refused. Nothing is retried. A container neither command finds has
    already gone or was not yet created, and one created afterwards is never
    started, because `docker run` starts it from the client, which is dead.
    It takes no CPU, but stays until removed by hand. A command that cannot
    be run at all is skipped, so the caller's budget error is the one
    reported.
    """
    for argv in ([docker, "kill", name], [docker, "rm", "-f", name]):
        try:
            done = subprocess.run(
                argv, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                timeout=KILL_TIMEOUT_S,
            )
        except (OSError, subprocess.TimeoutExpired):
            continue
        if done.returncode == 0:
            return


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
    def command(docker: str, tag: str, model: ModelManifest, io_dir: Path, name: str) -> list[str]:
        """Build the run command.

        Isolation is part of the benchmark, not an operational detail, so this
        is asserted by the test suite rather than trusted. Inputs and the
        request are mounted read-only; only the output directory is writable.
        The image's own working directory is preserved so a relative entrypoint
        resolves exactly as it did when the image was built. The container
        takes the caller's `name`, so a run cut off at its time budget can be
        found and killed.
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
            "--name", name,
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
        name = container_name(model)
        argv = self.command(docker, tag, model, io_dir, name)
        try:
            proc = subprocess.run(
                argv, capture_output=True, text=True, timeout=probe.max_runtime_s
            )
        except subprocess.TimeoutExpired:
            # The timeout ends the docker client, not the container, which
            # runs on in the daemon until it finishes by itself; --rm removes
            # it only then. On a shared host that load pushes the next cases
            # over their budgets too, so the container is killed by name
            # before the budget error is raised.
            kill_container(docker, name)
            raise RunnerError(
                f"{model.name}: container exceeded the time budget for {probe.id}"
            ) from None

        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-12:]
            raise RunnerError(
                f"{model.name}: container exited {proc.returncode}\n"
                + "\n".join("    " + line for line in tail)
            )

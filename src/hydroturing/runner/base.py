"""Runner interface.

Two implementations, chosen by the model manifest.

`subprocess` is reserved for the reference models that ship in this repo.
They are trusted, they are pure Python, and running them without Docker keeps
the probe acceptance gate fast enough to sit on every pull request.

`docker` is what every submitted model uses. The container never sees the
probe code, the criteria, or the tolerances. It sees /io/input and nothing
else, so a model cannot read the threshold it is being judged against.
"""

from __future__ import annotations

import time
from abc import ABC, abstractmethod
from pathlib import Path

from hydroturing.protocol import Case, RunResult, read_result, stage
from hydroturing.spec import ModelManifest, ProbeSpec


class RunnerError(RuntimeError):
    """The model failed to produce a result."""


class Runner(ABC):
    name = "abstract"

    def run(self, model: ModelManifest, probe: ProbeSpec, case: Case, io_dir: Path) -> RunResult:
        io_dir.mkdir(parents=True, exist_ok=True)
        request_path = stage(io_dir, case, probe, model)
        started = time.monotonic()
        self.invoke(model, probe, io_dir, request_path)
        elapsed = time.monotonic() - started
        return read_result(io_dir, case, probe, elapsed)

    @abstractmethod
    def invoke(self, model: ModelManifest, probe: ProbeSpec, io_dir: Path, request_path: Path) -> None:
        ...


def get_runner(model: ModelManifest) -> Runner:
    from hydroturing.runner.docker_runner import DockerRunner
    from hydroturing.runner.subprocess_runner import SubprocessRunner

    if model.runner == "subprocess":
        return SubprocessRunner()
    return DockerRunner()

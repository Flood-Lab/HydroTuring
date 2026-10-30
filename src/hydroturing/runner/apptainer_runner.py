"""Apptainer (Singularity) runner, for hosts that cannot run a Docker daemon.

A shared cluster gives its users neither a daemon nor root, so `--runner
apptainer` puts the same model behind the same isolation with what such a host
does have. The model's Dockerfile stays the only recipe: it is translated here,
stage by stage, into Apptainer definition files, so a model is never packaged
twice and the two runtimes cannot drift apart in what they install.

What the translation holds to Docker's meaning is what a Dockerfile in this
repository uses: multi-stage builds, a stage built on an earlier stage,
`COPY --from`, ARG, ENV and WORKDIR. Anything else raises rather than building
an image that quietly differs. USER is ignored, because the container runs as
the invoking user on either runtime.

What the runtime cannot enforce is said here rather than implied. Apptainer
needs cgroup delegation to cap memory and CPU, which a batch node does not
grant, so the manifest's `resources` are not applied per container: the
scheduler's limits on the job are the cap. The rest of the isolation holds and
is asserted by the test suite.
"""

from __future__ import annotations

import fcntl
import glob
import hashlib
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from hydroturing.runner.base import Runner, RunnerError
from hydroturing.spec import ModelManifest, ProbeSpec

# Longer than Docker's limit because a first build here pulls every layer and
# writes a compressed image to a network filesystem. It still exists so that a
# stuck build is reported, not left to the batch system's time limit.
BUILD_TIMEOUT_S = 7200

# Where COPY's sources wait inside a stage while %post, which runs in order,
# reaches the instruction that places them. Removed before the stage is sealed.
CONTEXT_MOUNT = "/.ht_context"

# Docker instructions with no effect on what a model computes.
IGNORED = {"LABEL", "EXPOSE", "STOPSIGNAL", "MAINTAINER", "CMD", "ENTRYPOINT", "USER", "VOLUME"}

_ARG_REF = re.compile(r"\$(?:\{(\w+)\}|(\w+))")
_ENV_PAIR = re.compile(r"""(\w+)=((?:"[^"]*"|'[^']*'|[^\s"'])*)""")


def images_dir() -> Path:
    """`HT_APPTAINER_IMAGES`, since a cluster's home is rarely large enough."""
    configured = os.environ.get("HT_APPTAINER_IMAGES")
    return Path(configured) if configured else Path.home() / ".hydroturing" / "sif"


def require_apptainer() -> str:
    """The binary, under either of the names it is installed as."""
    binary = shutil.which("apptainer") or shutil.which("singularity")
    if binary is None:
        raise RunnerError(
            "neither apptainer nor singularity is on PATH. On a cluster they "
            "are often installed only on compute nodes, so run from a batch job."
        )
    return binary


@dataclass
class Instruction:
    keyword: str
    value: str


@dataclass
class Stage:
    base: str
    name: str | None
    instructions: list[Instruction] = field(default_factory=list)
    # ARG defaults declared before the first FROM, usable in FROM lines.
    args: dict[str, str] = field(default_factory=dict)


def logical_lines(text: str) -> list[str]:
    """Dockerfile lines with continuations joined and comments dropped.

    A comment line inside a continuation is removed without ending it, which
    is how the Dockerfiles here annotate a long RUN. Line breaks survive inside
    a joined instruction only as a space, as they do for Docker's parser.
    """
    lines: list[str] = []
    current = ""
    for raw in text.splitlines():
        stripped = raw.strip()
        if stripped.startswith("#"):
            continue
        if not stripped and not current:
            continue
        if stripped.endswith("\\"):
            current += stripped[:-1].rstrip() + " "
            continue
        current += stripped
        if current.strip():
            lines.append(current.strip())
        current = ""
    if current.strip():
        lines.append(current.strip())
    return lines


def parse_dockerfile(text: str) -> list[Stage]:
    stages: list[Stage] = []
    global_args: dict[str, str] = {}
    for line in logical_lines(text):
        parts = line.split(None, 1)
        keyword = parts[0].upper()
        value = parts[1] if len(parts) > 1 else ""
        if keyword == "FROM":
            words = [w for w in value.split() if not w.startswith("--")]
            base = _ARG_REF.sub(
                lambda m: global_args.get(m.group(1) or m.group(2), ""), words[0]
            )
            name = words[2] if len(words) == 3 and words[1].upper() == "AS" else None
            stages.append(Stage(base=base, name=name, args=dict(global_args)))
        elif not stages:
            if keyword != "ARG":
                raise RunnerError(f"Dockerfile: {keyword} before the first FROM")
            key, _, default = value.partition("=")
            global_args[key.strip()] = default.strip().strip('"')
        elif keyword in IGNORED:
            continue
        elif keyword in {"RUN", "ENV", "ARG", "WORKDIR", "COPY", "ADD"}:
            stages[-1].instructions.append(Instruction(keyword, value))
        else:
            raise RunnerError(f"Dockerfile: {keyword} has no Apptainer translation")
    if not stages:
        raise RunnerError("Dockerfile: no FROM")
    return stages


def env_pairs(value: str) -> list[tuple[str, str]]:
    """`ENV A=1 B=2` and the older `ENV A 1`, with the value's quoting intact."""
    if "=" not in value.split(None, 1)[0]:
        key, _, rest = value.partition(" ")
        return [(key, shlex.quote(rest.strip()))]
    # Each value stays as the shell will read it, quotes included, so that
    # `PATH=/opt/conda/bin:$PATH` still expands and "a b" stays one value.
    return _ENV_PAIR.findall(value)


def split_copy(value: str) -> tuple[dict[str, str], list[str], str]:
    """COPY's flags, sources and destination, in shell or JSON form."""
    flags: dict[str, str] = {}
    words = value.split()
    while words and words[0].startswith("--"):
        key, _, val = words.pop(0)[2:].partition("=")
        flags[key] = val
    rest = " ".join(words)
    parts = json.loads(rest) if rest.startswith("[") else shlex.split(rest)
    if len(parts) < 2:
        raise RunnerError(f"Dockerfile: COPY needs a source and a destination: {value}")
    return flags, parts[:-1], parts[-1]


def context_digest(model: ModelManifest) -> str:
    """A digest of everything a build reads, so a changed adapter rebuilds.

    Docker's layer cache does this for free. The README and manifest are left
    out: a build reads neither, and editing prose should not cost a rebuild.
    """
    digest = hashlib.sha256()
    for path in sorted(p for p in Path(model.path).rglob("*") if p.is_file()):
        relative = path.relative_to(model.path).as_posix()
        if relative in {"README.md", "model.yaml"} or "__pycache__" in relative:
            continue
        digest.update(relative.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()[:12]


def image_path(model: ModelManifest) -> Path:
    return images_dir() / f"{model.name}_{model.version}_{context_digest(model)}.sif"


def lineage(stage: Stage, stages: list[Stage]) -> list[Stage]:
    """A stage and the earlier stages it is built on, oldest first."""
    by_name = {s.name: s for s in stages if s.name}
    chain = [stage]
    while chain[0].base in by_name:
        chain.insert(0, by_name[chain[0].base])
    return chain


def stage_workdir(stage: Stage, stages: list[Stage]) -> str:
    """The directory a stage ends in, which one built on it starts from."""
    for ancestor in reversed(lineage(stage, stages)):
        for instruction in reversed(ancestor.instructions):
            if instruction.keyword == "WORKDIR":
                return _expand(instruction.value, ancestor)
    return "/"


def final_workdir(stages: list[Stage]) -> str:
    """The directory a relative entrypoint resolves against.

    `apptainer exec` does not keep an image's working directory, and the
    manifests rely on it, so it is read from the recipe and passed as --pwd.
    """
    return stage_workdir(stages[-1], stages)


def _expand(value: str, stage: Stage) -> str:
    """Substitute the ARG and ENV values a stage has set, before a shell does."""
    known = dict(stage.args)
    for instruction in stage.instructions:
        if instruction.keyword == "ARG":
            key, _, default = instruction.value.partition("=")
            known[key.strip()] = default.strip().strip('"')
        elif instruction.keyword == "ENV":
            known.update(env_pairs(instruction.value))
    return _ARG_REF.sub(lambda m: known.get(m.group(1) or m.group(2), m.group(0)), value)


class ImageBuilder:
    """Turns one Dockerfile into one image, a stage at a time."""

    def __init__(self, model: ModelManifest, apptainer: str, build_dir: Path, log) -> None:
        self.model = model
        self.apptainer = apptainer
        self.build_dir = build_dir
        self.log = log
        self.stages = parse_dockerfile((Path(model.path) / "Dockerfile").read_text())
        self.images: dict[str, Path] = {}

    def needed(self) -> list[Stage]:
        """The final stage and whatever it reads from, in build order."""
        by_name = {s.name: s for s in self.stages if s.name}
        wanted = {id(self.stages[-1])}
        for stage in reversed(self.stages):
            if id(stage) not in wanted:
                continue
            if stage.base in by_name:
                wanted.add(id(by_name[stage.base]))
            for instruction in stage.instructions:
                if instruction.keyword == "COPY":
                    source = split_copy(instruction.value)[0].get("from")
                    if source in by_name:
                        wanted.add(id(by_name[source]))
        return [s for s in self.stages if id(s) in wanted]

    def build(self, destination: Path) -> None:
        stages = self.needed()
        for index, stage in enumerate(stages):
            target = destination if stage is stages[-1] else self.build_dir / f"stage{index}.sif"
            self.build_stage(stage, index, target)
            if stage.name:
                self.images[stage.name] = target

    def build_stage(self, stage: Stage, index: int, target: Path) -> None:
        definition = self.definition(stage, index)
        def_path = self.build_dir / f"stage{index}.def"
        def_path.write_text(definition)
        self.call(
            [self.apptainer, "build", "--fakeroot", "--force", str(target), str(def_path)],
            f"stage {stage.name or index}",
        )

    def definition(self, stage: Stage, index: int) -> str:
        files: list[str] = []
        post: list[str] = [
            # A base image's own environment, which %post does not start with.
            'for f in /.singularity.d/env/10-docker2singularity.sh '
            '/.singularity.d/env/90-environment.sh; do '
            '[ -f "$f" ] && . "$f"; done; true',
        ]
        environment: list[str] = []
        # A stage built on an earlier one inherits its ENV and directory. The
        # image's own environment file may not carry them through a second
        # definition's %environment, so they are written out again.
        ancestors = lineage(stage, self.stages)[:-1]
        workdir = stage_workdir(ancestors[-1], self.stages) if ancestors else "/"
        for ancestor in ancestors:
            for instruction in ancestor.instructions:
                if instruction.keyword == "ENV":
                    environment += [f"export {k}={v}" for k, v in env_pairs(instruction.value)]
        count = 0
        for instruction in stage.instructions:
            keyword, value = instruction.keyword, instruction.value
            if keyword == "ARG":
                key, _, default = value.partition("=")
                post.append(f"export {key.strip()}={shlex.quote(default.strip().strip(chr(34)))}")
            elif keyword == "ENV":
                for key, val in env_pairs(value):
                    post.append(f"export {key}={val}")
                    environment.append(f"export {key}={val}")
            elif keyword == "WORKDIR":
                workdir = _expand(value, stage)
                post.append(f"mkdir -p {shlex.quote(workdir)}")
            elif keyword == "RUN":
                command = value
                if command.startswith("["):
                    command = shlex.join(json.loads(command))
                # Each RUN in its own subshell, starting where WORKDIR says,
                # because a `cd` in one RUN does not reach the next in Docker.
                post.append(f"cd {shlex.quote(workdir)}\n( {command}\n)")
            elif keyword in {"COPY", "ADD"}:
                count = self.stage_copy(value, stage, index, count, workdir, files, post)
        post.append(f"rm -rf {CONTEXT_MOUNT}")
        lines: list[str] = []
        if stage.base in self.images:
            lines += ["Bootstrap: localimage", f"From: {self.images[stage.base]}"]
        else:
            if stage.base.lower() == "scratch":
                raise RunnerError("Dockerfile: FROM scratch has no Apptainer translation")
            lines += ["Bootstrap: docker", f"From: {stage.base}"]
        lines.append("")
        if files:
            lines += ["%files", *files, ""]
        lines += ["%post", *post, ""]
        if environment:
            lines += ["%environment", *environment, ""]
        return "\n".join(lines)

    def stage_copy(
        self, value: str, stage: Stage, index: int, count: int, workdir: str,
        files: list[str], post: list[str],
    ) -> int:
        flags, sources, dest = split_copy(value)
        dest = _expand(dest, stage)
        if not dest.startswith("/"):
            dest = os.path.join(workdir, dest) + ("/" if dest.endswith("/") else "")
        into_dir = dest.endswith("/") or len(sources) > 1
        for source in sources:
            source = _expand(source, stage)
            if flags.get("from"):
                archives = self.export_from_stage(flags["from"], source, index, count)
            else:
                archives = self.pack_from_context(source, index, count)
            for archive, is_dir, basename in archives:
                member = f"{CONTEXT_MOUNT}/{archive.name}"
                files.append(f"{archive} {member}")
                quoted = shlex.quote(dest)
                if is_dir or into_dir:
                    # Docker copies a directory's contents, not the directory.
                    post.append(f"mkdir -p {quoted} && tar --no-same-owner -xf {member} -C {quoted}")
                else:
                    # A single file onto a path: unpack beside it, then rename.
                    parent = os.path.dirname(dest) or "/"
                    unpack = f"mkdir -p {shlex.quote(parent)} && tar --no-same-owner -xf {member} -C {shlex.quote(parent)}"
                    if basename != os.path.basename(dest):
                        unpack += f" && mv -f {shlex.quote(os.path.join(parent, basename))} {quoted}"
                    post.append(unpack)
                count += 1
        return count

    def pack_from_context(self, source: str, index: int, count: int) -> list[tuple[Path, bool, str]]:
        root = Path(self.model.path)
        matches = sorted(glob.glob(str(root / source)))
        if not matches:
            raise RunnerError(f"{self.model.name}: COPY source {source} is not in the build context")
        archives = []
        for offset, match in enumerate(matches):
            path = Path(match)
            archive = self.build_dir / f"context{index}_{count}_{offset}.tar"
            if path.is_dir():
                argv = ["tar", "-C", str(path), "-cf", str(archive), "."]
            else:
                argv = ["tar", "-C", str(path.parent), "-cf", str(archive), path.name]
            self.call(argv, f"packing {source}")
            archives.append((archive, path.is_dir(), path.name))
        return archives

    def export_from_stage(self, stage_name: str, source: str, index: int, count: int) -> list[tuple[Path, bool, str]]:
        image = self.images.get(stage_name)
        if image is None:
            raise RunnerError(f"{self.model.name}: COPY --from={stage_name}: no such stage built")
        archive = self.build_dir / f"context{index}_{count}_0.tar"
        script = (
            f's={shlex.quote(source)}; '
            'if [ -d "$s" ]; then tar -C "$s" -cf /__out/a.tar . && echo dir; '
            'else tar -C "$(dirname "$s")" -cf /__out/a.tar "$(basename "$s")" && echo file; fi'
        )
        with tempfile.TemporaryDirectory(dir=self.build_dir) as out:
            done = self.call(
                [self.apptainer, "exec", "--fakeroot", "--containall", "--no-home",
                 "--bind", f"{out}:/__out", str(image), "sh", "-c", script],
                f"exporting {source} from {stage_name}",
            )
            shutil.move(os.path.join(out, "a.tar"), archive)
        is_dir = done.strip().splitlines()[-1].strip() == "dir"
        return [(archive, is_dir, os.path.basename(source.rstrip("/")))]

    def call(self, argv: list[str], what: str) -> str:
        self.log.write(f"\n$ {shlex.join(argv)}\n")
        self.log.flush()
        proc = subprocess.run(argv, stdout=subprocess.PIPE, stderr=self.log, text=True)
        self.log.write(proc.stdout or "")
        if proc.returncode != 0:
            self.log.flush()
            self.log.seek(0)
            tail = self.log.read().strip().splitlines()[-25:]
            raise RunnerError(
                f"{self.model.name}: {what} failed (exit {proc.returncode})\n"
                + "\n".join("    " + line for line in tail)
            )
        return proc.stdout or ""


def build(model: ModelManifest, apptainer: str | None = None) -> Path:
    """The model's image, built once.

    A lock keeps two jobs on one filesystem from building the same image, and
    the rename at the end means a finished file is always a whole image.
    """
    apptainer = apptainer or require_apptainer()
    dockerfile = Path(model.path) / "Dockerfile"
    if not dockerfile.exists():
        raise RunnerError(f"{model.name}: no Dockerfile at {dockerfile}")
    target = image_path(model)
    if target.exists():
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(target.with_suffix(".lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if target.exists():
            return target
        with tempfile.TemporaryDirectory(dir=target.parent, prefix=f".build-{model.name}-") as tmp:
            build_dir = Path(tmp)
            partial = build_dir / "image.sif"
            log_path = target.with_suffix(".build.log")
            with open(log_path, "w+") as log:
                ImageBuilder(model, apptainer, build_dir, log).build(partial)
            os.replace(partial, target)
    return target


class ApptainerRunner(Runner):
    name = "apptainer"

    def __init__(self) -> None:
        self._image: Path | None = None
        self._workdir: str | None = None

    @staticmethod
    def command(
        apptainer: str, image: Path | str, workdir: str, model: ModelManifest, io_dir: Path
    ) -> list[str]:
        """Build the run command: Docker's isolation, in Apptainer's terms.

        --containall leaves the host's home, working directory, /tmp and
        environment outside. --no-mount removes the paths the site's
        configuration would bind for everyone, without which a cluster's
        project and scratch filesystems, this repository's probes among them,
        would be readable from inside. --network none leaves only loopback. A
        SIF is read-only, and the image runs as the invoking user with no
        capabilities. /tmp is a size-limited tmpfs that vanishes with the run.
        """
        request_path = (io_dir / "request.json").resolve()
        input_dir = (io_dir / "input").resolve()
        output_dir = (io_dir / "output").resolve()
        return [
            apptainer, "exec",
            "--containall",
            "--no-home",
            "--no-mount", "bind-paths,home,cwd,hostfs",
            "--no-privs",
            "--net", "--network", "none",
            "--pwd", workdir,
            "--bind", f"{request_path}:/io/request.json:ro",
            "--bind", f"{input_dir}:/io/input:ro",
            "--bind", f"{output_dir}:/io/output:rw",
            str(image),
            # Apptainer refuses to set HOME through --env, and a library that
            # wants a cache directory has to find one that vanishes with the run.
            "env", "HOME=/tmp",
            *model.entrypoint, "--request", "/io/request.json",
        ]

    def invoke(self, model: ModelManifest, probe: ProbeSpec, io_dir: Path, request_path: Path) -> None:
        apptainer = require_apptainer()
        if self._image is None:
            self._image = build(model, apptainer=apptainer)
            self._workdir = final_workdir(
                parse_dockerfile((Path(model.path) / "Dockerfile").read_text())
            )
        argv = self.command(apptainer, self._image, self._workdir or "/", model, io_dir)
        # In its own session, so the whole process tree can be signalled:
        # with --containall the model is PID 1's child in its own namespace and
        # dies with it, but that takes the starter being killed first.
        proc = subprocess.Popen(
            argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            start_new_session=True,
        )
        try:
            stdout, stderr = proc.communicate(timeout=probe.max_runtime_s)
        except subprocess.TimeoutExpired:
            _kill_group(proc)
            raise RunnerError(
                f"{model.name}: container exceeded the time budget for {probe.id}"
            ) from None
        except BaseException:
            _kill_group(proc)
            raise
        if proc.returncode != 0:
            tail = (stderr or stdout or "").strip().splitlines()[-12:]
            raise RunnerError(
                f"{model.name}: container exited {proc.returncode}\n"
                + "\n".join("    " + line for line in tail)
            )


def _kill_group(proc: subprocess.Popen) -> None:
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass
    proc.communicate()

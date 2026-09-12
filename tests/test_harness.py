"""Tests for the harness itself, distinct from the probe acceptance gate.

The gate asks whether the benchmark discriminates. These ask whether the
machinery underneath it is sound: generators are deterministic, the window
arithmetic is not off by one, and the contract is enforced.
"""

from __future__ import annotations

import json
from dataclasses import replace

import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.harness import (
    IncompatibleError,
    build_case,
    run_model,
    run_probe,
    verify_adapter_contract,
)
from hydroturing.protocol import ProtocolError, read_result, stage
from hydroturing.scoring import ERROR, FAIL, INCOMPATIBLE, INCOMPLETE, NOT_SCORED, OK, PASS, VIOLATION
from hydroturing.spec import SpecError
from hydroturing.seeds import gate_seeds


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/catchment-closure")


def test_every_probe_declares_a_discriminating_baseline():
    """A probe with no must_fail baseline proves nothing."""
    for spec in registry.all_probes():
        assert spec.must_pass, f"{spec.id} declares no must_pass model"
        assert spec.must_fail, f"{spec.id} declares no must_fail model"


def test_generator_is_deterministic(probe):
    """Same seed, byte-identical forcing. Without this a failure cannot be
    reproduced, which makes a randomised benchmark unusable in a dispute."""
    first, static_a = build_case(probe, 7), None
    second = build_case(probe, 7)
    pd.testing.assert_frame_equal(first.forcing, second.forcing)
    assert first.static == second.static


def test_different_seeds_give_different_forcing(probe):
    a, b = build_case(probe, 7), build_case(probe, 8)
    assert not a.forcing["pr"].equals(b.forcing["pr"])


def test_generated_length_matches_the_spec(probe):
    case = build_case(probe, 1)
    assert len(case.forcing) == int(probe.period_years * 365) + probe.spinup_days


def test_storage_change_uses_the_step_before_the_window(probe):
    """The classic off-by-one: taking the initial storage from the first
    scored step instead of the last spinup step manufactures a residual."""
    from hydroturing.criteria.base import make_window
    from hydroturing.runner import get_runner
    import tempfile
    from pathlib import Path

    model = registry.find_model("reference_bucket")
    case = build_case(probe, 3)
    run = get_runner(model).run(model, probe, case, Path(tempfile.mkdtemp()))
    window = make_window(run, probe)
    assert window.state0.equals(run.table.iloc[probe.spinup_days - 1])
    assert len(window.table) == len(run.table) - probe.spinup_days


def test_exact_model_passes(probe):
    outcome = run_probe(registry.find_model("reference_bucket"), probe, [11, 12])
    assert outcome.verdict == PASS, outcome.failing


def test_silent_sink_is_caught_by_closure(probe):
    outcome = run_probe(registry.find_model("reference_leaky"), probe, [11])
    assert outcome.verdict == FAIL
    assert outcome.reason == VIOLATION
    assert "closure" in outcome.failing


def test_closure_by_construction_survives_closure_but_not_bounds(probe):
    """The point of the whole design. Randomising the forcing cannot catch a
    model that solves for storage as the residual: its closure is exact on
    every seed. Only the physical bounds on that invented storage catch it."""
    outcome = run_probe(registry.find_model("reference_cheater"), probe, [11, 99, 12345])
    closure = next(c for c in outcome.criteria if c.name == "closure")
    assert closure.passed, "cheater should pass closure; that is the whole problem"
    assert "state_bounds" in outcome.failing
    assert outcome.verdict == FAIL


def test_degenerate_model_is_caught(probe):
    outcome = run_probe(registry.find_model("reference_degenerate"), probe, [11])
    assert "non_degenerate" in outcome.failing


def test_missing_variables_report_incomplete_not_violation(probe):
    """A streamflow-only model has not violated conservation. It has failed to
    say enough to be checked, and the report must not conflate the two: the
    probe is not scored, so it is neither a pass nor a fail."""
    outcome = run_probe(registry.find_model("reference_streamflow_only"), probe, [11])
    assert outcome.verdict == NOT_SCORED
    assert outcome.reason == INCOMPLETE
    assert "evspsbl" in outcome.missing


def test_incompatible_timestep_is_reported_before_execution(probe):
    model = replace(registry.find_model("reference_bucket"), timesteps=("PT1H",))
    outcome = run_probe(model, probe, [11])
    assert outcome.verdict == NOT_SCORED
    assert outcome.reason == INCOMPATIBLE
    assert "timestep" in outcome.incompatible[0]


def _rolled_up(*outcomes):
    """The model verdict, reason and summary over probes with these outcomes."""
    from hydroturing.scoring import ModelReport, ProbeOutcome

    probes = [ProbeOutcome(f"mass/p{i}", "mass", v, r) for i, (v, r) in enumerate(outcomes)]
    report = ModelReport("m", "1", "0.1.0", probes)
    return report.verdict, report.reason, report.summary


def test_a_probe_that_cannot_be_put_to_the_model_counts_neither_way():
    """An N/A probe asked the model nothing. It must not fail a model that
    passes everything else, must not mask a violation or an error, and is not
    in the total the passes are counted out of."""
    assert _rolled_up((PASS, OK), (NOT_SCORED, INCOMPLETE)) == (PASS, OK, "1/1 probes passed")
    assert _rolled_up(
        (FAIL, VIOLATION), (PASS, OK), (NOT_SCORED, INCOMPLETE), (NOT_SCORED, INCOMPATIBLE)
    ) == (FAIL, VIOLATION, "1/2 probes passed")
    assert _rolled_up((FAIL, ERROR), (FAIL, VIOLATION), (NOT_SCORED, INCOMPLETE))[:2] == (FAIL, ERROR)


def test_a_model_no_probe_could_score_has_not_passed():
    """Reporting nothing checkable must not earn a PASS, nor a count out of nothing."""
    assert _rolled_up((NOT_SCORED, INCOMPATIBLE), (NOT_SCORED, INCOMPLETE)) == (
        NOT_SCORED, INCOMPLETE, "no probe could be scored"
    )
    assert _rolled_up((NOT_SCORED, INCOMPATIBLE))[:2] == (NOT_SCORED, INCOMPATIBLE)
    assert _rolled_up()[:2] == (FAIL, ERROR)


def test_run_exits_1_for_a_model_no_probe_could_score(capsys):
    """N/A is not a PASS, so CI keeps it in the scorecard beside the FAILs;
    only ERROR exits 2."""
    from hydroturing.cli import main

    assert main(["run", "--model", "reference_streamflow_only", "--probe", "mass/catchment-closure"]) == 1
    assert "N/A (INCOMPLETE)" in capsys.readouterr().out


def test_missing_forcing_is_reported_as_incompatible(probe):
    model = replace(
        registry.find_model("reference_bucket"),
        needs_forcing=("pr", "unavailable_driver"),
    )
    outcome = run_probe(model, probe, [11])
    assert outcome.verdict == NOT_SCORED
    assert outcome.reason == INCOMPATIBLE
    assert "unavailable_driver" in outcome.incompatible[0]


def test_paired_probe_requires_declared_perturbation_support(probe):
    paired = replace(probe, variants=("control", "perturbed"))
    model = replace(registry.find_model("reference_bucket"), supports_perturbation=False)
    outcome = run_probe(model, paired, [11])
    assert outcome.verdict == NOT_SCORED
    assert outcome.reason == INCOMPATIBLE
    assert "perturbation" in outcome.incompatible[0]


def test_a_seed_the_model_cannot_consume_stops_the_probe_before_any_run(probe, monkeypatch):
    """N/A has to mean the probe asked the model nothing. A case found
    unusable on a later seed must stop the probe before the model runs on any
    seed; otherwise failures already measured on the earlier seeds would be
    discarded, and a failing model could pass."""
    from hydroturing import harness

    leaky = registry.find_model("reference_leaky")
    seeds = gate_seeds(probe.id, 2)
    assert run_probe(leaky, probe, seeds).verdict == FAIL

    real_issues, real_runner = harness.compatibility_issues, harness.get_runner
    ran = []

    def issues(model, probe, case=None, **kwargs):
        if case is not None and case.seed == seeds[1]:
            return ["this seed's case cannot be consumed"]
        return real_issues(model, probe, case, **kwargs)

    class Counting:
        def __init__(self, inner):
            self.inner = inner

        def run(self, model, probe, case, io_dir):
            ran.append(case.seed)
            return self.inner.run(model, probe, case, io_dir)

    monkeypatch.setattr(harness, "compatibility_issues", issues)
    monkeypatch.setattr(harness, "get_runner", lambda model: Counting(real_runner(model)))
    outcome = run_probe(leaky, probe, seeds)
    assert (outcome.verdict, outcome.reason) == (NOT_SCORED, INCOMPATIBLE)
    assert ran == []


def test_an_exception_without_a_message_is_still_an_error(probe, monkeypatch):
    """A bare assert or an exhausted next() carries no message. It is still
    the machinery failing, so the reason is ERROR, not a FAIL with reason OK."""
    from hydroturing import harness

    class Broken:
        def run(self, model, probe, case, io_dir):
            raise AssertionError()

    monkeypatch.setattr(harness, "get_runner", lambda model: Broken())
    outcome = run_probe(registry.find_model("reference_bucket"), probe, [11])
    assert (outcome.verdict, outcome.reason, outcome.error) == (FAIL, ERROR, "AssertionError")


def test_run_exits_2_when_an_error_sits_beside_unscored_probes(monkeypatch):
    """Unscored probes must not hide an error: ERROR still exits 2."""
    from hydroturing import cli
    from hydroturing.scoring import ModelReport, ProbeOutcome

    report = ModelReport("m", "1", "0.1.0", [
        ProbeOutcome("mass/a", "mass", NOT_SCORED, INCOMPLETE),
        ProbeOutcome("mass/b", "mass", FAIL, ERROR, error="adapter crashed"),
    ])
    monkeypatch.setattr(cli, "run_model", lambda *args, **kwargs: report)
    assert cli.main(["run", "--model", "reference_bucket", "--probe", "mass/catchment-closure"]) == 2


def test_a_command_that_crashes_exits_2(monkeypatch, capsys):
    """Python exits 1 on an uncaught exception, and CI accepts exit 1 as a
    FAIL or an N/A. A crash nothing anticipated, such as a model.yaml that is
    not YAML, is the harness failing and must exit 2 with its traceback."""
    import yaml

    from hydroturing import cli

    def unreadable(name):
        raise yaml.YAMLError(f"{name}/model.yaml is not YAML")

    monkeypatch.setattr(cli.registry, "find_model", unreadable)
    assert cli.main(["verify-adapter", "--model", "reference_bucket"]) == 2
    assert cli.main(["run", "--model", "reference_bucket"]) == 2
    assert "reference_bucket/model.yaml is not YAML" in capsys.readouterr().err


def test_adapter_verification_runs_an_incomplete_model(probe):
    """A probe that is N/A (INCOMPLETE) must not skip the contract smoke test."""
    model = registry.find_model("reference_streamflow_only")
    result = verify_adapter_contract(model, probe, gate_seeds(probe.id, 1)[0])
    assert len(result.table) == result.case.n_steps
    assert {"time", "mrro", "dis"} <= set(result.table.columns)


def test_adapter_verification_does_not_run_a_probe_the_model_cannot_consume(probe, monkeypatch):
    """A model needing a forcing the probe does not generate cannot be put to
    it. That is the probe's N/A (INCOMPATIBLE), decided before the adapter is
    invoked, and must not surface as the adapter breaking the contract."""
    from hydroturing import harness

    def no_runner(model):
        raise AssertionError(f"the adapter of {model.name} was run")

    monkeypatch.setattr(harness, "get_runner", no_runner)
    model = replace(
        registry.find_model("reference_bucket"),
        needs_forcing=("pr", "unavailable_driver"),
    )
    with pytest.raises(IncompatibleError) as refused:
        verify_adapter_contract(model, probe, gate_seeds(probe.id, 1)[0])
    assert refused.value.issues == ["forcing does not provide unavailable_driver"]


def test_all_seeds_must_pass(probe):
    """A model cannot get through on a lucky draw."""
    outcome = run_probe(registry.find_model("reference_leaky"), probe, gate_seeds(probe.id, 5))
    closure = next(c for c in outcome.criteria if c.name == "closure")
    assert len(closure.per_seed) == 5
    assert all(entry["status"] == "fail" for entry in closure.per_seed)


def test_short_result_is_rejected(probe, tmp_path):
    """Contract enforcement: a truncated table must not be silently scored."""
    case = build_case(probe, 5)
    (tmp_path / "output").mkdir(parents=True)
    frame = pd.DataFrame({
        "time": case.forcing["time"][:10],
        "pr": 0.0, "evspsbl": 0.0, "mrro": 0.0,
        "mrso": 0.0, "snw": 0.0, "canopy": 0.0,
    })
    frame.to_csv(tmp_path / "output" / "result.csv", index=False)
    with pytest.raises(ProtocolError, match="rows"):
        read_result(tmp_path, case, probe, 0.0)


def _write_valid_result(path, case):
    frame = pd.DataFrame({
        "time": case.forcing["time"],
        "pr": case.forcing["pr"],
        "evspsbl": 0.0,
        "mrro": 0.0,
        "mrso": 0.0,
        "snw": 0.0,
        "canopy": 0.0,
    })
    (path / "output").mkdir(parents=True)
    frame.to_csv(path / "output" / "result.csv", index=False)
    (path / "output" / "run.json").write_text(
        json.dumps({"status": "ok", "n_steps": case.n_steps})
    )


def test_run_metadata_is_required(probe, tmp_path):
    case = build_case(probe, 5)
    _write_valid_result(tmp_path, case)
    (tmp_path / "output" / "run.json").unlink()
    with pytest.raises(ProtocolError, match="required output/run.json"):
        read_result(tmp_path, case, probe, 0.0)


def test_shifted_time_axis_is_rejected(probe, tmp_path):
    case = build_case(probe, 5)
    _write_valid_result(tmp_path, case)
    frame = pd.read_csv(tmp_path / "output" / "result.csv")
    frame.loc[20, "time"] = "2099-01-01"
    frame.to_csv(tmp_path / "output" / "result.csv", index=False)
    with pytest.raises(ProtocolError, match="time axis"):
        read_result(tmp_path, case, probe, 0.0)


def test_output_size_limit_is_enforced(probe, tmp_path):
    case = build_case(probe, 5)
    _write_valid_result(tmp_path, case)
    tiny_limit = replace(probe, max_output_mb=0.001)
    with pytest.raises(ProtocolError, match="exceeding"):
        read_result(tmp_path, case, tiny_limit, 0.0)


def test_request_hides_probe_identity_and_generator_seed(probe, tmp_path):
    case = build_case(probe, 123456)
    model = registry.find_model("reference_bucket")
    request_path = stage(tmp_path, case, probe, model)
    request = json.loads(request_path.read_text())

    assert probe.id not in request["case_id"]
    assert request["seed"] != case.seed
    assert "spinup_steps" not in request
    assert request["request"]["fluxes"] == list(model.emits_fluxes)
    assert request["request"]["states"] == list(model.emits_states)


# --- container isolation ----------------------------------------------------
# The daemon is not available in every dev environment, so the isolation flags
# are asserted structurally here and the real build-and-run happens in CI,
# where Docker exists. Isolation is part of the benchmark: a model that could
# read the probe definition could read the tolerance it is judged against.


def test_container_runs_with_hardened_read_only_inputs(tmp_path):
    from hydroturing.runner.docker_runner import DockerRunner

    model = registry.find_model("reference_bucket")
    argv = DockerRunner.command("docker", "img:1", model, tmp_path, "hydroturing-test")

    assert argv[:3] == ["docker", "run", "--rm"]
    assert "--network" in argv and argv[argv.index("--network") + 1] == "none"
    assert "--read-only" in argv
    assert argv[argv.index("--cap-drop") + 1] == "ALL"
    # Not root: with every capability dropped, root could not even write the
    # output directory the host user owns, and it should not be able to write
    # anything else of the host's either. A host without POSIX ids (Windows)
    # passes no --user; Docker Desktop's file sharing handles ownership there.
    import os
    if hasattr(os, "getuid"):
        assert argv[argv.index("--user") + 1] == f"{os.getuid()}:{os.getgid()}"
    else:
        assert "--user" not in argv
    assert argv[argv.index("--env") + 1] == "HOME=/tmp"
    assert argv[argv.index("--security-opt") + 1] == "no-new-privileges"
    assert "--pids-limit" in argv
    assert "--workdir" not in argv

    mounts = [argv[i + 1] for i, a in enumerate(argv) if a == "--mount"]
    assert len(mounts) == 3
    assert any("target=/io/request.json,readonly" in mount for mount in mounts)
    assert any("target=/io/input,readonly" in mount for mount in mounts)
    assert any("target=/io/output" in mount and "readonly" not in mount for mount in mounts)

    assert "--memory" in argv and "--cpus" in argv
    assert argv[-2:] == ["--request", "/io/request.json"]
    assert not any(str(registry.PROBES_DIR) in a for a in argv), (
        "the probes directory must never be mounted into a model container"
    )


def test_cpu_request_is_capped_at_the_host(monkeypatch, tmp_path):
    """Docker refuses a CPU quota above the host's count outright. A manifest
    asking for eight cores must still run on a four-core runner."""
    from dataclasses import replace

    from hydroturing.runner import docker_runner
    from hydroturing.runner.docker_runner import DockerRunner

    model = replace(registry.find_model("reference_bucket"), resources={"cpu": 8, "memory_gb": 8})
    monkeypatch.setattr(docker_runner.os, "cpu_count", lambda: 4)
    argv = DockerRunner.command("docker", "img:1", model, tmp_path, "hydroturing-test")
    assert argv[argv.index("--cpus") + 1] == "4"

    monkeypatch.setattr(docker_runner.os, "cpu_count", lambda: 16)
    argv = DockerRunner.command("docker", "img:1", model, tmp_path, "hydroturing-test")
    assert argv[argv.index("--cpus") + 1] == "8"


def test_missing_daemon_is_reported_as_such(monkeypatch):
    """A missing binary and a stopped daemon need different fixes, so they
    must not both surface as 'image build failed'."""
    import subprocess as sp

    from hydroturing.runner import docker_runner
    from hydroturing.runner.base import RunnerError

    monkeypatch.setattr(docker_runner.shutil, "which", lambda _: "/usr/bin/docker")
    monkeypatch.setattr(
        docker_runner.subprocess, "run",
        lambda *a, **k: sp.CompletedProcess(a, 1, "", "cannot connect"),
    )
    with pytest.raises(RunnerError, match="daemon is not reachable"):
        docker_runner.require_docker()


def fake_docker(tmp_path, script):
    """Write `script` as an executable `docker`, so build() runs without a daemon."""
    import os
    import shutil

    sh = shutil.which("sh")
    if os.name != "posix" or sh is None:
        pytest.skip("the fake docker is a POSIX shell script")
    docker = tmp_path / "docker"
    docker.write_text(f"#!{sh}\n{script}\n")
    docker.chmod(0o755)
    return str(docker)


def test_image_build_returns_when_docker_exits_not_when_its_output_closes(tmp_path):
    """On Docker Desktop, `docker build` starts `docker-credential-desktop
    get`, and that helper can be orphaned still holding the stderr it
    inherited. Reading the build's output through a pipe then waits for an
    end of file that never comes: `ht run` sat for seven minutes after the
    image was built. This docker exits 0 at once and leaves such a child."""
    import os
    import shlex
    import signal
    import threading
    import time

    from hydroturing.runner import docker_runner

    model = registry.find_model("reference_bucket")
    pidfile = tmp_path / "helper.pid"
    docker = fake_docker(tmp_path, f"sleep 30 &\necho $! > {shlex.quote(str(pidfile))}\nexit 0")

    def end_helper():
        # Once only, so a pid the system has since reused is never signalled.
        try:
            pid = int(pidfile.read_text())
            pidfile.unlink()
            os.kill(pid, signal.SIGTERM)
        except (OSError, ValueError):
            pass

    # A regression blocks until the helper exits. The watchdog ends it after
    # ten seconds, so the test fails on the assertion rather than hanging.
    watchdog = threading.Timer(10, end_helper)
    watchdog.start()
    started = time.monotonic()
    try:
        tag = docker_runner.build(model, docker=docker)
    finally:
        elapsed = time.monotonic() - started
        watchdog.cancel()
        watchdog.join()
        end_helper()

    assert tag == docker_runner.image_tag(model)
    assert elapsed < 5, f"build() returned after {elapsed:.1f} s, waiting on the orphan"


def test_image_build_that_never_exits_names_the_image_and_the_credential_helper(tmp_path):
    """A build can also stall inside docker, waiting on a credential helper
    that never answers. It ends at the timeout with a reason, not silently."""
    from hydroturing.runner import docker_runner
    from hydroturing.runner.base import RunnerError

    model = registry.find_model("reference_bucket")
    docker = fake_docker(tmp_path, "exec sleep 30")
    with pytest.raises(RunnerError, match="docker-credential-desktop") as raised:
        docker_runner.build(model, docker=docker, timeout=0.5)
    assert docker_runner.image_tag(model) in str(raised.value)


def test_failed_image_build_still_reports_the_end_of_its_log(tmp_path):
    """The build log is no longer read from a pipe; a failed build must
    still say why it failed."""
    from hydroturing.runner import docker_runner
    from hydroturing.runner.base import RunnerError

    model = registry.find_model("reference_bucket")
    docker = fake_docker(tmp_path, "echo 'step 1/1: pip install failed' >&2\nexit 1")
    with pytest.raises(RunnerError, match="image build failed\n    step 1/1: pip install failed"):
        docker_runner.build(model, docker=docker)


@pytest.mark.parametrize("kill_status", [0, 1], ids=["killed", "created-not-started"])
def test_container_past_its_time_budget_is_killed_not_left_running(monkeypatch, tmp_path, kill_status):
    """The timeout ends the docker client, not the container: after
    extreme-rain timed out for google_flood_forecast, its container was still
    running five minutes later and taking CPU from the cases behind it. A
    timed-out run kills the container by the name it ran under, and removes
    it by force when the kill fails, as on one created but never started."""
    import shlex

    from hydroturing.runner import docker_runner
    from hydroturing.runner.base import RunnerError

    calls = tmp_path / "calls"
    calls.mkdir()
    log = shlex.quote(str(calls))
    docker = fake_docker(tmp_path, f"""case "$1" in
  run) printf '%s\\n' "$@" > {log}/run; exec sleep 30 ;;
  kill) printf '%s\\n' "$@" > {log}/kill; exit {kill_status} ;;
  rm) printf '%s\\n' "$@" > {log}/rm ;;
esac""")
    monkeypatch.setattr(docker_runner, "require_docker", lambda: docker)

    model = registry.find_model("reference_bucket")
    # The fake must start and record its arguments inside the budget, so the
    # budget leaves a loaded CI runner room for that and nothing more.
    probe = replace(registry.find_probe("mass/catchment-closure"), max_runtime_s=2.0)
    with pytest.raises(RunnerError, match=f"{model.name}: container exceeded the time budget for {probe.id}"):
        docker_runner.DockerRunner().invoke(model, probe, tmp_path, tmp_path / "request.json")

    run = (calls / "run").read_text().splitlines()
    name = run[run.index("--name") + 1]
    assert name.startswith(f"hydroturing-{model.name}-")
    assert (calls / "kill").read_text().splitlines() == ["kill", name]
    if kill_status == 0:
        assert not (calls / "rm").exists()
    else:
        assert (calls / "rm").read_text().splitlines() == ["rm", "-f", name]


def test_submitted_model_cannot_request_host_subprocess_access(tmp_path):
    import shutil

    from hydroturing.spec import load_model

    target = tmp_path / "submitted_model"
    shutil.copytree(registry.MODELS_DIR / "_template", target)
    manifest = target / "model.yaml"
    manifest.write_text(
        manifest.read_text()
        .replace("name: _template", "name: submitted_model")
        .replace("runner: docker", "runner: subprocess")
    )

    with pytest.raises(SpecError, match="reserved for trusted reference"):
        load_model(target)


# --- scaffolding ------------------------------------------------------------
# Contributing must not begin with a blank page, and the directory a
# contributor gets must already validate. If the scaffold produces something
# that fails `ht validate`, the first thing a newcomer sees is an error they
# did not cause.


def test_draft_template_is_not_scaffoldable_unfilled(tmp_path):
    from hydroturing.scaffold import scaffold_probe, write_draft

    draft = write_draft(tmp_path / "draft.yaml")
    with pytest.raises(SpecError, match="still the template's example"):
        scaffold_probe(draft)


def test_scaffolded_probe_validates(tmp_path, monkeypatch):
    import yaml

    from hydroturing import scaffold
    from hydroturing.spec import load_probe

    probes_root = tmp_path / "probes"
    monkeypatch.setattr(scaffold, "PROBES_DIR", probes_root)

    draft = scaffold.write_draft(tmp_path / "draft.yaml")
    text = draft.read_text().replace("id: mass/my-probe", "id: mass/worked-example")
    text = text.replace("name: Your Name", "name: A Contributor")
    draft.write_text(text)

    target = scaffold.scaffold_probe(draft)
    assert target == probes_root / "mass" / "worked-example"
    assert {p.name for p in target.iterdir()} == {"probe.yaml", "generate.py", "README.md"}

    spec = load_probe(target)
    assert spec.id == "mass/worked-example"
    assert spec.must_pass and spec.must_fail

    # Guidance is for the contributor, not for the finished probe.
    body = (target / "probe.yaml").read_text()
    assert "#!" not in body
    assert "\n\n\n" not in body, "stripping guidance left a gap"

    # The generator skeleton must import and honour the length contract, so a
    # contributor's first `ht gate` fails on their physics rather than on the
    # scaffold's own arithmetic.
    from hydroturing.harness import load_generator

    frame, _static = load_generator(spec).generate(1)
    assert len(frame) == int(spec.period_years * 365) + spec.spinup_days
    assert "time" in frame.columns


def test_id_and_law_must_agree(tmp_path, monkeypatch):
    from hydroturing import scaffold

    monkeypatch.setattr(scaffold, "PROBES_DIR", tmp_path / "probes")
    draft = scaffold.write_draft(tmp_path / "draft.yaml")
    draft.write_text(
        draft.read_text()
        .replace("id: mass/my-probe", "id: energy/mismatch")
        .replace("name: Your Name", "name: A Contributor")
    )
    with pytest.raises(SpecError, match="they must agree"):
        scaffold.scaffold_probe(draft)


def test_every_probe_credits_its_authors():
    """Credit lands at the unit of contribution, which is the probe."""
    for spec in registry.all_probes():
        assert spec.authors, f"{spec.id} lists no authors"
        assert all(a.get("name") for a in spec.authors)


# --- templates --------------------------------------------------------------
# Each template is a promise that `ht init-probe --template <kind>` produces
# something a contributor can run. Nothing in CI exercises the templates
# otherwise, so a template that no longer scaffolds would be discovered by the
# first person to try it, which is the worst place to discover it.


def _scaffold_template(scaffold, tmp_path, kind, slug):
    draft = scaffold.write_draft(tmp_path / f"{kind}.yaml", kind=kind)
    text = draft.read_text()
    example = next(
        line.split(": ", 1)[1].strip()
        for line in text.splitlines()
        if line.startswith("id: ")
    )
    draft.write_text(
        text.replace(f"id: {example}", f"id: mass/{slug}")
        .replace("name: Your Name", "name: A Contributor")
    )
    return scaffold.scaffold_probe(draft), example


def _template_kinds():
    from hydroturing.scaffold import available_templates

    return sorted(available_templates())


def test_external_probe_root_is_discovered(tmp_path, monkeypatch):
    from hydroturing import scaffold

    private_root = tmp_path / "private-probes"
    monkeypatch.setattr(scaffold, "PROBES_DIR", private_root)
    target, _ = _scaffold_template(
        scaffold, tmp_path, "default", "private-evaluation-case"
    )

    probes = registry.all_probes([registry.PROBES_DIR, private_root])
    assert target in [probe.path for probe in probes]
    assert registry.find_probe(
        "mass/private-evaluation-case", [registry.PROBES_DIR, private_root]
    ).path == target


def test_unknown_criterion_is_rejected(tmp_path, monkeypatch):
    import yaml

    from hydroturing import scaffold
    from hydroturing.spec import load_probe

    monkeypatch.setattr(scaffold, "PROBES_DIR", tmp_path / "probes")
    target, _ = _scaffold_template(scaffold, tmp_path, "default", "unknown-criterion")
    spec_file = target / "probe.yaml"
    raw = yaml.safe_load(spec_file.read_text())
    original = next(iter(raw["criteria"][0]))
    raw["criteria"][0] = {"not_registered": {}}
    for model, criterion in raw["baselines"]["must_fail"].items():
        if criterion == original:
            raw["baselines"]["must_fail"][model] = "not_registered"
    spec_file.write_text(yaml.safe_dump(raw, sort_keys=False))

    with pytest.raises(SpecError, match="unknown criteria"):
        load_probe(target)


@pytest.mark.parametrize("kind", _template_kinds())
def test_every_template_scaffolds_into_a_valid_probe(kind, tmp_path, monkeypatch):
    from hydroturing import scaffold
    from hydroturing.harness import build_case
    from hydroturing.spec import load_probe

    monkeypatch.setattr(scaffold, "PROBES_DIR", tmp_path / "probes")
    target, _ = _scaffold_template(scaffold, tmp_path, kind, f"template-{kind}")

    spec = load_probe(target)
    assert "#!" not in (target / "probe.yaml").read_text()
    assert spec.must_pass and spec.must_fail

    # Every model the template names in its baselines has to exist, or the
    # contributor's first `ht gate` fails on a missing reference rather than on
    # their own physics.
    for name in (*spec.must_pass, *spec.must_fail):
        registry.find_model(name)

    # The generator that comes with the template has to honour the length
    # contract for every variant the probe declares.
    for variant in spec.variants or (None,):
        case = build_case(spec, 1, variant)
        assert len(case.forcing) == int(spec.period_years * 365) + spec.spinup_days


@pytest.mark.parametrize("kind", [k for k in _template_kinds() if k != "default"])
def test_shaped_templates_get_a_matching_generator(kind, tmp_path, monkeypatch):
    """A template whose criteria need labelled regimes or paired variants must
    scaffold a generator that produces them, not the plain skeleton."""
    from hydroturing import scaffold
    from hydroturing.harness import build_case
    from hydroturing.spec import load_probe

    monkeypatch.setattr(scaffold, "PROBES_DIR", tmp_path / "probes")
    target, _ = _scaffold_template(scaffold, tmp_path, kind, f"shaped-{kind}")
    spec = load_probe(target)

    names = {c.name for c in spec.criteria}
    if "regime_transfer" in names:
        case = build_case(spec, 1)
        labels = set(case.forcing["_regime"].unique())
        params = next(c.params for c in spec.criteria if c.name == "regime_transfer")
        assert {params["reference"], params["extrapolation"]} <= labels

    if spec.variants:
        cases = {v: build_case(spec, 1, v) for v in spec.variants}
        first, second = (cases[v].forcing for v in spec.variants[:2])
        assert not first.equals(second), (
            "the variants are identical, so the paired criterion compares "
            "a case with itself and passes for any model at all"
        )


def test_annotations_are_not_staged_for_the_model(tmp_path):
    """A `_regime` column tells the criteria which steps are the extrapolation.
    Handing it to the model would tell it which part it is being judged on."""
    import pandas as pd

    from hydroturing.protocol import Case, stage

    probe = registry.find_probe("mass/catchment-closure")
    model = registry.find_model("reference_bucket")
    case = build_case(probe, 5)
    case = Case(
        probe_id=case.probe_id,
        seed=case.seed,
        forcing=case.forcing.assign(_regime="anomaly"),
        static=case.static,
        spinup_steps=case.spinup_steps,
    )

    stage(tmp_path, case, probe, model)
    staged = pd.read_csv(tmp_path / "input" / "forcing.csv")
    assert "_regime" not in staged.columns
    assert "pr" in staged.columns


def test_paired_criteria_and_variants_must_agree(tmp_path, monkeypatch):
    """Declaring one without the other is a silent no-op at run time."""
    import yaml

    from hydroturing import scaffold
    from hydroturing.spec import load_probe

    monkeypatch.setattr(scaffold, "PROBES_DIR", tmp_path / "probes")
    target, _ = _scaffold_template(scaffold, tmp_path, "counterfactual", "pairing")

    spec_file = target / "probe.yaml"
    raw = yaml.safe_load(spec_file.read_text())
    del raw["case"]["variants"]
    spec_file.write_text(yaml.safe_dump(raw, sort_keys=False))
    with pytest.raises(SpecError, match="case.variants"):
        load_probe(target)


def test_regime_transfer_separates_the_stretches():
    """The point of the criterion: a leak confined to one stretch of the record
    is invisible to a whole-window budget and must not be invisible here."""
    import numpy as np
    import pandas as pd

    from hydroturing.criteria.base import Window, segments

    n = 400
    forcing = pd.DataFrame({
        "pr": np.full(n, 10.0),
        "_regime": ["ordinary"] * 300 + ["anomaly"] * 100,
    })
    blocks = segments(Window(forcing, forcing, forcing.iloc[0], 1.0))
    assert blocks == [("ordinary", 0, 300), ("anomaly", 300, 400)]


@pytest.mark.parametrize("kind", [k for k in _template_kinds() if k != "default"])
def test_every_template_discriminates_out_of_the_box(kind, tmp_path, monkeypatch):
    """The templates ship filled-in baselines, so they make a claim: scaffold
    this and the acceptance gate already separates the reference models. One
    seed rather than the full set, because this is checking that the template
    is coherent, not running the gate."""
    from hydroturing import scaffold
    from hydroturing.spec import load_probe

    monkeypatch.setattr(scaffold, "PROBES_DIR", tmp_path / "probes")
    target, _ = _scaffold_template(scaffold, tmp_path, kind, f"gate-{kind}")
    spec = load_probe(target)
    seeds = gate_seeds(spec.id, 1)

    for name in spec.must_pass:
        outcome = run_probe(registry.find_model(name), spec, seeds)
        assert outcome.verdict == PASS, (
            f"{kind}: {name} is an exact model and the template fails it on "
            f"{outcome.failing}; the template is wrong, not the model"
        )

    for name, expected in spec.must_fail.items():
        outcome = run_probe(registry.find_model(name), spec, seeds)
        assert expected in outcome.failing, (
            f"{kind}: {name} was expected to trip '{expected}', "
            f"tripped {outcome.failing or 'nothing'}"
        )


# --- report marks -----------------------------------------------------------
# A scored verdict is a bit, so the report should read as one at a glance,
# and N/A as neither value. The marks are decoration over the words, never a
# replacement for them: a log someone greps for FAIL has to keep finding it.


def _report(model_name):
    return run_model(registry.find_model(model_name), [registry.find_probe("mass/catchment-closure")], [11])


def test_text_report_marks_pass_and_fail(monkeypatch):
    from hydroturing import report

    monkeypatch.delenv("HT_ASCII", raising=False)
    monkeypatch.setattr(report, "use_emoji", lambda: True)

    good = report.to_text(_report("reference_bucket"))
    bad = report.to_text(_report("reference_leaky"))

    assert good.count(report.PASS_MARK) >= 6  # model, probe, five criteria
    assert report.FAIL_MARK not in good
    assert report.FAIL_MARK in bad and report.PASS_MARK in bad

    # The words survive alongside the marks, in both directions.
    assert "PASS" in good and "FAIL" in bad
    assert "VIOLATION" in bad


def test_ht_ascii_drops_the_marks_without_doubling_the_word(monkeypatch):
    """`FAIL FAIL` helps nobody. Where the marker would only repeat the word
    already on the line, it is dropped rather than printed twice."""
    from hydroturing import report

    monkeypatch.setenv("HT_ASCII", "1")
    text = report.to_text(_report("reference_leaky"))

    assert report.FAIL_MARK not in text and report.PASS_MARK not in text
    assert "FAIL FAIL" not in text
    assert text.startswith("reference_leaky v1.0.0  ->  FAIL (VIOLATION)")
    assert "  FAIL  mass/catchment-closure" in text


def test_an_unscored_probe_reads_as_neither_pass_nor_fail(monkeypatch):
    """N/A is not a FAIL, so neither its mark nor its word may read as one.
    The reason says why the probe was not scored, and the summary does not
    count the model out of a probe that asked it nothing."""
    from hydroturing import report

    monkeypatch.setattr(report, "use_emoji", lambda: True)
    unscored = _report("reference_streamflow_only")
    text = report.to_text(unscored)
    assert report.FAIL_MARK not in text and report.PASS_MARK not in text
    assert text.count(report.NOT_SCORED_MARK) == 2  # model and probe
    assert "N/A (INCOMPLETE)  [no probe could be scored]" in text
    assert f"{report.NOT_SCORED_MARK} **N/A** (INCOMPLETE)" in report.to_markdown(unscored)

    monkeypatch.setattr(report, "use_emoji", lambda: False)
    text = report.to_text(unscored)
    assert "FAIL" not in text
    assert "  N/A   mass/catchment-closure" in text

    rows = report.to_csv_rows(unscored, "2026-09-11")
    assert (rows[0]["verdict"], rows[0]["reason"]) == ("N/A", "INCOMPLETE")


def test_marks_are_dropped_when_the_stream_cannot_carry_them(monkeypatch):
    """An ASCII stdout must produce a plain report, not a UnicodeEncodeError
    halfway through one."""
    import io

    from hydroturing import report

    monkeypatch.delenv("HT_ASCII", raising=False)
    monkeypatch.setattr(report.sys, "stdout", io.TextIOWrapper(io.BytesIO(), encoding="ascii"))
    assert report.use_emoji() is False

    text = report.to_text(_report("reference_bucket"))
    text.encode("ascii")  # raises if a mark slipped through


def test_marked_columns_stay_aligned(monkeypatch):
    """PASS and FAIL are the same length and both marks are the same width, so
    a mixed report must not leave the criterion names in a ragged column."""
    from hydroturing import report

    monkeypatch.setattr(report, "use_emoji", lambda: True)
    text = report.to_text(_report("reference_leaky"))

    names = {c.name for c in _report("reference_leaky").probes[0].criteria}
    criterion_lines = [l for l in text.splitlines() if l.startswith(" " * 8)]
    assert len(criterion_lines) == len(names)
    columns = {
        line.index(name)
        for line in criterion_lines
        for name in names
        if name in line
    }
    assert len(columns) == 1, f"criterion names start at differing columns: {columns}"


def test_trusted_models_run_on_the_harness_interpreter():
    """Manifests say `python3` because a container has one; the host may
    not (Windows has `python` or `py`), and a virtual environment's
    interpreter is the one with the dependencies. A trusted in-repo model
    runs on whatever is running the harness."""
    import sys
    from hydroturing.runner.subprocess_runner import SubprocessRunner

    assert SubprocessRunner.resolve_entrypoint(["python3", "ht_adapter.py"]) == [sys.executable, "ht_adapter.py"]
    assert SubprocessRunner.resolve_entrypoint(["python", "x.py"])[0] == sys.executable
    assert SubprocessRunner.resolve_entrypoint(["./model.exe", "--go"]) == ["./model.exe", "--go"]

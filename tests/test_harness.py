"""Tests for the harness itself, distinct from the probe acceptance gate.

The gate asks whether the benchmark discriminates. These ask whether the
machinery underneath it is sound: generators are deterministic, the window
arithmetic is not off by one, and the contract is enforced.
"""

from __future__ import annotations

import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.harness import build_case, run_probe
from hydroturing.protocol import ProtocolError, read_result
from hydroturing.scoring import FAIL, INCOMPLETE, PASS, VIOLATION
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
    say enough to be checked, and the report must not conflate the two."""
    outcome = run_probe(registry.find_model("reference_streamflow_only"), probe, [11])
    assert outcome.verdict == FAIL
    assert outcome.reason == INCOMPLETE
    assert "evspsbl" in outcome.missing


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


# --- container isolation ----------------------------------------------------
# The daemon is not available in every dev environment, so the isolation flags
# are asserted structurally here and the real build-and-run happens in CI,
# where Docker exists. Isolation is part of the benchmark: a model that could
# read the probe definition could read the tolerance it is judged against.


def test_container_runs_with_no_network_and_one_mount(tmp_path):
    from hydroturing.runner.docker_runner import DockerRunner

    model = registry.find_model("reference_bucket")
    argv = DockerRunner.command("docker", "img:1", model, tmp_path)

    assert argv[:3] == ["docker", "run", "--rm"]
    assert "--network" in argv and argv[argv.index("--network") + 1] == "none"

    mounts = [argv[i + 1] for i, a in enumerate(argv) if a == "--mount"]
    assert len(mounts) == 1, "the container must see the case and nothing else"
    assert mounts[0] == f"type=bind,source={tmp_path.resolve()},target=/io"

    assert "--memory" in argv and "--cpus" in argv
    assert argv[-2:] == ["--request", "/io/request.json"]
    assert not any(str(registry.PROBES_DIR) in a for a in argv), (
        "the probes directory must never be mounted into a model container"
    )


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


# --- scaffolding ------------------------------------------------------------
# Contributing must not begin with a blank page, and the directory a
# contributor gets must already validate. If the scaffold produces something
# that fails `ht validate`, the first thing a newcomer sees is an error they
# did not cause.


def test_draft_template_is_not_scaffoldable_unfilled(tmp_path):
    from hydroturing.scaffold import scaffold_probe, write_draft

    draft = write_draft(tmp_path / "draft.yaml")
    with pytest.raises(SpecError, match="still 'mass/my-probe'"):
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

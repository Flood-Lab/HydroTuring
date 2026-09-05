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
    build_case,
    run_model,
    run_probe,
    verify_adapter_contract,
)
from hydroturing.protocol import ProtocolError, read_result, stage
from hydroturing.scoring import FAIL, INCOMPATIBLE, INCOMPLETE, PASS, VIOLATION
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


def test_incompatible_timestep_is_reported_before_execution(probe):
    model = replace(registry.find_model("reference_bucket"), timesteps=("PT1H",))
    outcome = run_probe(model, probe, [11])
    assert outcome.verdict == FAIL
    assert outcome.reason == INCOMPATIBLE
    assert "timestep" in outcome.incompatible[0]


def test_missing_forcing_is_reported_as_incompatible(probe):
    model = replace(
        registry.find_model("reference_bucket"),
        needs_forcing=("pr", "unavailable_driver"),
    )
    outcome = run_probe(model, probe, [11])
    assert outcome.reason == INCOMPATIBLE
    assert "unavailable_driver" in outcome.incompatible[0]


def test_paired_probe_requires_declared_perturbation_support(probe):
    paired = replace(probe, variants=("control", "perturbed"))
    model = replace(registry.find_model("reference_bucket"), supports_perturbation=False)
    outcome = run_probe(model, paired, [11])
    assert outcome.reason == INCOMPATIBLE
    assert "perturbation" in outcome.incompatible[0]


def test_adapter_verification_runs_an_incomplete_model(probe):
    """A scientific INCOMPLETE verdict must not skip the contract smoke test."""
    model = registry.find_model("reference_streamflow_only")
    result = verify_adapter_contract(model, probe, gate_seeds(probe.id, 1)[0])
    assert len(result.table) == result.case.n_steps
    assert {"time", "mrro", "dis"} <= set(result.table.columns)


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
    argv = DockerRunner.command("docker", "img:1", model, tmp_path)

    assert argv[:3] == ["docker", "run", "--rm"]
    assert "--network" in argv and argv[argv.index("--network") + 1] == "none"
    assert "--read-only" in argv
    assert argv[argv.index("--cap-drop") + 1] == "ALL"
    # Not root: with every capability dropped, root could not even write the
    # output directory the host user owns, and it should not be able to write
    # anything else of the host's either.
    import os
    assert argv[argv.index("--user") + 1] == f"{os.getuid()}:{os.getgid()}"
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
    argv = DockerRunner.command("docker", "img:1", model, tmp_path)
    assert argv[argv.index("--cpus") + 1] == "4"

    monkeypatch.setattr(docker_runner.os, "cpu_count", lambda: 16)
    argv = DockerRunner.command("docker", "img:1", model, tmp_path)
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
# The verdict is a bit, so the report should read as one at a glance. The
# marks are decoration over the words, never a replacement for them: a log
# someone greps for FAIL has to keep finding it.


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

"""Tests for the evaluation window.

A submitted model is scored on a flood event rather than the full record, so
that a heavy model fits its time budget. These check that the window lands
on the event the probe's own physics puts there, that it carries its spinup,
that the default is what the manifest and the command line say it is, and
that the probe still discriminates on it.
"""

from __future__ import annotations

import csv
from dataclasses import replace

import numpy as np
import pytest

from hydroturing import registry
from hydroturing.harness import (
    WindowBounds,
    WindowError,
    build_case,
    resolve_window_days,
    run_model,
    run_probe,
    select_window,
    verify_adapter_contract,
    window_case,
)
from hydroturing.protocol import Case
from hydroturing.report import append_csv, to_dict, to_text
from hydroturing.scoring import FAIL, PASS
from hydroturing.seeds import gate_seeds
from hydroturing.spec import DEFAULT_WINDOW_DAYS, FULL_WINDOW, SpecError, load_model


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/catchment-closure")


def _submitted(model_name="reference_bucket", **overrides):
    """A reference model dressed up as a submission, so the defaults for a
    submitted model apply without needing Docker."""
    model = registry.find_model(model_name)
    return replace(model, name="submitted_model", **overrides)


def test_reference_models_default_to_the_full_record(probe):
    for model in registry.all_models():
        if model.runner == "subprocess":
            assert resolve_window_days(model, probe) is None, model.name


def test_submitted_models_default_to_a_flood_event(probe):
    assert resolve_window_days(_submitted(), probe) == DEFAULT_WINDOW_DAYS["PT1D"]
    hourly = replace(probe, timestep="PT1H")
    assert resolve_window_days(_submitted(timesteps=("PT1H",)), hourly) == DEFAULT_WINDOW_DAYS["PT1H"]


def test_manifest_and_command_line_decide_the_window(probe):
    assert resolve_window_days(_submitted(window_days=45), probe) == 45
    assert resolve_window_days(_submitted(window_days=FULL_WINDOW), probe) is None
    # The command line wins over the manifest, in both directions.
    assert resolve_window_days(_submitted(window_days=45), probe, override=7) == 7
    assert resolve_window_days(_submitted(window_days=45), probe, override=FULL_WINDOW) is None
    assert resolve_window_days(registry.find_model("reference_bucket"), probe, override=30) == 30


def test_window_lands_on_the_flood_the_reference_model_produces(probe):
    """The event is where the exact physics puts the runoff peak, not where
    the most rain falls: in a catchment with a snowpack those are months
    apart."""
    import tempfile
    from pathlib import Path

    from hydroturing.runner import get_runner

    case = build_case(probe, gate_seeds(probe.id, 1)[0])
    bounds = select_window(case, probe, 30)
    offset, rows = bounds.rows(case.dt_days)
    start, stop = case.spinup_steps + offset, case.spinup_steps + offset + rows
    assert rows == 30

    reference = registry.find_model("reference_bucket")
    run = get_runner(reference).run(reference, probe, case, Path(tempfile.mkdtemp()))
    runoff = run.table["mrro"].to_numpy()
    peak = case.spinup_steps + int(np.argmax(runoff[case.spinup_steps :]))
    assert start <= peak < stop


def test_window_keeps_its_spinup_and_records_what_was_scored(probe):
    case = build_case(probe, 5)
    bounds = select_window(case, probe, 30)
    offset, rows = bounds.rows(case.dt_days)
    start, stop = case.spinup_steps + offset, case.spinup_steps + offset + rows
    cut = window_case(case, bounds)

    assert cut.n_steps == case.spinup_steps + 30
    assert cut.spinup_steps == case.spinup_steps
    assert cut.timestep == case.timestep
    assert cut.forcing["time"].iloc[0] == case.forcing["time"].iloc[start - case.spinup_steps]
    assert cut.window == {
        "days": 30,
        "rows": 30,
        "start": case.forcing["time"].iloc[start],
        "end": case.forcing["time"].iloc[stop - 1],
    }
    # The scored stretch is exactly the rows the selection named.
    scored = cut.forcing.iloc[cut.spinup_steps :].reset_index(drop=True)
    expected = case.forcing.iloc[start:stop].reset_index(drop=True)
    assert scored.equals(expected)


def test_window_longer_than_the_record_is_the_record(probe):
    case = build_case(probe, 5)
    cut = window_case(case, select_window(case, probe, 10_000))
    assert cut.n_steps == case.n_steps
    assert cut.window["rows"] == case.n_steps - case.spinup_steps


def test_window_that_drops_a_labelled_stretch_is_incompatible(probe):
    """A regime probe labels the stretch it scores. Cutting it away would
    evaluate the wrong thing, so the model is INCOMPATIBLE with the window
    rather than quietly scored on ordinary weather."""
    case = build_case(probe, 5)
    labels = np.array(["ordinary"] * case.n_steps, dtype=object)
    labels[-100:] = "anomaly"
    labelled = Case(
        probe_id=case.probe_id,
        seed=case.seed,
        forcing=case.forcing.assign(_regime=labels),
        static=case.static,
        spinup_steps=case.spinup_steps,
        timestep=case.timestep,
    )
    with pytest.raises(WindowError, match="anomaly"):
        window_case(labelled, WindowBounds(offset_days=10.0, days=30))

    # A window that keeps every labelled stretch can be scored.
    scored_rows = labelled.n_steps - labelled.spinup_steps
    cut = window_case(labelled, WindowBounds(offset_days=float(scored_rows - 200), days=200))
    assert set(cut.forcing["_regime"].iloc[cut.spinup_steps :]) == {"ordinary", "anomaly"}


def test_exact_model_passes_on_a_flood_window(probe):
    """The window must not fail the physics it is supposed to test. The
    runoff ratio is climatological and is reported rather than judged on a
    month; everything else must hold on the event."""
    seeds = gate_seeds(probe.id, 3)
    outcome = run_probe(_submitted(), probe, seeds)
    assert outcome.window_days == DEFAULT_WINDOW_DAYS["PT1D"]
    assert [w["seed"] for w in outcome.windows] == seeds
    assert all(w["rows"] == 30 for w in outcome.windows)
    assert outcome.verdict == PASS, outcome.failing

    degenerate = next(c for c in outcome.criteria if c.name == "non_degenerate")
    assert "not applied" in degenerate.diagnostics["runoff_ratio_check"]


def test_broken_models_are_still_caught_on_a_flood_window(probe):
    seeds = gate_seeds(probe.id, 2)
    for name, criterion in probe.must_fail.items():
        outcome = run_probe(_submitted(name), probe, seeds)
        assert outcome.verdict == FAIL
        assert criterion in outcome.failing, (name, outcome.failing)


def test_runoff_ratio_is_still_judged_on_the_full_record(probe):
    """Skipping the ratio on a short window must not disable it on the
    record the gate is defined on."""
    outcome = run_probe(registry.find_model("reference_degenerate"), probe, [11])
    degenerate = next(c for c in outcome.criteria if c.name == "non_degenerate")
    assert "runoff ratio" in degenerate.message
    assert "runoff_ratio_check" not in degenerate.diagnostics


def test_verify_adapter_uses_the_window(probe):
    model = _submitted("reference_streamflow_only")
    result = verify_adapter_contract(model, probe, gate_seeds(probe.id, 1)[0])
    assert result.case.window is not None
    assert result.case.n_steps == probe.spinup_days + DEFAULT_WINDOW_DAYS["PT1D"]
    assert len(result.table) == result.case.n_steps

    full = verify_adapter_contract(model, probe, gate_seeds(probe.id, 1)[0], window=FULL_WINDOW)
    assert full.case.window is None
    assert full.case.n_steps == int(probe.period_years * 365) + probe.spinup_days


def test_manifest_window_is_validated(tmp_path):
    import shutil

    target = tmp_path / "submitted_model"
    shutil.copytree(registry.MODELS_DIR / "_template", target)
    manifest = target / "model.yaml"
    text = manifest.read_text().replace("name: _template", "name: submitted_model")

    manifest.write_text(text + "\nwindow_days: 45\n")
    assert load_model(target).window_days == 45

    manifest.write_text(text + "\nwindow_days: full\n")
    assert load_model(target).window_days == FULL_WINDOW

    manifest.write_text(text + "\nwindow_days: 0\n")
    with pytest.raises(SpecError, match="window_days"):
        load_model(target)


def test_report_says_which_stretch_was_scored(probe):
    seeds = gate_seeds(probe.id, 1)
    report = run_model(_submitted(), [probe], seeds)
    payload = to_dict(report)["probes"][0]
    assert payload["window"]["days"] == DEFAULT_WINDOW_DAYS["PT1D"]
    assert payload["window"]["cases"][0]["seed"] == seeds[0]
    assert "window: 30-day flood event" in to_text(report)

    full = run_model(registry.find_model("reference_bucket"), [probe], seeds)
    assert to_dict(full)["probes"][0]["window"] is None
    assert "window:" not in to_text(full)


def test_csv_archive_appends_one_row_per_probe(probe, tmp_path):
    archive = tmp_path / "result.csv"
    seeds = gate_seeds(probe.id, 1)
    append_csv(run_model(registry.find_model("reference_bucket"), [probe], seeds), archive, "2026-09-03")
    append_csv(run_model(registry.find_model("reference_leaky"), [probe], seeds), archive, "2026-09-03")

    with open(archive, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert [r["model"] for r in rows] == ["reference_bucket", "reference_leaky"]
    assert rows[0]["verdict"] == "PASS" and rows[0]["window"] == "full record"
    assert rows[1]["reason"] == "VIOLATION" and rows[1]["detail"].startswith("closure:")
    assert rows[1]["seeds"] == str(seeds[0])
    # One header, however many appends.
    assert archive.read_text().count("run_date,") == 1


def test_cli_window_and_csv_flags(probe, tmp_path, capsys):
    from hydroturing.cli import main

    archive = tmp_path / "result.csv"
    code = main([
        "run", "--model", "reference_bucket", "--seed", "11",
        "--window", "30", "--csv", str(archive),
    ])
    assert code == 0
    assert "30-day flood event" in capsys.readouterr().out
    assert "30-day flood event" in archive.read_text()

    with pytest.raises(SystemExit):
        main(["run", "--model", "reference_bucket", "--window", "soon"])


def test_models_that_never_ran_are_labelled_so(probe):
    """An INCOMPLETE verdict is issued before the container starts. The
    archive must not say the model was scored on a record it never saw."""
    from hydroturing.report import to_csv_rows, window_label

    outcome = run_probe(_submitted("reference_streamflow_only"), probe, [11])
    assert outcome.reason == "INCOMPLETE"
    assert window_label(outcome) == "not run"

    report = run_model(_submitted("reference_streamflow_only"), [probe], [11])
    assert to_csv_rows(report, "2026-09-03")[0]["window"] == "not run"


def test_contract_check_can_be_archived(probe, tmp_path):
    """For a streamflow-only model the contract row is the only evidence it
    was built and run."""
    from hydroturing.cli import main

    archive = tmp_path / "result.csv"
    assert main(["verify-adapter", "--model", "reference_streamflow_only", "--csv", str(archive)]) == 0
    with open(archive, newline="") as fh:
        (row,) = list(csv.DictReader(fh))
    assert row["probe"] == "mass/catchment-closure (adapter contract)"
    assert row["verdict"] == "PASS" and row["reason"] == "OK"
    assert row["window"] == "full record"
    assert row["detail"].startswith("adapter contract OK: 4015 rows")
    assert row["seeds"] == str(gate_seeds(probe.id, 1)[0])

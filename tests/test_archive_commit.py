"""An archive row names the commit it ran at, and a stale archive refuses new rows."""

from __future__ import annotations

import csv
import subprocess

import pytest

from hydroturing import report
from hydroturing.report import CSV_COLUMNS, append_csv_rows, contract_row, source_commit


def row(**fields):
    return {name: fields.get(name, "") for name in CSV_COLUMNS}


def test_a_row_records_the_commit(monkeypatch):
    monkeypatch.setattr(report, "source_commit", lambda: "abc")
    assert contract_row("m", "1", "0", "docker", "p", 1, error="boom")["commit"] == "abc"


def test_the_commit_is_marked_when_the_tree_differs(monkeypatch):
    answers = {"rev-parse": "abc", "status": " M src/hydroturing/report.py"}

    def fake(command, **_):
        return subprocess.CompletedProcess(command, 0, answers[command[3]] + "\n", "")

    monkeypatch.setattr(report.subprocess, "run", fake)
    assert source_commit() == "abc+dirty"
    answers["status"] = ""
    assert source_commit() == "abc"


def test_outside_a_checkout_the_commit_is_empty(monkeypatch):
    def missing(*_, **__):
        raise FileNotFoundError("git")

    monkeypatch.setattr(report.subprocess, "run", missing)
    assert source_commit() == ""


def test_rows_append_under_the_header(tmp_path):
    path = tmp_path / "result.csv"
    append_csv_rows([row(model="a", commit="abc")], path)
    append_csv_rows([row(model="b", commit="def")], path)
    with open(path, newline="") as fh:
        assert [r["commit"] for r in csv.DictReader(fh)] == ["abc", "def"]


def test_an_archive_with_other_columns_refuses_rows(tmp_path):
    path = tmp_path / "result.csv"
    path.write_text("run_date,model\n2026-01-01,a\n")
    with pytest.raises(ValueError, match="expected"):
        append_csv_rows([row(model="b")], path)
    assert path.read_text() == "run_date,model\n2026-01-01,a\n"

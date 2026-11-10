"""The rate rescaling behind HYPE's hourly result, which runs without HYPE."""

import importlib.util
import os

import pytest

from hydroturing.spec import REPO_ROOT


@pytest.fixture
def hype(monkeypatch, tmp_path):
    (tmp_path / "par.txt").write_text(
        "!! a comment naming cevp 9\n"
        "cevp\t0.2\t0.4\n"
        "rrcs1\t0.3\n"
        "wcfc\t0.1\t0.2\n",
        encoding="latin-1",
    )
    monkeypatch.setitem(os.environ, "HYPE_TEMPLATE", str(tmp_path))
    spec = importlib.util.spec_from_file_location("hype_adapter", REPO_ROOT / "models" / "hype" / "ht_adapter.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _values(text, name):
    for line in text.splitlines():
        parts = line.split()
        if parts and parts[0] == name:
            return [float(v) for v in parts[1:]]
    raise AssertionError(name)


def test_a_daily_step_leaves_par_txt_alone(hype):
    text, changed = hype.par_txt(1.0)
    assert changed == []
    assert _values(text, "cevp") == [0.2, 0.4]
    assert _values(text, "rrcs1") == [0.3]


def test_an_hourly_step_scales_amounts_and_fractions_differently(hype):
    text, changed = hype.par_txt(1.0 / 24.0)
    assert _values(text, "cevp") == pytest.approx([0.2 / 24.0, 0.4 / 24.0])
    assert _values(text, "rrcs1") == pytest.approx([1.0 - 0.7 ** (1.0 / 24.0)])
    assert {c["parameter"] for c in changed} == {"cevp", "rrcs1"}


def test_an_hourly_step_keeps_comments_and_other_parameters(hype):
    text, _ = hype.par_txt(1.0 / 24.0)
    assert text.startswith("!! a comment naming cevp 9\n")
    assert _values(text, "wcfc") == [0.1, 0.2]

"""How LISFLOOD-FP's mass-file samples become rows."""

import importlib.util

import pytest

from hydroturing.spec import REPO_ROOT


@pytest.fixture
def adapter():
    spec = importlib.util.spec_from_file_location(
        "lisflood_fp_adapter", REPO_ROOT / "models" / "lisflood_fp" / "ht_adapter.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _mass_file(path, samples):
    # Time is column 0, Vol column 5, Qout column 8, as in the file's own header.
    lines = ["Time Tstep MinTstep NumTsteps Area Vol Qin Hds Qout"]
    lines += [f"{t} 0 0 0 0 {v} 0 0 {q}" for t, v, q in samples]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_a_sample_on_a_rows_end_belongs_to_that_row(adapter, tmp_path):
    samples = [(t, t * 10.0, 1.0) for t in range(0, 7201, 60)]
    _, volumes = adapter.read_mass(_mass_file(tmp_path / "m", samples), 2, 3600.0)
    assert volumes == [36000.0, 72000.0]


def test_outflow_is_the_interval_mean(adapter, tmp_path):
    samples = [(t, 0.0, t / 3600.0) for t in range(0, 3601, 60)]
    means, _ = adapter.read_mass(_mass_file(tmp_path / "m", samples), 1, 3600.0)
    assert means == pytest.approx([0.5])


def test_a_row_with_no_sample_is_refused(adapter, tmp_path):
    with pytest.raises(RuntimeError, match="no mass samples"):
        adapter.read_mass(_mass_file(tmp_path / "m", [(60.0, 0.0, 0.0)]), 2, 3600.0)

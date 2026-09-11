"""The /io contract between the harness and a model.

The boundary is a process boundary and the payload is files, because a
submitted model may be written in any language. The harness stages a case,
invokes the adapter once, and reads back a table.

    /io/
      request.json          <- harness writes
      input/forcing.csv     <- harness writes, read only
      input/static.json     <- harness writes, read only
      output/result.csv     -> adapter writes
      output/run.json       -> adapter writes

CSV is the required format. NetCDF is accepted when present but never
required, because demanding a netCDF stack inside every submitted container
would exclude exactly the language diversity this contract exists to allow.
"""

from __future__ import annotations

import hashlib
import json
import stat
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from hydroturing.spec import TIMESTEP_DAYS, ModelManifest, ProbeSpec, UNITS

REQUEST_FILE = "request.json"
FORCING_FILE = "input/forcing.csv"
STATIC_FILE = "input/static.json"
RESULT_CSV = "output/result.csv"
RESULT_NC = "output/result.nc"
RUN_FILE = "output/run.json"

TIME_COL = "time"


class ProtocolError(RuntimeError):
    """The adapter did not honour the contract."""


def _opaque_case_metadata(case: Case) -> tuple[str, int]:
    """Identifiers reproducible by the harness but useless for probe detection.

    The generator seed remains in the host-side report. The model gets a stable,
    derived seed for any stochastic inference it performs, without being handed
    the seed that selects the evaluation case.
    """
    # The variant suffix is dropped so that every variant of one seed hands
    # the model the same seed: a paired probe compares runs, and a stochastic
    # model must draw the same numbers in both or the comparison sees its
    # sampling noise instead of the perturbation.
    base_probe = case.probe_id.split("@", 1)[0]
    material = f"{base_probe}\0{case.seed}".encode()
    digest = hashlib.sha256(material).digest()
    case_id = f"case-{digest[:12].hex()}"
    model_seed = int.from_bytes(digest[12:16], "big") & 0x7FFFFFFF
    return case_id, model_seed


@dataclass
class Case:
    """One generated instance of a probe: the forcing, the attributes, the seed."""

    probe_id: str
    seed: int
    forcing: pd.DataFrame
    static: dict[str, Any]
    # Rows of spinup in front of the scored record, at this case's step.
    spinup_steps: int
    # The step this case runs at, as an ISO 8601 duration. A paired probe
    # may run its variants at different steps.
    timestep: str = "PT1D"
    # Set when the record was cut down to an evaluation window: the days
    # asked for and the first and last scored timestamps. The model is never
    # told; it sees a shorter forcing and nothing else.
    window: dict[str, Any] | None = None

    @property
    def n_steps(self) -> int:
        return len(self.forcing)

    @property
    def dt_days(self) -> float:
        return TIMESTEP_DAYS[self.timestep]

    def after_spinup(self, frame: pd.DataFrame) -> pd.DataFrame:
        return frame.iloc[self.spinup_steps :].reset_index(drop=True)


@dataclass
class RunResult:
    """What a model returned for one case, after contract validation."""

    case: Case
    table: pd.DataFrame
    meta: dict[str, Any]
    wall_seconds: float


def stage(io_dir: Path, case: Case, probe: ProbeSpec, model: ModelManifest) -> Path:
    """Write request.json and the read-only inputs. Returns the request path."""
    (io_dir / "input").mkdir(parents=True, exist_ok=True)
    (io_dir / "output").mkdir(parents=True, exist_ok=True)

    # Columns beginning with an underscore are the probe's own annotations —
    # which steps are the extrapolation, which site a block belongs to — and
    # they are stripped here. The criteria need them; the model must not have
    # them, or a regime probe would hand the model a label saying "this is the
    # part you are being tested on".
    visible = [c for c in case.forcing.columns if not c.startswith("_")]
    case.forcing[visible].to_csv(io_dir / FORCING_FILE, index=False)
    with open(io_dir / STATIC_FILE, "w") as fh:
        json.dump(case.static, fh, indent=2)

    case_id, model_seed = _opaque_case_metadata(case)
    request = {
        "case_id": case_id,
        "seed": model_seed,
        "timestep": case.timestep,
        "n_steps": case.n_steps,
        "request": {
            # Asking for every declared output keeps this part of the request
            # invariant across probes and prevents it identifying the criterion.
            "fluxes": list(model.emits_fluxes),
            "states": list(model.emits_states),
            "diagnostics": list(model.emits_diagnostics),
        },
        "input": {"forcing": FORCING_FILE, "static": STATIC_FILE},
        "output": {"table": RESULT_CSV, "run": RUN_FILE},
        "units": {v: UNITS[v] for v in model.emitted if v in UNITS},
        "notes": (
            "States are absolute storages, not tendencies. The harness "
            "differences them itself."
        ),
    }
    request_path = io_dir / REQUEST_FILE
    with open(request_path, "w") as fh:
        json.dump(request, fh, indent=2)
    return request_path


def _validate_output_files(io_dir: Path, probe: ProbeSpec) -> None:
    """Reject unsafe or oversized adapter output before parsing any of it."""
    output_dir = io_dir / "output"
    if not output_dir.exists():
        return

    allowed = {RESULT_CSV, RESULT_NC, RUN_FILE}
    total = 0
    for path in output_dir.iterdir():
        relative = f"output/{path.name}"
        if relative not in allowed:
            raise ProtocolError(f"adapter wrote unexpected output '{relative}'")
        info = path.lstat()
        if not stat.S_ISREG(info.st_mode):
            raise ProtocolError(f"adapter output '{relative}' is not a regular file")
        total += info.st_size

    limit = int(probe.max_output_mb * 1024 * 1024)
    if total > limit:
        raise ProtocolError(
            f"adapter output is {total / (1024 * 1024):.2f} MB, exceeding the "
            f"{probe.max_output_mb:g} MB limit"
        )


def _validate_time_axis(table: pd.DataFrame, case: Case) -> None:
    expected = case.forcing[TIME_COL].reset_index(drop=True)
    actual = table[TIME_COL].reset_index(drop=True)

    # Datetime-capable formats are compared as instants; otherwise the contract
    # falls back to exact textual equality.
    expected_dt = pd.to_datetime(expected, errors="coerce", utc=True)
    actual_dt = pd.to_datetime(actual, errors="coerce", utc=True)
    if not expected_dt.isna().any() and not actual_dt.isna().any():
        equal = expected_dt.equals(actual_dt)
        mismatch = expected_dt != actual_dt
    else:
        expected_text = expected.astype(str)
        actual_text = actual.astype(str)
        equal = expected_text.equals(actual_text)
        mismatch = expected_text != actual_text

    if not equal:
        first = int(mismatch.to_numpy().nonzero()[0][0])
        raise ProtocolError(
            f"result time axis differs from the forcing at row {first}: "
            f"got {actual.iloc[first]!r}, expected {expected.iloc[first]!r}"
        )


def read_result(io_dir: Path, case: Case, probe: ProbeSpec, wall_seconds: float) -> RunResult:
    """Read and validate what the adapter produced."""
    _validate_output_files(io_dir, probe)
    csv_path = io_dir / RESULT_CSV
    nc_path = io_dir / RESULT_NC

    if csv_path.exists():
        table = pd.read_csv(csv_path)
    elif nc_path.exists():
        try:
            import xarray as xr
        except ImportError:  # pragma: no cover - optional path
            raise ProtocolError(
                "adapter wrote result.nc but xarray is not installed; "
                "install hydroturing[netcdf] or write result.csv instead"
            ) from None
        table = xr.open_dataset(nc_path).to_dataframe().reset_index()
    else:
        raise ProtocolError(f"adapter wrote neither {RESULT_CSV} nor {RESULT_NC}")

    if TIME_COL not in table.columns:
        raise ProtocolError(f"result is missing the '{TIME_COL}' column")

    if len(table) != case.n_steps:
        raise ProtocolError(
            f"result has {len(table)} rows, expected {case.n_steps} "
            "(one row per forcing step, spinup included)"
        )

    _validate_time_axis(table, case)

    missing = [v for v in probe.required_vars if v not in table.columns]
    if missing:
        raise ProtocolError(f"result is missing requested variables: {missing}")

    for var in probe.required_vars:
        col = pd.to_numeric(table[var], errors="coerce")
        # Existing flux/state criteria handle infinities and retain the full
        # scorecard. Keep that behavior; new diagnostics require finite values
        # at the contract boundary, including interval-end temperatures.
        invalid = ~np.isfinite(col) if var in probe.requires_diagnostics else col.isna()
        if invalid.any():
            n_bad = int(invalid.sum())
            raise ProtocolError(f"column '{var}' has {n_bad} non-finite values")
        table[var] = col

    meta: dict[str, Any] = {}
    run_path = io_dir / RUN_FILE
    if not run_path.exists():
        raise ProtocolError(f"adapter did not write required {RUN_FILE}")
    try:
        with open(run_path) as fh:
            meta = json.load(fh)
    except json.JSONDecodeError:
        raise ProtocolError("run.json is not valid JSON") from None
    if not isinstance(meta, dict):
        raise ProtocolError("run.json must contain a JSON object")
    if meta.get("status") != "ok":
        raise ProtocolError(f"adapter reported status={meta.get('status')!r}, expected 'ok'")
    if "n_steps" in meta and meta["n_steps"] != case.n_steps:
        raise ProtocolError(
            f"run.json reports n_steps={meta['n_steps']!r}, expected {case.n_steps}"
        )

    return RunResult(case=case, table=table, meta=meta, wall_seconds=wall_seconds)

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

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from hydroturing.spec import ModelManifest, ProbeSpec, UNITS

REQUEST_FILE = "request.json"
FORCING_FILE = "input/forcing.csv"
STATIC_FILE = "input/static.json"
RESULT_CSV = "output/result.csv"
RESULT_NC = "output/result.nc"
RUN_FILE = "output/run.json"

TIME_COL = "time"


class ProtocolError(RuntimeError):
    """The adapter did not honour the contract."""


@dataclass
class Case:
    """One generated instance of a probe: the forcing, the attributes, the seed."""

    probe_id: str
    seed: int
    forcing: pd.DataFrame
    static: dict[str, Any]
    spinup_days: int

    @property
    def n_steps(self) -> int:
        return len(self.forcing)

    def after_spinup(self, frame: pd.DataFrame) -> pd.DataFrame:
        return frame.iloc[self.spinup_days :].reset_index(drop=True)


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

    request = {
        "case_id": f"{case.probe_id}#seed={case.seed}",
        "seed": case.seed,
        "timestep": probe.timestep,
        "n_steps": case.n_steps,
        "spinup_steps": case.spinup_days,
        "request": {
            "fluxes": list(probe.requires_fluxes),
            "states": list(probe.requires_states),
        },
        "input": {"forcing": FORCING_FILE, "static": STATIC_FILE},
        "output": {"table": RESULT_CSV, "run": RUN_FILE},
        "units": {v: UNITS[v] for v in probe.required_vars if v in UNITS},
        "notes": (
            "States are absolute storages, not tendencies. The harness "
            "differences them itself."
        ),
    }
    request_path = io_dir / REQUEST_FILE
    with open(request_path, "w") as fh:
        json.dump(request, fh, indent=2)
    return request_path


def read_result(io_dir: Path, case: Case, probe: ProbeSpec, wall_seconds: float) -> RunResult:
    """Read and validate what the adapter produced."""
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

    missing = [v for v in probe.required_vars if v not in table.columns]
    if missing:
        raise ProtocolError(f"result is missing requested variables: {missing}")

    for var in probe.required_vars:
        col = pd.to_numeric(table[var], errors="coerce")
        if col.isna().any():
            n_bad = int(col.isna().sum())
            raise ProtocolError(f"column '{var}' has {n_bad} non-finite values")
        table[var] = col

    meta: dict[str, Any] = {}
    run_path = io_dir / RUN_FILE
    if run_path.exists():
        try:
            with open(run_path) as fh:
                meta = json.load(fh)
        except json.JSONDecodeError:
            raise ProtocolError("run.json is not valid JSON") from None
        if meta.get("status") not in (None, "ok"):
            raise ProtocolError(f"adapter reported status={meta.get('status')!r}")

    return RunResult(case=case, table=table, meta=meta, wall_seconds=wall_seconds)

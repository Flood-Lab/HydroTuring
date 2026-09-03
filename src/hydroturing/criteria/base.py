"""Criterion registry and the shared evaluation window.

Every criterion is binary. A probe passes only when all of its criteria pass,
and a model passes HydroTuring only when all probes pass. Criteria still
report the quantity they measured, because a bare bit cannot show that one
model leaks six percent and another forty.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
import pandas as pd

from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

PASS = "pass"
FAIL = "fail"


@dataclass
class CriterionResult:
    name: str
    status: str
    message: str
    value: float | None = None
    threshold: float | None = None
    diagnostics: dict[str, Any] = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return self.status == PASS


CriterionFn = Callable[[RunResult, ProbeSpec, dict], CriterionResult]
CRITERIA: dict[str, CriterionFn] = {}


def criterion(name: str) -> Callable[[CriterionFn], CriterionFn]:
    def decorate(fn: CriterionFn) -> CriterionFn:
        CRITERIA[name] = fn
        return fn

    return decorate


def get(name: str) -> CriterionFn:
    if name not in CRITERIA:
        raise KeyError(
            f"unknown criterion '{name}'. Known: {sorted(CRITERIA)}"
        )
    return CRITERIA[name]


@dataclass
class Window:
    """The post-spinup evaluation window.

    Storage change is the difference between the state at the end of the
    window and the state at the last spinup step, so `state0` sits one row
    before the window starts. Getting this off by one is the classic way to
    manufacture a spurious residual, which is why it lives in one place.
    """

    forcing: pd.DataFrame
    table: pd.DataFrame
    state0: pd.Series
    dt_days: float

    def storage(self, variables: tuple[str, ...]) -> np.ndarray:
        present = [v for v in variables if v in self.table.columns]
        if not present:
            return np.zeros(len(self.table))
        return self.table[present].to_numpy().sum(axis=1)

    def storage_initial(self, variables: tuple[str, ...]) -> float:
        present = [v for v in variables if v in self.state0.index]
        return float(sum(float(self.state0[v]) for v in present))

    def volume(self, series: np.ndarray | pd.Series) -> np.ndarray:
        """Convert a rate (per day) to a per-step depth."""
        return np.asarray(series, dtype=float) * self.dt_days


def make_window(run: RunResult, probe: ProbeSpec) -> Window:
    case = run.case
    start = case.spinup_days
    if start >= len(run.table):
        raise ValueError("spinup consumes the entire record; nothing left to score")

    prior = max(start - 1, 0)
    return Window(
        forcing=case.forcing.iloc[start:].reset_index(drop=True),
        table=run.table.iloc[start:].reset_index(drop=True),
        state0=run.table.iloc[prior],
        dt_days=probe.dt_days,
    )

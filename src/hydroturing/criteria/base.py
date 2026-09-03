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
# A paired criterion is handed every variant of one seed instead of one run,
# because what it asserts is a relationship between them: that doubling the
# rain moves the budget, or that renaming the units does not.
PairedFn = Callable[[dict[str, RunResult], ProbeSpec, dict], CriterionResult]

CRITERIA: dict[str, CriterionFn | PairedFn] = {}
PAIRED: set[str] = set()


def criterion(name: str, paired: bool = False) -> Callable[[CriterionFn], CriterionFn]:
    def decorate(fn):
        CRITERIA[name] = fn
        if paired:
            PAIRED.add(name)
        return fn

    return decorate


def get(name: str):
    if name not in CRITERIA:
        raise KeyError(
            f"unknown criterion '{name}'. Known: {sorted(CRITERIA)}"
        )
    return CRITERIA[name]


def is_paired(name: str) -> bool:
    """Whether this criterion is scored across variants rather than one run.

    Unknown names answer False rather than raising, so that a probe naming a
    criterion that does not exist fails where that is actually diagnosed
    instead of here, in a consistency check about something else.
    """
    return name in PAIRED


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
    """The scored stretch of one run, at that run's own step.

    The step comes from the case rather than the probe, because a paired
    probe may run its variants at different steps and each run's per-step
    depths have to be integrated with its own dt.
    """
    case = run.case
    start = case.spinup_steps
    if start >= len(run.table):
        raise ValueError("spinup consumes the entire record; nothing left to score")

    prior = max(start - 1, 0)
    return Window(
        forcing=case.forcing.iloc[start:].reset_index(drop=True),
        table=run.table.iloc[start:].reset_index(drop=True),
        state0=run.table.iloc[prior],
        dt_days=case.dt_days,
    )


def segments(window: Window, column: str = "_regime") -> list[tuple[str, int, int]]:
    """Contiguous blocks of the scored window that share a forcing label.

    Returns (label, start, stop) with stop exclusive, in the order they occur.
    A regime that appears twice yields two blocks rather than one merged one,
    because storage carries across the record and a budget can only be closed
    over an interval that is actually continuous in time.
    """
    if column not in window.forcing.columns:
        raise ValueError(
            f"expected a '{column}' column in the forcing; the generator for a "
            "regime-aware probe has to label which steps are in range and "
            "which are the extrapolation"
        )

    labels = window.forcing[column].astype(str).to_numpy()
    if len(labels) == 0:
        return []

    blocks: list[tuple[str, int, int]] = []
    start = 0
    for i in range(1, len(labels)):
        if labels[i] != labels[start]:
            blocks.append((str(labels[start]), start, i))
            start = i
    blocks.append((str(labels[start]), start, len(labels)))
    return blocks


def storage_at(window: Window, states: tuple[str, ...], index: int) -> float:
    """Total reported storage one step before `index`, as closure measures it.

    Index 0 means the state carried in from spinup, which lives outside the
    window. Anywhere else it is the previous row. Same off-by-one that Window
    exists to centralise, applied to a block boundary instead of the start.
    """
    if index <= 0:
        return window.storage_initial(states)
    present = [v for v in states if v in window.table.columns]
    if not present:
        return 0.0
    return float(window.table.iloc[index - 1][present].sum())

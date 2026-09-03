"""Verdict roll-up.

Binary all the way up. A criterion passes or fails. A probe passes only when
every criterion passes on every seed. A model passes HydroTuring only when
every probe passes.

The reason a model failed is recorded separately from the verdict, because
VIOLATION and INCOMPLETE mean completely different things scientifically.
VIOLATION says the model reported its budget and the budget did not close.
INCOMPLETE says the model never reported enough to be checked at all, which
is where every streamflow-only model lands today. INCOMPATIBLE says the model
cannot consume the probe as declared. ERROR is reserved for adapter or harness
failures rather than scientific outcomes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

PASS = "PASS"
FAIL = "FAIL"

OK = "OK"
VIOLATION = "VIOLATION"
INCOMPLETE = "INCOMPLETE"
INCOMPATIBLE = "INCOMPATIBLE"
ERROR = "ERROR"

# Worst reason wins when rolling up.
_REASON_RANK = {OK: 0, VIOLATION: 1, INCOMPATIBLE: 2, INCOMPLETE: 3, ERROR: 4}


@dataclass
class CriterionOutcome:
    name: str
    status: str
    message: str
    value: float | None = None
    threshold: float | None = None
    worst_seed: int | None = None
    per_seed: list[dict[str, Any]] = field(default_factory=list)
    diagnostics: dict[str, Any] = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return self.status == "pass"


@dataclass
class ProbeOutcome:
    probe_id: str
    law: str
    verdict: str
    reason: str
    seeds: list[int] = field(default_factory=list)
    criteria: list[CriterionOutcome] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    incompatible: list[str] = field(default_factory=list)
    error: str | None = None
    flags: list[str] = field(default_factory=list)
    authors: list[dict[str, str]] = field(default_factory=list)
    # None when the full record was scored. Otherwise the days asked for and,
    # per seed, the stretch that was actually scored.
    window_days: int | None = None
    windows: list[dict[str, Any]] = field(default_factory=list)

    @property
    def failing(self) -> list[str]:
        return [c.name for c in self.criteria if not c.passed]


@dataclass
class ModelReport:
    model_name: str
    model_version: str
    suite_version: str
    probes: list[ProbeOutcome] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    runner: str = ""

    @property
    def verdict(self) -> str:
        return PASS if all(p.verdict == PASS for p in self.probes) and self.probes else FAIL

    @property
    def reason(self) -> str:
        if not self.probes:
            return ERROR
        return max((p.reason for p in self.probes), key=lambda r: _REASON_RANK[r])

    @property
    def summary(self) -> str:
        passed = sum(1 for p in self.probes if p.verdict == PASS)
        return f"{passed}/{len(self.probes)} probes passed"


def reason_for(
    failing: list[str],
    missing: list[str],
    error: str | None,
    incompatible: list[str] | None = None,
) -> str:
    if error:
        return ERROR
    if missing:
        return INCOMPLETE
    if incompatible:
        return INCOMPATIBLE
    if failing:
        return VIOLATION
    return OK

"""Verdict roll-up.

A criterion passes or fails. A probe passes only when every criterion passes
on every seed. A model passes HydroTuring only when every probe that could be
put to it passes.

The reason a model failed is recorded separately from the verdict, because
VIOLATION and INCOMPLETE mean completely different things scientifically.
VIOLATION says the model reported its budget and the budget did not close.
INCOMPLETE says the model never reported enough to be checked at all, which
is where every streamflow-only model lands on the budget probes. INCOMPATIBLE
says the model cannot consume the probe as declared. ERROR is reserved for
adapter or harness failures rather than scientific outcomes.

Neither INCOMPLETE nor INCOMPATIBLE is a failure. A probe stopped by either
never asked the model anything, so it is not scored: its verdict is N/A, and
the model's verdict is decided by the probes that did ask. A model that is N/A
on the energy probes and passes every other one passes. A model that no probe
could ask anything is N/A itself, since it has not earned a pass.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

PASS = "PASS"
FAIL = "FAIL"
NOT_SCORED = "N/A"

OK = "OK"
VIOLATION = "VIOLATION"
INCOMPLETE = "INCOMPLETE"
INCOMPATIBLE = "INCOMPATIBLE"
ERROR = "ERROR"

# Worst reason wins when rolling up. The scored probes (OK, VIOLATION, ERROR)
# are only ever compared with each other; INCOMPATIBLE and INCOMPLETE decide
# the reason only for a model that no probe could score.
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
    # The criteria the probe exists to score, which a pass is reported by.
    # Empty when the criteria never ran. See ProbeSpec.headline.
    headline: list[str] = field(default_factory=list)

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
    def scored(self) -> list[ProbeOutcome]:
        """The probes that asked the model something, and so decide its verdict."""
        return [p for p in self.probes if p.verdict != NOT_SCORED]

    @property
    def verdict(self) -> str:
        if not self.probes:
            return FAIL
        if not self.scored:
            return NOT_SCORED
        return PASS if all(p.verdict == PASS for p in self.scored) else FAIL

    @property
    def reason(self) -> str:
        if not self.probes:
            return ERROR
        deciding = self.scored or self.probes
        return max((p.reason for p in deciding), key=lambda r: _REASON_RANK[r])

    @property
    def summary(self) -> str:
        passed = sum(1 for p in self.probes if p.verdict == PASS)
        parts = [f"{passed}/{len(self.probes)} probes passed"]
        for reason in (INCOMPLETE, INCOMPATIBLE):
            count = sum(1 for p in self.probes if p.verdict == NOT_SCORED and p.reason == reason)
            if count:
                parts.append(f"{count} {reason}")
        return ", ".join(parts)


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

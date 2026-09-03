"""Report rendering. One machine-readable artifact, one human-readable table.

The verdict is a single bit, as designed. Everything under it stays
quantitative so that a paper can show progress before anyone crosses the
line, and so a failing model's author can see where the budget went.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from hydroturing.scoring import ModelReport, ProbeOutcome

TICK = {"PASS": "PASS", "FAIL": "FAIL"}


def to_dict(report: ModelReport) -> dict[str, Any]:
    return {
        "model": {"name": report.model_name, "version": report.model_version},
        "suite_version": report.suite_version,
        "verdict": report.verdict,
        "reason": report.reason,
        "summary": report.summary,
        "flags": report.flags,
        "probes": [_probe_dict(p) for p in report.probes],
    }


def _probe_dict(probe: ProbeOutcome) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "id": probe.probe_id,
        "law": probe.law,
        "verdict": probe.verdict,
        "reason": probe.reason,
        "seeds": probe.seeds,
    }
    if probe.missing:
        payload["missing"] = probe.missing
    if probe.error:
        payload["error"] = probe.error
    if probe.criteria:
        payload["criteria"] = {
            c.name: {
                "status": c.status,
                "message": c.message,
                "value": c.value,
                "threshold": c.threshold,
                "worst_seed": c.worst_seed,
                "per_seed": c.per_seed,
                "diagnostics": c.diagnostics,
            }
            for c in probe.criteria
        }
    return payload


def write_json(report: ModelReport, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(to_dict(report), fh, indent=2, default=float)
        fh.write("\n")
    return path


def to_markdown(report: ModelReport) -> str:
    lines = [
        f"### HydroTuring `{report.model_name}` v{report.model_version}",
        "",
        f"**{report.verdict}** ({report.reason}) &middot; {report.summary} "
        f"&middot; suite {report.suite_version}",
        "",
        "| Probe | Verdict | Reason | Detail |",
        "| --- | --- | --- | --- |",
    ]
    for probe in report.probes:
        detail = _detail(probe)
        lines.append(
            f"| `{probe.probe_id}` | {TICK[probe.verdict]} | {probe.reason} | {detail} |"
        )
    if report.flags:
        lines += ["", "Flags: " + ", ".join(f"`{f}`" for f in report.flags)]
    lines += [
        "",
        "A model passes HydroTuring only when every criterion of every probe passes.",
    ]
    return "\n".join(lines)


def _detail(probe: ProbeOutcome) -> str:
    if probe.error:
        return probe.error.splitlines()[0][:160]
    if probe.missing:
        return "does not report " + ", ".join(f"`{v}`" for v in probe.missing)
    failing = [c for c in probe.criteria if not c.passed]
    if not failing:
        closure = next((c for c in probe.criteria if c.name == "closure"), None)
        if closure and closure.value is not None:
            return f"residual {closure.value:.3%} of driver"
        return "all criteria pass"
    return "; ".join(f"**{c.name}**: {c.message}" for c in failing)[:400]


def to_text(report: ModelReport) -> str:
    lines = [
        f"{report.model_name} v{report.model_version}  ->  "
        f"{report.verdict} ({report.reason})  [{report.summary}]"
    ]
    for probe in report.probes:
        lines.append(f"  {probe.verdict:<4}  {probe.probe_id}")
        if probe.error:
            lines.append(f"          error: {probe.error.splitlines()[0]}")
        if probe.missing:
            lines.append(f"          missing: {', '.join(probe.missing)}")
        for c in probe.criteria:
            mark = "  ok " if c.passed else " FAIL"
            lines.append(f"      {mark}  {c.name:<18} {c.message}")
    return "\n".join(lines)

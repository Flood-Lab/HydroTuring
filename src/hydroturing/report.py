"""Report rendering. One machine-readable artifact, one human-readable table.

The verdict is a single bit, as designed. Everything under it stays
quantitative so that a paper can show progress before anyone crosses the
line, and so a failing model's author can see where the budget went.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

from hydroturing.scoring import ModelReport, ProbeOutcome

PASS_MARK, FAIL_MARK = "\u2705", "\u274c"
# Same width as each other, so columns line up whichever set is in use.
ASCII_PASS, ASCII_FAIL = "ok  ", "FAIL"


def use_emoji() -> bool:
    """Whether the marks can be printed where this output is going.

    Set HT_ASCII to force the plain words: a terminal that renders emoji at
    the wrong width, or a log someone is grepping, is a good enough reason.
    Otherwise the deciding question is whether the stream can encode them at
    all, because a UnicodeEncodeError halfway through a report is worse than
    a report with no pictures in it.
    """
    if os.environ.get("HT_ASCII"):
        return False
    encoding = getattr(sys.stdout, "encoding", None) or "ascii"
    try:
        PASS_MARK.encode(encoding)
    except (LookupError, UnicodeEncodeError):
        return False
    return True


def mark(passed: bool) -> str:
    """The pass or fail marker, for a column that carries no word of its own."""
    if use_emoji():
        return PASS_MARK if passed else FAIL_MARK
    return ASCII_PASS if passed else ASCII_FAIL


def verdict_label(verdict: str, pad: int = 0) -> str:
    """A verdict with its mark, where the mark adds something.

    Without emoji the marker and the word are the same word, and `FAIL FAIL`
    helps nobody, so the mark is dropped rather than doubled. `pad` widens the
    word, not the label: PASS and FAIL are the same length, so padding the word
    is what keeps the column after it straight in either alphabet.
    """
    word = f"{verdict:<{pad}}" if pad else verdict
    return f"{mark(verdict == 'PASS')} {word}" if use_emoji() else word


def prefix(passed: bool) -> str:
    """A leading mark for a line that already states the outcome in words."""
    return f"{mark(passed)} " if use_emoji() else ""


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
        # Credit travels with the result, not just with the source file.
        "authors": probe.authors,
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
        f"{prefix(report.verdict == 'PASS')}**{report.verdict}** ({report.reason}) "
        f"&middot; {report.summary} &middot; suite {report.suite_version}",
        "",
        "| Probe | Verdict | Reason | Detail |",
        "| --- | --- | --- | --- |",
    ]
    for probe in report.probes:
        detail = _detail(probe)
        lines.append(
            f"| `{probe.probe_id}` | {verdict_label(probe.verdict)} | "
            f"{probe.reason} | {detail} |"
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
    """The terminal report.

    The verdict keeps its word at model and probe level as well as its mark.
    A marker is faster to find; the word is what survives being pasted into an
    issue, and what someone greps a CI log for.
    """
    lines = [
        f"{report.model_name} v{report.model_version}  ->  "
        f"{verdict_label(report.verdict)} ({report.reason})  [{report.summary}]"
    ]
    for probe in report.probes:
        lines.append(f"  {verdict_label(probe.verdict, 4)}  {probe.probe_id}")
        if probe.error:
            lines.append(f"          error: {probe.error.splitlines()[0]}")
        if probe.missing:
            lines.append(f"          missing: {', '.join(probe.missing)}")
        for c in probe.criteria:
            lines.append(f"        {mark(c.passed)}  {c.name:<18} {c.message}")
    return "\n".join(lines)

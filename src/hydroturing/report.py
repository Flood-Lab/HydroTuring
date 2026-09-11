"""Report rendering. One machine-readable artifact, one human-readable table.

A scored verdict is a single bit, as designed, and a probe that could not be
put to the model is N/A rather than either value. Everything under the
verdict stays quantitative so that a paper can show progress before anyone
crosses the line, and so a failing model's author can see where the budget
went.
"""

from __future__ import annotations

import csv
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from hydroturing.scoring import (
    ERROR,
    FAIL,
    INCOMPATIBLE,
    NOT_SCORED,
    OK,
    PASS,
    ModelReport,
    ProbeOutcome,
)

PASS_MARK, FAIL_MARK, NOT_SCORED_MARK = "\u2705", "\u274c", "\u2796"
# Same width as each other, so columns line up whichever set is in use.
ASCII_PASS, ASCII_FAIL, ASCII_NOT_SCORED = "ok  ", "FAIL", "n/a "


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


def mark(passed: bool | None) -> str:
    """The pass or fail marker, for a column that carries no word of its own.

    None marks a probe that was never put to the model, which neither passed
    nor failed.
    """
    if passed is None:
        return NOT_SCORED_MARK if use_emoji() else ASCII_NOT_SCORED
    if use_emoji():
        return PASS_MARK if passed else FAIL_MARK
    return ASCII_PASS if passed else ASCII_FAIL


def _passed(verdict: str) -> bool | None:
    """A verdict as `mark` reads it, N/A being neither."""
    return None if verdict == NOT_SCORED else verdict == PASS


def verdict_label(verdict: str, pad: int = 0) -> str:
    """A verdict with its mark, where the mark adds something.

    Without emoji the marker and the word are the same word, and `FAIL FAIL`
    helps nobody, so the mark is dropped rather than doubled. `pad` widens the
    word, not the label: PASS and FAIL are the same length and N/A is one
    shorter, so padding the word is what keeps the column after it straight in
    either alphabet.
    """
    word = f"{verdict:<{pad}}" if pad else verdict
    return f"{mark(_passed(verdict))} {word}" if use_emoji() else word


def prefix(passed: bool | None) -> str:
    """A leading mark for a line that already states the outcome in words."""
    return f"{mark(passed)} " if use_emoji() else ""


def to_dict(report: ModelReport) -> dict[str, Any]:
    return {
        "model": {"name": report.model_name, "version": report.model_version},
        "suite_version": report.suite_version,
        "runner": report.runner,
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
        # Which stretch of the record was scored. Null is the full record.
        "window": (
            {"days": probe.window_days, "cases": probe.windows}
            if probe.window_days is not None
            else None
        ),
    }
    if probe.missing:
        payload["missing"] = probe.missing
    if probe.incompatible:
        payload["incompatible"] = probe.incompatible
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
        f"{prefix(_passed(report.verdict))}**{report.verdict}** ({report.reason}) "
        f"&middot; {report.summary} &middot; suite {report.suite_version}",
        "",
        "| Probe | Verdict | Reason | Window | Detail |",
        "| --- | --- | --- | --- | --- |",
    ]
    for probe in report.probes:
        detail = _detail(probe)
        lines.append(
            f"| `{probe.probe_id}` | {verdict_label(probe.verdict)} | "
            f"{probe.reason} | {window_label(probe)} | {detail} |"
        )
    if report.flags:
        lines += ["", "Flags: " + ", ".join(f"`{f}`" for f in report.flags)]
    lines += [
        "",
        (
            "A model passes HydroTuring only when every criterion of every probe that "
            "could be put to it passes. A probe it does not report enough for, or "
            "cannot consume, is N/A and counts neither way."
        ),
    ]
    return "\n".join(lines)


def window_label(probe: ProbeOutcome) -> str:
    """`full record`, or the N-day flood event the model was scored on.

    A model that never ran, because it cannot report what the probe needs or
    cannot consume the probe at all, was scored on nothing, and the label
    says so rather than implying a record it never saw.
    """
    if probe.missing or probe.incompatible:
        return "not run"
    if probe.window_days is None:
        return "full record"
    return f"{probe.window_days}-day flood event"


def _window_line(probe: ProbeOutcome) -> str:
    """The scored stretch of each seed, for the terminal report."""
    stretches = "; ".join(
        f"seed {w['seed']}: {w['start']} to {w['end']}" for w in probe.windows
    )
    return f"{window_label(probe)} after spinup" + (f" ({stretches})" if stretches else "")


def _detail(probe: ProbeOutcome, plain: bool = False) -> str:
    """One line on where the verdict came from. `plain` drops the markdown."""
    code = (lambda s: s) if plain else (lambda s: f"`{s}`")
    bold = (lambda s: s) if plain else (lambda s: f"**{s}**")
    if probe.error:
        return probe.error.splitlines()[0][:160]
    parts = []
    if probe.missing:
        parts.append("does not report " + ", ".join(code(v) for v in probe.missing))
    if probe.incompatible:
        parts.append("; ".join(probe.incompatible))
    failing = [c for c in probe.criteria if not c.passed]
    if failing:
        parts.append("; ".join(f"{bold(c.name)}: {c.message}" for c in failing))
    if parts:
        return "; ".join(parts)[:400]
    closure = next((c for c in probe.criteria if c.name == "closure"), None)
    if closure and closure.value is not None:
        return f"residual {closure.value:.3%} of driver"
    return "all criteria pass"


CSV_COLUMNS = [
    "run_date",
    "model",
    "version",
    "suite_version",
    "runner",
    "probe",
    "verdict",
    "reason",
    "window",
    "seeds",
    "detail",
]


def to_csv_rows(report: ModelReport, run_date: str | None = None) -> list[dict[str, str]]:
    """One row per probe, flat enough to live in a spreadsheet for years.

    The archive is a log rather than a leaderboard: every evaluation appends
    its rows, dated, so the same model can be seen before and after a fix and
    the suite can be seen growing around it.
    """
    run_date = run_date or datetime.now(timezone.utc).strftime("%Y-%m-%d")
    rows = []
    for probe in report.probes:
        rows.append({
            "run_date": run_date,
            "model": report.model_name,
            "version": report.model_version,
            "suite_version": report.suite_version,
            "runner": report.runner,
            "probe": probe.probe_id,
            "verdict": probe.verdict,
            "reason": probe.reason,
            "window": window_label(probe),
            "seeds": " ".join(str(s) for s in probe.seeds),
            "detail": _detail(probe, plain=True),
        })
    return rows


def contract_row(
    model_name: str,
    model_version: str,
    suite_version: str,
    runner: str,
    probe_id: str,
    seed: int,
    *,
    result=None,
    error: str | None = None,
    incompatible: list[str] | None = None,
    run_date: str | None = None,
) -> dict[str, str]:
    """The adapter contract check as one archive row.

    For a model that reports only discharge this is the only line saying
    that it was actually built and run on a budget probe: that probe is
    N/A (INCOMPLETE) before the container is ever started.

    A check on a probe the model cannot consume is N/A (INCOMPATIBLE), not
    an ERROR: the adapter was never invoked, and the row says why.
    """
    run_date = run_date or datetime.now(timezone.utc).strftime("%Y-%m-%d")
    if error is not None:
        verdict, reason, window, detail = FAIL, ERROR, "not run", error.splitlines()[0][:160]
    elif incompatible:
        verdict, reason, window = NOT_SCORED, INCOMPATIBLE, "not run"
        detail = "; ".join(incompatible)[:400]
    else:
        case = result.case
        window = (
            f"{case.window['days']}-day flood event" if case.window else "full record"
        )
        verdict, reason = PASS, OK
        detail = (
            f"adapter contract OK: {len(result.table)} rows in {result.wall_seconds:.1f}s, "
            f"columns {', '.join(c for c in result.table.columns if c != 'time')}"
        )
    return {
        "run_date": run_date,
        "model": model_name,
        "version": model_version,
        "suite_version": suite_version,
        "runner": runner,
        "probe": f"{probe_id} (adapter contract)",
        "verdict": verdict,
        "reason": reason,
        "window": window,
        "seeds": str(seed),
        "detail": detail,
    }


def append_csv_rows(rows: list[dict[str, str]], path: str | Path) -> Path:
    """Append rows to a CSV archive, writing the header if it is new."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    is_new = not path.exists() or path.stat().st_size == 0
    with open(path, "a", newline="") as fh:
        # The csv module ends rows in \r\n on every platform unless told
        # otherwise, and the archive is committed as LF text.
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, lineterminator="\n")
        if is_new:
            writer.writeheader()
        writer.writerows(rows)
    return path


def append_csv(report: ModelReport, path: str | Path, run_date: str | None = None) -> Path:
    """Append the report to a CSV archive, one row per probe."""
    return append_csv_rows(to_csv_rows(report, run_date), path)


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
        if probe.incompatible:
            lines.append(f"          incompatible: {'; '.join(probe.incompatible)}")
        if probe.window_days is not None:
            lines.append(f"          window: {_window_line(probe)}")
        for c in probe.criteria:
            lines.append(f"        {mark(c.passed)}  {c.name:<18} {c.message}")
    return "\n".join(lines)

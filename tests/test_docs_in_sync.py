"""The prose that describes the suite has to agree with the suite.

The harness discovers probes and models from the filesystem. The README, the
contributors file, the roadmap, the site and the archive are written by hand
and drift the moment a pull request merges without touching them, which is
exactly what happened on the first contributed probe. These tests turn the
checklist at the end of AGENTS.md into something the probe workflow enforces.
"""

from __future__ import annotations

import csv
import re

import pytest
import yaml

from hydroturing.spec import DIAG_VARS, FLUX_VARS, REPO_ROOT, STATE_VARS, load_model, load_probe

PROBES = {p.id: p for p in (load_probe(x) for x in REPO_ROOT.glob("probes/*/*/probe.yaml"))}
PROBE_IDS = sorted(PROBES)

# Probes merged before the proposal form asked for an affiliation. Each entry
# is a request outstanding with the author; remove it when the answer lands
# in probe.yaml. Nothing merged after this list was written may be added to it.
AFFILIATION_PENDING: set[str] = set()
EVALUATED_MODELS = sorted(
    load_model(p).name
    for p in REPO_ROOT.glob("models/*/model.yaml")
    if not p.parent.name.startswith(("_", "reference_"))
)


def read(name: str) -> str:
    return (REPO_ROOT / name).read_text(encoding="utf-8")


@pytest.mark.parametrize("probe_id", PROBE_IDS)
def test_readme_lists_the_probe(probe_id):
    assert f"[`{probe_id}`](probes/{probe_id})" in read("README.md"), (
        f"{probe_id} is missing from the README probes table"
    )


@pytest.mark.parametrize("probe_id", PROBE_IDS)
def test_contributors_credit_the_probe(probe_id):
    assert f"| `{probe_id}` |" in read("CONTRIBUTORS.md"), (
        f"{probe_id} has no row in CONTRIBUTORS.md; a merged probe earns co-authorship"
    )


@pytest.mark.parametrize("probe_id", PROBE_IDS)
def test_probe_authors_carry_an_affiliation(probe_id):
    """A merged probe is co-authorship, and a paper needs more than a handle."""
    if probe_id in AFFILIATION_PENDING:
        pytest.xfail(f"{probe_id}: affiliation requested from the author, not yet supplied")
    for author in PROBES[probe_id].authors:
        assert author.get("affiliation", "").strip(), (
            f"{probe_id}: author {author.get('name')!r} has no affiliation in probe.yaml; "
            "CONTRIBUTORS.md, CITATION.cff and the paper's author list are built from it"
        )


def _cff_authors() -> set[tuple[str, str]]:
    cff = yaml.safe_load(read("CITATION.cff"))
    return {
        (a.get("given-names", "").strip(), a.get("family-names", "").strip())
        for a in cff.get("authors", [])
    }


@pytest.mark.parametrize("probe_id", PROBE_IDS)
def test_citation_lists_every_probe_author(probe_id):
    """CITATION.cff is the software's author list; a probe author belongs on it."""
    listed = _cff_authors()
    for author in PROBES[probe_id].authors:
        given, _, family = author["name"].strip().rpartition(" ")
        assert (given, family) in listed, (
            f"{probe_id}: author {author['name']!r} is not in CITATION.cff "
            f"(expected given-names {given!r}, family-names {family!r})"
        )


@pytest.mark.parametrize("probe_id", PROBE_IDS)
def test_roadmap_marks_the_probe_merged(probe_id):
    assert f"### `{probe_id}` &middot; **merged**" in read("ROADMAP.md"), (
        f"{probe_id} is not marked merged in ROADMAP.md"
    )


@pytest.mark.parametrize("probe_id", PROBE_IDS)
def test_site_lists_the_probe_as_merged_in_every_language(probe_id):
    rows = re.findall(
        rf'<tr><td class="mono">{re.escape(probe_id)}</td>.*?</tr>', read("site/index.html")
    )
    assert len(rows) == 3, f"{probe_id} should appear once per language on the site, found {len(rows)}"
    for row in rows:
        assert 'class="tag merged"' in row, f"{probe_id} is not tagged merged on the site: {row[:120]}"


def _archived_standings(path=REPO_ROOT / "models" / "result.csv") -> dict[str, dict]:
    """Each model's standing at the version archived last.

    Per model: that version, the merged probes archived at it, and its passes
    and scored probes there. A probe that could not be put to a model is N/A
    and in neither number, so every model is counted out of its own total: a
    runoff-only model out of the probes that can ask about runoff, not out of
    the whole suite. Rows for probes no longer in the suite are not counted.
    """
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if r["probe"] in PROBES]
    version = {r["model"]: r["version"] for r in rows}
    latest = {
        (r["model"], r["probe"]): r["verdict"] for r in rows if r["version"] == version[r["model"]]
    }
    standings = {
        model: {"version": v, "probes": set(), "passed": 0, "scored": 0}
        for model, v in version.items()
    }
    for (model, probe), verdict in latest.items():
        standing = standings[model]
        standing["probes"].add(probe)
        if verdict != "N/A":
            standing["scored"] += 1
            standing["passed"] += verdict == "PASS"
    return standings


def _assert_counts_match_archive(
    where: str,
    listed: list[str],
    found: list[tuple[str, str, str]],
    standings: dict[str, dict] | None = None,
) -> None:
    """Every model row states its passes out of its scored probes, as archived.

    The standing is taken at the version archived last, which has to have been
    run on every merged probe: `ht run --probe` appends a single probe's row,
    and a newer version archived on one probe alone would otherwise count as
    1 of 1. A row may leave the count out only for a model no probe could
    score, which has nothing to count.
    """
    standings = _archived_standings() if standings is None else standings
    partial = [
        f"{model} {standings[model]['version']} ({len(standings[model]['probes'])} of {len(PROBE_IDS)})"
        for model in listed
        if model in standings and len(standings[model]["probes"]) < len(PROBE_IDS)
    ]
    assert not partial, (
        f"{where} states a standing for {partial}, but the newest archived version was not run "
        "on every probe; archive a full run before stating it"
    )
    counted = {model for model, _, _ in found}
    uncounted = [
        model for model in listed if model not in counted and standings.get(model, {}).get("scored")
    ]
    assert not uncounted, f"{where} gives no count for {uncounted}, though the archive scored them"
    wrong = []
    for model, passed, total in found:
        standing = standings.get(model)
        expected = (standing["passed"], standing["scored"]) if standing else None
        if expected != (int(passed), int(total)):
            archive = f"{expected[0]} of {expected[1]}" if expected else "none"
            wrong.append((model, f"{passed} of {total}", archive))
    assert not wrong, f"{where} counts disagree with models/result.csv (model, stated, archive): {wrong}"


def test_site_counts_each_model_out_of_the_probes_it_was_scored_on():
    rows = "".join(re.findall(r'"models\.rows":`(.*?)`', read("site/index.html"), re.DOTALL))
    listed = re.findall(r'<tr><td class="mono">([a-z0-9_]+)</td>', rows)
    found = re.findall(r'<tr><td class="mono">([a-z0-9_]+)</td>.*?<td>(\d+) / (\d+)</td>', rows)
    assert listed, "no model rows found in the site's models table"
    _assert_counts_match_archive("the site", listed, found)


def test_readme_counts_each_model_out_of_the_probes_it_was_scored_on():
    readme = read("README.md")
    listed = re.findall(r"^\| \[`([a-z0-9_]+)`\]\(models/", readme, re.MULTILINE)
    found = re.findall(
        r"^\| \[`([a-z0-9_]+)`\]\(models/[a-z0-9_]+\) \|.*?"
        r"\*\*(?:PASS|FAIL \([A-Z]+\)|N/A \([A-Z]+\))\*\*, (\d+) of (\d+)",
        readme,
        re.MULTILINE,
    )
    assert listed, "no model rows found in the README's models table"
    _assert_counts_match_archive("the README", listed, found)


def test_site_flowchart_shows_every_probe():
    html = read("site/index.html")
    svg = re.search(r'<svg class="infra".*?</svg>', html, re.DOTALL).group(0)
    labels = re.findall(r'data-i18n="(infra\.probe\.[a-z]+)"', svg)
    assert len(labels) == len(PROBE_IDS), (
        f"the flowchart shows {len(labels)} probes, the suite has {len(PROBE_IDS)}"
    )
    for key in labels:
        assert html.count(f'"{key}":') == 3, f"{key} lacks a translation in one of the three languages"


@pytest.mark.parametrize("variable", FLUX_VARS + STATE_VARS + DIAG_VARS)
def test_agents_doc_defines_every_variable(variable):
    assert f"| `{variable}` |" in read("AGENTS.md"), (
        f"{variable} is accepted by spec.py but the table in AGENTS.md does not define it"
    )


def test_archive_has_every_evaluated_model_on_every_probe():
    with open(REPO_ROOT / "models" / "result.csv", newline="", encoding="utf-8") as fh:
        archived = {(row["model"], row["probe"]) for row in csv.DictReader(fh)}
    missing = [
        (model, probe)
        for model in EVALUATED_MODELS
        for probe in PROBE_IDS
        if (model, probe) not in archived
    ]
    assert not missing, f"models/result.csv has no row for: {missing}"

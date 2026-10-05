#!/usr/bin/env python3
"""Write the per-law standings the site's two charts are drawn from.

The site shows the merged probes as a pie split by conservation law, and
each evaluated model as one bar split the same way: how many probes under
each law could score the model, and how many of those it passed. Both come
from the repository. A probe's law is the one its probe.yaml declares, and
a model's standing is the one tests/test_docs_in_sync.py checks the tables
against: the rows archived at the model's newest version in
models/result.csv, with N/A left out of both numbers.

The numbers live in a JSON block inside site/index.html, so the page stays
one self-contained file. Run this after archiving rows or merging a probe:

    python3 scripts/site_standings.py          # rewrite the block in place
    python3 scripts/site_standings.py --check  # exit 1 if the block is stale

Nothing here needs a dependency.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / "site" / "index.html"
LAWS = ("mass", "energy", "momentum")
BLOCK = re.compile(
    r'(<script type="application/json" id="ht-standings">)(.*?)(</script>)', re.DOTALL
)


LAW_LINE = re.compile(r"^law:\s*['\"]?(\w+)", re.MULTILINE)
ID_LINE = re.compile(r"^id:\s*['\"]?([\w./-]+)", re.MULTILINE)


def probe_laws() -> dict[str, str]:
    """Every merged probe's id, mapped to the law its probe.yaml declares.

    Read from the `law:` and `id:` lines rather than the directory, because a
    probe's id need not start with its law (`real/<slug>` is allowed).
    """
    laws = {}
    for path in ROOT.glob("probes/*/*/probe.yaml"):
        text = path.read_text(encoding="utf-8")
        law, probe_id = LAW_LINE.search(text), ID_LINE.search(text)
        if not law or not probe_id:
            raise SystemExit(f"{path} has no top-level law: or id: line")
        if law.group(1) not in LAWS:
            raise SystemExit(f"{path} declares law {law.group(1)!r}, not one of {LAWS}")
        laws[probe_id.group(1)] = law.group(1)
    return laws


def evaluated_models() -> list[str]:
    """The models the site lists: everything but the reference models and the template."""
    return sorted(
        path.parent.name
        for path in ROOT.glob("models/*/model.yaml")
        if not path.parent.name.startswith(("_", "reference_"))
    )


def standings() -> dict:
    laws = probe_laws()
    with open(ROOT / "models" / "result.csv", newline="", encoding="utf-8") as fh:
        rows = [row for row in csv.DictReader(fh) if row["probe"] in laws]
    version = {row["model"]: row["version"] for row in rows}
    latest = {
        (row["model"], row["probe"]): row["verdict"]
        for row in rows
        if row["version"] == version[row["model"]]
    }
    models = {}
    for model in evaluated_models():
        if model not in version:
            continue
        counts = {law: [0, 0] for law in LAWS}
        for (name, probe), verdict in latest.items():
            if name != model or verdict == "N/A":
                continue
            counts[laws[probe]][1] += 1
            counts[laws[probe]][0] += verdict == "PASS"
        models[model] = counts
    probes = {law: sum(1 for value in laws.values() if value == law) for law in LAWS}
    return {"probes": probes, "models": models}


def render(data: dict) -> str:
    """One model per line, so a re-archived model shows up as a one-line diff."""
    lines = [f'{{"probes":{json.dumps(data["probes"], separators=(",", ":"))},', '"models":{']
    items = list(data["models"].items())
    for i, (model, counts) in enumerate(items):
        comma = "," if i < len(items) - 1 else ""
        lines.append(f'"{model}":{json.dumps(counts, separators=(",", ":"))}{comma}')
    lines.append("}}")
    return "\n" + "\n".join(lines) + "\n"


def embedded(html: str) -> dict:
    match = BLOCK.search(html)
    if not match:
        raise SystemExit(f'{SITE} has no <script type="application/json" id="ht-standings"> block')
    return json.loads(match.group(2))


def main(argv: list[str]) -> int:
    html = SITE.read_text(encoding="utf-8")
    data = standings()
    if "--check" in argv:
        if embedded(html) != data:
            print("site/index.html standings are stale; run python3 scripts/site_standings.py")
            return 1
        print("site/index.html standings match models/result.csv")
        return 0
    embedded(html)  # fail loudly if the block is missing
    updated = BLOCK.sub(lambda m: m.group(1) + render(data) + m.group(3), html, count=1)
    SITE.write_text(updated, encoding="utf-8")
    print(f"wrote standings for {len(data['models'])} models and {sum(data['probes'].values())} probes")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

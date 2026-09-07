#!/usr/bin/env python3
"""Write the badge endpoints the README reads its contribution counts from.

Each file is the JSON that shields.io's `endpoint` badge consumes. The
counts come from the repository itself: a probe is a directory under
probes/ with a probe.yaml, a model is a directory under models/ whose
model.yaml runs in a container, i.e. a submitted model rather than a
reference or physical one the gate uses. Run by the pages workflow; nothing here needs a dependency.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def probes() -> int:
    return sum(1 for _ in ROOT.glob("probes/*/*/probe.yaml"))


def models() -> int:
    """Models evaluated against the suite: the submitted ones, which run in
    a container, and the physical models the gate requires every probe to
    pass. The deliberately broken reference models (reference_*) are there
    to trip criteria and are not counted; nor is the template."""
    count = 0
    for path in ROOT.glob("models/*/model.yaml"):
        name = path.parent.name
        if name.startswith(("_", "reference_")):
            continue
        count += 1
    return count


def badge(label: str, count: int) -> dict:
    return {
        "schemaVersion": 1,
        "label": label,
        "message": str(count),
        "color": "1f6f8b",
    }


def main(out_dir: str) -> int:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, count in (("probes", probes()), ("models", models())):
        (out / f"{name}.json").write_text(json.dumps(badge(f"{name} merged" if name == "probes" else f"{name} evaluated", count)))
        print(f"{name}: {count}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "site/badges"))

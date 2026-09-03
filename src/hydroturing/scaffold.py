"""Scaffolding: turn a filled-in draft into a probe directory.

Contributing should not begin with a blank page. The flow is: get a template,
fill in the fields, run one command, get a working directory that already
validates. What is left is the part only the contributor can do, which is the
physics.

    ht init-probe                       writes probe-draft.yaml to fill in
    ht init-probe --from probe-draft.yaml   creates probes/<law>/<slug>/
"""

from __future__ import annotations

import shutil
from pathlib import Path

import yaml

from hydroturing.spec import REPO_ROOT, SpecError

TEMPLATE_DIR = REPO_ROOT / "templates"
PROBES_DIR = REPO_ROOT / "probes"
MODELS_DIR = REPO_ROOT / "models"

GUIDANCE_PREFIX = "#!"


def strip_guidance(text: str) -> str:
    """Remove instructional lines, keep real comments.

    Lines starting with `#!` explain how to fill the template in and have no
    business in the finished probe. Ordinary `#` comments are the
    contributor's own and are kept.
    """
    kept = []
    for line in text.splitlines():
        if line.lstrip().startswith(GUIDANCE_PREFIX):
            continue
        # Removing a guidance block leaves a gap on either side of it. Collapse
        # runs of blank lines so the generated file reads as if hand-written.
        if not line.strip() and kept and not kept[-1].strip():
            continue
        kept.append(line)
    return "\n".join(kept).strip("\n") + "\n"


def write_draft(dest: Path | None = None) -> Path:
    dest = Path(dest) if dest else Path("probe-draft.yaml")
    if dest.exists():
        raise SpecError(f"{dest} already exists; pass -o to choose another path")
    shutil.copyfile(TEMPLATE_DIR / "probe.template.yaml", dest)
    return dest


def _render(template_name: str, **fields: object) -> str:
    return (TEMPLATE_DIR / template_name).read_text().format(**fields)


def scaffold_probe(draft_path: Path, force: bool = False) -> Path:
    draft_path = Path(draft_path)
    raw = draft_path.read_text()
    draft = yaml.safe_load(raw)

    if not isinstance(draft, dict):
        raise SpecError(f"{draft_path}: not a YAML mapping")

    for field in ("id", "law", "case"):
        if field not in draft:
            raise SpecError(f"{draft_path}: missing required field '{field}'")

    probe_id = str(draft["id"])
    if probe_id.count("/") != 1:
        raise SpecError(f"id must be '<law>/<slug>', got '{probe_id}'")
    law, slug = probe_id.split("/")

    if law != draft["law"] and law != "real":
        raise SpecError(
            f"id says '{law}' but law says '{draft['law']}'; they must agree"
        )
    if probe_id == "mass/my-probe":
        raise SpecError(
            "the template id is still 'mass/my-probe'. Fill in the draft before "
            "scaffolding, starting with id, title and authors."
        )

    target = PROBES_DIR / law / slug
    if target.exists() and not force:
        raise SpecError(f"{target} already exists; pass --force to overwrite")
    target.mkdir(parents=True, exist_ok=True)

    (target / "probe.yaml").write_text(strip_guidance(raw))

    case = draft["case"]
    generator = case.get("generator", "generate.py")
    (target / generator).write_text(
        _render(
            "generate.template.py",
            id=probe_id,
            period_years=case.get("period_years", 10),
            spinup_days=case.get("spinup_days", 365),
        )
    )
    (target / "README.md").write_text(_render("README.template.md", id=probe_id))
    return target


def scaffold_model(name: str, force: bool = False) -> Path:
    target = MODELS_DIR / name
    if target.exists() and not force:
        raise SpecError(f"{target} already exists; pass --force to overwrite")

    shutil.copytree(MODELS_DIR / "_template", target, dirs_exist_ok=force)
    for filename in ("model.yaml", "ht_adapter.py"):
        path = target / filename
        path.write_text(path.read_text().replace("_template", name))
    (target / "ht_adapter.py").chmod(0o755)
    return target

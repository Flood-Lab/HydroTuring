"""Scaffolding: turn a filled-in draft into a probe directory.

Contributing should not begin with a blank page. The flow is: get a template,
fill in the fields, run one command, get a working directory that already
validates. What is left is the part only the contributor can do, which is the
physics.

    ht init-probe                       writes probe-draft.yaml to fill in
    ht init-probe --template <kind>     the same, for a probe of that shape
    ht init-probe --from probe-draft.yaml   creates probes/<law>/<slug>/

A template is a pair: `probe.<kind>.template.yaml` and, where the shape needs
one, `generate.<kind>.template.py`. Which pair a draft came from is recorded
in a `#!` guidance line at the top of the draft, so the generator skeleton
matches the criteria without the contributor having to say so twice. Guidance
lines are stripped from the finished probe, so the marker leaves no trace.
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
TEMPLATE_MARKER = "#! template:"
DEFAULT_TEMPLATE = "default"


def available_templates() -> dict[str, Path]:
    """Every probe template in templates/, keyed by the name the CLI takes."""
    found = {DEFAULT_TEMPLATE: TEMPLATE_DIR / "probe.template.yaml"}
    for path in sorted(TEMPLATE_DIR.glob("probe.*.template.yaml")):
        found[path.name[len("probe."):-len(".template.yaml")]] = path
    return found


def template_kind(draft_text: str) -> str:
    """Which template a draft came from, read from its guidance marker."""
    for line in draft_text.splitlines():
        if line.startswith(TEMPLATE_MARKER):
            return line[len(TEMPLATE_MARKER):].strip()
        if line.strip() and not line.lstrip().startswith(GUIDANCE_PREFIX):
            break
    return DEFAULT_TEMPLATE


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


def write_draft(dest: Path | None = None, kind: str = DEFAULT_TEMPLATE) -> Path:
    templates = available_templates()
    if kind not in templates:
        raise SpecError(
            f"no template '{kind}'. Available: {', '.join(sorted(templates))}"
        )
    dest = Path(dest) if dest else Path("probe-draft.yaml")
    if dest.exists():
        raise SpecError(f"{dest} already exists; pass -o to choose another path")
    shutil.copyfile(templates[kind], dest)
    return dest


def _template_ids() -> set[str]:
    """The example ids the templates ship with.

    Scaffolding one of these unedited would create a probe named after the
    template rather than after the physics, and would collide with the next
    contributor to do the same.
    """
    import yaml as _yaml

    ids = set()
    for path in available_templates().values():
        try:
            draft = _yaml.safe_load(strip_guidance(path.read_text()))
        except _yaml.YAMLError:  # pragma: no cover - a broken template is caught by tests
            continue
        if isinstance(draft, dict) and "id" in draft:
            ids.add(str(draft["id"]))
    return ids


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
    if probe_id in _template_ids():
        raise SpecError(
            f"the id is still the template's example, '{probe_id}'. Fill in the "
            "draft before scaffolding, starting with id, title and authors."
        )

    target = PROBES_DIR / law / slug
    if target.exists() and not force:
        raise SpecError(f"{target} already exists; pass --force to overwrite")
    target.mkdir(parents=True, exist_ok=True)

    (target / "probe.yaml").write_text(strip_guidance(raw))

    case = draft["case"]
    generator = case.get("generator", "generate.py")

    # A shaped template gets its matching generator skeleton, so that a probe
    # whose criteria read a `_regime` column or compare variants starts from a
    # generator that produces them. Anything else gets the plain one.
    kind = template_kind(raw)
    shaped = TEMPLATE_DIR / f"generate.{kind}.template.py"
    generator_template = shaped.name if shaped.exists() else "generate.template.py"

    (target / generator).write_text(
        _render(
            generator_template,
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

"""Discovery of probes and models on disk."""

from __future__ import annotations

from pathlib import Path

from hydroturing.spec import REPO_ROOT, ModelManifest, ProbeSpec, load_model, load_probe

PROBES_DIR = REPO_ROOT / "probes"
MODELS_DIR = REPO_ROOT / "models"


def probe_paths(root: Path | None = None) -> list[Path]:
    root = root or PROBES_DIR
    return sorted(p.parent for p in root.glob("*/*/probe.yaml"))


def model_paths(root: Path | None = None) -> list[Path]:
    root = root or MODELS_DIR
    return sorted(
        p.parent for p in root.glob("*/model.yaml") if not p.parent.name.startswith("_")
    )


def all_probes() -> list[ProbeSpec]:
    return [load_probe(p) for p in probe_paths()]


def all_models() -> list[ModelManifest]:
    return [load_model(p) for p in model_paths()]


def find_probe(identifier: str) -> ProbeSpec:
    for probe in all_probes():
        if identifier in (probe.id, probe.slug, probe.path.name):
            return probe
    known = ", ".join(p.id for p in all_probes()) or "none"
    raise KeyError(f"no probe matching '{identifier}'. Known probes: {known}")


def find_model(name: str) -> ModelManifest:
    for model in all_models():
        if model.name == name:
            return model
    known = ", ".join(m.name for m in all_models()) or "none"
    raise KeyError(f"no model named '{name}'. Known models: {known}")

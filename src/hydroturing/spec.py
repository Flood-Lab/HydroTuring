"""Specifications for probes and models, loaded from YAML and schema-validated.

Deliberately dataclass-based rather than pydantic: the harness must run with
only numpy, pandas, PyYAML and jsonschema so that CI stays fast and a
contributor can run it without building a scientific Python stack.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import jsonschema
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_DIR = REPO_ROOT / "schemas"

# Canonical variable names. Fluxes are mm per timestep-day; states are mm.
# `dis` is the one exception and is m3 s-1, because that is what models report.
FLUX_VARS = ("pr", "evspsbl", "mrro", "dis")
STATE_VARS = ("mrso", "snw", "canopy")

UNITS = {
    "pr": "mm day-1",
    "evspsbl": "mm day-1",
    "mrro": "mm day-1",
    "dis": "m3 s-1",
    "mrso": "mm",
    "snw": "mm",
    "canopy": "mm",
}

TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0}


class SpecError(ValueError):
    """A probe or model manifest is malformed."""


def _load_schema(name: str) -> dict[str, Any]:
    with open(SCHEMA_DIR / name) as fh:
        return json.load(fh)


def _validate(payload: dict[str, Any], schema_name: str, source: Path) -> None:
    try:
        jsonschema.validate(payload, _load_schema(schema_name))
    except jsonschema.ValidationError as exc:
        loc = "/".join(str(p) for p in exc.absolute_path) or "<root>"
        raise SpecError(f"{source}: at {loc}: {exc.message}") from None


@dataclass(frozen=True)
class Criterion:
    """One pass/fail check. Every criterion is binary; there is no partial credit."""

    name: str
    params: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ProbeSpec:
    id: str
    title: str
    law: str
    track: str
    version: int
    requires_fluxes: tuple[str, ...]
    requires_states: tuple[str, ...]
    generator: str
    n_seeds: int
    timestep: str
    period_years: float
    spinup_days: int
    max_output_mb: float
    max_runtime_s: float
    criteria: tuple[Criterion, ...]
    must_pass: tuple[str, ...]
    must_fail: dict[str, str]
    provenance: str
    path: Path

    @property
    def required_vars(self) -> tuple[str, ...]:
        return self.requires_fluxes + self.requires_states

    @property
    def slug(self) -> str:
        return self.id.replace("/", "__")

    @property
    def generator_path(self) -> Path:
        return self.path / self.generator

    @property
    def dt_days(self) -> float:
        return TIMESTEP_DAYS[self.timestep]


@dataclass(frozen=True)
class ModelManifest:
    name: str
    version: str
    entrypoint: tuple[str, ...]
    timestep: str
    emits_fluxes: tuple[str, ...]
    emits_states: tuple[str, ...]
    runner: str
    needs_forcing: tuple[str, ...]
    supports_perturbation: bool
    resources: dict[str, Any]
    authors: tuple[str, ...]
    license: str
    description: str
    path: Path

    @property
    def emitted(self) -> tuple[str, ...]:
        return self.emits_fluxes + self.emits_states

    def missing_for(self, probe: ProbeSpec) -> list[str]:
        """Variables the probe needs that this model never reports.

        A non-empty result means the verdict is FAIL with reason INCOMPLETE:
        the model cannot demonstrate conservation because it never says
        enough to be checked.
        """
        return [v for v in probe.required_vars if v not in self.emitted]


def load_probe(path: str | Path) -> ProbeSpec:
    path = Path(path)
    spec_file = path / "probe.yaml" if path.is_dir() else path
    directory = spec_file.parent
    with open(spec_file) as fh:
        raw = yaml.safe_load(fh)
    if not isinstance(raw, dict):
        raise SpecError(f"{spec_file}: probe.yaml must be a mapping")
    _validate(raw, "probe.schema.json", spec_file)

    expected_id = f"{directory.parent.name}/{directory.name}"
    if raw["id"] != expected_id:
        raise SpecError(
            f"{spec_file}: id '{raw['id']}' does not match its location "
            f"(expected '{expected_id}')"
        )

    case = raw["case"]
    if not (directory / case["generator"]).exists():
        raise SpecError(f"{spec_file}: generator '{case['generator']}' not found")

    criteria = []
    for item in raw["criteria"]:
        (name, params), = item.items()
        criteria.append(Criterion(name=name, params=params or {}))

    names = [c.name for c in criteria]
    if len(names) != len(set(names)):
        raise SpecError(f"{spec_file}: duplicate criterion names in `criteria`")

    for model, criterion in raw["baselines"]["must_fail"].items():
        if criterion not in names:
            raise SpecError(
                f"{spec_file}: baselines.must_fail['{model}'] names criterion "
                f"'{criterion}', which this probe does not define"
            )

    requires = raw.get("requires", {})
    return ProbeSpec(
        id=raw["id"],
        title=raw["title"],
        law=raw["law"],
        track=raw["track"],
        version=raw["version"],
        requires_fluxes=tuple(requires.get("fluxes", [])),
        requires_states=tuple(requires.get("states", [])),
        generator=case["generator"],
        n_seeds=case["n_seeds"],
        timestep=case["timestep"],
        period_years=case["period_years"],
        spinup_days=case["spinup_days"],
        max_output_mb=case.get("max_output_mb", 5.0),
        max_runtime_s=case.get("max_runtime_s", 120.0),
        criteria=tuple(criteria),
        must_pass=tuple(raw["baselines"]["must_pass"]),
        must_fail=dict(raw["baselines"]["must_fail"]),
        provenance=raw.get("provenance", ""),
        path=directory,
    )


def load_model(path: str | Path) -> ModelManifest:
    path = Path(path)
    spec_file = path / "model.yaml" if path.is_dir() else path
    directory = spec_file.parent
    with open(spec_file) as fh:
        raw = yaml.safe_load(fh)
    if not isinstance(raw, dict):
        raise SpecError(f"{spec_file}: model.yaml must be a mapping")
    _validate(raw, "model.schema.json", spec_file)

    if raw["name"] != directory.name:
        raise SpecError(
            f"{spec_file}: name '{raw['name']}' does not match directory "
            f"'{directory.name}'"
        )

    unknown = [v for v in raw["emits"]["fluxes"] if v not in FLUX_VARS]
    unknown += [v for v in raw["emits"]["states"] if v not in STATE_VARS]
    if unknown:
        raise SpecError(f"{spec_file}: unknown variables in emits: {unknown}")

    return ModelManifest(
        name=raw["name"],
        version=str(raw["version"]),
        entrypoint=tuple(raw["entrypoint"]),
        timestep=raw["timestep"],
        emits_fluxes=tuple(raw["emits"]["fluxes"]),
        emits_states=tuple(raw["emits"]["states"]),
        runner=raw.get("runner", "docker"),
        needs_forcing=tuple(raw.get("needs_forcing", [])),
        supports_perturbation=bool(raw.get("supports", {}).get("perturbation", False)),
        resources=raw.get("resources", {}),
        authors=tuple(raw.get("authors", [])),
        license=raw.get("license", ""),
        description=raw.get("description", ""),
        path=directory,
    )

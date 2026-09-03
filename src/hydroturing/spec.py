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

# Supported timesteps as ISO 8601 durations, and their length in days. Fluxes
# are always rates in mm per day whatever the step, so a per-step depth is
# the rate times the step length.
TIMESTEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}

# How much of the scored record a submitted model is evaluated on when its
# manifest does not say. A heavy model is tested on the largest flood event
# of the generated record, spinup included, rather than the full period: a
# month of daily output or a week of hourly output is enough to see the
# event and short enough to fit the container's time budget. Reference
# models always see the full record, because the acceptance gate runs on it.
DEFAULT_WINDOW_DAYS = {"PT1D": 30, "PT1H": 7, "PT15M": 7, "PT5M": 7, "PT1M": 7}
FULL_WINDOW = "full"

# These repository-owned baselines are the only code allowed to bypass the
# container boundary. A submitted manifest cannot opt itself into host access.
TRUSTED_SUBPROCESS_MODELS = {
    "reference_bucket",
    "reference_calendar",
    "reference_cheater",
    "reference_degenerate",
    "reference_in_sample",
    "reference_leaky",
    "reference_streamflow_only",
    "reference_fixed_step",
}


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
    authors: tuple[dict[str, str], ...]
    citation: str
    requires_fluxes: tuple[str, ...]
    requires_states: tuple[str, ...]
    generator: str
    n_seeds: int
    timestep: str
    period_years: float
    spinup_days: int
    max_output_mb: float
    max_runtime_s: float
    variants: tuple[str, ...]
    criteria: tuple[Criterion, ...]
    must_pass: tuple[str, ...]
    must_fail: dict[str, str]
    provenance: str
    path: Path
    # Length of the scored record in days. `period_years` is the usual way
    # to say it; a short sub-daily probe says `period_days` instead.
    period_days: float = 0.0
    # A paired probe may run its variants at different steps, which is how a
    # resolution transform is expressed. Absent variants use `timestep`.
    variant_timesteps: dict[str, str] = field(default_factory=dict)

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

    def timestep_for(self, variant: str | None) -> str:
        """The step a variant runs at; the probe's own step unless declared."""
        if variant is None:
            return self.timestep
        return self.variant_timesteps.get(variant, self.timestep)

    @property
    def timesteps(self) -> tuple[str, ...]:
        """Every step this probe runs a model at, control first."""
        steps = [self.timestep]
        for variant in self.variants:
            step = self.timestep_for(variant)
            if step not in steps:
                steps.append(step)
        return tuple(steps)

    def spinup_steps_for(self, variant: str | None = None) -> int:
        return int(round(self.spinup_days / TIMESTEP_DAYS[self.timestep_for(variant)]))

    def n_steps_for(self, variant: str | None = None) -> int:
        """Rows a generator must produce for a variant: spinup plus the period."""
        dt = TIMESTEP_DAYS[self.timestep_for(variant)]
        return int(round(self.period_days / dt)) + self.spinup_steps_for(variant)

    @property
    def control(self) -> str | None:
        """The variant every single-run criterion is scored against.

        None for an ordinary probe, which has one case per seed. For a paired
        probe the first declared variant is the control, and the others are
        the counterfactuals it is compared with.
        """
        return self.variants[0] if self.variants else None


@dataclass(frozen=True)
class ModelManifest:
    name: str
    version: str
    entrypoint: tuple[str, ...]
    # Every step the model can be run at, native step first. A model that
    # only works at one resolution lists one, and is INCOMPATIBLE with any
    # probe that needs another: that is a finding, not a failure to run.
    timesteps: tuple[str, ...]
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
    # Days of the scored record the model is evaluated on: a positive integer,
    # FULL_WINDOW for the whole record, or None to take the default for the
    # kind of model (see DEFAULT_WINDOW_DAYS).
    window_days: int | str | None = None

    @property
    def timestep(self) -> str:
        """The native step, which is what a single-step probe is matched on."""
        return self.timesteps[0]

    def supports_timestep(self, timestep: str) -> bool:
        return timestep in self.timesteps

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

    # JSON Schema can validate the shape of a criterion declaration, but the
    # executable registry is the authority on which criterion names exist.
    from hydroturing import criteria as criteria_mod  # noqa: PLC0415

    unknown = [name for name in names if name not in criteria_mod.CRITERIA]
    if unknown:
        raise SpecError(
            f"{spec_file}: unknown criteria {unknown}; known criteria are "
            f"{sorted(criteria_mod.CRITERIA)}"
        )

    for model, criterion in raw["baselines"]["must_fail"].items():
        if criterion not in names:
            raise SpecError(
                f"{spec_file}: baselines.must_fail['{model}'] names criterion "
                f"'{criterion}', which this probe does not define"
            )

    variants = tuple(case.get("variants", []))
    _check_variants(spec_file, criteria, variants)

    if "period_days" in case:
        period_days = float(case["period_days"])
    else:
        period_days = float(case["period_years"]) * 365.0
    period_years = float(case.get("period_years", period_days / 365.0))

    variant_timesteps = dict(case.get("timesteps", {}))
    unknown_variants = [v for v in variant_timesteps if v not in variants]
    if unknown_variants:
        raise SpecError(
            f"{spec_file}: case.timesteps names {unknown_variants}, which are not "
            f"in case.variants {list(variants)}"
        )
    if variants and variant_timesteps.get(variants[0], case["timestep"]) != case["timestep"]:
        raise SpecError(
            f"{spec_file}: the control variant '{variants[0]}' must run at "
            f"case.timestep ({case['timestep']})"
        )

    requires = raw.get("requires", {})
    return ProbeSpec(
        id=raw["id"],
        title=raw["title"],
        law=raw["law"],
        track=raw["track"],
        version=raw["version"],
        authors=tuple(raw["authors"]),
        citation=raw.get("citation", ""),
        requires_fluxes=tuple(requires.get("fluxes", [])),
        requires_states=tuple(requires.get("states", [])),
        generator=case["generator"],
        n_seeds=case["n_seeds"],
        timestep=case["timestep"],
        period_years=period_years,
        spinup_days=case["spinup_days"],
        period_days=period_days,
        variant_timesteps=variant_timesteps,
        max_output_mb=case.get("max_output_mb", 5.0),
        max_runtime_s=case.get("max_runtime_s", 120.0),
        variants=tuple(case.get("variants", [])),
        criteria=tuple(criteria),
        must_pass=tuple(raw["baselines"]["must_pass"]),
        must_fail=dict(raw["baselines"]["must_fail"]),
        provenance=raw.get("provenance", ""),
        path=directory,
    )


def _check_variants(spec_file: Path, criteria: list[Criterion], variants: tuple[str, ...]) -> None:
    """Paired criteria and `case.variants` have to agree.

    A paired criterion compares the model's answer across two runs of the same
    seed. Declaring one without variants asks for a comparison with nothing to
    compare against; declaring variants without one runs the model twice and
    then ignores the second answer. Both are silent no-ops at run time, which
    is why they are errors here.

    Imported inside the function: the criteria package imports this module, so
    the registry is only populated once this one has finished loading.
    """
    from hydroturing.criteria.base import is_paired  # noqa: PLC0415

    paired = [c.name for c in criteria if is_paired(c.name)]
    if paired and not variants:
        raise SpecError(
            f"{spec_file}: criteria {paired} compare runs of the same seed, "
            "so `case.variants` must name the cases to compare"
        )
    if variants and not paired:
        raise SpecError(
            f"{spec_file}: `case.variants` runs the model {len(variants)} times "
            "per seed, but no criterion compares the results"
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

    runner = raw.get("runner", "docker")
    if runner == "subprocess" and raw["name"] not in TRUSTED_SUBPROCESS_MODELS:
        raise SpecError(
            f"{spec_file}: runner 'subprocess' is reserved for trusted reference "
            "models; submitted models must use runner 'docker'"
        )

    unknown = [v for v in raw["emits"]["fluxes"] if v not in FLUX_VARS]
    unknown += [v for v in raw["emits"]["states"] if v not in STATE_VARS]
    if unknown:
        raise SpecError(f"{spec_file}: unknown variables in emits: {unknown}")

    declared = raw["timestep"]
    timesteps = tuple(declared) if isinstance(declared, list) else (declared,)
    if len(set(timesteps)) != len(timesteps):
        raise SpecError(f"{spec_file}: timestep lists a step twice")

    return ModelManifest(
        name=raw["name"],
        version=str(raw["version"]),
        entrypoint=tuple(raw["entrypoint"]),
        timesteps=timesteps,
        emits_fluxes=tuple(raw["emits"]["fluxes"]),
        emits_states=tuple(raw["emits"]["states"]),
        runner=runner,
        needs_forcing=tuple(raw.get("needs_forcing", [])),
        supports_perturbation=bool(raw.get("supports", {}).get("perturbation", False)),
        resources=raw.get("resources", {}),
        authors=tuple(raw.get("authors", [])),
        license=raw.get("license", ""),
        description=raw.get("description", ""),
        path=directory,
        window_days=raw.get("window_days"),
    )

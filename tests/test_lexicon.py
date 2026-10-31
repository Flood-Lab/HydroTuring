"""Everything that names a variable has to agree with schemas/lexicon.yaml."""

import ast
import re

import pytest
import yaml

from hydroturing import harness, registry, spec

LEXICON = yaml.safe_load((spec.SCHEMA_DIR / "lexicon.yaml").read_text())
OUTPUT, FORCING, STATIC = LEXICON["output"], LEXICON["forcing"], LEXICON["static"]
MODELS = registry.all_models()
PROBES = registry.all_probes()


def _agents_units(text: str) -> str:
    """AGENTS.md writes units for people, the lexicon as spec.py does."""
    for prose, unit in (("mm/day", "mm day-1"), ("W/m2", "W m-2"), ("m3/s", "m3 s-1")):
        text = text.replace(prose, unit)
    return text


def test_every_output_says_what_it_means_over_a_row():
    for name, entry in OUTPUT.items():
        assert entry["over_the_row"] in {"mean", "end", "instant"}, name
        if entry["group"] == "state":
            assert entry["over_the_row"] == "end", name


def test_agents_table_agrees_with_the_lexicon():
    agents = (spec.REPO_ROOT / "AGENTS.md").read_text()
    rows = re.findall(r"^\| `([^`]+)` \|.*\| ([^|]+?) \|$", agents, re.M)
    assert rows
    for name, units in rows:
        entry = OUTPUT.get(name, FORCING.get(name))
        assert entry is not None, f"{name} is in AGENTS.md's table but not the lexicon"
        assert _agents_units(units) == entry["units"], name
    assert {name for name, _ in rows} >= set(OUTPUT)


@pytest.mark.parametrize("model", MODELS, ids=lambda m: m.name)
def test_manifest_uses_lexicon_names(model):
    assert set(model.needs_forcing) | set(model.uses_forcing) <= set(FORCING)
    assert set(model.needs_static) | set(model.uses_static) <= set(STATIC)
    assert set(model.emitted) | set(model.emits_routing) <= set(OUTPUT)


@pytest.mark.parametrize("probe", PROBES, ids=lambda p: p.id)
def test_probe_uses_lexicon_names(probe):
    assert set(probe.requires_forcing) <= set(FORCING)
    assert set(probe.requires_static) <= set(STATIC)
    named = (*probe.requires_fluxes, *probe.requires_states, *probe.requires_diagnostics, *probe.requires_routing)
    assert set(named) <= set(OUTPUT)


@pytest.mark.parametrize("probe", PROBES, ids=lambda p: p.id)
def test_generated_case_uses_lexicon_names(probe):
    case = harness.build_case(probe, 1, probe.variants[0] if probe.variants else None)
    visible = {c for c in case.forcing.columns if not c.startswith("_")} - {"time"}
    assert visible <= set(FORCING)
    assert set(case.static) <= set(STATIC)


def _static_reads(source: str) -> set[str]:
    """Names read as static["x"], static.get("x") or "x" in static."""
    found = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Subscript) and getattr(node.value, "id", "") == "static":
            if isinstance(node.slice, ast.Constant):
                found.add(node.slice.value)
        elif (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "get"
              and getattr(node.func.value, "id", "") == "static" and isinstance(node.args[0], ast.Constant)):
            found.add(node.args[0].value)
        elif (isinstance(node, ast.Compare) and getattr(node.comparators[0], "id", "") == "static"
              and isinstance(node.left, ast.Constant)):
            found.add(node.left.value)
    return found


@pytest.mark.parametrize("model", [m for m in MODELS if (m.path / "ht_adapter.py").exists()], ids=lambda m: m.name)
def test_adapter_reads_only_lexicon_static_names(model):
    # A misspelt optional name is the dangerous one: static.get() would quietly
    # take the fallback on every case and nothing else would notice.
    assert _static_reads((model.path / "ht_adapter.py").read_text()) <= set(STATIC)


@pytest.mark.parametrize("model", [m for m in MODELS if not m.name.startswith("reference_")], ids=lambda m: m.name)
def test_crosswalk_covers_what_the_model_emits(model):
    crosswalk = LEXICON["crosswalk"][model.name]
    assert set(crosswalk) == (set(model.emitted) | set(model.emits_routing)) - {"pr"}
    package = "\n".join(
        path.read_text(errors="ignore")
        for path in model.path.rglob("*")
        if path.is_file() and path.name != "README.md" and path.suffix != ".csv" and "__pycache__" not in path.parts
    )
    for native in crosswalk.values():
        for token in re.findall(r"[A-Za-z_][A-Za-z_0-9]*", native):
            assert token in package, f"{model.name}'s crosswalk names {token}, which its package never mentions"


def _dicts_holding_status(tree: ast.AST) -> list[ast.Dict]:
    return [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Dict) and any(isinstance(k, ast.Constant) and k.value == "status" for k in node.keys)
    ]


@pytest.mark.parametrize("model", [m for m in MODELS if (m.path / "ht_adapter.py").exists()], ids=lambda m: m.name)
def test_zero_stores_and_not_reported_sit_at_the_top_of_run_json(model):
    tree = ast.parse((model.path / "ht_adapter.py").read_text())
    top = {k.value for d in _dicts_holding_status(tree) for k in d.keys if isinstance(k, ast.Constant)}
    named = {
        k.value for node in ast.walk(tree) if isinstance(node, ast.Dict) for k in node.keys
        if isinstance(k, ast.Constant) and k.value in ("zero_stores", "not_reported")
    }
    assert named <= top

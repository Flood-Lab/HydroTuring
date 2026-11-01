#!/usr/bin/env python3
"""Reproduce native MODFLOW evidence for the confined-aquifer storage envelope.

The benchmark adapter and executable are unchanged. An observer copies native
files after MODFLOW finishes and before its temporary directory is removed.
Default-profile result.csv must be byte-identical to a standard adapter run.
FloPy is an optional dependency needed only to read those native binary files.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from hydroturing import registry
from hydroturing.criteria import CRITERIA
from hydroturing.harness import build_case
from hydroturing.protocol import RunResult, read_result
from hydroturing.runner.docker_runner import DockerRunner, image_tag
from hydroturing.seeds import gate_seeds

PROBE_ID = "mass/aquifer-equilibrium-bounds"
NATIVE_HEAD_AUDIT_LIMIT_M = 1.0e-7
NATIVE_MAPPING_AUDIT_LIMIT_MM = 1.0e-6

# This observer only copies completed native output. It neither supplies a
# residual-derived boundary flux nor modifies the native integration.
CAPTURE = r'''
import argparse
from pathlib import Path
import shutil
import sys

sys.path.insert(0, "/model")
import ht_adapter as adapter
import flopy

parser = argparse.ArgumentParser()
parser.add_argument("--request", required=True)
parser.add_argument("--profile", default="standard", choices=("standard", "tight_solver", "grid20"))
args = parser.parse_args()
destination = Path(args.request).parent / "output" / "native"
destination.mkdir(parents=True, exist_ok=True)

if args.profile == "grid20":
    adapter.NROW = adapter.NCOL = 20
    adapter.RIVER_ROW = adapter.RIVER_COL = 10
if args.profile == "tight_solver":
    original_ims = flopy.mf6.ModflowIms
    def tight_ims(*a, **kw):
        kw["outer_dvclose"] = kw["inner_dvclose"] = 1e-10
        return original_ims(*a, **kw)
    flopy.mf6.ModflowIms = tight_ims

original_run = flopy.mf6.MFSimulation.run_simulation
def observed_run(simulation, *a, **kw):
    result = original_run(simulation, *a, **kw)
    if result[0]:
        source = Path(simulation.simulation_data.mfpath.get_sim_path())
        for path in source.iterdir():
            if path.suffix in (".hds", ".cbc", ".lst", ".sto", ".ims", ".ic", ".riv"):
                shutil.copy2(path, destination / path.name)
    return result
flopy.mf6.MFSimulation.run_simulation = observed_run
sys.argv = [sys.argv[0], "--request", args.request]
adapter.main()
'''


def compare_standard_output(standard: Path, captured: Path) -> bool:
    """Reject a diagnostic observer that changes the adapter's reported values."""
    if standard.read_bytes() != captured.read_bytes():
        raise ValueError("diagnostic capture changes the standard result.csv")
    return True


class NativeCaptureRunner(DockerRunner):
    """Extract audit-only files, then apply the unchanged output contract.

    The container has the normal isolated mounts. Only the observer writes
    extra files; the host moves those outside /io before the ordinary parser
    validates result.csv and run.json. A separate standard run checks that
    the observer did not change any reported value.
    """

    def invoke(self, model, probe, io_dir, request_path):
        super().invoke(model, probe, io_dir, request_path)
        source = io_dir / "output/native"
        destination = io_dir.with_name(io_dir.name + "-native")
        if destination.exists():
            raise ValueError(f"native audit destination already exists: {destination}")
        source.rename(destination)


def native_metrics(forcing, table, heads, sto_sy, static, *,
                   native_initial_head=None, native_stage=None, native_bed=None,
                   native_storativity=None) -> dict:
    """Check native head, connected-bed and storage-law evidence independently."""
    stage = forcing["sw_stage_m"].to_numpy(dtype=float)
    inventory = table["gw"].to_numpy(dtype=float)
    heads = np.asarray(heads, dtype=float).reshape(len(stage), -1)
    sto_sy = np.asarray(sto_sy, dtype=float).reshape(len(stage), -1)
    if not all(np.isfinite(a).all() for a in (stage, inventory, heads, sto_sy)):
        raise ValueError("native audit contains non-finite data")
    h0 = float(static["aquifer_initial_head_m"])
    top = float(static["aquifer_top_m"])
    bottom = float(static["aquifer_bottom_m"])
    bed_offset = float(static["river_bottom_offset_m"])
    native_initial_head = h0 if native_initial_head is None else native_initial_head
    native_stage = stage if native_stage is None else np.asarray(native_stage, dtype=float)
    native_bed = stage - bed_offset if native_bed is None else np.asarray(native_bed, dtype=float)
    if (not np.isfinite(native_initial_head)
            or not np.isfinite(native_stage).all() or not np.isfinite(native_bed).all()
            or native_stage.shape != stage.shape or native_bed.shape != stage.shape):
        raise ValueError("native serialized inputs must be finite and aligned")
    area = float(static["area_km2"]) * 1.0e6
    sigma = float(static["aquifer_storage_coefficient"])
    native_storativity = sigma if native_storativity is None else native_storativity
    if not np.isfinite(native_storativity) or native_storativity <= 0.0:
        raise ValueError("native storativity must be finite and positive")
    initial_inventory = float(static["aquifer_specific_yield"]) * (h0 - bottom) * 1000.0
    with np.errstate(over="ignore", invalid="ignore"):
        clearance = heads - native_bed[:, None]
        confinement = heads - top
        lower, upper = min(h0, float(stage.min())), max(h0, float(stage.max()))
        requested_excess = max(float(np.max(lower - heads)), float(np.max(heads - upper)), 0.0)
        native_lower = min(native_initial_head, float(native_stage.min()))
        native_upper = max(native_initial_head, float(native_stage.max()))
        envelope_excess = max(float(np.max(native_lower - heads)), float(np.max(heads - native_upper)), 0.0)
        input_roundoff = max(abs(native_initial_head - h0), float(np.max(np.abs(native_stage - stage))))
        sy_gross = float(np.abs(sto_sy).sum() / area * 1000.0)
        predicted_change = 1000.0 * native_storativity * (heads.mean(axis=1) - native_initial_head)
        mapping_error = float(np.max(np.abs(inventory - initial_inventory - predicted_change)))
        requested_mapping_error = float(np.max(np.abs(
            inventory - initial_inventory - 1000.0 * sigma * (heads.mean(axis=1) - h0))))
    values = [float(heads.min()), float(heads.max()), float(confinement.min()),
              float(clearance.min()), envelope_excess, sy_gross, mapping_error,
              requested_excess, input_roundoff, requested_mapping_error,
              abs(native_storativity - sigma)]
    if not np.isfinite(values).all():
        raise ValueError("native audit arithmetic overflows")
    if float(confinement.min()) <= 0.0:
        raise ValueError("native heads leave the confined regime")
    if float(clearance.min()) < 1.0:
        raise ValueError("native riverbed clearance is below 1 m")
    # These are audit precision checks, separate from the probe's storage allowance.
    # 1e-7 m is ten times the standard IMS head tolerance; the 1e-6 mm
    # mapping check also covers double-precision output and reductions.
    if envelope_excess > NATIVE_HEAD_AUDIT_LIMIT_M:
        raise ValueError("native heads leave the prescribed head envelope")
    if sy_gross > 1.0e-8:
        raise ValueError("native STO-SY is active in a supposedly confined case")
    if mapping_error > NATIVE_MAPPING_AUDIT_LIMIT_MM:
        raise ValueError("native head change disagrees with reported confined storage")
    return dict(zip((
        "minimum_head_m", "maximum_head_m", "minimum_head_above_top_m",
        "minimum_riverbed_clearance_m", "native_head_envelope_excess_m",
        "sto_sy_gross_mm", "storage_head_mapping_max_error_mm",
        "requested_head_envelope_excess_m", "input_head_roundoff_m",
        "requested_storage_head_mapping_max_error_mm", "input_storativity_roundoff",
    ), values))


def read_native_inputs(directory: Path, forcing, static) -> dict:
    """Verify FloPy's serialized heads rather than assuming text is lossless.

    The pinned adapter writes one uniform IC constant at eight decimal places,
    and one RIV record per day in eight-place scientific notation. Reject any
    other value; the rounded inputs explain precision, never changed physics.
    """
    initial_text = (directory / "ht_gw.ic").read_text()
    constants = re.findall(r"\bCONSTANT\s+([^\s]+)", initial_text, re.IGNORECASE)
    if len(constants) != 1:
        raise ValueError("native serialized initial head is not one uniform constant")
    initial = float(constants[0])
    sto = (directory / "ht_gw.sto").read_text()
    ss = re.findall(r"\bss\s+CONSTANT\s+([^\s]+)", sto, re.IGNORECASE)
    if "STORAGECOEFFICIENT" not in sto.upper() or len(ss) != 1:
        raise ValueError("native serialized storage is not uniform confined storativity")
    storativity = float(ss[0])
    river_rows = []
    for line in (directory / "ht_gw.riv").read_text().splitlines():
        fields = line.split()
        if len(fields) == 6 and all(f.isdigit() for f in fields[:3]):
            river_rows.append([float(f) for f in fields[3:]])
    values = np.asarray(river_rows)
    stage = forcing["sw_stage_m"].to_numpy(dtype=float)
    expected_stage = np.array([float(f"{s:.8E}") for s in stage])
    expected_bed = np.array([float(f"{s-static['river_bottom_offset_m']:.8E}") for s in stage])
    expected_c = float(f"{static['river_conductance_m2_per_day']:.8E}")
    if (initial != float(f"{static['aquifer_initial_head_m']:.8f}")
            or storativity != float(f"{static['aquifer_storage_coefficient']:.8f}")
            or values.shape != (len(stage), 3) or not np.isfinite(values).all()
            or not np.array_equal(values[:, 0], expected_stage)
            or not np.array_equal(values[:, 2], expected_bed)
            or not np.all(values[:, 1] == expected_c)):
        raise ValueError("native serialized river or initial heads differ from prescribed inputs")
    return {"native_initial_head": initial, "native_stage": values[:, 0],
            "native_bed": values[:, 2], "native_storativity": storativity}


def read_native(directory: Path):
    import flopy
    headfile = flopy.utils.HeadFile(directory / "ht_gw.hds", precision="double")
    budget = flopy.utils.CellBudgetFile(directory / "ht_gw.cbc", precision="double")
    try:
        heads = headfile.get_alldata()
        sy_records = budget.get_data(text="STO-SY")
        if len(sy_records) != len(heads):
            raise ValueError("native STO-SY records are missing or misaligned")
        sy = np.stack(sy_records)
        times = headfile.get_times()
        if not np.allclose(times, np.arange(1, len(heads) + 1), rtol=0.0, atol=1.0e-10):
            raise ValueError("native head records are not daily interval-end values")
        return heads, sy
    finally:
        headfile.close()
        budget.close()


def native_case(seed, variant, profile, workdir, standard_workdir, runners):
    probe = registry.find_probe(PROBE_ID)
    model = registry.find_model("modflow6")
    case = build_case(probe, seed, variant)
    root = workdir / str(seed) / variant
    captured = root / f"captured-{profile}"
    captured_input = captured / "input"
    captured_input.mkdir(parents=True, exist_ok=True)
    (captured_input / "capture_modflow.py").write_text(CAPTURE)
    audit_model = replace(model, entrypoint=(
        "python3", "/io/input/capture_modflow.py", "--profile", profile,
    ))
    run = runners["capture"].run(audit_model, probe, case, captured)
    output = captured / "output" / "result.csv"
    matched = False
    if profile == "standard":
        if standard_workdir:
            standard = standard_workdir / f"modflow6__{probe.slug}__{seed}__{variant}"
        else:
            standard = root / "standard"
        if not (standard / "output/result.csv").exists():
            runners["standard"].run(model, probe, case, standard)
        matched = compare_standard_output(standard / "output/result.csv", output)
    native_directory = captured.with_name(captured.name + "-native")
    heads, sy = read_native(native_directory)
    serialized = read_native_inputs(native_directory, case.forcing, case.static)
    metrics = native_metrics(case.forcing, run.table, heads, sy, case.static, **serialized)
    for item in probe.criteria:
        if item.name != "aquifer_equilibrium_bounds":
            result = CRITERIA[item.name](run, probe, item.params)
            if not result.passed:
                raise ValueError(f"native {profile}/{variant}: {result.message}")
    static = case.static
    b = 1000 * static["aquifer_storage_coefficient"]
    target = float(case.forcing["sw_stage_m"].iloc[case.spinup_steps])
    bound = b * (target - static["aquifer_initial_head_m"])
    initial = float(run.table["gw"].iloc[case.spinup_steps - 1])
    trace = (run.table["gw"].iloc[case.spinup_steps:].to_numpy(dtype=float) - initial).tolist()
    return {
        "seed": seed, "variant": variant, "profile": profile,
        "standard_output_byte_identical": matched if profile == "standard" else None,
        "result_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "wall_seconds": run.wall_seconds, "head_change_m": target - static["aquifer_initial_head_m"],
        "signed_storage_bound_mm": bound, "storage_change_mm": trace,
        **metrics,
    }


def reference_sweep(seeds):
    """Check all declared criteria and isolated faults on predetermined seeds."""
    probe = registry.find_probe(PROBE_ID)
    path = ROOT / "models/reference_confined_aquifer/ht_adapter.py"
    spec = importlib.util.spec_from_file_location("confined_reference_audit", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    controls = {
        "reference_confined_aquifer": {},
        "reference_confined_aquifer_implicit": {"method": "backward_euler"},
        "reference_amplified_river_head": {"head_gain": 1.4},
        "reference_inflated_aquifer_storage": {"storage_gain": 1.4},
    }
    models = {name: registry.find_model(name) for name in controls}
    results = {name: {"seeds": 0, "minimum_bound_ratio": None, "maximum_bound_ratio": 0.0}
               for name in controls}
    for seed in seeds:
        cases = {v: build_case(probe, seed, v) for v in probe.variants}
        for name, kwargs in controls.items():
            model = models[name]
            runs = {
                v: RunResult(case, pd.DataFrame(module.simulate(
                    case.forcing.to_dict("records"), case.static, **kwargs,
                )), {}, 0.0, model=model)
                for v, case in cases.items()
            }
            failures = set()
            headline = None
            for item in probe.criteria:
                if item.name == "aquifer_equilibrium_bounds":
                    result = CRITERIA[item.name](runs, probe, item.params)
                    headline = result.value
                    if not result.passed:
                        failures.add(item.name)
                else:
                    for run in runs.values():
                        if not CRITERIA[item.name](run, probe, item.params).passed:
                            failures.add(item.name)
            expected = set() if name in probe.must_pass else {"aquifer_equilibrium_bounds"}
            if failures != expected:
                raise ValueError(f"reference sweep {seed}/{name}: {failures}, expected {expected}")
            entry = results[name]
            entry["seeds"] += 1
            entry["minimum_bound_ratio"] = min(headline, entry["minimum_bound_ratio"]) if entry["minimum_bound_ratio"] is not None else headline
            entry["maximum_bound_ratio"] = max(headline, entry["maximum_bound_ratio"])
    return results


def precision_budget(params):
    """Conservative paired precision budget over the generator's fixed domain.

    STO constants have eight decimal places; RIV heads at these elevations
    have at most 5e-7 m rounding per run. Two runs enter every paired bound.
    The native audit separately enforces the head and storage-mapping limits.
    This is conditional on that writer, domain and those audit checks.
    """
    sigma_min, head_change_min, bound_min, bound_max = 0.002, 0.5, 1.0, 15.0
    relative = float(params["rel_tol"])
    absolute = float(params["abs_tol_mm"])
    if not np.isfinite((relative, absolute)).all() or relative < 0 or absolute <= 0:
        raise ValueError("precision allowances must be finite and nonnegative")
    sigma_roundoff = 5.0e-9 / sigma_min
    input_budget = sigma_roundoff + (1 + sigma_roundoff) * 1.0e-6 / head_change_min
    head_budget = (1 + sigma_roundoff) * 2 * NATIVE_HEAD_AUDIT_LIMIT_M / head_change_min
    mapping_budget = 2 * NATIVE_MAPPING_AUDIT_LIMIT_MM
    combined = input_budget + head_budget + mapping_budget / bound_min
    minimum_fraction = max(relative, absolute / bound_max)
    return {
        "minimum_storativity": sigma_min, "minimum_head_change_m": head_change_min,
        "storage_bound_range_mm": [bound_min, bound_max],
        "per_run_stage_roundoff_m": 5.0e-7, "per_run_storativity_roundoff": 5.0e-9,
        "per_run_head_audit_limit_m": NATIVE_HEAD_AUDIT_LIMIT_M,
        "per_run_mapping_audit_limit_mm": NATIVE_MAPPING_AUDIT_LIMIT_MM,
        "input_rounding_relative_budget": input_budget,
        "head_audit_relative_budget": head_budget,
        "paired_mapping_absolute_budget_mm": mapping_budget,
        "combined_relative_budget": combined,
        "minimum_allowance_to_budget_ratio": minimum_fraction / combined,
    }


def score_native_groups(records, workdir):
    """Re-read actual native outputs and judge each full paired experiment."""
    probe, model = registry.find_probe(PROBE_ID), registry.find_model("modflow6")
    item = next(c for c in probe.criteria if c.name == "aquifer_equilibrium_bounds")
    results = []
    for seed, profile in sorted({(r["seed"], r["profile"]) for r in records}):
        variants = {r["variant"] for r in records if (r["seed"], r["profile"]) == (seed, profile)}
        if variants != set(probe.variants):
            raise ValueError(f"native {seed}/{profile}: incomplete paired experiment")
        runs = {
            v: read_result(workdir / str(seed) / v / f"captured-{profile}",
                           build_case(probe, seed, v), probe, 0.0, model=model)
            for v in probe.variants
        }
        result = CRITERIA[item.name](runs, probe, item.params)
        if not result.passed:
            raise ValueError(f"native {seed}/{profile}: {result.message}")
        results.append({
            "seed": seed, "profile": profile, "passed": result.passed, "ratio": result.value,
            "maximum_envelope_excess_mm": max(
                d["max_excess_mm"] for d in result.diagnostics["variants"].values()),
            "diagnostics": result.diagnostics,
        })
    return results


def fault_magnitude_study(seeds, gains=(1.0, 1.00005, 1.0005, 1.001, 1.01, 1.4)):
    """Compare policies on conservative faults, without fitting a trajectory."""
    probe = registry.find_probe(PROBE_ID)
    item = next(c for c in probe.criteria if c.name == "aquifer_equilibrium_bounds")
    policies = {
        "previous": dict(item.params, rel_tol=0.05, abs_tol_mm=0.001,
                         initial_abs_tol_mm=0.001, neutral_abs_tol_mm=0.001),
        "calibrated": dict(item.params),
    }
    path = ROOT / "models/reference_confined_aquifer/ht_adapter.py"
    spec = importlib.util.spec_from_file_location("confined_fault_study", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    model = registry.find_model("reference_confined_aquifer")
    controls = {kind: [dict(gain=g, ordinary_criteria_passed=0, **{
        policy: {"passed": 0, "failed": 0, "minimum_ratio": None, "maximum_ratio": 0.0}
        for policy in policies
    }) for g in gains] for kind in ("head", "storage")}
    for seed in seeds:
        cases = {v: build_case(probe, seed, v) for v in probe.variants}
        for kind, entries in controls.items():
            for entry in entries:
                runs = {v: RunResult(case, pd.DataFrame(module.simulate(
                    case.forcing.to_dict("records"), case.static,
                    **{f"{kind}_gain": entry["gain"]},
                )), {}, 0.0, model=model) for v, case in cases.items()}
                for criterion in probe.criteria:
                    if criterion.name != item.name and not all(
                        CRITERIA[criterion.name](run, probe, criterion.params).passed
                        for run in runs.values()
                    ):
                        raise ValueError(f"fault study {seed}/{kind}: accounting or sign failure")
                entry["ordinary_criteria_passed"] += 1
                for policy, params in policies.items():
                    result = CRITERIA[item.name](runs, probe, params)
                    if result.value is None or not np.isfinite(result.value):
                        raise ValueError("fault study produced an invalid criterion ratio")
                    stats = entry[policy]
                    stats["passed" if result.passed else "failed"] += 1
                    stats["minimum_ratio"] = min(result.value, stats["minimum_ratio"]) if stats["minimum_ratio"] is not None else result.value
                    stats["maximum_ratio"] = max(result.value, stats["maximum_ratio"])
    return {"seeds": list(seeds), "policies": policies, "controls": controls}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--standard-workdir", type=Path)
    parser.add_argument("--report", type=Path, default=ROOT / "results/aquifer-equilibrium-bounds/validation.json")
    parser.add_argument("--audit-seeds", type=int, nargs="*", default=[7, 42, 13, 101, 2026])
    parser.add_argument("--skip-sensitivity", action="store_true")
    args = parser.parse_args(argv)
    probe = registry.find_probe(PROBE_ID)
    model = registry.find_model("modflow6")
    seeds = list(dict.fromkeys([*gate_seeds(PROBE_ID, probe.n_seeds), *args.audit_seeds]))
    runners = {"standard": DockerRunner(), "capture": NativeCaptureRunner()}
    records = []
    for seed in seeds:
        for variant in probe.variants:
            record = native_case(seed, variant, "standard", args.workdir, args.standard_workdir, runners)
            records.append(record)
            print(f"native {seed}/{variant}: confined, connected, STO-SY=0; standard output unchanged", flush=True)
    if not args.skip_sensitivity:
        for profile in ("tight_solver", "grid20"):
            for variant in probe.variants:
                record = native_case(seeds[0], variant, profile, args.workdir, None, runners)
                records.append(record)
                print(f"sensitivity {profile}/{variant}: native assumptions and budget verified", flush=True)
    reference_seeds = list(dict.fromkeys([*range(100), *gate_seeds(PROBE_ID, probe.n_seeds)]))
    refs = reference_sweep(reference_seeds)
    criterion_params = next(c.params for c in probe.criteria if c.name == "aquifer_equilibrium_bounds")
    native_results = score_native_groups(records, args.workdir)
    study = fault_magnitude_study(reference_seeds)
    image = json.loads(subprocess.check_output(["docker", "image", "inspect", image_tag(model)], text=True))[0]
    report = {
        "probe": PROBE_ID, "model": model.name, "model_version": model.version,
        "modflow_version": "6.7.0", "flopy_version": "3.11.0",
        "docker_image_id": image["Id"], "docker_architecture": image["Architecture"],
        "criterion_params": criterion_params,
        "precision_budget": precision_budget(criterion_params),
        "native_records": records, "native_probe_results": native_results,
        "reference_sweep": refs, "fault_magnitude_study": study,
        "audit_notes": [
            "Native heads and STO-SY independently check confinement, river connection and the storage law.",
            "The standard diagnostic observer must preserve result.csv bytes.",
            "Serialized IC/RIV inputs are checked against the pinned writer's rounding; both native and requested head-envelope errors are retained.",
            "Sensitivity changes only IMS head tolerance or the spatial grid; daily reporting is preserved.",
            "A grid-sensitivity pass is a bounds check, not a claim of identical transient hydrographs or grid convergence.",
            "Frozen output can pass this necessary storage bound and is explicitly reported as unexercised.",
        ],
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"wrote {args.report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

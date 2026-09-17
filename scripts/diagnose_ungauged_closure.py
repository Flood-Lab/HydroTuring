"""Annual budgets and sampled attributes; diagnostics never replace ht verdicts."""
from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

from hydroturing.criteria.closure import closure
from hydroturing.harness import build_case, load_generator
from hydroturing.protocol import read_result
from hydroturing.registry import find_probe


def diagnose(run, probe):
    if run.case.timestep != "PT1D":
        raise ValueError("this diagnostic is defined for daily cases")
    start = run.case.spinup_steps
    if start < 1 or (len(run.table) - start) % 365:
        raise ValueError("annual blocks need a preceding state and complete 365-day years")
    params = dict(next(c.params for c in probe.criteria if c.name == "closure"))
    full = closure(run, probe, params)
    years = []
    for lo in range(start, len(run.table), 365):
        hi = lo + 365
        case = replace(run.case, forcing=run.case.forcing.iloc[:hi], spinup_steps=lo)
        result = closure(replace(run, case=case, table=run.table.iloc[:hi]), probe, params)
        years.append(dict(start=str(run.table.time.iloc[lo]), end=str(run.table.time.iloc[hi-1]),
                          residual_mm=result.diagnostics["cumulative_residual"],
                          precipitation_mm=result.diagnostics["denominator_total"],
                          storage_change_mm=result.diagnostics["storage_change"],
                          relative_residual=result.value,
                          exceeds_full_window_tolerance=not result.passed))
    generator = load_generator(probe)
    seed = run.case.seed
    return dict(diagnostic_only=True, seed=seed,
                attribute_category=generator.ATTRIBUTES[seed % len(generator.ATTRIBUTES)][0],
                climate=generator.climate_for_seed(seed)[0], static=run.case.static,
                full_window_closure=dict(status=full.status, relative_residual=full.value),
                years=years,
                interpretation="Annual residuals do not change the official ht verdict. "
                               "Full-window closure does not imply annual closure.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", type=Path, help="one completed ht case directory")
    parser.add_argument("--seed", type=int, required=True, help="generator seed from the ht report")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    probe = find_probe("mass/ungauged-basin-closure")
    case = build_case(probe, args.seed)
    root = args.case_dir.resolve()
    # Reject a different seed/version/window instead of silently pairing its
    # model output with freshly generated but different forcing or attributes.
    if (root / "input/forcing.csv").read_text() != case.forcing.to_csv(index=False):
        raise ValueError("stored forcing differs from the generated case")
    if json.loads((root / "input/static.json").read_text()) != case.static:
        raise ValueError("stored attributes differ from the generated case")
    request = json.loads((root / "request.json").read_text())
    if request["n_steps"] != case.n_steps or request["timestep"] != case.timestep:
        raise ValueError("stored request has a different length or timestep")
    run = read_result(root, case, probe, wall_seconds=0.0)
    args.output.write_text(json.dumps(diagnose(run, probe), indent=2) + "\n")


if __name__ == "__main__":
    main()

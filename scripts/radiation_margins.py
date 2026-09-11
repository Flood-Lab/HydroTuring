#!/usr/bin/env python3
"""Reproduce the radiation probe's margin table through real adapters.

Run five gate seeds, twenty additional seeds, seed 36 near the emissivity
ceiling, and gate seeds at both emissivity bounds. Check output parity and
report per-case maximum/minimum residual-to-tolerance ratios and day/night
violation counts. Each negative control's maximum ratio must reach 1.3;
the criterion's threshold remains 1. This extended sweep is separate from
the fast regression tests.

    python3 scripts/radiation_margins.py [output.json]
"""

from __future__ import annotations

import json
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.criteria.radiation import STEFAN_BOLTZMANN
from hydroturing.harness import build_case
from hydroturing.runner import get_runner
from hydroturing.seeds import gate_seeds

PROBE = registry.find_probe("energy/radiation-consistency")
POSITIVE = "reference_radiative"
NEGATIVES = ("reference_air_emitter", "reference_no_reflection")
CRITERION = get("radiative_identity")
PARAMS = dict(PROBE.criteria[0].params)
KELVIN = 273.15
MARGIN = 1.3


def run(name, case):
    # Read back what the model declares, as verify-adapter does, so that the
    # coupled reference can be run on a probe that asks for more than it has.
    model = registry.find_model(name)
    asks = replace(
        PROBE, requires_fluxes=model.emits_fluxes, requires_states=model.emits_states,
        requires_diagnostics=model.emits_diagnostics,
    )
    with tempfile.TemporaryDirectory(prefix="ht-margins-") as workdir:
        return get_runner(model).run(model, asks, case, Path(workdir))


def score(result):
    r = CRITERION(result, PROBE, PARAMS)
    case = result.case
    forcing = case.forcing.iloc[case.spinup_steps:].reset_index(drop=True)
    hours = pd.to_datetime(forcing["time"]).dt.hour
    day = ((hours >= 6) & (hours < 18)).to_numpy()
    # The criterion reports the worst step; the split needs every step, so
    # the residual is recomputed here with the criterion's own constants.
    scored = result.table.iloc[case.spinup_steps:].reset_index(drop=True)
    eps = case.static[PARAMS["emissivity"]]
    sigma = PARAMS.get("sigma", STEFAN_BOLTZMANN)
    expected = eps * sigma * scored["ts"] ** 4 + (1.0 - eps) * forcing["rlds"]
    slack = ((scored["rlus"] - expected).abs()
             / np.maximum(PARAMS["rel_tol"] * scored["rlus"].abs(), PARAMS["abs_floor"])).to_numpy()
    violating = slack > 1.0
    return {
        "status": r.status,
        "worst_slack": r.diagnostics.get("worst_slack", float("nan")),
        # Distinguish the weakest single-step margin from the worst step.
        "min_slack": float(slack.min()),
        "violating": int(violating.sum()),
        "violating_day": int((violating & day).sum()),
        "violating_night": int((violating & ~day).sum()),
    }


def sweep(label, seeds, eps=None):
    rows = []
    for seed in seeds:
        case = build_case(PROBE, seed)
        if eps is not None:
            case.static["eps"] = eps
        tables = {}
        row = {"label": label, "seed": seed, "eps": case.static["eps"]}
        for name in (POSITIVE, *NEGATIVES):
            result = run(name, case)
            tables[name] = result.table
            row[name] = score(result)
        coupled = run("reference_coupled", case).table

        positive = tables[POSITIVE]
        row["positive_matches_coupled"] = all(
            np.array_equal(positive[c].to_numpy(), coupled[c].to_numpy()) for c in coupled.columns
        )
        for name in NEGATIVES:
            same = [
                c for c in positive.columns
                if np.array_equal(positive[c].to_numpy(), tables[name][c].to_numpy())
            ]
            row[f"{name}_differs_only_in_rlus"] = set(positive.columns) - set(same) == {"rlus"}
        ts = positive["ts"].to_numpy()
        row["ts_min_k"], row["ts_max_k"] = float(ts.min()), float(ts.max())
        row["ts_frozen_steps"] = int((ts <= KELVIN).sum())
        row["ok"] = (
            row["positive_matches_coupled"]
            and all(row[f"{n}_differs_only_in_rlus"] for n in NEGATIVES)
            and row["ts_frozen_steps"] == 0
            and row[POSITIVE]["status"] == "pass"
            and all(row[n]["status"] == "fail" and row[n]["worst_slack"] >= MARGIN for n in NEGATIVES)
        )
        rows.append(row)
    return rows


def main() -> int:
    gate = gate_seeds(PROBE.id, PROBE.n_seeds)
    groups = [
        ("gate", gate, None),
        ("additional", list(range(20)), None),
        # Seed 36 draws the highest emissivity among seeds 0-49, 0.9891: a
        # drawn case near the ceiling. Forced 0.99 also checks the endpoint.
        ("hi-eps", [36], None),
        ("eps=0.99", gate, 0.99),
        ("eps=0.95", gate, 0.95),
    ]
    rows = [row for label, seeds, eps in groups for row in sweep(label, seeds, eps)]
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(json.dumps(rows, indent=1))

    print(f"{'group':10} {'seed':>10} {'eps':>6} | {'positive':>9} | "
          f"{'air_emitter':>11} {'viol':>4} {'day/night':>9} | "
          f"{'no_reflection':>13} {'viol':>4} {'day/night':>9} |")
    for r in rows:
        p, a, n = r[POSITIVE], r[NEGATIVES[0]], r[NEGATIVES[1]]
        print(f"{r['label']:10} {r['seed']:>10} {r['eps']:>6.4f} | {p['worst_slack']:>9.1e} | "
              f"{a['worst_slack']:>11.2f} {a['violating']:>4} {a['violating_day']:>4}/{a['violating_night']:<4} | "
              f"{n['worst_slack']:>13.2f} {n['violating']:>4} {n['violating_day']:>4}/{n['violating_night']:<4} | "
              f"{'ok' if r['ok'] else 'PROBLEM'}")

    for label, _, _ in groups:
        sub = [r for r in rows if r["label"] == label]
        print(f"\n{label}: {len(sub)} seed(s); positive worst step "
              f"{max(r[POSITIVE]['worst_slack'] for r in sub):.1e}; "
              f"min over seeds of the worst step: air_emitter "
              f"{min(r[NEGATIVES[0]]['worst_slack'] for r in sub):.2f}, no_reflection "
              f"{min(r[NEGATIVES[1]]['worst_slack'] for r in sub):.2f}; "
              f"smallest single step: air_emitter "
              f"{min(r[NEGATIVES[0]]['min_slack'] for r in sub):.2f}, no_reflection "
              f"{min(r[NEGATIVES[1]]['min_slack'] for r in sub):.2f}; "
              f"skin {min(r['ts_min_k'] for r in sub):.1f} to {max(r['ts_max_k'] for r in sub):.1f} K")
    return 0 if all(r["ok"] for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())

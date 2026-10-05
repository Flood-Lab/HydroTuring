from __future__ import annotations
import numpy as np
from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.criteria.response import pick
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

@criterion("datum_flux_invariance", paired=True)
def datum_flux_invariance(runs: dict[str, RunResult], probe: ProbeSpec, params: dict) -> CriterionResult:
    """Compare exchange trajectories after a pure vertical-datum translation.

    The absolute difference is summed before any sign cancellation.  The
    allowance is a relative share of the control's gross exchange, with a
    fixed floor for output precision.  Gross activity is checked separately:
    a model that reports zero (or a constant bookkeeping flux) must not pass
    merely because both datum variants agree.  ``exchange_directions`` is
    paired with this criterion in the probe to require a real two-way stage
    response in every variant.
    """
    control_name = str(params.get("control", "control"))
    transformed = params.get("transformed", "transformed")
    if isinstance(transformed, str):
        transformed = [transformed]
    transformed = [str(x) for x in transformed]
    if control_name not in runs:
        raise ValueError(f"datum_flux_invariance needs control variant {control_name!r}")
    control = make_window(pick(runs, params, "control", control_name), probe)
    exchange = str(params.get("exchange", "gw_sw_exchange"))
    if exchange not in control.table.columns:
        raise ValueError(f"datum_flux_invariance needs '{exchange}' in the model result")
    q0 = control.volume(control.table[exchange])
    if not np.isfinite(q0).all():
        return CriterionResult(name="datum_flux_invariance", status=FAIL, threshold=1.0, message="control exchange contains non-finite values")
    gross = float(np.abs(q0).sum())
    minimum = float(params.get("minimum_activity_mm", 0.1))
    activity_ok = gross >= minimum
    rtol = float(params.get("rtol", 1e-6))
    floor = float(params.get("abs_tol_mm", 1e-5))
    allowed = max(rtol * gross, floor)
    diagnostics = {"control_gross_exchange_mm": gross, "allowed_residual_mm": allowed, "variants": {}}
    worst_name, worst = None, 0.0
    all_ok = activity_ok
    for name in transformed:
        if name not in runs:
            raise ValueError(f"datum_flux_invariance needs variant {name!r}")
        other = make_window(runs[name], probe)
        if len(other.table) != len(control.table):
            raise ValueError("datum variants have different scored window lengths")
        if exchange not in other.table.columns:
            raise ValueError(f"datum_flux_invariance needs '{exchange}' in variant {name!r}")
        q = other.volume(other.table[exchange])
        if not np.isfinite(q).all():
            all_ok = False
            residual = float("inf")
        else:
            residual = float(np.abs(q - q0).sum())
            all_ok = all_ok and residual <= allowed
        diagnostics["variants"][name] = {"residual_mm": residual}
        if residual > worst:
            worst_name, worst = name, residual
    ok = all_ok
    message = (f"datum shift leaves {exchange} invariant (worst residual {worst:.3e} mm)" if ok else
               (f"datum shift changes {exchange} by {worst:.6g} mm in {worst_name!r} (limit {allowed:.6g} mm)" if activity_ok else
                f"exchange activity is only {gross:.6g} mm (minimum {minimum:g} mm)"))
    diagnostics["activity_ok"] = activity_ok
    return CriterionResult(name="datum_flux_invariance", status=PASS if ok else FAIL,
                           value=max(worst / max(allowed, np.finfo(float).tiny), minimum / max(gross, np.finfo(float).tiny) if not activity_ok else 0.0),
                           threshold=1.0, message=message, diagnostics=diagnostics)

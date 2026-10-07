"""Regression tests for momentum/wave-celerity-bounds."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get, is_paired
from hydroturing.harness import build_case, compatibility_issues, run_probe
from hydroturing.protocol import Case, RunResult
from hydroturing.seeds import gate_seeds
from hydroturing.spec import Criterion, ProbeSpec

ROOT = Path(__file__).resolve().parents[1]
STATES = ("low", "medium", "high")
STATE_Q = (8.0, 16.0, 32.0)
EVENT_COLUMNS = {
    "low": "_pulse_low",
    "medium": "_pulse_medium",
    "high": "_pulse_high",
}


def _probe() -> ProbeSpec:
    return ProbeSpec(
        id="momentum/wave-celerity-bounds",
        title="test",
        law="momentum",
        track="synthetic",
        version=1,
        authors=({"name": "test", "affiliation": "test"},),
        citation="",
        requires_fluxes=("dis",),
        requires_states=(),
        requires_forcing=("q_in",),
        requires_static=(
            "width_m", "cross_section_shape", "slope", "manning_n",
            "reach_length_m",
        ),
        generator="generate.py",
        n_seeds=1,
        timestep="PT1H",
        period_years=24.0 / 365.0,
        period_days=24.0,
        spinup_days=2,
        max_output_mb=2.0,
        max_runtime_s=300.0,
        variants=("short", "long"),
        criteria=(Criterion("wave_celerity_bounds", {}),),
        must_pass=("reference_saint_venant",),
        must_fail={"reference_fixed_celerity": "wave_celerity_bounds"},
        provenance="synthetic",
        path=ROOT / "probes" / "momentum" / "wave-celerity-bounds",
    )


def _params(**overrides) -> dict:
    params = {
        "states": list(STATES),
        "event_columns": dict(EVENT_COLUMNS),
        "relative_tolerance": 0.05,
        "base_flow_relative_tolerance": 0.01,
        "baseline_hours": 12,
        "response_hours": 72,
        "max_tail_fraction": 0.01,
        "timing_floor_hours": 0.05,
        "ordering_margin_fraction": 0.0,
    }
    params.update(overrides)
    return params


def _normal_depth(q: float, width: float, slope: float, n: float) -> float:
    def f(h):
        a = width * h
        r = a / (width + 2*h)
        return a * r ** (2/3) * slope ** 0.5 / n

    lo, hi = 1e-6, 1.0
    while f(hi) < q:
        hi *= 2
    for _ in range(80):
        mid = (lo + hi) / 2
        if f(mid) < q:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _c_kin(q: float, width=75.0, slope=0.0012, n=0.033) -> float:
    h = _normal_depth(q, width, slope, n)
    a = width * h
    p = width + 2*h
    return q * (5/(3*a) - 4/(3*width*p))


def _synthetic_runs(celerities: tuple[float, float, float] | None = None) -> dict[str, RunResult]:
    if celerities is None:
        celerities = tuple(_c_kin(q) for q in STATE_Q)

    probe = _probe()
    spin = 48
    block = 96
    n = spin + 3 * block
    time = pd.date_range("2001-01-01", periods=n, freq="h")
    forcing = pd.DataFrame({
        "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "pr": np.zeros(n),
        "tas": np.full(n, 15.0),
        "pet": np.zeros(n),
        **{column: np.zeros(n) for column in EVENT_COLUMNS.values()},
    })

    base_by_row = np.full(n, STATE_Q[0], dtype=float)
    pulse_centres: dict[str, float] = {}
    for i, (state, q) in enumerate(zip(STATES, STATE_Q)):
        start = spin + i * block
        stop = start + block
        base_by_row[start:stop] = q
        pulse_start = start + 24
        forcing.loc[pulse_start:pulse_start + 5, EVENT_COLUMNS[state]] = 0.05 * q
        pulse_centres[state] = pulse_start + 3.0

    q_in = base_by_row.copy()
    for state, q in zip(STATES, STATE_Q):
        start = int(round(pulse_centres[state] - 3.0))
        q_in[start:start + 6] += 0.05 * q
    forcing["q_in"] = q_in

    runs = {}
    for side, length in (("short", 4000.0), ("long", 20000.0)):
        dis = base_by_row.copy()
        centres = np.arange(n, dtype=float) + 0.5
        for state, q, celerity in zip(STATES, STATE_Q, celerities):
            travel_h = length / celerity / 3600.0
            response = np.exp(
                -0.5 * ((centres - (pulse_centres[state] + travel_h)) / 1.2) ** 2
            )
            dis += 0.2 * q * response

        case = Case(
            probe_id=f"{probe.id}@{side}",
            seed=1,
            forcing=forcing.copy(),
            static={
                "area_km2": 100.0,
                "width_m": 75.0,
                "cross_section_shape": "rectangular",
                "bed_elevation_m": 50.0,
                "slope": 0.0012,
                "manning_n": 0.033,
                "reach_length_m": length,
            },
            spinup_steps=spin,
            timestep="PT1H",
        )
        table = pd.DataFrame({"time": forcing["time"], "dis": dis})
        runs[side] = RunResult(case, table, {"status": "ok"}, 0.0)
    return runs



def _hayami_kernel(hours: np.ndarray, length_m: float, celerity: float, diffusivity: float) -> np.ndarray:
    """Discrete Hayami impulse response for a constant-parameter diffusion wave."""
    seconds = np.asarray(hours, dtype=float) * 3600.0
    kernel = np.zeros_like(seconds)
    positive = seconds > 0.0
    t = seconds[positive]
    kernel[positive] = (
        length_m
        / (2.0 * np.sqrt(np.pi * diffusivity * t**3))
        * np.exp(-((length_m - celerity * t) ** 2) / (4.0 * diffusivity * t))
    )
    total = float(kernel.sum())
    assert total > 0.0
    return kernel / total


def _hayami_synthetic_runs() -> dict[str, RunResult]:
    """Pair reaches with the right celerity but substantial diffusive broadening."""
    probe = _probe()
    spin = 48
    block = 96
    n_rows = spin + 3 * block
    time = pd.date_range("2001-01-01", periods=n_rows, freq="h")
    forcing = pd.DataFrame({
        "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "pr": np.zeros(n_rows),
        "tas": np.full(n_rows, 15.0),
        "pet": np.zeros(n_rows),
        **{column: np.zeros(n_rows) for column in EVENT_COLUMNS.values()},
    })

    base = np.full(n_rows, STATE_Q[0], dtype=float)
    pulse_starts = {}
    for i, (state, q) in enumerate(zip(STATES, STATE_Q)):
        start = spin + i * block
        stop = start + block
        base[start:stop] = q
        pulse_start = start + 24
        pulse_starts[state] = pulse_start
        forcing.loc[pulse_start:pulse_start + 5, EVENT_COLUMNS[state]] = 0.05 * q
    q_in = base.copy()
    for state, q in zip(STATES, STATE_Q):
        pulse_start = pulse_starts[state]
        q_in[pulse_start:pulse_start + 6] += 0.05 * q
    forcing["q_in"] = q_in

    width = 75.0
    slope = 0.0012
    manning_n = 0.033
    runs = {}
    for side, length in (("short", 4000.0), ("long", 20000.0)):
        discharge = base.copy()
        for state, q in zip(STATES, STATE_Q):
            celerity = _c_kin(q, width=width, slope=slope, n=manning_n)
            # Five times the low-inertia Hayami value Q/(2 B S) makes the
            # hydrograph visibly broader while leaving the first moment at L/c.
            diffusivity = 5.0 * q / (2.0 * width * slope)
            kernel = _hayami_kernel(
                np.arange(60, dtype=float) + 0.5,
                length,
                celerity,
                diffusivity,
            )
            routed_pulse = np.convolve(np.ones(6, dtype=float), kernel)
            start = pulse_starts[state]
            stop = min(start + len(routed_pulse), n_rows)
            discharge[start:stop] += 0.05 * q * routed_pulse[:stop - start]

        case = Case(
            probe_id=f"{probe.id}@{side}",
            seed=1,
            forcing=forcing.copy(),
            static={
                "area_km2": 100.0,
                "width_m": width,
                "cross_section_shape": "rectangular",
                "bed_elevation_m": 50.0,
                "slope": slope,
                "manning_n": manning_n,
                "reach_length_m": length,
            },
            spinup_steps=spin,
            timestep="PT1H",
        )
        table = pd.DataFrame({"time": forcing["time"], "dis": discharge})
        runs[side] = RunResult(case, table, {"status": "ok"}, 0.0)
    return runs


def test_centroid_pair_recovers_state_dependent_celerity_under_hayami_diffusion():
    result = get("wave_celerity_bounds")(
        _hayami_synthetic_runs(),
        _probe(),
        _params(),
    )
    assert result.passed, result.message
    for state, q in zip(STATES, STATE_Q):
        diag = result.diagnostics["states"][state]
        assert abs(diag["relative_residual"]) < 0.05
        expected_diffusivity = 5.0 * q / (2.0 * 75.0 * 0.0012)
        inferred_diffusivity = diag["paired_diffusivity_m2_s"]
        assert inferred_diffusivity is not None
        assert abs(inferred_diffusivity / expected_diffusivity - 1.0) < 0.20



def test_hayami_diffusion_biases_peak_lag_but_not_first_moment():
    """Peak timing is shape-sensitive; the first temporal moment identifies translation."""
    peak_errors = []
    for q in STATE_Q:
        celerity = _c_kin(q)
        diffusivity = 5.0 * q / (2.0 * 75.0 * 0.0012)
        peak_h = []
        centroid_h = []
        for length in (4000.0, 20000.0):
            kernel = _hayami_kernel(
                np.arange(60, dtype=float) + 0.5,
                length,
                celerity,
                diffusivity,
            )
            response = np.convolve(np.ones(6, dtype=float), kernel)
            centres = np.arange(len(response), dtype=float) + 0.5
            peak_h.append(float(centres[int(np.argmax(response))]))
            centroid_h.append(float(np.dot(centres, response) / response.sum()))

        centroid_c = 16000.0 / ((centroid_h[1] - centroid_h[0]) * 3600.0)
        assert abs(centroid_c / celerity - 1.0) < 0.05

        peak_dt_h = peak_h[1] - peak_h[0]
        assert peak_dt_h > 0.0
        peak_c = 16000.0 / (peak_dt_h * 3600.0)
        peak_errors.append(abs(peak_c / celerity - 1.0))

    # At least one state is badly biased if a peak index is used as the clock.
    assert max(peak_errors) > 0.20


def test_generated_manning_characteristic_speed_has_resolved_state_ordering():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    for seed in gate_seeds(probe.id, probe.n_seeds):
        case = build_case(probe, seed, "short")
        width = float(case.static["width_m"])
        slope = float(case.static["slope"])
        manning_n = float(case.static["manning_n"])
        celerities = [
            _c_kin(q, width=width, slope=slope, n=manning_n)
            for q in STATE_Q
        ]
        for left, right in zip(celerities, celerities[1:]):
            assert right > left
            assert right - left > 0.05 * max(left, right)


def test_criterion_is_registered_as_paired():
    assert is_paired("wave_celerity_bounds")


def test_non_rectangular_geometry_is_rejected():
    runs = _synthetic_runs()
    for side in ("short", "long"):
        run = runs[side]
        static = dict(run.case.static)
        static["cross_section_shape"] = "trapezoidal"
        runs[side] = RunResult(
            Case(
                probe_id=run.case.probe_id,
                seed=run.case.seed,
                forcing=run.case.forcing,
                static=static,
                spinup_steps=run.case.spinup_steps,
                timestep=run.case.timestep,
            ),
            run.table,
            run.meta,
            run.wall_seconds,
        )
    with pytest.raises(ValueError, match="rectangular"):
        get("wave_celerity_bounds")(runs, _probe(), _params())


def test_probe_uses_direct_river_inflow_contract():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    assert probe.requires_forcing == ("q_in",)
    seed = gate_seeds(probe.id, 1)[0]
    short = build_case(probe, seed, "short")
    long = build_case(probe, seed, "long")
    assert np.all(short.forcing["pr"].to_numpy(float) == 0.0)
    pd.testing.assert_series_equal(short.forcing["q_in"], long.forcing["q_in"])
    assert set(short.forcing["q_in"].unique()) >= set(STATE_Q)


def test_wave_probe_declared_inputs_are_real_model_opt_ins():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    case = build_case(probe, gate_seeds(probe.id, 1)[0], "short")
    for name in ("reference_saint_venant", "reference_fixed_celerity", "wflow_sbm"):
        model = registry.find_model(name)
        assert "q_in" in model.uses_forcing
        for key in ("width_m", "cross_section_shape", "slope", "manning_n", "reach_length_m"):
            assert key in model.needs_static + model.uses_static
        assert compatibility_issues(model, probe, case) == []


def test_registered_probe_and_references_are_compatible():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    case = build_case(probe, gate_seeds(probe.id, 1)[0], "short")
    for name in ("reference_saint_venant", "reference_fixed_celerity"):
        assert compatibility_issues(registry.find_model(name), probe, case) == []



def test_generator_pairs_change_only_reach_length():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    seed = gate_seeds(probe.id, 1)[0]
    short = build_case(probe, seed, "short")
    long = build_case(probe, seed, "long")
    pd.testing.assert_frame_equal(short.forcing, long.forcing, check_exact=True)
    s, l = dict(short.static), dict(long.static)
    assert l.pop("reach_length_m") > s.pop("reach_length_m")
    assert s == l




def test_pair_rejects_any_static_change_besides_reach_length():
    runs = _synthetic_runs()
    long = runs["long"]
    static = dict(long.case.static)
    static["width_m"] = float(static["width_m"]) * 1.01
    runs["long"] = RunResult(
        Case(
            probe_id=long.case.probe_id,
            seed=long.case.seed,
            forcing=long.case.forcing,
            static=static,
            spinup_steps=long.case.spinup_steps,
            timestep=long.case.timestep,
        ),
        long.table,
        long.meta,
        long.wall_seconds,
    )
    with pytest.raises(ValueError, match="differ in static 'width_m'"):
        get("wave_celerity_bounds")(runs, _probe(), _params())


def test_pair_rejects_a_shifted_base_state():
    runs = _synthetic_runs()
    long = runs["long"]
    table = long.table.copy()
    table["dis"] = table["dis"] * 1.02
    runs["long"] = RunResult(long.case, table, long.meta, long.wall_seconds)
    result = get("wave_celerity_bounds")(runs, _probe(), _params())
    assert not result.passed
    assert "does not match prescribed q_in" in result.message


def test_model_cannot_choose_its_own_base_q_to_make_theory_self_consistent():
    # Under the old implementation, scaling the model output to half the
    # prescribed steady flow and routing at c_kin(Q/2) could make both the
    # measured speed and the theory target agree, because c_kin was evaluated
    # at the model's own base discharge.  The operating point is now q_in, so
    # that self-consistent but non-conservative construction must fail.
    wrong_c = tuple(_c_kin(0.5 * q) for q in STATE_Q)
    runs = _synthetic_runs(wrong_c)
    for side in ("short", "long"):
        run = runs[side]
        table = run.table.copy()
        table["dis"] = 0.5 * table["dis"]
        runs[side] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get("wave_celerity_bounds")(runs, _probe(), _params())
    assert not result.passed
    assert "prescribed q_in" in result.message


def test_generator_has_three_isolated_pulses_and_response_tail():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    case = build_case(probe, gate_seeds(probe.id, 1)[0], "short")
    scored = case.after_spinup(case.forcing)
    starts = []
    for state in STATES:
        marker = EVENT_COLUMNS[state]
        indices = np.flatnonzero(scored[marker].to_numpy(float) > 0)
        assert len(indices) == 6
        assert np.all(np.diff(indices) == 1)
        starts.append(int(indices[0]))
        assert int(indices[-1]) + 72 < len(scored)
    assert np.diff(starts).tolist() == [8 * 24, 8 * 24]



def test_generator_keeps_the_linearization_subcritical_and_small():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    gravity = 9.80665
    for seed in gate_seeds(probe.id, probe.n_seeds):
        case = build_case(probe, seed, "short")
        width = float(case.static["width_m"])
        slope = float(case.static["slope"])
        manning_n = float(case.static["manning_n"])

        # The highest base state is the worst Froude number in this design.
        depth = _normal_depth(STATE_Q[-1], width, slope, manning_n)
        velocity = STATE_Q[-1] / (width * depth)
        froude = velocity / np.sqrt(gravity * depth)
        assert froude < 0.5

        scored = case.after_spinup(case.forcing)
        for state, q in zip(STATES, STATE_Q):
            marker = EVENT_COLUMNS[state]
            pulse_rows = scored[marker].to_numpy(float) > 0.0
            assert pulse_rows.sum() == 6
            pulse_q = scored.loc[pulse_rows, marker].to_numpy(float)
            np.testing.assert_allclose(
                pulse_q / q,
                0.05,
                rtol=0.0,
                atol=1.0e-12,
            )
            np.testing.assert_allclose(
                scored.loc[pulse_rows, "q_in"].to_numpy(float),
                1.05 * q,
                rtol=0.0,
                atol=1.0e-12,
            )


def test_generator_stays_inside_wide_channel_allowance():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    for seed in gate_seeds(probe.id, probe.n_seeds):
        case = build_case(probe, seed, "short")
        width = float(case.static["width_m"])
        slope = float(case.static["slope"])
        n = float(case.static["manning_n"])
        q = 32.0
        exact = _c_kin(q, width=width, slope=slope, n=n)
        wide_depth = (q * n / (width * slope ** 0.5)) ** (3.0 / 5.0)
        wide = (5.0 / 3.0) * q / (width * wide_depth)
        assert abs(wide - exact) / exact < 0.05


def test_analytic_measurement_passes_when_celerity_matches_manning():
    result = get("wave_celerity_bounds")(
        _synthetic_runs(),
        _probe(),
        _params(),
    )
    assert result.passed, result.message


def _muskingum_synthetic_runs(x: float = 0.3) -> dict[str, RunResult]:
    """Route each six-hour pulse through one Muskingum subreach."""
    runs = _synthetic_runs()
    for side, length in (("short", 4000.0), ("long", 20000.0)):
        run = runs[side]
        forcing = run.case.forcing
        discharge = np.empty(len(forcing), dtype=float)
        # Rebuild the three state blocks from their prescribed steady states.
        discharge[:] = STATE_Q[0]
        discharge[48 + 96:48 + 2 * 96] = STATE_Q[1]
        discharge[48 + 2 * 96:] = STATE_Q[2]

        for i, (state, q) in enumerate(zip(STATES, STATE_Q)):
            block_start = 48 + i * 96
            block_stop = block_start + 96
            pulse_start = block_start + 24
            inflow = np.zeros(96, dtype=float)
            local_start = pulse_start - block_start
            inflow[local_start:local_start + 6] = 0.05 * q

            k_hours = length / _c_kin(q) / 3600.0
            denominator = 2.0 * k_hours * (1.0 - x) + 1.0
            c0 = (1.0 - 2.0 * k_hours * x) / denominator
            c1 = (1.0 + 2.0 * k_hours * x) / denominator
            c2 = (2.0 * k_hours * (1.0 - x) - 1.0) / denominator
            routed = np.zeros(96, dtype=float)
            for j in range(95):
                routed[j + 1] = (
                    c0 * inflow[j + 1] + c1 * inflow[j] + c2 * routed[j]
                )
            discharge[block_start:block_stop] += routed

        runs[side] = RunResult(
            run.case,
            pd.DataFrame({"time": forcing["time"], "dis": discharge}),
            run.meta,
            run.wall_seconds,
        )
    return runs


def test_signed_first_moment_accepts_dipping_muskingum_response():
    # With X=0.3 and one coarse subreach the long-reach C0 coefficient is
    # negative, so the response initially dips below base.  Its signed first
    # moment is nevertheless the correct K=L/c and must not be clipped away.
    runs = _muskingum_synthetic_runs(x=0.3)
    assert any(
        np.any(
            runs["long"].table["dis"].to_numpy(float)[48 + i * 96:48 + (i + 1) * 96]
            < q - 1.0e-12
        )
        for i, q in enumerate(STATE_Q)
    )
    result = get("wave_celerity_bounds")(runs, _probe(), _params())
    assert result.passed, result.message
    for state in STATES:
        assert abs(result.diagnostics["states"][state]["relative_residual"]) < 0.01


def test_truncated_slow_tail_is_rejected_before_centroid_timing():
    runs = _synthetic_runs()
    run = runs["long"]
    table = run.table.copy()
    pulse_start = 48 + 24
    hours = np.arange(80, dtype=float)
    table.loc[pulse_start:pulse_start + 79, "dis"] += 1.6 * np.exp(-hours / 80.0)
    runs["long"] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get("wave_celerity_bounds")(runs, _probe(), _params())
    assert not result.passed
    assert "response is not contained in the response window" in result.message
    assert result.diagnostics["response_tail_fraction"] > 0.01


def test_tail_volume_rejects_response_even_when_last_sample_returns_to_base():
    runs = _synthetic_runs()
    run = runs["long"]
    table = run.table.copy()
    pulse_start = 48 + 24
    response_stop = pulse_start + 72

    # A pointwise edge check would accept this because the final sample is
    # exactly back at base, even though the preceding five hours still carry a
    # material amount of response volume.
    table.loc[response_stop - 6:response_stop - 2, "dis"] += 0.2
    table.loc[response_stop - 1, "dis"] = STATE_Q[0]
    runs["long"] = RunResult(run.case, table, run.meta, run.wall_seconds)

    result = get("wave_celerity_bounds")(runs, _probe(), _params())
    assert not result.passed
    assert "response is not contained in the response window" in result.message
    assert result.diagnostics["response_tail_edge_m3s"] == pytest.approx(0.0)
    assert result.diagnostics["response_tail_fraction"] > 0.01


def test_unresolved_timing_diagnostics_are_strict_json():
    runs = _synthetic_runs((1.0, 1.0, 1.0))
    long = runs["long"]
    short = runs["short"]
    # Copy the short response onto the long reach so the paired centroid
    # difference is exactly zero and the criterion takes the unresolved branch.
    runs["long"] = RunResult(long.case, short.table.copy(), long.meta, long.wall_seconds)
    result = get("wave_celerity_bounds")(runs, _probe(), _params(relative_tolerance=10.0))
    assert not result.passed
    assert "propagation is not downstream/resolved" in result.message
    json.dumps(asdict(result), allow_nan=False)
    assert "c_obs_m_s" not in result.diagnostics
    assert "relative_residual" not in result.diagnostics


def test_zero_response_is_a_scientific_failure_not_an_exception():
    runs = _synthetic_runs()
    run = runs["long"]
    table = run.table.copy()
    flat = np.full(len(table), STATE_Q[0], dtype=float)
    flat[48:144] = STATE_Q[0]
    flat[144:240] = STATE_Q[1]
    flat[240:336] = STATE_Q[2]
    table["dis"] = flat
    runs["long"] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get("wave_celerity_bounds")(runs, _probe(), _params())
    assert not result.passed
    assert result.diagnostics["variant"] == "long"
    assert "no measurable positive response" in result.message


def test_nonfinite_discharge_is_a_scientific_failure_not_an_exception():
    runs = _synthetic_runs()
    run = runs["short"]
    table = run.table.copy()
    table.loc[100, "dis"] = np.nan
    runs["short"] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get("wave_celerity_bounds")(runs, _probe(), _params())
    assert not result.passed
    assert result.diagnostics == {"variant": "short", "nonfinite_count": 1}


def test_state_response_catches_fixed_celerity_when_magnitude_check_is_neutralized():
    result = get("wave_celerity_bounds")(
        _synthetic_runs((1.0, 1.0, 1.0)),
        _probe(),
        _params(relative_tolerance=10.0, ordering_margin_fraction=0.05),
    )
    assert not result.passed
    assert "low->medium" in result.message or "medium->high" in result.message


def test_registered_fixed_celerity_control_trips_gate():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    outcome = run_probe(
        registry.find_model("reference_fixed_celerity"),
        probe,
        gate_seeds(probe.id, 1),
    )
    assert outcome.verdict == "FAIL"
    assert "wave_celerity_bounds" in outcome.failing


def test_saint_venant_daily_uniform_flow_contract_is_preserved():
    model = registry.find_model("reference_saint_venant")
    assert model.supports_timestep("PT1D")
    assert model.supports_timestep("PT1H")

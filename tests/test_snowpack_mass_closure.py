"""Unit tests for the snowpack mass-closure probe and shared criteria.

These tests pin the behaviour introduced by mass/snowpack-mass-closure.

They cover:

- the year-aligned two-cycle probe calendar;
- regime-wise closure preventing opposite residuals from cancelling;
- storage continuity at the start of the scored window and across regimes;
- optional versus required sinks;
- explicit storage scoping;
- the reviewed 5% closure threshold;
- timestep-aware precipitation integration in snowpack_response;
- rejection of incomplete or mislabelled snow cycles.
"""

from __future__ import annotations

import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria.base import FAIL, PASS
from hydroturing.criteria.closure import closure
from hydroturing.criteria.degeneracy import snowpack_response
from hydroturing.harness import build_case
from hydroturing.protocol import Case, RunResult


PROBE_ID = "mass/snowpack-mass-closure"


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe(PROBE_ID)


@pytest.fixture(scope="module")
def closure_params(probe):
    return dict(next(c.params for c in probe.criteria if c.name == "closure"))


@pytest.fixture(scope="module")
def response_params(probe):
    return dict(
        next(c.params for c in probe.criteria if c.name == "snowpack_response")
    )


def _run(
    *,
    timestep: str,
    regimes: list[str],
    precipitation: list[float],
    snow: list[float],
    snm: list[float] | None = None,
    initial_snow: float = 0.0,
    extra_table: dict[str, list[float]] | None = None,
) -> RunResult:
    """Build a small synthetic run with one unscored spinup row.

    The calendar date is intentionally arbitrary here. These tests exercise
    shared criterion arithmetic, not season-dependent snow physics.

    ``precipitation``, ``snow`` and ``snm`` contain scored values only.
    ``extra_table`` contains complete arrays including the spinup row.
    """
    n = len(regimes)
    if len(precipitation) != n or len(snow) != n:
        raise ValueError(
            "regimes, precipitation and snow must have equal length"
        )
    if snm is not None and len(snm) != n:
        raise ValueError("snm must have the same length as regimes")

    frequencies = {
        "PT1D": "D",
        "PT1H": "h",
    }
    if timestep not in frequencies:
        raise ValueError(f"unsupported test timestep {timestep}")

    time = pd.date_range(
        "2000-01-01",
        periods=n + 1,
        freq=frequencies[timestep],
    )

    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "pr": [0.0, *precipitation],
            "tas": [0.0] * (n + 1),
            "pet": [0.0] * (n + 1),
            "_regime": ["spinup", *regimes],
        }
    )

    table_data = {
        "time": forcing["time"].tolist(),
        "snw": [initial_snow, *snow],
    }

    if snm is not None:
        table_data["snm"] = [0.0, *snm]

    if extra_table is not None:
        for name, values in extra_table.items():
            if len(values) != n + 1:
                raise ValueError(
                    f"extra table column '{name}' must include the spinup row"
                )
            table_data[name] = values

    case = Case(
        probe_id=PROBE_ID,
        seed=0,
        forcing=forcing,
        static={},
        spinup_steps=1,
        timestep=timestep,
    )

    return RunResult(
        case=case,
        table=pd.DataFrame(table_data),
        meta={},
        wall_seconds=0.0,
    )


def test_probe_case_uses_year_aligned_two_cycle_calendar(probe):
    """The real probe must keep both snow cycles at comparable seasons.

    One dry spinup day precedes two 360-day
    accumulation-storage-melt cycles. Each stage lasts 120 days.
    """
    case = build_case(probe, seed=0)

    assert case.spinup_steps == 1
    assert len(case.forcing) == 721

    time = pd.to_datetime(case.forcing["time"])

    assert time.iloc[0] == pd.Timestamp("1999-08-31")
    assert time.iloc[1] == pd.Timestamp("1999-09-01")
    assert time.iloc[-1] == pd.Timestamp("2001-08-20")

    assert case.forcing["_regime"].iloc[0] == "spinup"
    assert case.forcing["pr"].iloc[0] == pytest.approx(0.0)

    scored = case.forcing.iloc[case.spinup_steps :].reset_index(drop=True)

    expected_regimes = (
        ["accumulation"] * 120
        + ["storage"] * 120
        + ["melt"] * 120
        + ["accumulation"] * 120
        + ["storage"] * 120
        + ["melt"] * 120
    )

    assert scored["_regime"].tolist() == expected_regimes


def test_segment_residuals_do_not_cancel(probe, closure_params):
    """Opposite errors in different regimes must not cancel.

    Over the whole record the +10 mm and -10 mm residuals sum to zero.
    Segment scoring must instead accumulate their absolute magnitudes.
    """
    run = _run(
        timestep="PT1D",
        regimes=["accumulation", "melt"],
        precipitation=[10.0, 0.0],
        snow=[0.0, 0.0],
        snm=[0.0, 10.0],
    )

    segmented = closure(run, probe, closure_params)

    whole_record_params = dict(closure_params)
    whole_record_params.pop("segment_column")
    whole_record = closure(run, probe, whole_record_params)

    assert whole_record.status == PASS
    assert whole_record.value == pytest.approx(0.0)

    assert segmented.status == FAIL
    assert segmented.value == pytest.approx(2.0)

    residuals = segmented.diagnostics["segment_residuals"]
    assert len(residuals) == 2

    assert residuals[0]["label"] == "accumulation"
    assert residuals[0]["residual"] == pytest.approx(10.0)

    assert residuals[1]["label"] == "melt"
    assert residuals[1]["residual"] == pytest.approx(-10.0)


def test_first_segment_uses_spinup_storage(probe, closure_params):
    """The first scored block starts from the final spinup storage.

    The run starts with 5 mm of snow, receives another 10 mm, and ends the
    first scored step with 15 mm. Closure therefore requires the initial
    5 mm carried from spinup to be used as the segment-start state.
    """
    run = _run(
        timestep="PT1D",
        regimes=["accumulation"],
        precipitation=[10.0],
        snow=[15.0],
        snm=[0.0],
        initial_snow=5.0,
    )

    result = closure(run, probe, closure_params)

    assert result.status == PASS
    assert result.value == pytest.approx(0.0)

    residuals = result.diagnostics["segment_residuals"]
    assert len(residuals) == 1
    assert residuals[0]["residual"] == pytest.approx(0.0)


def test_segment_storage_starts_from_previous_row(probe, closure_params):
    """A later regime starts from storage at the end of the preceding row.

    The first block stores 10 mm. The second block releases exactly those
    10 mm. Both blocks therefore close individually.
    """
    run = _run(
        timestep="PT1D",
        regimes=["accumulation", "melt"],
        precipitation=[10.0, 0.0],
        snow=[10.0, 0.0],
        snm=[0.0, 10.0],
    )

    result = closure(run, probe, closure_params)

    assert result.status == PASS
    assert result.value == pytest.approx(0.0)

    residuals = result.diagnostics["segment_residuals"]
    assert len(residuals) == 2
    assert residuals[0]["residual"] == pytest.approx(0.0)
    assert residuals[1]["residual"] == pytest.approx(0.0)


def test_missing_optional_sublimation_is_allowed(probe, closure_params):
    """The snow budget may omit sbl because it is declared optional."""
    assert "sbl" in closure_params["optional_sinks"]

    run = _run(
        timestep="PT1D",
        regimes=["accumulation"],
        precipitation=[10.0],
        snow=[0.0],
        snm=[10.0],
    )

    result = closure(run, probe, closure_params)

    assert result.status == PASS
    assert result.value == pytest.approx(0.0)


def test_missing_required_sink_raises(probe, closure_params):
    """A sink that is not optional must still be present."""
    assert "snm" in closure_params["sinks"]
    assert "snm" not in closure_params["optional_sinks"]

    run = _run(
        timestep="PT1D",
        regimes=["accumulation"],
        precipitation=[10.0],
        snow=[0.0],
        snm=None,
    )

    with pytest.raises(ValueError, match="closure needs 'snm'"):
        closure(run, probe, closure_params)


def test_explicit_states_ignore_unselected_stores(probe, closure_params):
    """Explicit snow states must keep unrelated stores out of the budget.

    The soil store gains 50 mm only to verify that a reported state outside
    the explicitly selected control volume does not enter snowpack closure.
    """
    run = _run(
        timestep="PT1D",
        regimes=["accumulation"],
        precipitation=[10.0],
        snow=[0.0],
        snm=[10.0],
        extra_table={
            "mrso": [0.0, 50.0],
        },
    )

    explicit = closure(run, probe, closure_params)

    inferred_params = dict(closure_params)
    inferred_params.pop("states")
    inferred = closure(run, probe, inferred_params)

    assert explicit.status == PASS
    assert explicit.value == pytest.approx(0.0)

    assert inferred.status == FAIL
    assert inferred.value > closure_params["threshold"]


@pytest.mark.parametrize(
    ("snm", "expected"),
    [
        (96.0, PASS),
        (94.0, FAIL),
    ],
)
def test_closure_resolves_leaks_around_the_five_percent_limit(
    probe,
    closure_params,
    snm,
    expected,
):
    """A 4% residual passes while a 6% residual fails the 5% rule."""
    run = _run(
        timestep="PT1D",
        regimes=["accumulation"],
        precipitation=[100.0],
        snow=[0.0],
        snm=[snm],
    )

    result = closure(run, probe, closure_params)

    assert result.status == expected


def test_snowpack_response_is_timestep_invariant(probe, response_params):
    """The same 100 mm accumulation must score identically at PT1D and PT1H.

    At PT1H the precipitation rate is 2400 mm/day for one hour, which
    integrates to the same 100 mm depth as 100 mm/day for one day.
    """
    daily = _run(
        timestep="PT1D",
        regimes=["accumulation", "storage", "melt"],
        precipitation=[100.0, 0.0, 0.0],
        snow=[100.0, 100.0, 0.0],
    )

    hourly = _run(
        timestep="PT1H",
        regimes=["accumulation", "storage", "melt"],
        precipitation=[2400.0, 0.0, 0.0],
        snow=[100.0, 100.0, 0.0],
    )

    daily_result = snowpack_response(
        daily,
        probe,
        response_params,
    )
    hourly_result = snowpack_response(
        hourly,
        probe,
        response_params,
    )

    assert daily_result.status == PASS
    assert hourly_result.status == PASS

    assert daily_result.diagnostics[
        "cycle_1_peak_fraction"
    ] == pytest.approx(1.0)

    assert hourly_result.diagnostics[
        "cycle_1_peak_fraction"
    ] == pytest.approx(1.0)

    assert daily_result.diagnostics[
        "cycle_1_melt_fraction"
    ] == pytest.approx(0.0)

    assert hourly_result.diagnostics[
        "cycle_1_melt_fraction"
    ] == pytest.approx(0.0)


def test_snowpack_response_rejects_incomplete_cycle(
    probe,
    response_params,
):
    """An accumulation-storage sequence without melt is not scoreable."""
    run = _run(
        timestep="PT1D",
        regimes=["accumulation", "storage"],
        precipitation=[100.0, 0.0],
        snow=[100.0, 100.0],
    )

    with pytest.raises(
        ValueError,
        match="complete accumulation-storage-melt cycles",
    ):
        snowpack_response(run, probe, response_params)


def test_snowpack_response_rejects_mislabelled_cycle(
    probe,
    response_params,
):
    """A complete three-block record must still use the declared order."""
    run = _run(
        timestep="PT1D",
        regimes=["accumulation", "melt", "storage"],
        precipitation=[100.0, 0.0, 0.0],
        snow=[100.0, 0.0, 0.0],
    )

    with pytest.raises(
        ValueError,
        match="snowpack_response expects regimes",
    ):
        snowpack_response(run, probe, response_params)

def test_bypass_control_is_caught_by_closure_alone(probe):
    """`reference_snow_bypass` must trip `closure` and nothing else.

    The model exists to show `closure` catching water that leaves the
    snowpack while the whole-catchment budget still closes, which is what
    separates it from `reference_leaky`. If it also trips
    `snowpack_response`, a reader cannot tell which criterion the probe is
    demonstrating, and a regression in `closure` alone would leave the model
    failing for the wrong reason.

    `ht gate` cannot catch that: `cli.py` asks only that the declared
    criterion be *among* those tripped, so the gate stays green however many
    extra criteria a control picks up.
    """
    from hydroturing.harness import run_probe

    model = registry.find_model("reference_snow_bypass")
    outcome = run_probe(model, probe, seeds=[0])

    assert outcome.verdict == "FAIL"
    assert set(outcome.failing) == {"closure"}


def test_bypass_control_clears_the_peak_fraction_floor(probe, response_params):
    """The control's retained pack must sit above the floor, with margin.

    `reference_snow_bypass` keeps `1 - BYPASS_FRACTION` of every snowfall by
    construction, so its peak fraction is a constant rather than a
    distribution. Raising `min_peak_fraction` past it, or raising
    `BYPASS_FRACTION`, is what would silently cost the control its isolation.
    """
    from hydroturing.harness import run_probe

    model = registry.find_model("reference_snow_bypass")
    outcome = run_probe(model, probe, seeds=[0])

    response = next(c for c in outcome.criteria if c.name == "snowpack_response")
    peaks = [
        v for k, v in response.diagnostics.items() if k.endswith("peak_fraction")
    ]

    assert peaks, "snowpack_response reported no peak fraction"
    floor = float(response_params["min_peak_fraction"])
    for peak in peaks:
        assert peak > floor, f"peak fraction {peak} is at or below the floor {floor}"
    assert min(peaks) - floor >= 0.04

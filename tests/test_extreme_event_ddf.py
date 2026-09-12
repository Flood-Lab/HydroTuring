"""Independent numerical references for rainfall DDF calibration and reports.

Sample L-moments use the published lmoments3 usage example:
https://lmoments3.readthedocs.io/stable/usage.html
The GEV comparison parameters were calculated with Hosking's PELGEV rational
approximation, independently of the generator's numerical root solver:
https://raw.githubusercontent.com/cran/lmom/master/src/lmoments.f
Support and quantiles use Hosking's k convention, the negative of the usual
extreme-value shape xi. These tests introduce no statistical dependency.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.harness import load_generator


@pytest.fixture(scope="module")
def generator():
    return load_generator(registry.find_probe("mass/extreme-event-closure"))


PUBLISHED_SAMPLE = [2.0, 3.0, 4.0, 2.4, 5.5, 1.2, 5.4, 2.2, 7.1, 1.3, 1.5]
GUMBEL_TAU3 = 0.169925001442312363
EULER_GAMMA = 0.577215664901532861


def fit_parameters(k=0.0, location=10.0, scale=2.0):
    return {"location": location, "scale": scale, "shape_k": k}


def test_sample_l_moments_match_published_software_fixture(generator):
    assert generator.sample_l_moments(PUBLISHED_SAMPLE) == pytest.approx(
        (3.2363636363636363, 1.1418181818181818, 0.27388535031847133),
        abs=1e-13,
    )


def test_unbiased_sample_l_moments_have_hand_calculated_values(generator):
    sample = np.array([4.0, 1.0, 3.0, 2.0])
    assert generator.sample_l_moments(sample) == pytest.approx((2.5, 5 / 6, 0.0), abs=1e-14)
    np.testing.assert_array_equal(sample, [4, 1, 3, 2])


def test_gev_fit_matches_independent_hosking_rational_approximation(generator):
    fit = generator.fit_gev_lmom(PUBLISHED_SAMPLE)
    # Hosking's approximation and a converged numerical root are not identical.
    assert fit["location"] == pytest.approx(2.1792883907597185, abs=2e-6)
    assert fit["scale"] == pytest.approx(1.3956403711604544, abs=2e-6)
    assert fit["shape_k"] == pytest.approx(-0.15556094560377193, abs=2e-6)
    assert fit["l1"] == pytest.approx(3.2363636363636363)
    assert fit["l2"] == pytest.approx(1.1418181818181818)
    assert fit["tau3"] == pytest.approx(0.27388535031847133)
    assert fit["n"] == len(PUBLISHED_SAMPLE)
    json.dumps(fit, allow_nan=False)


def test_gev_fit_reaches_the_gumbel_limit(generator):
    # For [0, 1, x], the unbiased sample L-skewness is (x - 2) / x.
    # Choose x to make it the known population Gumbel L-skewness.
    x = 2 / (1 - GUMBEL_TAU3)
    fit = generator.fit_gev_lmom([0.0, 1.0, x])
    expected_scale = (x / 3) / math.log(2)
    expected_location = (1 + x) / 3 - EULER_GAMMA * expected_scale
    assert fit["shape_k"] == pytest.approx(0.0, abs=1e-8)
    assert fit["scale"] == pytest.approx(expected_scale, abs=1e-7)
    assert fit["location"] == pytest.approx(expected_location, abs=1e-7)


@pytest.mark.parametrize("sample", [[], [1], [1, 2], [1, 1, 1], [1, 2, np.nan],
                                         [1, 2, np.inf], [[1, 2], [3, 4]]])
def test_gev_fit_rejects_undefined_sample_l_moments(generator, sample):
    with pytest.raises(ValueError):
        generator.fit_gev_lmom(sample)


@pytest.mark.parametrize("period, expected", [
    (2, 10.733025841163329),
    (10, 14.50073465462489),
    (100, 19.20029845355316),
    (500, 22.427214528174922),
])
def test_gumbel_return_levels_match_closed_form_numbers(generator, period, expected):
    fit = fit_parameters()
    assert generator.gev_return_level(fit, period) == pytest.approx(expected, abs=1e-12)
    assert generator.gev_return_period(fit, expected) == pytest.approx(period, rel=1e-12)


@pytest.mark.parametrize("k, expected", [
    (0.5, 10 + 4 * (1 - math.sqrt(math.log(2)))),
    (-0.5, 10 + 4 * (1 / math.sqrt(math.log(2)) - 1)),
])
def test_nonzero_shape_medians_match_square_root_closed_forms(generator, k, expected):
    fit = fit_parameters(k=k)
    assert generator.gev_return_level(fit, 2) == pytest.approx(expected, abs=1e-12)
    assert generator.gev_return_period(fit, expected) == pytest.approx(2, abs=1e-12)


@pytest.mark.parametrize("k", [-0.5, 0.0, 0.5])
def test_location_has_the_same_probability_for_every_shape(generator, k):
    # F(location) = exp(-1), independently of shape or scale.
    assert generator.gev_return_period(fit_parameters(k=k), 10) == pytest.approx(
        1.5819767068693265, abs=1e-13,
    )


@pytest.mark.parametrize("depth", [14.0, 15.0])
def test_positive_hosking_shape_has_a_finite_upper_endpoint(generator, depth):
    assert math.isinf(generator.gev_return_period(fit_parameters(k=0.5), depth))


@pytest.mark.parametrize("depth", [6.0, 5.0])
def test_negative_hosking_shape_has_a_finite_lower_endpoint(generator, depth):
    assert generator.gev_return_period(fit_parameters(k=-0.5), depth) == 1.0


def test_extreme_probabilities_do_not_round_the_cdf_to_one(generator):
    # 1 - 1e-18 rounds to 1 in float64; the direct survival calculation must
    # still resolve this finite return period without forming that subtraction.
    fit = fit_parameters()
    expected_depth = 92.89306334778564
    assert generator.gev_return_level(fit, 1e18) == pytest.approx(expected_depth, abs=1e-11)
    assert generator.gev_return_period(fit, expected_depth) == pytest.approx(1e18, rel=1e-12)


def test_gumbel_probability_overflow_limits_are_explicit(generator):
    fit = fit_parameters(location=0, scale=1)
    assert generator.gev_return_period(fit, -1000.0) == 1.0
    assert math.isinf(generator.gev_return_period(fit, 1000.0))


@pytest.mark.parametrize("period", [0, 1, -10, np.nan, np.inf])
def test_invalid_return_period_is_rejected(generator, period):
    with pytest.raises(ValueError):
        generator.gev_return_level(fit_parameters(), period)


def test_annual_maxima_assign_cross_boundary_windows_by_their_end(generator):
    # Two prefix days + two five-day years. A 40, 50 mm storm ends year 1.
    # Its 40 + 50 + 0 window ends in year 2 and must not be discarded.
    rain = np.array([10, 20, 0, 0, 0, 40, 50, 0, 0, 0, 0, 0], dtype=float)
    result = generator.annual_maxima(rain, years=2, durations=(1, 3), days_per_year=5)
    assert set(result) == {1, 3}
    np.testing.assert_array_equal(result[1], [50, 0])
    np.testing.assert_array_equal(result[3], [90, 90])


def test_prefix_contributes_to_rolling_totals_but_not_annual_daily_maxima(generator):
    result = generator.annual_maxima(
        np.array([100, 200, 1, 0, 0, 0, 0]), years=1, durations=(1, 3), days_per_year=5,
    )
    np.testing.assert_array_equal(result[1], [1])
    np.testing.assert_array_equal(result[3], [301])


def test_annual_rolling_totals_may_span_dry_days_and_separate_events(generator):
    result = generator.annual_maxima(
        np.array([0, 0, 20, 0, 30, 0, 0]), years=1, durations=(1, 3), days_per_year=5,
    )
    np.testing.assert_array_equal(result[1], [30])
    np.testing.assert_array_equal(result[3], [50])


@pytest.mark.parametrize("rain", [np.zeros(6), np.zeros(8),
                                       np.array([0, 0, 1, 2, -1, 0, 0]),
                                       np.array([0, 0, 1, 2, np.nan, 0, 0]),
                                       np.array([0, 0, 1, 2, np.inf, 0, 0])])
def test_annual_maxima_reject_wrong_length_or_invalid_rain(generator, rain):
    with pytest.raises(ValueError):
        generator.annual_maxima(rain, years=1, durations=(1, 3), days_per_year=5)


@pytest.mark.parametrize("durations", [(), (0,), (-1,), (1, 2.5)])
def test_annual_maxima_reject_invalid_durations(generator, durations):
    with pytest.raises(ValueError):
        generator.annual_maxima(np.zeros(7), years=1, durations=durations, days_per_year=5)


def test_calibration_has_one_maximum_per_year_and_serializes_reproducibly(generator):
    calibration = generator.calibrate_ddf(20260912)
    assert calibration == generator.calibrate_ddf(20260912)
    assert calibration["years"] == 100
    assert calibration["days_per_year"] == 365
    assert calibration["prefix_days"] == 6
    assert calibration["random_stream"] == 1
    assert calibration["durations_days"] == [1, 3, 7]
    for duration in ("1", "3", "7"):
        maxima = np.array(calibration["annual_maxima_mm"][duration])
        assert maxima.shape == (100,)
        assert np.isfinite(maxima).all() and (maxima > 0).all()
        assert calibration["fits"][duration]["n"] == 100
        levels = calibration["return_levels_mm"][duration]
        assert set(levels) == {"2", "5", "10", "20", "50", "100", "200", "500"}
        assert np.all(np.diff([levels[str(t)] for t in (2, 5, 10, 20, 50, 100, 200, 500)]) > 0)
    json.dumps(calibration, allow_nan=False)
    other = generator.calibrate_ddf(20260913)
    assert calibration["annual_maxima_mm"] != other["annual_maxima_mm"]


def test_changing_calibration_preserves_baseline_weather_history_and_rainfall_mass(generator, monkeypatch):
    baseline, baseline_static = generator.generate_baseline(20260912)
    forcing, static = generator.generate(20260912)
    calibrate = generator.calibrate_ddf
    monkeypatch.setattr(generator, "calibrate_ddf", lambda seed: calibrate(seed, years=50))
    alternate, alternate_static = generator.generate(20260912)
    alternate_baseline, alternate_baseline_static = generator.generate_baseline(20260912)
    assert baseline.to_csv(index=False) == alternate_baseline.to_csv(index=False)
    assert static == alternate_static == baseline_static == alternate_baseline_static
    # Calibration now controls how many source events are combined. Their
    # rainfall timing may change, while the weather source and total mass do not.
    for frame in (forcing, alternate):
        for column in ("time", "tas", "pet"):
            np.testing.assert_array_equal(frame[column], baseline[column])
        selection = frame.attrs["rainfall_diagnostics"]["selection"]
        start, stop = selection["start_row"], selection["stop_row"]
        np.testing.assert_array_equal(frame.pr.iloc[:start], baseline.pr.iloc[:start])
        np.testing.assert_array_equal(frame.pr.iloc[stop:], baseline.pr.iloc[stop:])
        assert frame.pr.sum() == pytest.approx(baseline.pr.sum(), rel=0, abs=1e-8)
        assert frame.pr.iloc[start:stop].sum() == pytest.approx(baseline.pr.iloc[start:stop].sum(), rel=0, abs=1e-8)
    assert forcing.attrs["rainfall_diagnostics"]["calibration"]["years"] == 100
    assert alternate.attrs["rainfall_diagnostics"]["calibration"]["years"] == 50
    assert forcing.attrs["rainfall_diagnostics"]["selection"] == alternate.attrs["rainfall_diagnostics"]["selection"]


def test_generated_diagnostics_identify_median_wet_year_and_are_json_safe(generator):
    forcing, _ = generator.generate(20260912)
    baseline, _ = generator.generate_baseline(20260912)
    report = forcing.attrs["rainfall_diagnostics"]
    assert report["seed"] == 20260912
    assert "group_size" not in report
    assert report["group_count"] == forcing.loc[forcing["_event_id"] > 0, "_event_id"].nunique()
    assert report["group_count"] == len(report["groups"])
    annual_rain = baseline.pr.to_numpy()[365:].reshape(20, 365).sum(axis=1)
    selection = report["selection"]
    # Both central ranks are equidistant for an even sample: take the earlier
    # chronological block rather than let floating-point subtraction choose.
    ordered = np.argsort(annual_rain)
    selected_index = min(ordered[9:11])
    start, stop = 365 * (selected_index + 1), 365 * (selected_index + 2)
    assert selection["year_number"] == selected_index + 1
    assert (selection["start_row"], selection["stop_row"]) == (start, stop)
    assert selection["baseline_precip_mm"] == pytest.approx(annual_rain[selected_index])
    assert selection["median_precip_mm"] == pytest.approx(np.median(annual_rain))
    assert selection["scored_precip_mm"] == pytest.approx(annual_rain.sum())
    assert selection["selected_fraction_of_scored_precip"] == pytest.approx(
        annual_rain[selected_index] / annual_rain.sum(),
    )
    np.testing.assert_allclose(selection["annual_precip_mm"], annual_rain)
    assert selection["start_time"] == baseline.time.iloc[start]
    assert selection["end_time"] == baseline.time.iloc[stop - 1]
    assert isinstance(selection["selection"], str) and selection["selection"]
    assert "final_year" not in report
    assert [row["duration_days"] for row in report["modified_year"]] == [1, 3, 7]
    assert report["modified_year"] == generator.describe_modified_year(
        forcing, report["calibration"], start, stop,
    )
    assert json.loads(json.dumps(report, allow_nan=False)) == report


def test_constructed_group_depths_and_strict_threshold_diagnostics_agree(generator):
    forcing, _ = generator.generate(20260912)
    report = forcing.attrs["rainfall_diagnostics"]
    construction = report["construction"]
    assert construction["target_return_period_years"] == 100
    assert report["calibration"]["years"] == 100
    assert construction["inter_event_dry_days"] == 1
    anchor = construction["anchor_row"]
    assert forcing.time.iloc[anchor] == construction["anchor_time"]
    assert pd.Timestamp(construction["anchor_time"]).strftime("%m-%d") == "05-01"
    start, stop = report["selection"]["start_row"], report["selection"]["stop_row"]
    assert isinstance(construction["placement"], str) and construction["placement"]
    assert isinstance(construction["window_convention"], str) and construction["window_convention"]
    thresholds = construction["thresholds_mm"]
    assert thresholds == {
        str(d): report["calibration"]["return_levels_mm"][str(d)]["100"] for d in (1, 3, 7)
    }
    assert len(report["events"]) == 3 * len(report["groups"])
    event_rows = {(row["event_id"], row["duration_days"]): row for row in report["events"]}
    assert len(event_rows) == len(report["events"])
    for index, group in enumerate(report["groups"], start=1):
        assert group["event_id"] == index
        assert group["source_event_count"] == len(group["source_events"])
        assert group["source_event_count"] >= 1
        assert start <= group["start"] < group["stop"] <= stop
        assert group["original_start"] == group["source_events"][0][0]
        assert group["original_stop"] - group["original_start"] == group["stop"] - group["start"]
        if index == 1:
            assert group["start"] == construction["anchor_row"]
            assert forcing.time.iloc[group["start"]] == construction["anchor_time"]
        else:
            previous = report["groups"][index - 2]
            assert group["start"] == previous["stop"] + 1
            assert forcing.pr.iloc[previous["stop"]] == 0
        isolated_rain = forcing.pr.iloc[group["start"]:group["stop"]].to_numpy()
        for duration in (1, 3, 7):
            # For nonnegative isolated rainfall, adding zeros outside the
            # event yields its total if shorter than d, otherwise the largest
            # full d-day sum. Neighboring rain is deliberately excluded.
            depth = max(
                isolated_rain[i:i + duration].sum()
                for i in range(max(1, len(isolated_rain) - duration + 1))
            )
            assert group["depths_mm"][str(duration)] == pytest.approx(depth, rel=0, abs=1e-6)
            row = event_rows[(index, duration)]
            assert row["source_event_count"] == group["source_event_count"]
            assert row["event_duration_days"] == len(isolated_rain)
            assert row["start_time"] == forcing.time.iloc[group["start"]]
            assert row["end_time"] == forcing.time.iloc[group["stop"] - 1]
            assert row["precip_mm"] == pytest.approx(isolated_rain.sum(), rel=0, abs=1e-6)
            assert row["max_depth_mm"] == group["depths_mm"][str(duration)]
            assert row["threshold_mm"] == thresholds[str(duration)]
            assert row["exceeds_threshold"] is (group["depths_mm"][str(duration)] > thresholds[str(duration)])
            if row["exceeds_threshold"]:
                assert row["extrapolates_record"] is True
                assert row["return_period_years"] is None or row["return_period_years"] > 100
            assert row["target_exceeded"] is group["target_exceeded"]
            assert row["stop_reason"] == group["stop_reason"]
        exceeded = [d for d in (1, 3, 7) if group["depths_mm"][str(d)] > thresholds[str(d)]]
        assert group["exceeded_durations_days"] == exceeded
        assert group["target_exceeded"] is bool(exceeded)
        assert group["stop_reason"] == ("threshold_exceeded" if exceeded else "end_of_year")
    assert report["target_exceeded_group_count"] == sum(g["target_exceeded"] for g in report["groups"])


def test_constructed_event_report_keeps_below_target_groups_and_zero_extends_short_events(generator):
    forcing = empty_forcing(generator)
    forcing.loc[3650, "pr"] = 500
    forcing.loc[3652:3653, "pr"] = [20, 50]
    forcing.loc[3655, "pr"] = 1000
    group = {
        "event_id": 1, "start": 3652, "stop": 3654, "source_event_count": 2,
        "target_exceeded": False, "stop_reason": "end_of_year",
    }
    # Set the Gumbel Q100 threshold at 70 mm, exactly the event's 3/7-day
    # zero-extended total, using the independent standard Gumbel quantile.
    q100 = -math.log(-math.log1p(-1 / 100))
    calibration = description_calibration(fit=fit_parameters(location=0, scale=70 / q100))
    calibration["return_levels_mm"] = {str(d): {"100": 70.0} for d in (1, 3, 7)}
    rows = generator.describe_constructed_events(forcing, [group], calibration)
    assert [row["duration_days"] for row in rows] == [1, 3, 7]
    assert [row["max_depth_mm"] for row in rows] == [50, 70, 70]
    for row in rows:
        assert row["event_duration_days"] == 2
        assert row["precip_mm"] == 70
        assert row["exceeds_threshold"] is False  # Equality does not trigger construction.
        assert row["target_exceeded"] is False
        assert row["stop_reason"] == "end_of_year"
    assert rows[1]["return_period_years"] == pytest.approx(100)
    assert rows[2]["return_period_years"] == pytest.approx(100)
    actual_year = generator.describe_modified_year(forcing, calibration, 3650, 4015)
    assert all(row["max_depth_mm"] >= 1000 for row in actual_year)
    json.dumps(rows, allow_nan=False)


def test_packing_preserves_isolated_event_depths_but_can_raise_actual_multiday_maxima(generator):
    original = empty_forcing(generator)
    original.loc[3660:3661, "pr"] = [20, 50]
    original.loc[3680:3681, "pr"] = [50, 20]
    groups = [{
        "event_id": event_id, "start": start, "stop": start + 2,
        "source_events": ((start, start + 2),), "source_event_count": 1,
        "depths_mm": {"1": 50.0, "3": 70.0, "7": 70.0},
        "exceeded_durations_days": [1], "target_exceeded": True,
        "stop_reason": "threshold_exceeded",
    } for event_id, start in ((1, 3660), (2, 3680))]
    thresholds = {1: 40.0, 3: 70.0, 7: 70.0}
    q100 = -math.log(-math.log1p(-1 / 100))
    calibration = description_calibration()
    calibration["fits"] = {
        str(d): fit_parameters(location=0, scale=depth / q100) for d, depth in thresholds.items()
    }
    calibration["return_levels_mm"] = {str(d): {"100": depth} for d, depth in thresholds.items()}
    packed_rain, packed_groups = generator.pack_event_groups(original.pr.to_numpy(), groups, gap_days=1)
    packed = original.copy(deep=True)
    packed["pr"] = packed_rain
    assert [(group["start"], group["stop"]) for group in packed_groups] == [(3660, 3662), (3663, 3665)]
    assert [(group["original_start"], group["original_stop"]) for group in packed_groups] == [(3660, 3662), (3680, 3682)]
    assert [group["source_events"] for group in packed_groups] == [group["source_events"] for group in groups]
    assert packed.pr.sum() == original.pr.sum() == 140
    before_events = generator.describe_constructed_events(original, groups, calibration)
    after_events = generator.describe_constructed_events(packed, packed_groups, calibration)
    for before, after in zip(before_events, after_events, strict=True):
        for field in ("event_id", "source_event_count", "event_duration_days", "duration_days",
                      "precip_mm", "max_depth_mm", "threshold_mm", "exceeds_threshold",
                      "return_period_years", "target_exceeded"):
            assert after[field] == before[field]
    assert before_events[3]["start_time"] == original.time.iloc[3680]
    assert after_events[3]["start_time"] == packed.time.iloc[3663]
    assert after_events[3]["end_time"] == packed.time.iloc[3664]
    before_year = generator.describe_modified_year(original, calibration, 3650, 4015)
    after_year = generator.describe_modified_year(packed, calibration, 3650, 4015)
    assert [row["max_depth_mm"] for row in before_year] == [50, 70, 70]
    assert [row["max_depth_mm"] for row in after_year] == [50, 100, 140]


def description_calibration(durations=(1, 3, 7), fit=None):
    """A fixed calibration suffices to check window selection independently."""
    return {
        "years": 100,
        "durations_days": list(durations),
        "fits": {str(d): dict(fit or fit_parameters()) for d in durations},
        "annual_maxima_mm": {str(d): [1000.0] for d in durations},
    }


def empty_forcing(generator):
    return pd.DataFrame({
        "time": pd.date_range("2000-01-01", periods=generator.N_STEPS).strftime("%Y-%m-%d"),
        "pr": np.zeros(generator.N_STEPS),
    })


@pytest.mark.parametrize("year_number", [1, 5, 20])
def test_modified_year_window_can_use_prefix_and_durations_can_choose_different_storms(generator, year_number):
    forcing = empty_forcing(generator)
    first, stop = 365 * year_number, 365 * (year_number + 1)
    forcing.loc[first - 10, "pr"] = 500  # Earlier maxima are outside every eligible window.
    forcing.loc[first - 6, "pr"] = 7
    forcing.loc[first - 1:first, "pr"] = [100, 1]
    forcing.loc[first + 30, "pr"] = 50
    if stop < len(forcing):
        forcing.loc[stop, "pr"] = 1000  # A window ending here belongs to the following year.
        forcing.loc[stop + 30, "pr"] = 2000
    daily, three_day, seven_day = generator.describe_modified_year(
        forcing, description_calibration(), first, stop,
    )
    assert daily["max_depth_mm"] == 50
    assert daily["start_time"] == daily["end_time"] == forcing.time.iloc[first + 30]
    for record, duration, depth in ((three_day, 3, 101), (seven_day, 7, 108)):
        assert record["max_depth_mm"] == depth
        assert record["start_time"] == forcing.time.iloc[first - duration + 1]
        assert record["end_time"] == forcing.time.iloc[first]


@pytest.mark.parametrize("year_number", [1, 5, 20])
def test_modified_year_maximum_includes_last_day_and_dry_gaps_but_excludes_stop(generator, year_number):
    forcing = empty_forcing(generator)
    start, stop = 365 * year_number, 365 * (year_number + 1)
    forcing.loc[stop - 3:stop - 1, "pr"] = [20, 0, 30]
    if stop < len(forcing):
        forcing.loc[stop, "pr"] = 1000
    daily, three_day, seven_day = generator.describe_modified_year(
        forcing, description_calibration(), start, stop,
    )
    assert [record["max_depth_mm"] for record in (daily, three_day, seven_day)] == [30, 50, 50]
    assert all(record["end_time"] == forcing.time.iloc[stop - 1] for record in (daily, three_day, seven_day))


@pytest.mark.parametrize("fit, depth, status", [
    (fit_parameters(k=0.5), 15.0, "above_fitted_upper_endpoint"),
    (fit_parameters(location=0, scale=1), 1000.0, "tail_not_representable"),
])
def test_unrepresentable_return_period_is_null_with_explanatory_status(generator, fit, depth, status):
    forcing = empty_forcing(generator)
    start, stop = 5 * 365, 6 * 365
    forcing.loc[stop - 1, "pr"] = depth
    report = generator.describe_modified_year(forcing, description_calibration((1,), fit), start, stop)[0]
    assert report["return_period_years"] is None
    assert report["return_period_status"] == status
    assert report["extrapolates_record"] is True
    json.dumps(report, allow_nan=False)


def test_return_period_extrapolation_and_observed_maximum_flags_are_distinct(generator):
    forcing = empty_forcing(generator)
    # Gumbel(10, 2) has this closed-form 200-year depth, below our fixture's
    # observed maximum. Neither flag should be inferred from the other.
    start, stop = 5 * 365, 6 * 365
    forcing.loc[stop - 1, "pr"] = 10 - 2 * math.log(-math.log1p(-1 / 200))
    report = generator.describe_modified_year(forcing, description_calibration((1,)), start, stop)[0]
    assert report["return_period_years"] == pytest.approx(200)
    assert report["return_period_status"] == "estimated"
    assert report["extrapolates_record"] is True
    assert report["exceeds_calibration_max"] is False


def test_modified_year_description_requires_the_full_record(generator):
    with pytest.raises(ValueError):
        generator.describe_modified_year(empty_forcing(generator).iloc[:-1], description_calibration(), 365, 730)


@pytest.mark.parametrize("start, stop", [
    (0, 365),             # Spinup is not one of the twenty scored years.
    (366, 731),           # A full-length window still has to align with the year blocks.
    (365, 729),
    (365, 731),
    (3650, 4016),
    (7665, 8030),
    (-365, 0),
])
def test_modified_year_description_rejects_nonannual_or_out_of_record_bounds(generator, start, stop):
    with pytest.raises(ValueError):
        generator.describe_modified_year(empty_forcing(generator), description_calibration(), start, stop)

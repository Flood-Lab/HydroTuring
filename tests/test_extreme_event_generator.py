"""Greedy Q100 storms packed from May 1 in the median-wet scored year."""

from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.harness import (
    WindowBounds, build_case, load_generator, resolve_window_days, window_case,
)
from hydroturing.protocol import stage


HIGH_THRESHOLDS = {1: 1000.0, 3: 1000.0, 7: 1000.0}


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe("mass/extreme-event-closure")


@pytest.fixture(scope="module")
def generator(probe):
    return load_generator(probe)


def test_a_single_dry_day_separates_events(generator):
    assert generator.rainfall_events(np.array([1, 2, 0, 3, 0, 4, 5])) == [
        (0, 2), (3, 4), (5, 7)
    ]
    assert generator.rainfall_events(np.array([])) == []


@pytest.mark.parametrize("rain,expected", [
    ([6], {1: 6, 3: 6, 7: 6}),
    ([2, 4], {1: 4, 3: 6, 7: 6}),
    ([1, 2, 3, 4, 5, 6, 7, 8], {1: 8, 3: 21, 7: 35}),
])
def test_duration_depths_use_only_the_event_with_zero_padding(generator, rain, expected):
    assert generator.event_duration_maxima(np.array(rain)) == expected


@pytest.mark.parametrize("trigger,limit", [(1, 1.5), (3, 4.5), (7, 10.5)])
def test_any_one_duration_can_close_a_group(generator, trigger, limit):
    # One seven-day source has depths 1, 3, 7; two aligned sources have 2, 6, 14.
    rain = np.array([0] + [1] * 7 + [0] + [1] * 7 + [0] + [1] * 7 + [0], dtype=float)
    thresholds = {**HIGH_THRESHOLDS, trigger: limit}
    result, groups = generator.overlap_event_groups(rain, start=0, thresholds=thresholds)
    expected = rain.copy()
    expected[1:8] = 2.0
    expected[9:16] = 0.0
    np.testing.assert_array_equal(result, expected)
    assert len(groups) == 2
    first, remainder = groups
    assert first["source_events"] == ((1, 8), (9, 16))
    assert first["source_event_count"] == 2
    assert first["exceeded_durations_days"] == [trigger]
    assert first["depths_mm"] == {"1": 2, "3": 6, "7": 14}
    assert first["target_exceeded"]
    assert first["stop_reason"] == "threshold_exceeded"
    assert remainder["source_events"] == ((17, 24),)
    assert not remainder["target_exceeded"]
    assert remainder["stop_reason"] == "end_of_year"


def test_equal_threshold_keeps_adding_then_starts_at_the_next_unused_event(generator):
    rain = np.array([0, 3, 0, 2, 0, 1, 0, 9, 0], dtype=float)
    result, groups = generator.overlap_event_groups(
        rain, start=0, thresholds={**HIGH_THRESHOLDS, 1: 5},
    )
    np.testing.assert_array_equal(result, [0, 6, 0, 0, 0, 0, 0, 9, 0])
    assert [group["source_events"] for group in groups] == [
        ((1, 2), (3, 4), (5, 6)), ((7, 8),),
    ]
    assert [group["source_event_count"] for group in groups] == [3, 1]
    assert all(group["target_exceeded"] for group in groups)


def test_one_micro_mm_above_the_threshold_is_not_hidden_by_a_tolerance(generator):
    rain = np.array([0, 1, 0, 0.000001, 0])
    result, groups = generator.overlap_event_groups(
        rain, start=0, thresholds={**HIGH_THRESHOLDS, 1: 1},
    )
    assert groups[0]["source_event_count"] == 2
    assert groups[0]["target_exceeded"]
    assert groups[0]["depths_mm"]["1"] == 1.000001
    np.testing.assert_array_equal(result, [0, 1.000001, 0, 0, 0])


def test_an_already_exceeding_first_event_closes_as_a_singleton(generator):
    rain = np.array([0, 6, 0, 2, 0, 2, 0, 1, 0], dtype=float)
    result, groups = generator.overlap_event_groups(
        rain, start=0, thresholds={**HIGH_THRESHOLDS, 1: 5},
    )
    np.testing.assert_array_equal(result, [0, 6, 0, 5, 0, 0, 0, 0, 0])
    assert groups[0]["source_events"] == ((1, 2),)
    assert groups[0]["target_exceeded"]
    assert groups[1]["source_event_count"] == 3
    assert groups[1]["depths_mm"]["1"] == 5
    assert not groups[1]["target_exceeded"]
    assert groups[1]["stop_reason"] == "end_of_year"


def test_added_sources_are_consumed_whole_even_if_the_first_day_crosses(generator):
    rain = np.array([0, 2, 0, 4, 5, 6, 0, 1, 0], dtype=float)
    before = rain.copy()
    result, groups = generator.overlap_event_groups(
        rain, start=0, thresholds={**HIGH_THRESHOLDS, 1: 5},
    )
    np.testing.assert_array_equal(result, [0, 6, 5, 6, 0, 0, 0, 1, 0])
    np.testing.assert_array_equal(rain, before)
    assert not np.shares_memory(result, rain)
    assert result.sum() == rain.sum()
    assert groups[0]["source_events"] == ((1, 2), (3, 6))
    assert (groups[0]["start"], groups[0]["stop"]) == (1, 4)
    assert groups[1]["source_events"] == ((7, 8),)


def test_threshold_windows_never_borrow_rain_from_a_neighboring_event(generator):
    rain = np.array([0, 6, 0, 6, 0, 6, 0], dtype=float)
    result, groups = generator.overlap_event_groups(
        rain, start=0, thresholds={1: 100, 3: 10, 7: 10},
    )
    # A full-record 3-day window already contains 12 mm, but the first isolated
    # event has only 6 mm. It must consume the second event before exceeding.
    np.testing.assert_array_equal(result, [0, 12, 0, 0, 0, 6, 0])
    assert groups[0]["source_event_count"] == 2
    assert groups[0]["depths_mm"] == {"1": 12, "3": 12, "7": 12}
    assert groups[0]["exceeded_durations_days"] == [3, 7]
    assert not groups[1]["target_exceeded"]


@pytest.mark.parametrize("count", [1, 2, 4])
def test_all_remaining_events_form_one_unmet_group_including_a_singleton(generator, count):
    rain = np.zeros(2 * count + 1)
    rain[1::2] = np.arange(1, count + 1)
    result, groups = generator.overlap_event_groups(rain, start=0, thresholds=HIGH_THRESHOLDS)
    expected = np.zeros_like(rain)
    expected[1] = {1: 1, 2: 3, 4: 10}[count]
    np.testing.assert_array_equal(result, expected)
    assert len(groups) == 1
    assert groups[0]["source_event_count"] == count
    assert not groups[0]["target_exceeded"]
    assert groups[0]["exceeded_durations_days"] == []
    assert groups[0]["stop_reason"] == "end_of_year"
    assert result.sum() == rain.sum()


def test_crossing_and_unfinished_events_are_not_cut_or_grouped(generator):
    rain = np.array([0, 1, 2, 0, 3, 0, 4, 0, 5, 6], dtype=float)
    result, groups = generator.overlap_event_groups(rain, start=2, thresholds=HIGH_THRESHOLDS)
    np.testing.assert_array_equal(result, [0, 1, 2, 0, 7, 0, 0, 0, 5, 6])
    assert [group["source_events"] for group in groups] == [((4, 5), (6, 7))]


def test_event_at_the_window_start_is_complete_after_an_observed_dry_day(generator):
    rain = np.array([0, 1, 0, 2, 0], dtype=float)
    result, groups = generator.overlap_event_groups(rain, start=1, thresholds=HIGH_THRESHOLDS)
    np.testing.assert_array_equal(result, [0, 3, 0, 0, 0])
    assert groups[0]["start"] == 1


def test_event_crossing_an_interior_stop_is_unchanged(generator):
    rain = np.array([0, 1, 0, 2, 0, 3, 4, 0], dtype=float)
    result, groups = generator.overlap_event_groups(
        rain, start=0, stop=6, thresholds=HIGH_THRESHOLDS,
    )
    np.testing.assert_array_equal(result, [0, 3, 0, 0, 0, 3, 4, 0])
    assert len(groups) == 1


def test_record_edges_have_no_unobserved_dry_neighbors(generator):
    rain = np.array([9, 0, 1, 0, 2, 0, 8], dtype=float)
    result, groups = generator.overlap_event_groups(rain, start=0, thresholds=HIGH_THRESHOLDS)
    np.testing.assert_array_equal(result, [9, 0, 3, 0, 0, 0, 8])
    assert len(groups) == 1


@pytest.mark.parametrize("rain", [[], [0, 0, 0], [2, 0, 3]])
def test_no_complete_event_returns_unchanged_rain(generator, rain):
    result, groups = generator.overlap_event_groups(
        np.array(rain), start=0, thresholds=HIGH_THRESHOLDS,
    )
    np.testing.assert_array_equal(result, rain)
    assert groups == []


@pytest.mark.parametrize("rain", [[0, -1], [0, np.nan], [0, np.inf], [[1, 2]]])
def test_invalid_rain_is_a_generator_error(generator, rain):
    with pytest.raises(ValueError):
        generator.overlap_event_groups(np.array(rain), start=0, thresholds=HIGH_THRESHOLDS)


@pytest.mark.parametrize("thresholds", [
    {}, {1: 0}, {1: -1}, {1: np.nan}, {1: np.inf}, {0: 10}, {1.5: 10},
])
def test_thresholds_require_positive_finite_depths_and_integer_durations(generator, thresholds):
    with pytest.raises(ValueError):
        generator.overlap_event_groups(np.array([0, 1, 0]), start=0, thresholds=thresholds)


@pytest.mark.parametrize("start,stop", [(-1, 3), (0, 4), (2, 1)])
def test_invalid_window_is_a_generator_error(generator, start, stop):
    with pytest.raises(ValueError, match="bounds"):
        generator.overlap_event_groups(
            np.array([0, 1, 0]), start=start, stop=stop, thresholds=HIGH_THRESHOLDS,
        )


def test_packing_clears_old_locations_before_placing_consecutive_storms(generator):
    # New group locations overlap old locations, including the second group's
    # own source span. Clearing after writing would erase some of these depths.
    rain = np.array([
        9, 0, 6, 0, 0, 0, 7, 8, 9, 10, 0, 0, 11, 12, 0, 0, 0, 0, 1, 0, 0, 5,
    ], dtype=float)
    overlapped, original_groups = generator.overlap_event_groups(
        rain, start=0, thresholds={**HIGH_THRESHOLDS, 1: 5},
    )
    before_rain = overlapped.copy()
    before_groups = deepcopy(original_groups)
    packed, groups = generator.pack_event_groups(overlapped, original_groups)
    np.testing.assert_array_equal(packed, [
        9, 0, 6, 0, 7, 8, 9, 10, 0, 11, 12, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 5,
    ])
    np.testing.assert_array_equal(overlapped, before_rain)
    assert original_groups == before_groups
    assert not np.shares_memory(packed, overlapped)
    assert groups is not original_groups
    assert packed.sum() == rain.sum()
    assert [group["start"] for group in groups] == [2, 4, 9, 12]
    assert [group["stop"] for group in groups] == [3, 8, 11, 13]
    assert groups[0]["start"] == original_groups[0]["start"]
    assert not groups[-1]["target_exceeded"]
    for original, group in zip(original_groups, groups):
        assert group["original_start"] == original["start"]
        assert group["original_stop"] == original["stop"]
        for key, value in original.items():
            if key not in {"start", "stop"}:
                assert group[key] == value, key
        np.testing.assert_array_equal(
            packed[group["start"]:group["stop"]],
            overlapped[original["start"]:original["stop"]],
        )
    for previous, following in zip(groups[:-1], groups[1:]):
        assert following["start"] == previous["stop"] + 1
        assert packed[previous["stop"]] == 0


def test_packing_preserves_events_excluded_by_window_boundaries(generator):
    rain = np.array([0, 8, 9, 0, 6, 0, 0, 7, 0, 2, 3, 0], dtype=float)
    overlapped, groups = generator.overlap_event_groups(
        rain, start=2, stop=10, thresholds={**HIGH_THRESHOLDS, 1: 5},
    )
    packed, groups = generator.pack_event_groups(overlapped, groups)
    np.testing.assert_array_equal(packed, [0, 8, 9, 0, 6, 0, 7, 0, 0, 2, 3, 0])
    assert [(g["start"], g["stop"]) for g in groups] == [(4, 5), (6, 7)]
    assert packed.sum() == rain.sum()


def test_packing_a_single_group_preserves_its_anchor_and_shape(generator):
    rain = np.array([9, 0, 0, 6, 7, 0, 8], dtype=float)
    overlapped, groups = generator.overlap_event_groups(rain, start=0, thresholds=HIGH_THRESHOLDS)
    packed, groups = generator.pack_event_groups(overlapped, groups)
    np.testing.assert_array_equal(packed, rain)
    assert len(groups) == 1
    assert (groups[0]["start"], groups[0]["stop"]) == (3, 5)
    assert (groups[0]["original_start"], groups[0]["original_stop"]) == (3, 5)


@pytest.mark.parametrize("anchor,expected_spans", [
    (2, [(2, 4), (5, 6), (7, 8)]),
    (12, [(12, 14), (15, 16), (17, 18)]),
])
def test_explicit_anchor_moves_whole_sequence_earlier_or_later(generator, anchor, expected_spans):
    rain = np.zeros(24)
    rain[0], rain[-1] = 9, 8  # Ungrouped record-edge rain stays untouched.
    rain[5:7], rain[11], rain[17] = [6, 7], 8, 1
    overlapped, original_groups = generator.overlap_event_groups(
        rain, start=0, thresholds={**HIGH_THRESHOLDS, 1: 5},
    )
    before_groups = deepcopy(original_groups)
    packed, groups = generator.pack_event_groups(
        overlapped, original_groups, anchor=anchor, window=(2, 22),
    )
    expected = np.zeros_like(rain)
    expected[0], expected[-1] = 9, 8
    for (a, b), hyetograph in zip(expected_spans, ([6, 7], [8], [1])):
        expected[a:b] = hyetograph
    # Both placements reuse an original location: every old location must be
    # cleared before any new group is written, without erasing moved rainfall.
    np.testing.assert_array_equal(packed, expected)
    np.testing.assert_array_equal(overlapped, rain)
    assert original_groups == before_groups
    assert not np.shares_memory(packed, overlapped)
    assert packed.sum() == rain.sum()
    assert [(g["start"], g["stop"]) for g in groups] == expected_spans
    assert not groups[-1]["target_exceeded"]
    for original, placed in zip(original_groups, groups):
        assert (placed["original_start"], placed["original_stop"]) == (
            original["start"], original["stop"],
        )
        for key in original.keys() - {"start", "stop"}:
            assert placed[key] == original[key], key


@pytest.mark.parametrize("anchor", [1, 9])
def test_explicit_anchor_can_use_the_full_record_interior_by_default(generator, anchor):
    rain = np.zeros(12)
    rain[5:7] = [6, 7]
    overlapped, groups = generator.overlap_event_groups(rain, start=0, thresholds=HIGH_THRESHOLDS)
    packed, groups = generator.pack_event_groups(overlapped, groups, anchor=anchor)
    expected = np.zeros_like(rain)
    expected[anchor:anchor + 2] = [6, 7]
    np.testing.assert_array_equal(packed, expected)
    assert (groups[0]["start"], groups[0]["stop"]) == (anchor, anchor + 2)


def test_packed_wet_span_can_exactly_fill_window_with_dry_neighbors_outside(generator):
    rain = np.zeros(12)
    rain[5:7] = [6, 7]
    overlapped, groups = generator.overlap_event_groups(rain, start=0, thresholds=HIGH_THRESHOLDS)
    packed, groups = generator.pack_event_groups(overlapped, groups, anchor=2, window=(2, 4))
    np.testing.assert_array_equal(packed, [0, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0])
    assert (groups[0]["start"], groups[0]["stop"]) == (2, 4)


@pytest.mark.parametrize("anchor,window", [
    (True, None), (1.5, None), (-1, None),
    (0, None), (10, None),  # A wet step at a record edge lacks an observed dry neighbor.
    (1, (2, 10)), (9, (2, 10)),  # The complete placed span must fit the chosen window.
    (2, (-1, 12)), (2, (0, 13)), (2, (8, 7)),
    (2, (2.5, 10)), (2, (True, 10)),
])
def test_invalid_anchor_or_window_is_rejected(generator, anchor, window):
    rain = np.zeros(12)
    rain[5:7] = [6, 7]
    overlapped, groups = generator.overlap_event_groups(rain, start=0, thresholds=HIGH_THRESHOLDS)
    with pytest.raises(ValueError):
        generator.pack_event_groups(overlapped, groups, anchor=anchor, window=window)


@pytest.mark.parametrize("unrelated_wet_row", [7, 8, 9, 10])
def test_moving_a_group_cannot_overwrite_or_touch_ungrouped_rain(generator, unrelated_wet_row):
    rain = np.zeros(14)
    rain[2:4] = [6, 7]
    rain[unrelated_wet_row] = 9
    overlapped, groups = generator.overlap_event_groups(
        rain, start=1, stop=5, thresholds=HIGH_THRESHOLDS,
    )
    with pytest.raises(ValueError):
        generator.pack_event_groups(overlapped, groups, anchor=8)
    np.testing.assert_array_equal(overlapped, rain)


@pytest.mark.parametrize("rain", [[], [1, 0, 2]])
def test_packing_no_groups_returns_an_unchanged_copy(generator, rain):
    source = np.array(rain, dtype=float)
    packed, groups = generator.pack_event_groups(source, [])
    np.testing.assert_array_equal(packed, source)
    assert not np.shares_memory(packed, source)
    assert groups == []


def test_explicit_two_day_gap_preserves_two_zero_rows(generator):
    rain = np.array([0, 6, 0, 0, 0, 7, 0, 0, 0, 0, 1, 0], dtype=float)
    overlapped, groups = generator.overlap_event_groups(
        rain, start=0, thresholds={**HIGH_THRESHOLDS, 1: 5},
    )
    packed, groups = generator.pack_event_groups(overlapped, groups, gap_days=2)
    np.testing.assert_array_equal(packed, [0, 6, 0, 0, 7, 0, 0, 1, 0, 0, 0, 0])
    assert [group["start"] for group in groups] == [1, 4, 7]


@pytest.mark.parametrize("gap", [0, -1, 1.5, True])
def test_packing_requires_a_positive_integer_dry_gap(generator, gap):
    rain = np.array([0, 6, 0], dtype=float)
    overlapped, groups = generator.overlap_event_groups(rain, start=0, thresholds=HIGH_THRESHOLDS)
    with pytest.raises(ValueError):
        generator.pack_event_groups(overlapped, groups, gap_days=gap)


def test_baseline_preserves_existing_weather_formulas_on_the_record_substream(generator, monkeypatch):
    original = load_generator(registry.find_probe("mass/catchment-closure"))
    monkeypatch.setattr(original, "N_STEPS", generator.N_STEPS)
    # default_rng accepts a Generator: the original climate code can consume
    # the new record stream without depending on the former unsplit seed.
    expected, expected_static = original.generate(generator._rng(20260912, 0))
    actual, actual_static = generator.generate_baseline(20260912)
    pd.testing.assert_frame_equal(actual, expected)
    assert actual_static == expected_static


def test_record_and_calibration_streams_are_reproducible_and_independent(generator):
    expected = generator._rng(42, 0).random(16)
    calibration = generator._rng(42, 1)
    assert not np.array_equal(expected, calibration.random(16))
    calibration.random(10000)
    np.testing.assert_array_equal(generator._rng(42, 0).random(16), expected)


def annual_rain_fixture(generator, totals):
    """Give each scored block its supplied depth, excluding a wet spinup."""
    rain = np.zeros(generator.N_STEPS)
    rain[0] = 1_000_000
    for i, total in enumerate(totals):
        rain[365 + i * 365 + 1] = total
    return rain


@pytest.mark.parametrize("totals,expected_year,expected_rank", [
    (list(range(1, 21)), 10, 10),
    (list(range(20, 0, -1)), 10, 11),
    ([7] * 20, 1, 1),
])
def test_median_wet_selection_uses_central_pair_and_earliest_tie(
    generator, totals, expected_year, expected_rank,
):
    rain = annual_rain_fixture(generator, totals)
    before = rain.copy()
    selected = generator.select_median_wet_year(rain)
    np.testing.assert_array_equal(rain, before)
    assert selected["year_number"] == expected_year
    assert selected["annual_rank_ascending"] == expected_rank
    assert selected["median_precip_mm"] == np.median(totals)
    assert selected["baseline_precip_mm"] == totals[expected_year - 1]
    assert selected["annual_precip_mm"] == totals
    assert selected["scored_precip_mm"] == sum(totals)
    assert selected["selected_fraction_of_scored_precip"] == totals[expected_year - 1] / sum(totals)
    assert selected["distance_from_median_mm"] == abs(totals[expected_year - 1] - np.median(totals))
    assert (selected["start_row"], selected["stop_row"]) == (365 * expected_year, 365 * (expected_year + 1))
    rain[:365] = 0
    assert generator.select_median_wet_year(rain) == selected


def test_six_decimal_annual_ties_are_exact_even_when_the_earlier_year_is_wetter(generator):
    totals = [11.000001, 10] + list(range(1, 10)) + list(range(12, 21))
    rain = annual_rain_fixture(generator, totals)
    rain[366], rain[368] = 0.1, 10.900001
    selected = generator.select_median_wet_year(rain)
    assert selected["year_number"] == 1
    assert selected["annual_rank_ascending"] == 11
    assert selected["baseline_precip_mm"] == 11.000001
    assert selected["median_precip_mm"] == 10.5000005
    assert selected["distance_from_median_mm"] == 0.5000005


def test_all_dry_scored_years_have_no_defined_precipitation_fraction(generator):
    selected = generator.select_median_wet_year(annual_rain_fixture(generator, [0] * 20))
    assert selected["year_number"] == 1
    assert selected["scored_precip_mm"] == selected["median_precip_mm"] == 0
    assert selected["selected_fraction_of_scored_precip"] is None
    json.dumps(selected, allow_nan=False)


def test_median_wet_year_fraction_is_reported_not_assumed_to_be_five_percent(generator):
    selected = generator.select_median_wet_year(
        annual_rain_fixture(generator, [0] * 9 + [10, 11] + [12] * 9),
    )
    assert selected["year_number"] == 10
    assert selected["selected_fraction_of_scored_precip"] == 10 / 129
    assert selected["selected_fraction_of_scored_precip"] > 0.05


@pytest.mark.parametrize("length", [4015, 7664, 7666])
def test_year_selection_rejects_missing_or_partial_scored_blocks(generator, length):
    with pytest.raises(ValueError, match="full 7665-row"):
        generator.select_median_wet_year(np.zeros(length))


def test_generate_modifies_the_median_year_and_preserves_both_sides(generator, monkeypatch):
    baseline, static = generator.generate_baseline(42)
    baseline["pr"] = annual_rain_fixture(generator, list(range(20)))
    baseline.loc[[3651, 3653, 3655], "pr"] = [2, 3, 4]
    monkeypatch.setattr(
        generator, "generate_baseline", lambda seed: (baseline.copy(deep=True), dict(static)),
    )
    frame, actual_static = generator.generate(42)
    selected = frame.attrs["rainfall_diagnostics"]["selection"]
    assert selected["year_number"] == 10
    assert (selected["start_row"], selected["stop_row"]) == (3650, 4015)
    assert selected["median_precip_mm"] == 9.5
    pd.testing.assert_series_equal(frame["pr"].iloc[:3650], baseline["pr"].iloc[:3650])
    pd.testing.assert_series_equal(frame["pr"].iloc[4015:], baseline["pr"].iloc[4015:])
    assert frame.loc[3773, "time"] == "2010-05-01"
    assert frame.loc[3773, "pr"] == 9
    assert frame.loc[[3651, 3653, 3655], "pr"].eq(0).all()
    assert frame["pr"].sum() == baseline["pr"].sum()
    assert actual_static == static
    pd.testing.assert_frame_equal(frame[["time", "tas", "pet"]], baseline[["time", "tas", "pet"]])


@pytest.mark.parametrize("record_start,expected_date,expected_anchor", [
    ("2000-01-01", "2010-05-01", 3773),
    ("2001-01-01", "2011-05-01", 3772),
])
def test_may_first_anchor_uses_actual_calendar_dates_not_a_fixed_row_offset(
    generator, monkeypatch, record_start, expected_date, expected_anchor,
):
    baseline, static = generator.generate_baseline(42)
    baseline["time"] = pd.date_range(record_start, periods=len(baseline)).strftime("%Y-%m-%d")
    baseline["pr"] = annual_rain_fixture(generator, list(range(1, 21)))
    baseline.loc[[3651, 3653, 3655], "pr"] = [2, 3, 5]
    monkeypatch.setattr(
        generator, "generate_baseline", lambda seed: (baseline.copy(deep=True), dict(static)),
    )
    frame, _ = generator.generate(42)
    construction = frame.attrs["rainfall_diagnostics"]["construction"]
    group = frame.attrs["rainfall_diagnostics"]["groups"][0]
    assert construction["anchor_row"] == expected_anchor
    assert construction["anchor_time"] == expected_date
    assert group["start"] == expected_anchor
    assert frame.loc[expected_anchor, "time"] == expected_date
    assert frame.loc[expected_anchor, "pr"] == 10
    assert expected_anchor != 3650 + 120
    assert frame.loc[[3651, 3653, 3655], "pr"].eq(0).all()
    assert frame["pr"].sum() == baseline["pr"].sum()


@pytest.mark.parametrize("seed", [7, 42, 20260912])
def test_generated_groups_use_own_seed_q100_and_preserve_weather_and_volume(generator, seed):
    baseline, baseline_static = generator.generate_baseline(seed)
    frame, static = generator.generate(seed)
    repeated, repeated_static = generator.generate(seed)
    assert len(frame) == 7665
    assert generator.PERIOD_YEARS == 20
    assert generator.TARGET_RETURN_PERIOD == 100
    assert generator.INTER_EVENT_DRY_DAYS == 1
    assert generator.STORM_START_MONTH_DAY == (5, 1)
    assert frame.to_csv(index=False) == repeated.to_csv(index=False)
    diagnostics = frame.attrs["rainfall_diagnostics"]
    assert json.dumps(diagnostics, allow_nan=False, sort_keys=True) == json.dumps(
        repeated.attrs["rainfall_diagnostics"], allow_nan=False, sort_keys=True,
    )
    assert static == baseline_static == repeated_static
    pd.testing.assert_frame_equal(frame[["time", "tas", "pet"]], baseline[["time", "tas", "pet"]])
    selected = diagnostics["selection"]
    assert 1 <= selected["year_number"] <= 20
    start, stop = selected["start_row"], selected["stop_row"]
    assert (start, stop) == (365 * selected["year_number"], 365 * (selected["year_number"] + 1))
    assert pd.Timestamp(selected["start_time"]) == pd.Timestamp(baseline["time"].iloc[start])
    assert pd.Timestamp(selected["end_time"]) == pd.Timestamp(baseline["time"].iloc[stop - 1])
    assert "modified_year" in diagnostics
    assert diagnostics["calibration"]["years"] == 100
    construction = diagnostics["construction"]
    assert construction["target_return_period_years"] == 100
    assert construction["inter_event_dry_days"] == 1
    anchor = construction["anchor_row"]
    assert start <= anchor < stop
    assert construction["anchor_time"] == baseline["time"].iloc[anchor]
    assert pd.Timestamp(construction["anchor_time"]).strftime("%m-%d") == "05-01"
    thresholds = {int(d): value for d, value in construction["thresholds_mm"].items()}
    assert set(thresholds) == {1, 3, 7}
    for d, value in thresholds.items():
        assert value == diagnostics["calibration"]["return_levels_mm"][str(d)]["100"]
    pd.testing.assert_series_equal(frame["pr"].iloc[:start], baseline["pr"].iloc[:start])
    pd.testing.assert_series_equal(frame["pr"].iloc[stop:], baseline["pr"].iloc[stop:])
    for a, b in [(0, len(frame)), (start, stop)]:
        assert frame["pr"].iloc[a:b].sum() == pytest.approx(
            baseline["pr"].iloc[a:b].sum(), abs=1e-9, rel=0,
        )
    assert list(frame.columns) == ["time", "pr", "tas", "pet"]
    assert np.isfinite(frame[["pr", "tas", "pet"]].to_numpy()).all()
    assert (frame[["pr", "pet"]].to_numpy() >= 0).all()

    overlapped, original_groups = generator.overlap_event_groups(
        baseline["pr"].to_numpy(), start=start, stop=stop, thresholds=thresholds,
    )
    expected_rain, groups = generator.pack_event_groups(
        overlapped, original_groups, anchor=anchor, window=(start, stop),
    )
    np.testing.assert_array_equal(frame["pr"], np.round(expected_rain, 6))
    assert diagnostics["group_count"] == len(groups)
    assert diagnostics["target_exceeded_group_count"] == sum(g["target_exceeded"] for g in groups)
    assert len(diagnostics["groups"]) == len(groups)
    assert groups[0]["start"] == anchor
    for original, packed, reported in zip(original_groups, groups, diagnostics["groups"]):
        assert packed["original_start"] == reported["original_start"] == original["start"]
        assert packed["original_stop"] == reported["original_stop"] == original["stop"]
        assert reported["start"] == packed["start"]
        assert reported["stop"] == packed["stop"]
        assert reported["source_events"] == [list(source) for source in original["source_events"]]
        for key in ("event_id", "source_event_count", "depths_mm", "target_exceeded",
                    "exceeded_durations_days", "stop_reason"):
            assert original[key] == packed[key] == reported[key], key
    for previous, following in zip(groups[:-1], groups[1:]):
        assert following["start"] == previous["stop"] + 1
        assert frame["pr"].iloc[previous["stop"]] == 0
    eligible = [
        (a, b) for a, b in generator.rainfall_events(baseline["pr"].to_numpy())
        if start <= a and b <= stop and a > 0 and b < len(baseline)
    ]
    assert [source for group in groups for source in group["source_events"]] == eligible
    summary = generator.event_summary(frame)
    assert not summary.empty
    assert len(summary) == len(groups)
    assert set(summary["event_id"]) == {group["event_id"] for group in groups}
    all_runs = generator.rainfall_events(frame["pr"].to_numpy())
    for group in groups:
        first, last = group["start"], group["stop"] - 1
        event = frame.iloc[first:last + 1]
        assert start <= first <= last < stop
        assert (first, last + 1) in all_runs
        assert frame["pr"].iloc[first - 1] == frame["pr"].iloc[last + 1] == 0
        assert (event["pr"] > 0).all()


def test_summary_reports_daily_depth_and_intensity_of_all_groups(generator):
    frame = pd.DataFrame({
        "time": ["2000-01-01", "2000-01-02", "2000-01-03"],
        "pr": [6.0, 2.0, 0.0],
    })
    frame.attrs["rainfall_diagnostics"] = {"groups": [{"event_id": 1, "start": 0, "stop": 2}]}
    assert generator.event_summary(frame).to_dict("records") == [{
        "event_id": 1, "start": "2000-01-01", "end": "2000-01-02",
        "duration_days": 2, "precip_mm": 8.0,
        "mean_mm_per_day": 4.0, "peak_mm_per_day": 6.0,
    }]
    frame.attrs["rainfall_diagnostics"]["groups"] = []
    assert generator.event_summary(frame).empty


def test_annotations_and_ddf_construction_metadata_stay_off_model_inputs(probe, tmp_path):
    case = build_case(probe, seed=20260912)
    stage(tmp_path, case, probe, registry.find_model("reference_bucket"))
    staged = pd.read_csv(tmp_path / "input" / "forcing.csv")
    assert list(staged.columns) == ["time", "pr", "tas", "pet"]
    assert len(staged) == case.n_steps == 7665
    assert list(case.forcing.columns) == ["time", "pr", "tas", "pet"]
    assert "rainfall_diagnostics" in case.forcing.attrs
    request = json.loads((tmp_path / "request.json").read_text())
    assert "rainfall_diagnostics" not in request
    static = json.loads((tmp_path / "input" / "static.json").read_text())
    assert static == case.static
    assert "rainfall_diagnostics" not in static
    assert "thresholds_mm" not in static


def test_a_short_submitted_model_window_cannot_discard_the_experiment(probe):
    submitted = replace(registry.find_model("reference_bucket"), name="submitted_model", window_days=30)
    assert resolve_window_days(submitted, probe) == 7300


def test_diagnostic_attributes_do_not_create_hidden_window_requirements(probe):
    case = build_case(probe, seed=20260912)
    # Exercise slicing directly: only the explicit probe minimum should
    # require the full experiment, not diagnostic event IDs in forcing.
    sliced = window_case(case, WindowBounds(offset_days=0, days=30))
    assert sliced.n_steps == case.spinup_steps + 30
    assert list(sliced.forcing.columns) == ["time", "pr", "tas", "pet"]

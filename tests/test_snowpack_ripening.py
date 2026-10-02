"""energy/snowpack-ripening: water may not leave a pack that still holds cold content.

The synthetic cases below are one accumulation stage and one dry melt block.
Rates are in mm/day and states in mm at every step, as the contract has them,
so the same physical pack can be written at a daily or an hourly step.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get
from hydroturing.criteria.base import FAIL, PASS, CriterionIncompatibleError
from hydroturing.harness import compatibility_issues, run_probe
from hydroturing.protocol import Case, RunResult
from hydroturing.scoring import FAIL as PROBE_FAIL, PASS as PROBE_PASS
from hydroturing.seeds import gate_seeds

C_ICE = 2100.0
PROBE_ID = "energy/snowpack-ripening"
ACCUMULATION, MELT = 10, 40
PACK = 500.0            # mm, above the 300 mm floor
OPENING_K = 15.0        # the pack opens 15 K below freezing


def _params(**overrides):
    params = {
        "threshold": 0.05,
        "min_peak_pack_mm": 300.0,
        "cold_content_tolerance_j": 1.0,
        "min_opening_depression_k": 2.0,
        "sensible_heat_coefficient": 10.0,
        "min_opening_ice_share": 0.5,
        "pack_temperature_margin_k": 2.0,
        "max_ground_melt_mm_per_day": 0.3,
    }
    params.update(overrides)
    return params


def _pack(ripe_on=12, drain=15.0, leak=0.0, leak_days=None, opening_k=OPENING_K,
          steps_per_day=1):
    """A pack that warms steadily, ripens at the end of day `ripe_on`, then drains.

    Built at its own step, so the warming is spread evenly across the day rather
    than arriving in one step. `leak` is water leaving while the pack is still
    cold, in mm/day, for the first `leak_days` days of the block (all its cold
    days when None).
    """
    k = steps_per_day
    n_acc, n_melt = ACCUMULATION * k, MELT * k
    c0 = C_ICE * PACK * opening_k
    csnow = np.full(n_acc + n_melt, c0)
    snm = np.zeros(n_acc + n_melt)
    snw = np.full(n_acc + n_melt, PACK)
    for j in range(n_melt):
        t = n_acc + j
        csnow[t] = max(0.0, c0 * (1.0 - (j + 1) / (ripe_on * k)))
        if csnow[t] > 0.0:
            snm[t] = leak if leak_days is None or j < leak_days * k else 0.0
        else:
            snm[t] = drain
        snw[t] = snw[t - 1] - snm[t] / k
    tas = np.concatenate([np.full(n_acc, -15.0), np.linspace(-15.0, 6.0, n_melt)])
    return {
        "snm": snm, "snw": snw, "csnow": csnow, "tas": tas,
        "rn": np.full(n_acc + n_melt, 60.0),
        "regime": np.array(["accumulation"] * n_acc + ["melt"] * n_melt),
        "pr": np.zeros(n_acc + n_melt),
    }


def _run(series, timestep="PT1D"):
    n = len(series["snm"])
    freq = "D" if timestep == "PT1D" else "h"
    time = pd.date_range("2001-01-01", periods=n, freq=freq)
    forcing = pd.DataFrame({
        "time": time, "pr": series["pr"], "tas": series["tas"], "rn": series["rn"],
        "_regime": series["regime"],
    })
    case = Case(probe_id=PROBE_ID, seed=1, forcing=forcing, static={},
                spinup_steps=0, timestep=timestep)
    table = pd.DataFrame({
        "time": time, "snm": series["snm"], "snw": series["snw"], "csnow": series["csnow"],
    })
    if "lwsnl" in series:
        table["lwsnl"] = series["lwsnl"]
    return RunResult(case, table, {}, 0.0)


def _score(series, timestep="PT1D", **overrides):
    return get("snowpack_ripening")(_run(series, timestep), None, _params(**overrides))


# --- the rule -------------------------------------------------------------

def test_a_pack_that_ripens_before_it_drains_passes():
    result = _score(_pack())
    assert result.status == PASS
    assert result.value == 0.0
    assert result.diagnostics["blocks"][0]["cold_steps"] == 11


def test_water_leaving_a_cold_pack_fails():
    # 3 mm/day through the eleven cold days is 33 mm, 6.6% of the 500 mm pack.
    result = _score(_pack(leak=3.0))
    assert result.status == FAIL
    assert result.value == pytest.approx(33.0 / PACK)
    assert "drained from a pack that had not ripened" in result.message


def test_ground_melt_inside_the_allowance_passes():
    # Snow-17's own daygm: 0.1 mm/day for eleven cold days, 0.22% of the pack.
    result = _score(_pack(leak=0.1))
    assert result.status == PASS
    assert result.value == pytest.approx(1.1 / PACK)


def test_the_share_does_not_depend_on_the_step():
    # snm is a rate; summing it without the step would inflate an hourly record
    # 24-fold, the fault the snowpack-mass-closure review caught in #102. The leak
    # sits on the first five days, well inside the cold stretch, because the step
    # in which the pack ripens is excused whole: a day at a daily step, an hour at
    # an hourly one, which is the one place the two may differ.
    daily = _score(_pack(leak=6.0, leak_days=5), "PT1D")
    hourly = _score(_pack(leak=6.0, leak_days=5, steps_per_day=24), "PT1H")
    assert daily.value == pytest.approx(30.0 / PACK)
    assert hourly.value == pytest.approx(daily.value)
    assert daily.status == hourly.status == FAIL


# --- what stops a model emptying the check ---------------------------------

def test_a_pack_reported_ripe_on_the_opening_row_fails():
    series = _pack(leak=3.0)
    series["csnow"][:] = 0.0
    result = _score(series)
    assert result.status == FAIL
    assert result.diagnostics["outcome"] == "opens_ripe"


def test_a_pack_opening_just_below_the_depression_fails_and_just_above_passes():
    assert _score(_pack(opening_k=1.9)).diagnostics.get("outcome") == "opens_ripe"
    assert _score(_pack(opening_k=2.1)).status == PASS


def test_a_pack_reported_liquid_on_the_opening_row_fails():
    # The degree-day leak with no cold content, and the whole pack declared
    # liquid on the row before the block: zero ice would otherwise skip the
    # opening check, and nothing that drains would be counted.
    series = _pack(leak=3.0)
    series["csnow"][:] = 0.0
    series["lwsnl"] = np.zeros_like(series["snw"])
    series["lwsnl"][ACCUMULATION - 1] = series["snw"][ACCUMULATION - 1]
    result = _score(series)
    assert result.status == FAIL
    assert result.diagnostics["outcome"] == "opens_without_ice"


def test_a_token_pack_on_the_opening_row_cannot_carry_a_token_cold_content():
    # snw dips to 1 mm on that one row with a cold content that reads -20 C over
    # it; afterwards the pack is back and reported ripe while it leaks.
    series = _pack(leak=3.0)
    series["csnow"][:] = 0.0
    series["snw"][ACCUMULATION - 1] = 1.0
    series["csnow"][ACCUMULATION - 1] = C_ICE * 1.0 * 20.0
    assert _score(series).diagnostics["outcome"] == "opens_without_ice"


def test_cold_content_cannot_fall_faster_than_energy_arrives():
    # The whole deficit gone in one day: 15.75 MJ/m2 against about 5 MJ/m2 of
    # radiation and sensible heat from air 15 K warmer than a -30 C mean pack.
    series = _pack(ripe_on=1)
    result = _score(series)
    assert result.status == FAIL
    assert result.diagnostics["outcome"] == "ripens_too_fast"


def _written_off(series, row):
    """Cold content reported zero from `row` on, while the pack keeps leaking."""
    series["csnow"][row:] = 0.0
    return _score(series)


@pytest.mark.parametrize("construction", ["liquid_one_row", "snw_dip", "liquid_throughout"])
def test_a_reported_fall_in_ice_does_not_write_the_cold_content_off(construction):
    # Each makes the pack's reported ice fall on the block's first row, so a
    # credit for "ice that left" would excuse the whole deficit vanishing there.
    # Ice that melts is at 0 C and takes no cold content with it.
    series = _pack(leak=3.0)
    t = ACCUMULATION
    if construction == "liquid_one_row":
        series["lwsnl"] = np.zeros_like(series["snw"])
        series["lwsnl"][t] = series["snw"][t]
    elif construction == "snw_dip":
        series["snw"][t] = 1.0
    else:
        series["lwsnl"] = np.zeros_like(series["snw"])
        series["lwsnl"][t:] = series["snw"][t:]
    assert _written_off(series, t).diagnostics["outcome"] == "ripens_too_fast"


def test_a_pack_reported_liquid_cannot_read_as_impossibly_cold():
    # One row with almost no ice and its cold content kept makes the pack read
    # millions of kelvin below freezing on the next step, and the sensible heat
    # that follows would cover any fall. Floored at the coldest air less 2 K,
    # the supply cannot.
    series = _pack(leak=3.0)
    t = ACCUMULATION
    series["lwsnl"] = np.zeros_like(series["snw"])
    series["lwsnl"][t] = series["snw"][t] - 1e-7
    assert _written_off(series, t + 1).diagnostics["outcome"] == "ripens_too_fast"


def test_negative_cold_content_and_negative_outflow_fail():
    series = _pack()
    series["csnow"][ACCUMULATION + 2] = -5.0e5
    assert _score(series).diagnostics["outcome"] == "negative_cold_content"
    series = _pack(leak=3.0)
    series["snm"][ACCUMULATION + 3] = -40.0
    assert _score(series).diagnostics["outcome"] == "negative_outflow"


# --- what the criterion refuses -------------------------------------------

def test_a_thin_pack_is_refused_rather_than_scored():
    series = _pack(leak=0.1)
    series["snw"] = series["snw"] * 0.2          # 100 mm: ground melt would read as a leak
    with pytest.raises(CriterionIncompatibleError, match="300 mm"):
        _score(series)


def test_a_thin_pack_that_leaks_beyond_ground_melt_is_still_failed():
    # The same pack at a fifth of the size, leaking 3 mm/day while cold: 33 mm
    # against the 3.3 mm that 0.3 mm/day of ground melt over eleven cold days
    # could explain. The floor exists for ground melt, not for this.
    series = _pack(leak=3.0)
    series["snw"] = series["snw"] * 0.2
    series["csnow"] = series["csnow"] * 0.2
    result = _score(series)
    assert result.status == FAIL
    assert result.value == pytest.approx(33.0 / (PACK * 0.2))


def test_a_wet_block_is_refused():
    series = _pack()
    series["pr"][ACCUMULATION + 5] = 4.0
    with pytest.raises(ValueError, match="has to be dry"):
        _score(series)


def test_unknown_and_non_finite_parameters_are_rejected():
    with pytest.raises(ValueError, match="unknown parameters"):
        _score(_pack(), min_peak_pack=300.0)
    with pytest.raises(ValueError, match="finite, non-negative"):
        _score(_pack(), threshold=float("nan"))


# --- the probe --------------------------------------------------------------

def test_the_probe_requires_no_forcing_so_temperature_index_models_are_scored():
    probe = registry.find_probe(PROBE_ID)
    assert probe.requires_forcing == ()
    assert compatibility_issues(registry.find_model("sacsma_snow17"), probe) == []


@pytest.mark.parametrize(
    ("model_name", "verdict", "outcome"),
    [
        ("reference_snow_energy", PROBE_PASS, None),
        ("sacsma_snow17", PROBE_PASS, None),
        ("reference_degree_day", PROBE_FAIL, None),
        ("reference_always_ripe", PROBE_FAIL, "opens_ripe"),
    ],
)
def test_each_baseline_is_decided_by_the_rule_it_exists_for(model_name, verdict, outcome, tmp_path):
    probe = registry.find_probe(PROBE_ID)
    report = run_probe(registry.find_model(model_name), probe, gate_seeds(PROBE_ID, 1), workdir=tmp_path)
    assert report.verdict == verdict
    result = next(c for c in report.criteria if c.name == "snowpack_ripening")
    # The degree-day pack is caught by the rule itself, on water leaving it while
    # cold; the always-ripe pack by the opening check, having emptied the rule.
    assert (result.diagnostics or {}).get("outcome") == outcome

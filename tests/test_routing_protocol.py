"""Reach-indexed routing declarations and /io tables obey the contract."""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest
import yaml

from hydroturing import registry
from hydroturing.harness import verify_adapter_contract
from hydroturing.protocol import Case, ProtocolError, read_result, stage
from hydroturing.spec import SpecError, load_model, load_probe


@pytest.fixture
def base_probe():
    return registry.find_probe("mass/catchment-closure")


@pytest.fixture
def routing_probe(base_probe):
    return replace(
        base_probe,
        requires_fluxes=(),
        requires_states=(),
        requires_diagnostics=(),
        requires_routing=("q_in", "q_out", "channel_storage"),
        requires_static=("routing_network",),
        variants=(),
    )


@pytest.fixture
def case():
    return Case(
        probe_id="mass/routing-contract-test",
        seed=11,
        forcing=pd.DataFrame(
            {
                "time": pd.date_range("2000-01-01", periods=3, freq="D"),
                "pr": [0.0, 1.0, 0.0],
            }
        ),
        static={"routing_network": {"reaches": ["A", "B"]}},
        spinup_steps=1,
    )


def _write_probe(path, base_probe, routing):
    path.mkdir(parents=True)
    raw = yaml.safe_load((base_probe.path / "probe.yaml").read_text())
    raw["id"] = "mass/routing-contract-test"
    raw["requires"] = {"routing": routing, "static": ["routing_network"]}
    (path / "generate.py").write_text("# Only manifest loading is tested.\n")
    (path / "probe.yaml").write_text(yaml.safe_dump(raw))


def _write_model(path, routing):
    path.mkdir()
    manifest = {
        "name": path.name,
        "version": "1",
        "entrypoint": ["python", "adapter.py"],
        "timestep": "PT1D",
        "emits": {"fluxes": [], "states": [], "routing": routing},
    }
    (path / "model.yaml").write_text(yaml.safe_dump(manifest))


def _routing_table(case):
    rows = []
    for time in case.forcing["time"]:
        rows.extend(
            [
                {
                    "time": time.strftime("%Y-%m-%d"),
                    "reach_id": "A",
                    "q_in": 2.0,
                    "q_out": 1.0,
                    "channel_storage": 86400.0,
                },
                {
                    "time": time.strftime("%Y-%m-%d"),
                    "reach_id": "B",
                    "q_in": 3.0,
                    "q_out": 2.0,
                    "channel_storage": 86400.0,
                },
            ]
        )
    return pd.DataFrame(rows)


def _write_outputs(path, case, routing=None):
    output = path / "output"
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"time": case.forcing["time"]}).to_csv(
        output / "result.csv", index=False
    )
    if routing is not None:
        routing.to_csv(output / "routing.csv", index=False)
    (output / "run.json").write_text(json.dumps({"status": "ok"}))


def test_routing_declarations_load_and_old_defaults_remain_empty(
    base_probe, tmp_path
):
    probe_dir = tmp_path / "mass" / "routing-contract-test"
    variables = ["q_in", "q_out", "channel_storage"]
    _write_probe(probe_dir, base_probe, variables)
    probe = load_probe(probe_dir)

    model_dir = tmp_path / "routing_model"
    _write_model(model_dir, variables)
    model = load_model(model_dir)

    assert probe.requires_routing == model.emits_routing == tuple(variables)
    assert model.missing_for(probe) == []
    assert base_probe.requires_routing == ()
    assert registry.find_model("reference_bucket").emits_routing == ()



def test_routing_probe_must_require_network(base_probe, tmp_path):
    path = tmp_path / "mass" / "routing-contract-test"
    _write_probe(path, base_probe, ["q_out"])
    raw = yaml.safe_load((path / "probe.yaml").read_text())
    raw["requires"].pop("static")
    (path / "probe.yaml").write_text(yaml.safe_dump(raw))

    with pytest.raises(SpecError, match="routing_network"):
        load_probe(path)


@pytest.mark.parametrize("kind", ["probe", "model"])
def test_unknown_routing_variable_is_rejected(base_probe, tmp_path, kind):
    if kind == "probe":
        path = tmp_path / "mass" / "routing-contract-test"
        _write_probe(path, base_probe, ["not_a_routing_variable"])
        loader = load_probe
    else:
        path = tmp_path / "routing_model"
        _write_model(path, ["not_a_routing_variable"])
        loader = load_model

    with pytest.raises(SpecError, match="unknown variables"):
        loader(path)


def test_stage_requests_routing_path_and_units(case, routing_probe, tmp_path):
    model = replace(
        registry.find_model("reference_bucket"),
        emits_routing=("q_in", "q_out", "channel_storage"),
    )
    request = json.loads(stage(tmp_path, case, routing_probe, model).read_text())

    assert request["request"]["routing"] == [
        "q_in",
        "q_out",
        "channel_storage",
    ]
    assert request["output"]["routing"] == "output/routing.csv"
    assert request["units"]["q_in"] == "m3 s-1"
    assert request["units"]["q_out"] == "m3 s-1"
    assert request["units"]["channel_storage"] == "m3"


def test_stage_omits_routing_for_nonrouting_probe(
    case, base_probe, tmp_path
):
    model = replace(
        registry.find_model("reference_bucket"),
        emits_routing=("q_in",),
    )
    request = json.loads(
        stage(tmp_path, case, base_probe, model).read_text()
    )

    assert request["request"]["routing"] == []
    assert "q_in" not in request["units"]


def test_valid_routing_table_roundtrips(case, routing_probe, tmp_path):
    expected = _routing_table(case)
    _write_outputs(tmp_path, case, expected)

    result = read_result(tmp_path, case, routing_probe, 0.0)

    assert result.routing is not None
    pd.testing.assert_frame_equal(result.routing, expected)


def test_routing_output_requires_declared_reaches(
    case, routing_probe, tmp_path
):
    case = replace(case, static={})
    _write_outputs(tmp_path, case, _routing_table(case))

    with pytest.raises(ProtocolError, match="routing_network.reaches"):
        read_result(tmp_path, case, routing_probe, 0.0)


@pytest.mark.parametrize(
    ("reaches", "message"),
    [
        ([], "non-empty"),
        (["A", "A"], "duplicates"),
        ([{"id": "A"}, {"id": "B"}], "string IDs"),
    ],
)
def test_invalid_declared_reach_ids_are_rejected(
    case, routing_probe, tmp_path, reaches, message
):
    case = replace(
        case,
        static={"routing_network": {"reaches": reaches}},
    )
    _write_outputs(tmp_path, case, _routing_table(case))

    with pytest.raises(ProtocolError, match=message):
        read_result(tmp_path, case, routing_probe, 0.0)


def test_leading_zero_reach_ids_are_preserved(
    case, routing_probe, tmp_path
):
    case = replace(
        case,
        static={"routing_network": {"reaches": ["007", "012"]}},
    )
    table = _routing_table(case)
    table["reach_id"] = table["reach_id"].replace(
        {"A": "007", "B": "012"}
    )
    _write_outputs(tmp_path, case, table)

    result = read_result(tmp_path, case, routing_probe, 0.0)

    assert result.routing is not None
    assert set(result.routing["reach_id"]) == {"007", "012"}


def test_missing_required_routing_file_is_rejected(case, routing_probe, tmp_path):
    _write_outputs(tmp_path, case)
    with pytest.raises(ProtocolError, match="required output/routing.csv"):
        read_result(tmp_path, case, routing_probe, 0.0)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda table: table.drop(columns="q_out"), "missing columns"),
        (
            lambda table: pd.concat([table, table.iloc[[0]]], ignore_index=True),
            "duplicate",
        ),
        (lambda table: table.loc[table["reach_id"] != "B"], "reach IDs"),
    ],
)
def test_malformed_routing_shape_is_rejected(
    case, routing_probe, tmp_path, mutation, message
):
    _write_outputs(tmp_path, case, mutation(_routing_table(case)))
    with pytest.raises(ProtocolError, match=message):
        read_result(tmp_path, case, routing_probe, 0.0)


def test_routing_time_axis_must_match_forcing(case, routing_probe, tmp_path):
    table = _routing_table(case)
    table.loc[table["reach_id"] == "A", "time"] = pd.date_range(
        "2001-01-01", periods=case.n_steps, freq="D"
    ).strftime("%Y-%m-%d")
    _write_outputs(tmp_path, case, table)

    with pytest.raises(ProtocolError, match="time axis"):
        read_result(tmp_path, case, routing_probe, 0.0)


@pytest.mark.parametrize("bad_value", [np.nan, np.inf, "not-a-number"])
def test_routing_values_must_be_finite(
    case, routing_probe, tmp_path, bad_value
):
    table = _routing_table(case)
    table["q_in"] = table["q_in"].astype(object)
    table.loc[0, "q_in"] = bad_value
    _write_outputs(tmp_path, case, table)

    with pytest.raises(ProtocolError, match="non-finite"):
        read_result(tmp_path, case, routing_probe, 0.0)


def test_verify_adapter_skips_routing_on_nonrouting_probe(
    base_probe, tmp_path
):
    model = replace(
        registry.find_model("reference_bucket"),
        emits_routing=("q_in",),
    )

    result = verify_adapter_contract(
        model, base_probe, 11, workdir=tmp_path
    )

    assert result.routing is None


def test_verify_adapter_checks_routing_on_network_probe(
    routing_probe, tmp_path, monkeypatch
):
    import hydroturing.harness as harness

    original_build_case = harness.build_case

    def build_case_with_network(*args, **kwargs):
        built = original_build_case(*args, **kwargs)
        static = dict(built.static)
        static["routing_network"] = {"reaches": ["A"]}
        return replace(built, static=static)

    monkeypatch.setattr(
        harness,
        "build_case",
        build_case_with_network,
    )
    model = replace(
        registry.find_model("reference_bucket"),
        emits_routing=("q_in",),
    )

    with pytest.raises(ProtocolError, match="routing.csv"):
        verify_adapter_contract(
            model, routing_probe, 11, workdir=tmp_path
        )

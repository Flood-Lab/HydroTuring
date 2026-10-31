"""A model may run a shorter solver step inside a row, under the manifest's record of it."""

import json
from dataclasses import replace

import pytest
import yaml

from hydroturing import harness, protocol, registry, spec

MANIFEST = {
    "name": "stepper",
    "version": "1",
    "entrypoint": ["python3", "ht_adapter.py"],
    "emits": {"fluxes": ["pr", "mrro"], "states": []},
}


def load(tmp_path, timestep, substep):
    path = tmp_path / MANIFEST["name"] / "model.yaml"
    path.parent.mkdir()
    path.write_text(yaml.safe_dump({**MANIFEST, "timestep": timestep, "substep": substep}))
    return spec.load_model(path)


def test_a_coarser_row_may_name_a_shorter_step(tmp_path):
    assert load(tmp_path, ["PT1H", "PT1D"], {"PT1D": "PT1H"}).substeps == {"PT1D": "PT1H"}


@pytest.mark.parametrize(
    ("timestep", "substep", "complaint"),
    [
        (["PT1D", "PT1H"], {"PT1D": "PT1H"}, "native step"),
        (["PT1H"], {"PT1D": "PT1H"}, "does not list"),
        (["PT1H", "PT1D"], {"PT1D": "PT1D"}, "evenly"),
        (["PT1M", "PT1H"], {"PT1H": "PT1D"}, "evenly"),
        (["PT1M", "PT5M"], {"PT5M": "PT15M"}, "evenly"),
    ],
)
def test_a_substep_that_does_not_fit_is_refused(tmp_path, timestep, substep, complaint):
    with pytest.raises(spec.SpecError, match=complaint):
        load(tmp_path, timestep, substep)


def _request(model, probe, tmp_path):
    variant = harness.select_variants(model, probe)[0]
    path = protocol.stage(tmp_path, harness.build_case(probe, 1, variant), probe, model)
    return json.loads(path.read_text())


def test_the_request_carries_it_only_on_that_row_step(tmp_path):
    model = replace(registry.find_model("reference_bucket"), timesteps=("PT1H", "PT1D"), substeps={"PT1D": "PT1H"})
    daily = _request(model, registry.find_probe("mass/catchment-closure"), tmp_path / "day")
    assert daily["timestep"] == "PT1D" and daily["substep"] == "PT1H"
    hourly = _request(model, registry.find_probe("energy/radiation-consistency"), tmp_path / "hour")
    assert hourly["timestep"] == "PT1H" and "substep" not in hourly
    bucket = _request(registry.find_model("reference_bucket"), registry.find_probe("mass/catchment-closure"), tmp_path / "plain")
    assert "substep" not in bucket


def test_a_resolution_probe_never_serves_the_substepped_row():
    probe = registry.find_probe("mass/resolution-invariance")
    model = replace(
        registry.find_model("reference_bucket"), timesteps=("PT1H", "PT1M", "PT1D"), substeps={"PT1D": "PT1H"}
    )
    served = {probe.timestep_for(v) for v in harness.select_variants(model, probe)}
    assert served == {"PT1H", "PT1M"}

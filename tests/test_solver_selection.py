import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from patchsim.core.model import CompartmentalModel, NetworkModel
from patchsim.core.model_runner import Model
from patchsim.core.simulation import (
    MODEL_TEMPLATE_CONFIGS,
    get_config_schema,
    get_run_settings,
    load_config,
    run_simulation,
    setup_simulation,
)
from patchsim.utils.geo import generate_contacts


def _write_project(tmp_path: Path) -> Path:
    patch_file = tmp_path / "patches.csv"
    seed_file = tmp_path / "seeds.csv"
    network_file = tmp_path / "network.csv"
    pd.DataFrame({"patch": ["A", "B"], "population": [100, 200]}).to_csv(patch_file, index=False)
    pd.DataFrame({"patch": ["A", "B"], "S": [99, 200], "I": [1, 0], "R": [0, 0]}).to_csv(seed_file, index=False)
    pd.DataFrame(
        {
            "day": [0, 0, 0, 0],
            "source": ["A", "A", "B", "B"],
            "target": ["A", "B", "A", "B"],
            "weight": [0.8, 0.2, 0.2, 0.8],
        }
    ).to_csv(network_file, index=False)
    config = {
        "ModelName": "solver-test",
        "PatchFile": str(patch_file),
        "SeedFile": str(seed_file),
        "NetworkFile": str(network_file),
        "OutputDir": str(tmp_path / "output"),
        "TMax": 5,
        "compartments": ["S", "I", "R"],
        "Parameters": {"beta": 0.08, "gamma": 0.1},
        "Transitions": {"S -> I": "beta", "I -> R": "gamma * I"},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path


def _run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(["patchsim", *args], capture_output=True, text=True)


def test_schema_declares_solver_and_time_step_defaults():
    properties = get_config_schema()["properties"]
    assert properties["Solver"] == {"type": "string", "enum": ["ode", "discrete"], "default": "ode"}
    assert properties["TimeStep"] == {
        "type": "number",
        "exclusiveMinimum": 0,
        "default": 1.0,
    }


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("Solver", "rk45", "Solver"),
        ("TMax", True, "TMax"),
        ("TMax", 1.5, "TMax"),
        ("TMax", 0, "TMax"),
        ("TimeStep", True, "TimeStep"),
        ("TimeStep", 0, "TimeStep"),
        ("TimeStep", float("nan"), "TimeStep"),
        ("TimeStep", float("inf"), "TimeStep"),
    ],
)
def test_run_settings_reject_invalid_values(field, value, message):
    config = {"TMax": 5, field: value}
    with pytest.raises(ValueError, match=message):
        get_run_settings(config)


@pytest.mark.parametrize("command", ["validate", "run"])
def test_cli_rejects_run_settings_before_creating_outputs(tmp_path, command):
    config_path = _write_project(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["TimeStep"] = 0
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    result = _run_cli(command, "-c", str(config_path))

    assert result.returncode != 0
    assert "TimeStep" in result.stderr
    assert not Path(config["OutputDir"]).exists()


def test_absent_and_explicit_ode_config_are_equivalent(tmp_path):
    config_path = _write_project(tmp_path)

    default_run = _run_cli("run", "-c", str(config_path), "--json")
    assert default_run.returncode == 0, default_run.stderr
    default_summary = json.loads(default_run.stdout)
    default_frame = pd.read_csv(default_summary["csv_path"])

    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["Solver"] = "ode"
    config["TimeStep"] = 1.0
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    explicit_run = _run_cli("run", "-c", str(config_path), "--json")
    assert explicit_run.returncode == 0, explicit_run.stderr
    explicit_summary = json.loads(explicit_run.stdout)
    explicit_frame = pd.read_csv(explicit_summary["csv_path"])

    assert default_summary["solver"] == explicit_summary["solver"] == "ode"
    assert default_summary["time_step"] == explicit_summary["time_step"] == 1.0
    assert default_summary["csv_path"] == explicit_summary["csv_path"]
    assert default_summary["plot_path"] == explicit_summary["plot_path"]
    pd.testing.assert_frame_equal(default_frame, explicit_frame)


def test_discrete_cli_reports_method_and_time_grid(tmp_path):
    config_path = _write_project(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config.update({"Solver": "discrete", "TimeStep": 0.5, "TMax": 5})
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    validation = _run_cli("validate", "-c", str(config_path), "--json")
    run = _run_cli("run", "-c", str(config_path), "--json")

    assert validation.returncode == 0, validation.stderr
    assert json.loads(validation.stdout)["solver"] == "discrete"
    assert json.loads(validation.stdout)["time_step"] == 0.5
    assert run.returncode == 0, run.stderr
    summary = json.loads(run.stdout)
    assert summary["solver"] == "discrete"
    assert summary["time_step"] == 0.5
    assert summary["csv_path"].endswith("_discrete.csv")
    assert summary["plot_path"].endswith("_discrete.png")
    frame = pd.read_csv(summary["csv_path"])
    assert frame["time"].tolist() == [0.0, 0.5, 1.0, 1.5, 2.0]
    assert frame.iloc[0][["S_0", "I_0", "R_0"]].tolist() == [99.0, 1.0, 0.0]
    log_text = next((Path(config["OutputDir"]) / "logs").glob("*.log")).read_text(encoding="utf-8")
    assert "Solver: discrete" in log_text
    assert "TimeStep: 0.5" in log_text


@pytest.mark.parametrize("template_name", sorted(MODEL_TEMPLATE_CONFIGS))
def test_shared_ode_derivatives_match_current_runner(template_name):
    template = MODEL_TEMPLATE_CONFIGS[template_name]
    transitions = [{"transition": key.replace(" ", ""), "rate": rate} for key, rate in template["Transitions"].items()]
    base = CompartmentalModel(template["compartments"], template["Parameters"], transitions)
    network = NetworkModel(base, 2, [[0.8, 0.2], [0.2, 0.8]])
    y0 = {f"{compartment}_{patch}": 0.0 for patch in range(2) for compartment in template["compartments"]}
    y0["S_0"] = 99.0
    y0["I_0"] = 1.0
    y0["S_1"] = 100.0
    times = np.arange(20, dtype=float)

    current = Model(network, template["compartments"]).solve(y0, times)
    _times, shared = network.simulate_ode(y0, times)

    for compartment in current:
        np.testing.assert_allclose(current[compartment], shared[compartment], rtol=1e-12, atol=1e-12)


def test_patch_parameters_affect_both_solver_paths():
    base = CompartmentalModel(
        ["I", "R"],
        {"gamma": 0.1},
        [{"transition": "I->R", "rate": "gamma * I"}],
    )
    network = NetworkModel(base, 2, [[0.0, 0.0], [0.0, 0.0]])
    network.patch_names = ["A", "B"]
    network.patch_parameters = {"A": {"gamma": 0.1}, "B": {"gamma": 0.2}}
    y0 = {"I_0": 100.0, "R_0": 0.0, "I_1": 100.0, "R_1": 0.0}
    times = [0.0, 0.5]

    discrete = network.simulate_discrete(y0, times)
    _times, ode = network.simulate_ode(y0, times)

    assert discrete["I_0"][-1] == pytest.approx(95.0)
    assert discrete["I_1"][-1] == pytest.approx(90.0)
    assert ode["I_1"][-1] < ode["I_0"][-1]


@pytest.mark.parametrize("solver", ["ode", "discrete"])
def test_both_solvers_consume_generated_contact_network(tmp_path, solver):
    regions = tmp_path / "regions.csv"
    network_file = tmp_path / "contacts.csv"
    regions.write_text("patch,lat,lon\nA,0,0\nB,0,1\n", encoding="utf-8")
    generate_contacts(
        regions,
        network_file,
        id_column="patch",
        kernel="distance",
        decay=2.0,
        min_distance_km=0.001,
        normalize="row",
        self_share=0.8,
    )
    config_path = _write_project(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config.update(
        {
            "NetworkFile": str(network_file),
            "OutputDir": str(tmp_path / f"output-{solver}"),
            "Solver": solver,
            "TimeStep": 0.5,
        }
    )
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    loaded = load_config(str(config_path))
    net, y0, patches, num_patches = setup_simulation(loaded)
    summary = run_simulation(loaded, loaded["ModelName"], net, y0, patches, num_patches)

    assert Path(summary["csv_path"]).is_file()
    assert summary["solver"] == solver

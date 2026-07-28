import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from patchsim.core.model import CompartmentalModel, NetworkModel
from patchsim.core.simulation import load_config, setup_simulation


def _write_grouped_project(tmp_path: Path) -> Path:
    pd.DataFrame(
        {
            "patch": ["A", "B"],
            "population": [200, 100],
        }
    ).to_csv(tmp_path / "patches.csv", index=False)
    pd.DataFrame(
        {
            "patch": ["A", "A", "B", "B"],
            "group": ["low", "high", "high", "low"],
            "population": [100, 100, 0, 100],
        }
    ).to_csv(tmp_path / "groups.csv", index=False)
    pd.DataFrame(
        {
            "patch": ["A", "A", "B", "B"],
            "group": ["low", "high", "low", "high"],
            "S": [99, 100, 100, 0],
            "I": [1, 0, 0, 0],
            "R": [0, 0, 0, 0],
        }
    ).to_csv(tmp_path / "seeds.csv", index=False)
    pd.DataFrame(
        {
            "day": [0, 0, 0, 0],
            "source": ["A", "A", "B", "B"],
            "target": ["A", "B", "A", "B"],
            "weight": [0.8, 0.2, 0.2, 0.8],
        }
    ).to_csv(tmp_path / "network.csv", index=False)
    pd.DataFrame(
        {
            "focal_group": ["low", "low", "high", "high"],
            "contributor_group": ["low", "high", "low", "high"],
            "weight": [2.0, 1.0, 3.0, 4.0],
        }
    ).to_csv(tmp_path / "interactions.csv", index=False)

    config = {
        "ModelName": "groups",
        "PatchFile": "patches.csv",
        "GroupFile": "groups.csv",
        "InteractionFile": "interactions.csv",
        "InteractionUnits": "contacts/person/day",
        "SeedFile": "seeds.csv",
        "NetworkFile": "network.csv",
        "OutputDir": "output",
        "TMax": 3,
        "Solver": "ode",
        "TimeStep": 1,
        "compartments": ["S", "I", "R"],
        "Parameters": {"beta": 0.08, "gamma": 0.1},
        "Transitions": {"S -> I": "beta", "I -> R": "gamma * I"},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path


def test_setup_loads_generic_groups_and_diagnostics(tmp_path):
    config_path = _write_grouped_project(tmp_path)
    config = load_config(str(config_path))

    net, y0, patches, num_patches = setup_simulation(config)

    assert patches == ["A", "B"]
    assert num_patches == 2
    assert net.groups == ["low", "high"]
    assert y0["I_0_0"] == 1
    assert y0["S_1_0"] == 100
    assert y0["S_1_1"] == 0
    np.testing.assert_allclose(net.interaction, [[2.0, 1.0], [3.0, 4.0]])
    np.testing.assert_allclose(net.compute_force_of_infection(y0), [[0.016, 0.024], [0.004, 0.006]])
    assert net.interaction_diagnostics["units"] == "contacts/person/day"
    assert (
        net.interaction_diagnostics["sha256"]
        == hashlib.sha256((tmp_path / "interactions.csv").read_bytes()).hexdigest()
    )
    assert net.interaction_diagnostics["reciprocity"] == "diagnostic_only"


def test_force_of_infection_uses_focal_rows_and_zero_population_contributes_zero():
    base = CompartmentalModel(
        ["S", "I", "R"],
        {"beta": 1.0, "gamma": 0.1},
        [{"transition": "S->I", "rate": "beta"}, {"transition": "I->R", "rate": "gamma * I"}],
    )
    model = NetworkModel(
        base,
        num_patches=1,
        network_matrix=[[0.0]],
        groups=["low", "high", "absent"],
        interaction_matrix=[
            [2.0, 1.0, 100.0],
            [3.0, 4.0, 100.0],
            [0.0, 0.0, 0.0],
        ],
    )
    state = {
        "S_0_0": 90.0,
        "I_0_0": 10.0,
        "R_0_0": 0.0,
        "S_0_1": 80.0,
        "I_0_1": 20.0,
        "R_0_1": 0.0,
        "S_0_2": 0.0,
        "I_0_2": 0.0,
        "R_0_2": 0.0,
    }

    force = model.compute_force_of_infection(state)

    np.testing.assert_allclose(force, [[0.4, 1.1, 0.0]])
    derivatives = model.compute_derivatives(state)
    assert derivatives["S_0_0"] == pytest.approx(-36.0)
    assert derivatives["S_0_1"] == pytest.approx(-88.0)
    assert all(np.isfinite(value) for value in derivatives.values())


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (
            lambda path: pd.DataFrame(
                {
                    "patch": ["A", "A", "B"],
                    "group": ["low", "high", "low"],
                    "population": [100, 100, 100],
                }
            ).to_csv(path / "groups.csv", index=False),
            "groups for patch 'B' do not match",
        ),
        (
            lambda path: pd.DataFrame(
                {
                    "focal_group": ["low", "high"],
                    "contributor_group": ["low", "low"],
                    "weight": [1, 1],
                }
            ).to_csv(path / "interactions.csv", index=False),
            "must mention every group in contributor_group",
        ),
        (
            lambda path: pd.DataFrame(
                {
                    "focal_group": ["low", "low", "high", "other"],
                    "contributor_group": ["low", "high", "low", "high"],
                    "weight": [1, 1, 1, 1],
                }
            ).to_csv(path / "interactions.csv", index=False),
            "unknown groups",
        ),
    ],
)
def test_group_input_coverage_is_validated(tmp_path, mutate, message):
    config_path = _write_grouped_project(tmp_path)
    mutate(tmp_path)

    with pytest.raises(ValueError, match=message):
        setup_simulation(load_config(str(config_path)))


def test_group_fields_are_required_together(tmp_path):
    config_path = _write_grouped_project(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    del config["InteractionUnits"]
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(ValueError, match="require.*together"):
        setup_simulation(load_config(str(config_path)))


@pytest.mark.parametrize("solver", ["ode", "discrete"])
def test_grouped_cli_validate_and_run(tmp_path, solver):
    config_path = _write_grouped_project(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["Solver"] = solver
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    validate = subprocess.run(
        ["patchsim", "validate", "-c", str(config_path), "--json"],
        capture_output=True,
        text=True,
    )
    assert validate.returncode == 0, validate.stderr
    validation = json.loads(validate.stdout)
    assert validation["groups"] == ["low", "high"]
    assert validation["num_groups"] == 2
    assert validation["interaction"]["units"] == "contacts/person/day"

    run = subprocess.run(
        ["patchsim", "run", "-c", str(config_path), "--json"],
        capture_output=True,
        text=True,
    )
    assert run.returncode == 0, run.stderr
    summary = json.loads(run.stdout)
    assert summary["groups"] == ["low", "high"]
    output = pd.read_csv(summary["csv_path"])
    assert output.columns.tolist() == [
        "time",
        "S_0_0",
        "I_0_0",
        "R_0_0",
        "S_0_1",
        "I_0_1",
        "R_0_1",
        "S_1_0",
        "I_1_0",
        "R_1_0",
        "S_1_1",
        "I_1_1",
        "R_1_1",
    ]
    for patch_idx, group_totals in enumerate(([100, 100], [100, 0])):
        for group_idx, expected_total in enumerate(group_totals):
            columns = [f"{compartment}_{patch_idx}_{group_idx}" for compartment in ("S", "I", "R")]
            np.testing.assert_allclose(output[columns].sum(axis=1), expected_total, atol=1e-8)
    assert Path(summary["plot_path"]).is_file()

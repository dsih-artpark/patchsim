from pathlib import Path

import pandas as pd
import pytest
import yaml

from patchsim.core.simulation import load_config, setup_simulation


def test_setup_simulation_returns_objects(tmp_data_dir):
    cfg = load_config(tmp_data_dir["config"])
    net, y0, patches, num_patches = setup_simulation(cfg)
    # basic assertions:
    assert num_patches == 2
    assert isinstance(y0, dict)
    # check that y0 keys include S_0 and I_1 etc.
    assert "S_0" in y0 and "I_1" in y0
    assert patches == ["A", "B"]


def test_setup_simulation_population_check(tmp_data_dir):
    # mutate seed so it doesn't sum to population to trigger error
    seed_df = pd.read_csv(tmp_data_dir["seed_csv"])
    seed_df.loc[0, "S"] = 500  # break conservation for PatchA
    seed_df.to_csv(tmp_data_dir["seed_csv"], index=False)
    cfg = load_config(tmp_data_dir["config"])
    with pytest.raises(ValueError):
        setup_simulation(cfg)


def test_setup_simulation_rejects_list_transition_format(tmp_data_dir):
    with open(tmp_data_dir["config"], "r") as f:
        cfg = yaml.safe_load(f)

    cfg["Transitions"] = [
        {"from": "S", "to": "I", "rate": "beta"},
        {"from": "I", "to": "R", "rate": "gamma * I"},
    ]

    with open(tmp_data_dir["config"], "w") as f:
        yaml.safe_dump(cfg, f)

    cfg = load_config(tmp_data_dir["config"])
    with pytest.raises(ValueError):
        setup_simulation(cfg)


def test_load_config_resolves_paths_relative_to_config(tmp_path):
    cfg_dir = tmp_path / "proj"
    (cfg_dir / "data" / "patch").mkdir(parents=True)
    (cfg_dir / "data" / "seeds").mkdir(parents=True)
    (cfg_dir / "output").mkdir(parents=True)

    patch_csv = cfg_dir / "data" / "patch" / "patch.csv"
    seed_csv = cfg_dir / "data" / "seeds" / "seed.csv"

    pd.DataFrame({"patch": ["A"], "Population": [100]}).to_csv(patch_csv, index=False)
    pd.DataFrame({"patch": ["A"], "S": [99], "I": [1], "R": [0]}).to_csv(seed_csv, index=False)

    config_path = cfg_dir / "config.yaml"
    with open(config_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(
            {
                "PatchFile": "data/patch/patch.csv",
                "SeedFile": "data/seeds/seed.csv",
                "OutputDir": "output",
                "TMax": 5,
                "compartments": ["S", "I", "R"],
                "Parameters": {"beta": 0.2, "gamma": 0.1},
                "Transitions": {"S -> I": "beta", "I -> R": "gamma * I"},
            },
            f,
        )

    cfg = load_config(str(config_path))
    assert Path(cfg["PatchFile"]).is_absolute()
    assert Path(cfg["SeedFile"]).is_absolute()
    assert Path(cfg["OutputDir"]).is_absolute()


def test_setup_simulation_rejects_unknown_transition_identifiers(tmp_data_dir):
    with open(tmp_data_dir["config"], "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    cfg["Transitions"] = {"S -> I": "betax", "I -> R": "gamma * I"}

    with open(tmp_data_dir["config"], "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f)

    cfg = load_config(tmp_data_dir["config"])
    with pytest.raises(ValueError, match="unknown names"):
        setup_simulation(cfg)


def test_setup_simulation_rejects_compartment_mismatch(tmp_data_dir):
    with open(tmp_data_dir["config"], "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    cfg["compartments"] = ["S", "I", "R", "E"]

    with open(tmp_data_dir["config"], "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f)

    cfg = load_config(tmp_data_dir["config"])
    with pytest.raises(ValueError, match="Compartment mismatch"):
        setup_simulation(cfg)

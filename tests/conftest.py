import yaml
import pandas as pd
import numpy as np
import os
import pytest

@pytest.fixture
def tmp_data_dir(tmp_path):
    """Create small patch, seed, network CSVs and a config YAML in a temp folder."""
    folder = tmp_path / "data"
    patch_folder = folder / "patch"
    net_folder = folder / "networks"
    seed_folder = folder / "seeds"
    patch_folder.mkdir(parents=True)
    net_folder.mkdir()
    seed_folder.mkdir()

    # sample patch CSV
    patch_df = pd.DataFrame({
        "patch": ["A", "B"],
        "Population": [1000, 800]
    })
    patch_csv = patch_folder / "sample-patch.csv"
    patch_df.to_csv(patch_csv, index=False)

    # sample seed CSV (S, I, R columns)
    seed_df = pd.DataFrame({
        "patch": ["A", "B"],
        "S": [990, 795],
        "I": [10, 5],
        "R": [0, 0]
    })
    seed_csv = seed_folder / "sample-seed.csv"
    seed_df.to_csv(seed_csv, index=False)

    # sample network CSV (day, source, target, weight)
    net_df = pd.DataFrame({
        "day": [0, 0],
        "source": ["A", "B"],
        "target": ["B", "A"],
        "weight": [0.1, 0.1]
    })
    net_csv = net_folder / "sample-net.csv"
    net_df.to_csv(net_csv, index=False)

    # sample yaml config with inline dicts (matching new format)
    config = {
        "PatchFile": str(patch_csv),
        "NetworkFile": str(net_csv),
        "SeedFile": str(seed_csv),
        "OutputDir": str(tmp_path / "output"),
        "TMax": 5,
        "compartments": ["S", "I", "R"],
        "PatchParameters": [
            {"patch": "A", "parameters": {"beta": 0.5, "gamma": 0.1}},
            {"patch": "B", "parameters": {"beta": 0.3, "gamma": 0.08}}
        ],
        "Transitions": [
            {"from": "S", "to": "I", "rate": "beta * S * I / (S + I + R)"},
            {"from": "I", "to": "R", "rate": "gamma * I"}
        ]
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f)

    return {
        "root": tmp_path,
        "patch_csv": str(patch_csv),
        "seed_csv": str(seed_csv),
        "net_csv": str(net_csv),
        "config": str(config_path),
    }

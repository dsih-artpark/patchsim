"""
Utility functions for loading and parsing data files.
"""

import logging
import random
from datetime import datetime

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger("PatchSimLogger")


def read_config(config_path: str) -> dict[str, str]:
    """Read and parse a YAML configuration file.

    Args:
        config_path: Path to the YAML configuration file

    Returns:
        Dictionary containing configuration parameters
    """
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def read_patch_population(file_path: str) -> list[dict[str, int]]:
    """Read and parse patch population data from a CSV file.

    Args:
        file_path: Path to the CSV file containing patch population data

    Returns:
        List of dictionaries containing patch population data
    """
    df = pd.read_csv(file_path)
    return df.to_dict("records")


def read_network(file_path: str) -> list[dict[str, float]]:
    """Read and parse the network file for patch connectivity.

    Args:
        file_path: Path to the CSV file containing network data

    Returns:
        List of dictionaries containing network connections
    """
    df = pd.read_csv(file_path)
    return df.to_dict("records")


def read_seeding_infection(
    file_path: str,
    start_date: datetime,
    end_date: datetime,
) -> list[dict[str, str]]:
    """Read and parse the seeding infection data.

    Args:
        file_path: Path to the CSV file containing seeding data
        start_date: Start date of the simulation period
        end_date: End date of the simulation period

    Returns:
        List of dictionaries containing seeding data
    """
    df = pd.read_csv(file_path)
    df["date"] = pd.to_datetime(df["date"])
    mask = (df["date"] >= start_date) & (df["date"] <= end_date)
    return df[mask].to_dict("records")


def apply_seeding_infections(
    patches: list[dict[str, int]],
    seeds: list[dict[str, str]],
    current_date: datetime,
) -> list[dict[str, int]]:
    """Apply seeding infections to patches based on the current date.

    Args:
        patches: List of dictionaries containing patch data
        seeds: List of dictionaries containing seeding data
        current_date: Current simulation date

    Returns:
        Updated list of patch data with seeding infections applied
    """
    for seed in seeds:
        if pd.to_datetime(seed["date"]) == current_date:
            patch_idx = int(seed["patch"]) - 1
            patches[patch_idx]["infected"] += int(seed["count"])
    return patches


def set_random_seed(*, seed: int) -> None:
    """Set the random seed for reproducibility.

    Args:
        seed (int): Random seed value.
    """
    np.random.seed(seed)
    random.seed(seed)
    logger.info(f"Random seed set to {seed}.")

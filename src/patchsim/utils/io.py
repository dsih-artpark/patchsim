import yaml
import csv
import numpy as np
import random
from datetime import datetime
from patchsim.utils.logger import setup_logger
from typing import List, Dict

logger = setup_logger()

def read_config(config_path: str) -> Dict[str, str]:
    """
    Read and parse a YAML configuration file.

    Args:
        config_path (str): Path to the configuration file.

    Returns:
        dict: Parsed configuration parameters.
    """
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
            logger.info(f"Configuration file '{config_path}' read successfully.")
            return config
    except Exception as e:
        logger.error(f"Failed to read configuration file '{config_path}': {e}")
        raise

def read_patch_population(file_path: str) -> List[Dict[str, int]]:
    """
    Read and parse patch population data from a CSV file.

    Args:
        file_path (str): Path to the CSV file.

    Returns:
        list: List of patches with population data.
    """
    patches = []
    try:
        with open(file_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                patches.append({
                    "region": row["region"],
                    "N0": int(row["population"])
                })
        logger.info(f"Patch population file '{file_path}' read successfully.")
    except Exception as e:
        logger.error(f"Failed to read patch population file '{file_path}': {e}")
        raise
    return patches

def read_network(file_path: str) -> List[Dict[str, float]]:
    """
    Read and parse the network file for patch connectivity.

    Args:
        file_path (str): Path to the network CSV file.

    Returns:
        list: List of edges with weights.
    """
    network = []
    try:
        with open(file_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if "source" not in row or "target" not in row or "weight" not in row:
                    logger.warning(f"Invalid format in network file '{file_path}'.")
                    continue
                network.append({
                    "source": row["source"],
                    "target": row["target"],
                    "weight": float(row["weight"]),
                })
        logger.info(f"Network file '{file_path}' read successfully.")
    except Exception as e:
        logger.error(f"Failed to read network file '{file_path}': {e}")
        raise
    return network

def read_seeding_infection(
    file_path: str, *, start_date: datetime, end_date: datetime
) -> List[Dict[str, str]]:
    """
    Read and parse the seeding infection data.

    Args:
        file_path (str): Path to the seeding infection CSV file.
        start_date (datetime): Start date from the configuration.
        end_date (datetime): End date from the configuration.

    Returns:
        list: List of seeding infection data with patch, date, and seed count.
    """
    seeds = []
    try:
        with open(file_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                seed_date = datetime.strptime(row["date"], "%Y-%m-%d")
                if seed_date < start_date or seed_date > end_date:
                    logger.warning(
                        f"Seeding date {seed_date} for region {row['region']} is outside the simulation period ({start_date} to {end_date})."
                    )
                seeds.append({
                    "region": row["region"],
                    "date": seed_date,
                    "seed_count": int(row["seed_count"]),
                })
        logger.info(f"Seeding infection file '{file_path}' read successfully.")
    except Exception as e:
        logger.error(f"Failed to read seeding infection file '{file_path}': {e}")
        raise
    return seeds

def apply_seeding_infections(
    patches: List[Dict[str, int]], seeds: List[Dict[str, str]], *, current_date: datetime
) -> List[Dict[str, int]]:
    """
    Apply seeding infections to patches based on the current date.

    Args:
        patches (list): List of patches with population data.
        seeds (list): List of seeding infection data with patch, date, and seed count.
        current_date (datetime): The current date for the simulation.

    Returns:
        list: Updated patches with seeding infections applied.
    """
    for seed in seeds:
        if seed["date"] <= current_date:  # Only apply seeds for dates up to the current date
            for patch in patches:
                if patch["region"] == seed["region"]:
                    # Cap infections at the total population
                    seed_count = min(seed["seed_count"], patch["N0"] - patch.get("I0", 0))
                    patch["I0"] = patch.get("I0", 0) + seed_count
                    patch["S0"] = max(0, patch["N0"] - patch["I0"])  # Update susceptibles
    logger.info(f"Seeding infections applied for date {current_date}.")
    return patches

def set_random_seed(*, seed: int) -> None:
    """
    Set the random seed for reproducibility.

    Args:
        seed (int): Random seed value.
    """
    np.random.seed(seed)
    random.seed(seed)
    logger.info(f"Random seed set to {seed}.")
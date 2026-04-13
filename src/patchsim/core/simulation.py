"""
This module is reserved for reusable simulation utilities or wrappers.
ODE and discrete simulation logic is implemented in core/model.py.
"""

import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from patchsim.core.model import CompartmentalModel, NetworkModel
from patchsim.core.model_runner import Model
from patchsim.utils.logger import setup_logger

EPSILON = 1e-6


def load_config(config_path: str) -> dict[str, Any]:
    """Load and validate configuration file."""
    cfg_path = Path(config_path).expanduser().resolve()
    with open(cfg_path, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    if not isinstance(config, dict):
        raise ValueError("Configuration must be a YAML mapping at the top level")

    # Validate required fields
    required_fields = ["PatchFile", "SeedFile", "OutputDir", "Transitions", "TMax"]
    for field in required_fields:
        if field not in config:
            raise ValueError(f"Missing required field '{field}' in config")

    # Resolve relative paths against the config file directory.
    cfg_dir = cfg_path.parent
    for key in ["PatchFile", "SeedFile", "NetworkFile", "OutputDir"]:
        val = config.get(key)
        if isinstance(val, str) and val.strip():
            p = Path(val).expanduser()
            if not p.is_absolute():
                config[key] = str((cfg_dir / p).resolve())

    return config


def setup_simulation(config: dict[str, Any]) -> tuple[NetworkModel, dict[str, float], list, int]:
    """Set up the simulation model and initial conditions."""
    # Load patch data
    patch_df = pd.read_csv(config["PatchFile"])
    patches = patch_df["patch"].tolist()
    populations = patch_df.set_index("patch")["Population"].to_dict()

    for p, pop in populations.items():
        if pop <= 0:
            raise ValueError(f"Population for patch {p} must be positive")

    # Load seed data
    seed_df = pd.read_csv(config["SeedFile"])
    seed_compartments = [col for col in seed_df.columns if col != "patch"]

    configured_compartments = config.get("compartments", config.get("Compartments"))
    if configured_compartments is None:
        compartments = seed_compartments
    else:
        if not isinstance(configured_compartments, list) or not configured_compartments:
            raise ValueError("'compartments' must be a non-empty list when provided")
        compartments = [str(c) for c in configured_compartments]

        missing_in_seed = sorted(set(compartments) - set(seed_compartments))
        extra_in_seed = sorted(set(seed_compartments) - set(compartments))
        if missing_in_seed or extra_in_seed:
            raise ValueError(
                "Compartment mismatch between config 'compartments' and SeedFile columns. "
                f"Missing in SeedFile: {missing_in_seed}; Extra in SeedFile: {extra_in_seed}"
            )

    for _, row in seed_df.iterrows():
        patch = row["patch"]
        if patch not in populations:
            raise ValueError(f"SeedFile contains unknown patch '{patch}' not present in PatchFile")
        total = sum(row[c] for c in compartments)
        if not all(row[c] >= 0 for c in compartments):
            raise ValueError(f"Seed values must be non-negative for patch {patch}")
        if abs(total - populations[patch]) >= EPSILON:
            raise ValueError(f"Seed values do not sum to population for patch {patch}")

    # Set up network
    num_patches = len(patches)
    if "NetworkFile" not in config or config["NetworkFile"] is None:
        # Multi-patch model with no network: use zero matrix
        network_matrix = np.zeros((num_patches, num_patches))
    else:
        # Multi-patch model
        net_df = pd.read_csv(config["NetworkFile"])
        net_df = net_df[net_df["day"] == 0]
        patch_idx = {p: i for i, p in enumerate(patches)}
        network_matrix = np.zeros((num_patches, num_patches))

        for _, row in net_df.iterrows():
            source = row["source"].strip('"')
            target = row["target"].strip('"')
            if source not in patch_idx:
                raise ValueError(f"NetworkFile contains unknown source patch '{source}'")
            if target not in patch_idx:
                raise ValueError(f"NetworkFile contains unknown target patch '{target}'")
            i = patch_idx[source]
            j = patch_idx[target]
            if row["weight"] < 0:
                raise ValueError(f"Network weight must be non-negative between {row['source']} and {row['target']}")
            network_matrix[i, j] = row["weight"]

    # Set up model
    global_params = config.get("Parameters", {})

    # Collect per-patch parameters if provided (needed for transition-name validation too)
    patch_params: dict[str, dict[str, Any]] = {}
    if "PatchParameters" in config:
        for entry in config["PatchParameters"]:
            patch_name = entry["patch"]
            patch_params[patch_name] = entry.get("parameters", {})

    transitions_cfg = config.get("Transitions", {})
    # Transitions must be provided as arrow-map syntax in config, e.g.:
    # Transitions: {S -> I: beta, I -> R: gamma * I}
    if not isinstance(transitions_cfg, dict) or not transitions_cfg:
        raise ValueError("'Transitions' must be a non-empty mapping in arrow syntax, e.g. {S -> I: 'beta'}.")

    transitions: list[dict[str, Any]] = []
    patch_param_names = set()
    for per_patch in patch_params.values():
        patch_param_names |= set(per_patch.keys())
    allowed_names = set(compartments) | set(global_params.keys()) | patch_param_names
    for k, v in transitions_cfg.items():
        parts = [p.strip() for p in str(k).split("->")]
        if len(parts) != 2 or not all(parts):
            raise ValueError(f"Invalid transition key '{k}'. Use 'S -> I' format.")
        source, target = parts[0], parts[1]
        if source not in compartments or target not in compartments:
            raise ValueError(
                f"Transition '{k}' references unknown compartments. Known compartments: {sorted(compartments)}"
            )

        if isinstance(v, str):
            identifiers = set(re.findall(r"[A-Za-z_]\w*", v))
            python_keywords = {"and", "or", "not", "True", "False", "None"}
            unknown_identifiers = sorted(identifiers - allowed_names - python_keywords)
            if unknown_identifiers:
                raise ValueError(
                    f"Transition '{k}' uses unknown names in expression '{v}': {unknown_identifiers}. "
                    f"Allowed names are compartments + Parameters keys: {sorted(allowed_names)}"
                )

        transitions.append({"transition": f"{source}->{target}", "rate": v})

    # Validate that all configured patches exist in PatchFile
    unknown_patches = set(patch_params) - set(patches)
    if unknown_patches:
        raise ValueError(f"PatchParameters contains unknown patches: {sorted(unknown_patches)}")

    # Ensure every patch has a full parameter set: global + per-patch override
    for p in patches:
        patch_params[p] = {**global_params, **patch_params.get(p, {})}

    # Initialize the base model (will hold default/global transitions)
    base_model = CompartmentalModel(compartments=compartments, parameters=global_params, transitions=transitions)

    # Prepare initial conditions
    patch_idx = {p: i for i, p in enumerate(patches)}
    y0 = {}
    for _, row in seed_df.iterrows():
        for c in compartments:
            y0[f"{c}_{patch_idx[row['patch']]}"] = row[c]

    # Create network model
    net = NetworkModel(base_model=base_model, num_patches=num_patches, network_matrix=network_matrix)

    # Attach per-patch parameters to the network model
    net.patch_parameters = patch_params
    net.patch_names = patches

    return net, y0, patches, num_patches


def run_simulation(
    config: dict[str, Any], model_name: str, net: NetworkModel, y0: dict[str, float], patches: list, num_patches: int
) -> None:
    """Run the simulation and save results."""
    # Create output directories
    for subdir in ["plots", "runs"]:
        dir_path = os.path.join(config["OutputDir"], subdir)
        os.makedirs(dir_path, exist_ok=True)

    plots_dir = os.path.join(config["OutputDir"], "plots")
    runs_dir = os.path.join(config["OutputDir"], "runs")

    # Set up logger
    logger = setup_logger(model_name, config, num_patches, patches, net.base_model)

    # Validate and construct time range
    t_max = config.get("TMax")
    if not isinstance(t_max, int) or t_max <= 0:
        raise ValueError("'TMax' must be a positive integer.")
    t_range = np.arange(t_max, dtype=float)

    # Run simulation
    model = Model(net, compartments=list(net.base_model.compartments))
    out_ode = model.solve(y0, t_range)

    # Save results
    out_df = pd.DataFrame(out_ode)
    out_df["time"] = t_range
    cols = ["time"] + [c for c in out_df.columns if c != "time"]
    out_df = out_df[cols]

    csv_path = os.path.join(runs_dir, f"all_patches_{model_name}_ode.csv")
    out_df.to_csv(csv_path, index=False)
    logger.info(f"Saved simulation output to {csv_path}")

    model.visualize(t_range, out_ode, patches, plots_dir, model_name)

    logger.info(f"Saved all patch subplots to {plots_dir}/patch_timeseries_{model_name}_ode.png")

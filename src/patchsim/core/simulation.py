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


MODEL_TEMPLATE_CONFIGS: dict[str, dict[str, Any]] = {
    "sir": {
        "compartments": ["S", "I", "R"],
        "Parameters": {"beta": 0.08, "gamma": 0.1},
        "Transitions": {"S -> I": "beta", "I -> R": "gamma * I"},
    },
    "seir": {
        "compartments": ["S", "E", "I", "R"],
        "Parameters": {"beta": 0.08, "sigma": 0.2, "gamma": 0.1},
        "Transitions": {"S -> E": "beta", "E -> I": "sigma * E", "I -> R": "gamma * I"},
    },
    "sirs": {
        "compartments": ["S", "I", "R"],
        "Parameters": {"beta": 0.08, "gamma": 0.1, "waning": 0.02},
        "Transitions": {"S -> I": "beta", "I -> R": "gamma * I", "R -> S": "waning * R"},
    },
    "sis": {
        "compartments": ["S", "I"],
        "Parameters": {"beta": 0.08, "gamma": 0.1},
        "Transitions": {"S -> I": "beta", "I -> S": "gamma * I"},
    },
}


def get_config_schema() -> dict[str, Any]:
    """Return the JSON Schema for PatchSim configuration files."""
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://dsih-artpark.github.io/patchsim/config.schema.json",
        "title": "PatchSim configuration",
        "type": "object",
        "required": ["PatchFile", "SeedFile", "OutputDir", "Transitions", "TMax"],
        "properties": {
            "PatchFile": {"type": "string"},
            "SeedFile": {"type": "string"},
            "NetworkFile": {"type": ["string", "null"]},
            "OutputDir": {"type": "string"},
            "ModelName": {"type": "string"},
            "TMax": {"type": "integer", "minimum": 1},
            "Tolerance": {"type": ["number", "string"]},
            "MaxIter": {"type": "integer", "minimum": 1},
            "StartDate": {"type": ["string", "null"]},
            "EndDate": {"type": ["string", "null"]},
            "Logging": {"type": ["boolean", "string"]},
            "compartments": {
                "type": "array",
                "items": {"type": "string"},
                "minItems": 1,
            },
            "Compartments": {
                "type": "array",
                "items": {"type": "string"},
                "minItems": 1,
            },
            "Parameters": {
                "type": "object",
                "additionalProperties": {"type": ["number", "integer", "string", "boolean"]},
            },
            "PatchParameters": {
                "type": "array",
                "items": {
                    "type": "object",
                    "required": ["patch"],
                    "properties": {
                        "patch": {"type": "string"},
                        "parameters": {"type": "object"},
                    },
                    "additionalProperties": True,
                },
            },
            "Transitions": {
                "type": "object",
                "minProperties": 1,
                "additionalProperties": {"type": "string"},
            },
        },
        "additionalProperties": True,
    }


def get_model_catalog() -> list[dict[str, str]]:
    """Return built-in model references and YAML templates for the CLI."""
    return [{"name": name, "kind": "yaml-template"} for name in sorted(MODEL_TEMPLATE_CONFIGS)]


def get_available_template_names() -> list[str]:
    """Return the built-in starter template names."""
    return sorted(MODEL_TEMPLATE_CONFIGS)


def get_init_template_config(template_name: str, project_name: str) -> dict[str, Any]:
    """Return a starter project config for the requested template."""
    try:
        template = MODEL_TEMPLATE_CONFIGS[template_name]
    except KeyError as e:
        raise ValueError(
            f"Unknown template '{template_name}'. Available templates: {sorted(MODEL_TEMPLATE_CONFIGS)}"
        ) from e

    config: dict[str, Any] = {
        "PatchFile": "data/patch/patch-population.csv",
        "NetworkFile": "data/networks/network-static.csv",
        "SeedFile": "data/seeds/seed-initial.csv",
        "Logging": False,
        "ModelName": project_name,
        "TMax": 60,
        "Tolerance": 1e-8,
        "MaxIter": 10000,
        "StartDate": "2020-01-01",
        "EndDate": "2022-12-31",
        "OutputDir": f"output/{project_name}",
    }
    config.update(template)
    return config


def load_config(config_path: str) -> dict[str, Any]:
    """Load and validate configuration file."""
    cfg_path = Path(config_path).expanduser().resolve()
    with open(cfg_path, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    if not isinstance(config, dict):
        raise ValueError(
            f"Configuration error in {cfg_path}: expected YAML mapping (key-value pairs), "
            "but got a list or scalar value instead. \n"
            "Ensure the config file has the format: key1: value1\\n  key2: value2"
        )

    # Validate required fields
    required_fields = ["PatchFile", "SeedFile", "OutputDir", "Transitions", "TMax"]
    for field in required_fields:
        if field not in config:
            available = ", ".join(sorted(config.keys()))
            raise ValueError(
                f"Configuration error: missing required field '{field}'.\n"
                f"Available fields in config: {available}\n"
                f"Please add '{field}' to your config file."
            )

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
    # Load patch data (accept either case for the 'patch'/'population' columns)
    patch_df = pd.read_csv(config["PatchFile"])
    patch_col = next((c for c in patch_df.columns if c.lower() == "patch"), None)
    pop_col = next((c for c in patch_df.columns if c.lower() == "population"), None)
    if patch_col is None or pop_col is None:
        raise ValueError(
            f"PatchFile ({config['PatchFile']}) must have 'patch' and 'population' columns.\n"
            f"Found columns: {list(patch_df.columns)}"
        )
    patches = patch_df[patch_col].tolist()
    populations = patch_df.set_index(patch_col)[pop_col].to_dict()

    for p, pop in populations.items():
        if pop <= 0:
            raise ValueError(
                f"Invalid population in {config['PatchFile']}: patch '{p}' has population {pop}.\n"
                "Population must be a positive number (> 0).\n"
                "Please correct the population value in your patch file."
            )

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
                f"Compartment mismatch between config and SeedFile ({config['SeedFile']}).\n"
                f"Config compartments: {sorted(compartments)}\n"
                f"SeedFile columns: {seed_compartments}\n"
                f"Missing in SeedFile: {missing_in_seed}\n"
                f"Extra in SeedFile: {extra_in_seed}\n"
                "Please ensure the SeedFile has columns matching your config compartments."
            )

    for _, row in seed_df.iterrows():
        patch = row["patch"]
        if patch not in populations:
            raise ValueError(
                f"Unknown patch '{patch}' in SeedFile ({config['SeedFile']}).\n"
                f"Known patches from PatchFile: {sorted(populations.keys())}\n"
                "Please ensure all patches in SeedFile match PatchFile."
            )
        total = sum(row[c] for c in compartments)
        if not all(row[c] >= 0 for c in compartments):
            neg_comps = [c for c in compartments if row[c] < 0]
            raise ValueError(
                f"Invalid seed data for patch '{patch}': negative values found.\n"
                f"Compartments with negative values: {neg_comps}\n"
                "All seed values must be non-negative."
            )
        if abs(total - populations[patch]) >= EPSILON:
            raise ValueError(
                f"Seed mismatch for patch '{patch}': seed sum ({total}) != population ({populations[patch]}).\n"
                f"Seed compartments: {dict((c, row[c]) for c in compartments)}\n"
                "Ensure seed values sum exactly to the patch population."
            )

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
                raise ValueError(
                    f"Unknown source patch '{source}' in NetworkFile ({config['NetworkFile']}).\n"
                    f"Known patches: {sorted(patch_idx.keys())}\n"
                    "Please ensure all patches in NetworkFile match PatchFile."
                )
            if target not in patch_idx:
                raise ValueError(
                    f"Unknown target patch '{target}' in NetworkFile ({config['NetworkFile']}).\n"
                    f"Known patches: {sorted(patch_idx.keys())}\n"
                    "Please ensure all patches in NetworkFile match PatchFile."
                )
            i = patch_idx[source]
            j = patch_idx[target]
            if row["weight"] < 0:
                raise ValueError(
                    f"Invalid network weight in NetworkFile ({config['NetworkFile']}): "
                    f"weight={row['weight']} from '{source}' to '{target}'.\n"
                    "Network weights must be non-negative. Please correct your network file."
                )
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
            raise ValueError(
                f"Invalid transition key '{k}'.\n"
                "Use arrow format: 'source -> target' (e.g., 'S -> I')\n"
                "Available compartments: {}".format(sorted(compartments))
            )
        source, target = parts[0], parts[1]
        if source not in compartments or target not in compartments:
            bad_comps = [p for p in [source, target] if p not in compartments]
            raise ValueError(
                f"Transition '{k}' uses unknown compartments: {bad_comps}\n"
                f"Known compartments: {sorted(compartments)}\n"
                "Please correct the transition definition."
            )

        if isinstance(v, str):
            identifiers = set(re.findall(r"[A-Za-z_]\w*", v))
            python_keywords = {"and", "or", "not", "True", "False", "None"}
            unknown_identifiers = sorted(identifiers - allowed_names - python_keywords)
            if unknown_identifiers:
                raise ValueError(
                    f"Transition '{k}' uses undefined names: {unknown_identifiers}\n"
                    f"Expression: '{v}'\n"
                    f"Defined names (compartments + parameters): {sorted(allowed_names)}\n"
                    "Please check your transition expression for typos or add missing parameters."
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
) -> dict[str, Any]:
    """Run the simulation and save results.

    Returns:
        A summary dictionary describing the generated artifacts.
    """
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
        raise ValueError(
            f"Invalid 'TMax' value: {t_max}\n"
            "'TMax' must be a positive integer (number of time steps).\n"
            "Example: TMax: 100"
        )
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

    plot_path = os.path.join(plots_dir, f"patch_timeseries_{model_name}_ode.png")
    logger.info(f"Saved all patch subplots to {plot_path}")

    return {
        "model_name": model_name,
        "output_dir": config["OutputDir"],
        "csv_path": csv_path,
        "plot_path": plot_path,
        "num_patches": num_patches,
        "patches": patches,
        "t_max": t_max,
    }

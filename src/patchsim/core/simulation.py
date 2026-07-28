"""
This module is reserved for reusable simulation utilities or wrappers.
ODE and discrete simulation logic is implemented in core/model.py.
"""

import hashlib
import os
import re
from numbers import Real
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from patchsim.core.model import CompartmentalModel, NetworkModel
from patchsim.utils.logger import setup_logger
from patchsim.utils.viz import plot_patch_subplots

EPSILON = 1e-6
DEFAULT_SOLVER = "ode"
DEFAULT_TIME_STEP = 1.0
SOLVERS = ("ode", "discrete")


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
            "GroupFile": {"type": "string"},
            "InteractionFile": {"type": "string"},
            "InteractionUnits": {"type": "string", "minLength": 1},
            "OutputDir": {"type": "string"},
            "ModelName": {"type": "string"},
            "TMax": {"type": "integer", "minimum": 1},
            "Solver": {"type": "string", "enum": list(SOLVERS), "default": DEFAULT_SOLVER},
            "TimeStep": {
                "type": "number",
                "exclusiveMinimum": 0,
                "default": DEFAULT_TIME_STEP,
            },
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
        "dependentRequired": {
            "GroupFile": ["InteractionFile", "InteractionUnits"],
            "InteractionFile": ["GroupFile", "InteractionUnits"],
            "InteractionUnits": ["GroupFile", "InteractionFile"],
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
        "Solver": DEFAULT_SOLVER,
        "TimeStep": DEFAULT_TIME_STEP,
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

    solver, _t_max, time_step = get_run_settings(config)
    config["Solver"] = solver
    config["TimeStep"] = time_step

    # Resolve relative paths against the config file directory.
    cfg_dir = cfg_path.parent
    for key in ["PatchFile", "SeedFile", "NetworkFile", "GroupFile", "InteractionFile", "OutputDir"]:
        val = config.get(key)
        if isinstance(val, str) and val.strip():
            p = Path(val).expanduser()
            if not p.is_absolute():
                config[key] = str((cfg_dir / p).resolve())

    return config


def get_run_settings(config: dict[str, Any]) -> tuple[str, int, float]:
    """Validate and return the solver, reporting-point count, and grid interval."""
    solver = config.get("Solver", DEFAULT_SOLVER)
    if not isinstance(solver, str) or solver not in SOLVERS:
        raise ValueError(f"'Solver' must be one of {list(SOLVERS)}; received {solver!r}")

    t_max = config.get("TMax")
    if isinstance(t_max, bool) or not isinstance(t_max, int) or t_max <= 0:
        raise ValueError(f"'TMax' must be a positive integer count of reporting points; received {t_max!r}")

    time_step = config.get("TimeStep", DEFAULT_TIME_STEP)
    if isinstance(time_step, bool) or not isinstance(time_step, Real):
        raise ValueError(f"'TimeStep' must be a finite positive number; received {time_step!r}")
    time_step = float(time_step)
    if not np.isfinite(time_step) or time_step <= 0:
        raise ValueError(f"'TimeStep' must be a finite positive number; received {time_step!r}")
    return solver, t_max, time_step


def _numeric_values(frame: pd.DataFrame, columns: list[str], source: str) -> None:
    """Convert selected columns to finite numbers in place."""
    for column in columns:
        try:
            frame[column] = pd.to_numeric(frame[column], errors="raise")
        except (TypeError, ValueError) as e:
            raise ValueError(f"{source} column '{column}' must contain only numbers.") from e
        if not np.all(np.isfinite(frame[column].to_numpy(dtype=float))):
            raise ValueError(f"{source} column '{column}' must contain only finite numbers.")


def _load_groups(
    config: dict[str, Any], patches: list[str], populations: dict[str, float]
) -> tuple[list[str], dict[tuple[str, str], float]]:
    """Load the optional patch-by-group population table."""
    fields = ("GroupFile", "InteractionFile", "InteractionUnits")
    present = [field for field in fields if config.get(field) is not None]
    if present and len(present) != len(fields):
        missing = [field for field in fields if config.get(field) is None]
        raise ValueError(f"Grouped simulations require {list(fields)} together; missing {missing}.")
    if not present:
        return [], {}

    units = config["InteractionUnits"]
    if not isinstance(units, str) or not units.strip():
        raise ValueError("'InteractionUnits' must be a non-empty unit description.")

    path = config["GroupFile"]
    header = pd.read_csv(path, nrows=0).columns
    patch_col = next((c for c in header if c.lower() == "patch"), None)
    group_col = next((c for c in header if c.lower() == "group"), None)
    pop_col = next((c for c in header if c.lower() == "population"), None)
    if patch_col is None or group_col is None or pop_col is None:
        raise ValueError(f"GroupFile ({path}) must have 'patch', 'group', and 'population' columns.")

    frame = pd.read_csv(
        path,
        dtype={patch_col: "string", group_col: "string"},
        keep_default_na=False,
    ).rename(columns={patch_col: "patch", group_col: "group", pop_col: "population"})
    _numeric_values(frame, ["population"], f"GroupFile ({path})")
    for column in ("patch", "group"):
        values = frame[column].astype(str)
        if any(not value or value != value.strip() for value in values):
            raise ValueError(f"GroupFile ({path}) contains an empty or whitespace-padded {column} identifier.")
        frame[column] = values
    if (frame["population"] < 0).any():
        raise ValueError(f"GroupFile ({path}) populations must be non-negative.")
    if frame.duplicated(["patch", "group"]).any():
        raise ValueError(f"GroupFile ({path}) contains duplicate patch/group pairs.")

    patch_set = set(frame["patch"])
    expected_patch_set = set(patches)
    if patch_set != expected_patch_set:
        raise ValueError(
            f"GroupFile ({path}) patch identifiers do not match PatchFile "
            f"(missing={sorted(expected_patch_set - patch_set)}, extra={sorted(patch_set - expected_patch_set)})."
        )

    first_patch_groups = frame.loc[frame["patch"] == patches[0], "group"].tolist()
    if not first_patch_groups:
        raise ValueError(f"GroupFile ({path}) has no groups for first patch '{patches[0]}'.")
    groups = list(first_patch_groups)
    expected_groups = set(groups)
    rows_by_patch = {patch: rows for patch, rows in frame.groupby("patch", sort=False)}
    for patch in patches:
        patch_rows = rows_by_patch[patch]
        actual = set(patch_rows["group"])
        if actual != expected_groups:
            raise ValueError(
                f"GroupFile ({path}) groups for patch '{patch}' do not match the first patch "
                f"(missing={sorted(expected_groups - actual)}, extra={sorted(actual - expected_groups)})."
            )
        total = float(patch_rows["population"].sum())
        if abs(total - populations[patch]) >= EPSILON:
            raise ValueError(
                f"Group populations for patch '{patch}' sum to {total}, not PatchFile population {populations[patch]}."
            )

    group_populations = {
        (row.patch, row.group): float(row.population)
        for row in frame[["patch", "group", "population"]].itertuples(index=False)
    }
    return groups, group_populations


def _load_interactions(
    config: dict[str, Any],
    groups: list[str],
    patches: list[str],
    group_populations: dict[tuple[str, str], float],
    network_matrix: np.ndarray,
) -> tuple[np.ndarray | None, dict[str, Any] | None]:
    """Load a shared group interaction matrix and return validation diagnostics."""
    if not groups:
        return None, None

    path = config["InteractionFile"]
    required = ["focal_group", "contributor_group", "weight"]
    header = pd.read_csv(path, nrows=0).columns.tolist()
    missing = [column for column in required if column not in header]
    if missing:
        raise ValueError(f"InteractionFile ({path}) is missing columns: {missing}.")

    frame = pd.read_csv(
        path,
        dtype={"focal_group": "string", "contributor_group": "string"},
        keep_default_na=False,
    )
    _numeric_values(frame, ["weight"], f"InteractionFile ({path})")
    for column in ("focal_group", "contributor_group"):
        values = frame[column].astype(str)
        if any(not value or value != value.strip() for value in values):
            raise ValueError(f"InteractionFile ({path}) contains an empty or whitespace-padded {column}.")
        frame[column] = values
    if (frame["weight"] < 0).any():
        raise ValueError(f"InteractionFile ({path}) weights must be non-negative.")
    if frame.duplicated(["focal_group", "contributor_group"]).any():
        raise ValueError(f"InteractionFile ({path}) contains duplicate group pairs.")

    expected_groups = set(groups)
    for column in ("focal_group", "contributor_group"):
        actual = set(frame[column])
        unknown = actual - expected_groups
        if unknown:
            raise ValueError(f"InteractionFile ({path}) contains unknown groups in {column}: {sorted(unknown)}.")
        absent = expected_groups - actual
        if absent:
            raise ValueError(
                f"InteractionFile ({path}) must mention every group in {column}; missing {sorted(absent)}."
            )

    group_idx = {group: idx for idx, group in enumerate(groups)}
    matrix = np.zeros((len(groups), len(groups)), dtype=float)
    for row in frame.itertuples(index=False):
        matrix[group_idx[row.focal_group], group_idx[row.contributor_group]] = float(row.weight)

    max_reciprocity_residual = 0.0
    for patch in patches:
        for focal in groups:
            for contributor in groups:
                i = group_idx[focal]
                j = group_idx[contributor]
                forward = group_populations[(patch, focal)] * matrix[i, j]
                reverse = group_populations[(patch, contributor)] * matrix[j, i]
                denominator = max(forward, reverse)
                residual = abs(forward - reverse) / denominator if denominator > 0 else 0.0
                max_reciprocity_residual = max(max_reciprocity_residual, residual)

    interaction_row_sums = matrix.sum(axis=1)
    effective_spatial_matrix = np.ones((1, 1), dtype=float) if len(patches) == 1 else network_matrix
    spatial_row_sums = effective_spatial_matrix.sum(axis=1)
    diagnostics = {
        "units": config["InteractionUnits"].strip(),
        "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        "interaction_row_sum": {
            "min": float(interaction_row_sums.min()),
            "max": float(interaction_row_sums.max()),
        },
        "spatial_row_sum": {
            "min": float(spatial_row_sums.min()),
            "max": float(spatial_row_sums.max()),
        },
        "max_local_reciprocity_residual": float(max_reciprocity_residual),
        "reciprocity": "diagnostic_only",
    }
    return matrix, diagnostics


def setup_simulation(config: dict[str, Any]) -> tuple[NetworkModel, dict[str, float], list, int]:
    """Set up the simulation model and initial conditions."""
    # Load patch data (accept either case for the 'patch'/'population' columns)
    patch_columns = pd.read_csv(config["PatchFile"], nrows=0).columns
    patch_col = next((c for c in patch_columns if c.lower() == "patch"), None)
    pop_col = next((c for c in patch_columns if c.lower() == "population"), None)
    if patch_col is None or pop_col is None:
        raise ValueError(
            f"PatchFile ({config['PatchFile']}) must have 'patch' and 'population' columns.\n"
            f"Found columns: {list(patch_columns)}"
        )
    patch_df = pd.read_csv(config["PatchFile"], converters={patch_col: str})
    patches = patch_df[patch_col].tolist()
    if not patches:
        raise ValueError(f"PatchFile ({config['PatchFile']}) must contain at least one patch.")
    if any(not patch.strip() for patch in patches):
        raise ValueError(f"PatchFile ({config['PatchFile']}) contains an empty patch identifier.")
    if len(set(patches)) != len(patches):
        raise ValueError(f"PatchFile ({config['PatchFile']}) contains duplicate patch identifiers.")
    populations = patch_df.set_index(patch_col)[pop_col].to_dict()

    for p, pop in populations.items():
        if not isinstance(pop, Real) or not np.isfinite(pop) or pop <= 0:
            raise ValueError(
                f"Invalid population in {config['PatchFile']}: patch '{p}' has population {pop}.\n"
                "Population must be a positive number (> 0).\n"
                "Please correct the population value in your patch file."
            )

    groups, group_populations = _load_groups(config, patches, populations)

    # Load seed data
    seed_converters = {"patch": str}
    if groups:
        seed_converters["group"] = str
    seed_df = pd.read_csv(config["SeedFile"], converters=seed_converters, keep_default_na=False)
    if groups and "group" not in seed_df.columns:
        raise ValueError(f"SeedFile ({config['SeedFile']}) must include a 'group' column for grouped simulations.")
    if not groups and "group" in seed_df.columns:
        raise ValueError("SeedFile includes 'group', but GroupFile and InteractionFile are not configured.")
    identifier_columns = {"patch", "group"} if groups else {"patch"}
    seed_compartments = [col for col in seed_df.columns if col not in identifier_columns]
    _numeric_values(seed_df, seed_compartments, f"SeedFile ({config['SeedFile']})")

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

    if groups:
        if seed_df.duplicated(["patch", "group"]).any():
            raise ValueError(f"SeedFile ({config['SeedFile']}) contains duplicate patch/group pairs.")
        expected_pairs = set(group_populations)
        actual_pairs = set(zip(seed_df["patch"], seed_df["group"], strict=True))
        if actual_pairs != expected_pairs:
            raise ValueError(
                f"SeedFile ({config['SeedFile']}) patch/group coverage is incomplete "
                f"(missing={sorted(expected_pairs - actual_pairs)}, extra={sorted(actual_pairs - expected_pairs)})."
            )

    for _, row in seed_df.iterrows():
        patch = row["patch"]
        if patch not in populations:
            raise ValueError(
                f"Unknown patch '{patch}' in SeedFile ({config['SeedFile']}).\n"
                f"Known patches from PatchFile: {sorted(populations.keys())}\n"
                "Please ensure all patches in SeedFile match PatchFile."
            )
        group = row["group"] if groups else None
        if groups and group not in groups:
            raise ValueError(f"Unknown group '{group}' in SeedFile ({config['SeedFile']}).")
        total = sum(row[c] for c in compartments)
        if not all(row[c] >= 0 for c in compartments):
            neg_comps = [c for c in compartments if row[c] < 0]
            raise ValueError(
                f"Invalid seed data for patch '{patch}': negative values found.\n"
                f"Compartments with negative values: {neg_comps}\n"
                "All seed values must be non-negative."
            )
        expected_population = group_populations[(patch, group)] if groups else populations[patch]
        if abs(total - expected_population) >= EPSILON:
            stratum = f", group '{group}'" if groups else ""
            raise ValueError(
                f"Seed mismatch for patch '{patch}'{stratum}: "
                f"seed sum ({total}) != population ({expected_population}).\n"
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
        net_df = pd.read_csv(
            config["NetworkFile"],
            converters={"source": str, "target": str},
        )
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

    interaction_matrix, interaction_diagnostics = _load_interactions(
        config,
        groups,
        patches,
        group_populations,
        network_matrix,
    )

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

    # Create network model
    net = NetworkModel(
        base_model=base_model,
        num_patches=num_patches,
        network_matrix=network_matrix,
        groups=groups,
        interaction_matrix=interaction_matrix,
    )

    # Prepare initial conditions
    patch_idx = {p: i for i, p in enumerate(patches)}
    group_idx = {group: i for i, group in enumerate(groups)}
    y0 = {}
    for _, row in seed_df.iterrows():
        for c in compartments:
            group_position = group_idx[row["group"]] if groups else 0
            y0[net.state_key(c, patch_idx[row["patch"]], group_position)] = row[c]

    # Attach per-patch parameters to the network model
    net.patch_parameters = patch_params
    net.patch_names = patches
    net.interaction_diagnostics = interaction_diagnostics

    return net, y0, patches, num_patches


def run_simulation(
    config: dict[str, Any], model_name: str, net: NetworkModel, y0: dict[str, float], patches: list, num_patches: int
) -> dict[str, Any]:
    """Run the simulation and save results.

    Returns:
        A summary dictionary describing the generated artifacts.
    """
    solver, t_max, time_step = get_run_settings(config)
    config["Solver"] = solver
    config["TimeStep"] = time_step

    # Create output directories
    for subdir in ["plots", "runs"]:
        dir_path = os.path.join(config["OutputDir"], subdir)
        os.makedirs(dir_path, exist_ok=True)

    plots_dir = os.path.join(config["OutputDir"], "plots")
    runs_dir = os.path.join(config["OutputDir"], "runs")

    # Set up logger
    logger = setup_logger(model_name, config, num_patches, patches, net.base_model)

    t_range = np.arange(t_max, dtype=float) * time_step

    # Run simulation
    if solver == "ode":
        _times, results = net.simulate_ode(y0, t_range)
    else:
        results = net.simulate_discrete(y0, t_range)

    # Save results
    out_df = pd.DataFrame(results)
    out_df["time"] = t_range
    cols = ["time"] + [c for c in out_df.columns if c != "time"]
    out_df = out_df[cols]

    csv_path = os.path.join(runs_dir, f"all_patches_{model_name}_{solver}.csv")
    out_df.to_csv(csv_path, index=False)
    logger.info(f"Saved simulation output to {csv_path}")

    plot_patch_subplots(
        t_range,
        results,
        patches,
        plots_dir,
        model_name,
        compartments=list(net.base_model.compartments),
        groups=net.groups,
        solver=solver,
    )

    plot_path = os.path.join(plots_dir, f"patch_timeseries_{model_name}_{solver}.png")
    logger.info(f"Saved all patch subplots to {plot_path}")

    summary = {
        "model_name": model_name,
        "output_dir": config["OutputDir"],
        "csv_path": csv_path,
        "plot_path": plot_path,
        "num_patches": num_patches,
        "patches": patches,
        "t_max": t_max,
        "solver": solver,
        "time_step": time_step,
    }
    if net.groups:
        summary.update({"num_groups": net.num_groups, "groups": net.groups})
    return summary

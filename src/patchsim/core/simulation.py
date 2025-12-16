"""
This module is reserved for reusable simulation utilities or wrappers.
ODE and discrete simulation logic is implemented in core/model.py.
"""

import os
from typing import Any
import numpy as np
import pandas as pd
import yaml
from scipy.integrate import odeint
from patchsim.core.model import CompartmentalModel, NetworkModel
from patchsim.utils.logger import setup_logger
from patchsim.utils.viz import plot_patch_subplots

EPSILON = 1e-6

def load_config(config_path: str) -> dict[str, Any]:
    """Load and validate configuration file."""
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Validate required fields
    required_fields = ['PatchFile', 'SeedFile', 'OutputDir']
    for field in required_fields:
        if field not in config:
            raise ValueError(f"Missing required field '{field}' in config")

    return config


def setup_simulation(config: dict[str, Any]) -> tuple[NetworkModel, dict[str, float], list, int]:
    """Set up the simulation model and initial conditions."""
    # Load patch data
    patch_df = pd.read_csv(config['PatchFile'])
    patches = patch_df['patch'].tolist()
    populations = patch_df.set_index('patch')['Population'].to_dict()

    for p, pop in populations.items():
        if pop <= 0:
            raise ValueError(f"Population for patch {p} must be positive")

    # Load seed data
    seed_df = pd.read_csv(config['SeedFile'])
    compartments = [col for col in seed_df.columns if col != 'patch']

    for _, row in seed_df.iterrows():
        patch = row['patch']
        total = sum(row[c] for c in compartments)
        if not all(row[c] >= 0 for c in compartments):
            raise ValueError(f"Seed values must be non-negative for patch {patch}")
        if abs(total - populations[patch]) >= EPSILON:
            raise ValueError(f"Seed values do not sum to population for patch {patch}")

    # Set up network
    num_patches = len(patches)
    if 'NetworkFile' not in config or config['NetworkFile'] is None:
        # Single patch model
        network_matrix = [[0]]
        num_patches = 1
    else:
        # Multi-patch model
        net_df = pd.read_csv(config['NetworkFile'])
        net_df = net_df[net_df['day'] == 0]
        patch_idx = {p: i for i, p in enumerate(patches)}
        network_matrix = np.zeros((num_patches, num_patches))

        for _, row in net_df.iterrows():
            i = patch_idx[row['source'].strip('"')]
            j = patch_idx[row['target'].strip('"')]
            if row['weight'] < 0:
                raise ValueError(f"Network weight must be non-negative between {row['source']} and {row['target']}")
            network_matrix[i, j] = row['weight']

    # Set up model
    global_params = config.get('Parameters', {})
    transitions = config.get('Transitions', [])
    if not transitions:
        raise ValueError("No transitions defined in config")

    # Collect per-patch parameters if provided
    patch_params = {}
    if 'PatchParameters' in config:
        for entry in config['PatchParameters']:
            patch_name = entry['patch']
            patch_params[patch_name] = entry.get('parameters', {})

    # Validate patch parameter entries
    for p in patches:
        if p not in patch_params:
            # fallback: use global params if not defined for this patch
            patch_params[p] = global_params.copy()

    # Initialize the base model (will hold default/global transitions)
    base_model = CompartmentalModel(
        compartments=compartments,
        parameters=global_params,
        transitions=transitions
    )

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
    
    return net, y0, patches, num_patches



def run_simulation(
    config: dict[str, Any],
    model_name: str,
    net: NetworkModel,
    y0: dict[str, float],
    patches: list,
    num_patches: int
) -> None:
    """Run the simulation and save results."""
    # Create output directories
    for subdir in ['plots', 'runs']:
        dir_path = os.path.join(config['OutputDir'], subdir)
        os.makedirs(dir_path, exist_ok=True)

    plots_dir = os.path.join(config['OutputDir'], 'plots')
    runs_dir = os.path.join(config['OutputDir'], 'runs')

    # Set up logger
    logger = setup_logger(model_name, config, num_patches, patches, net.base_model)

    # Run simulation
    t_range = np.linspace(0, config['TMax']-1, int(config['TMax']))
    _, out_ode = net.simulate_ode(y0, t_range, odeint)

    # Save results
    out_df = pd.DataFrame(out_ode)
    out_df['time'] = t_range
    cols = ['time'] + [c for c in out_df.columns if c != 'time']
    out_df = out_df[cols]

    csv_path = os.path.join(runs_dir, f"all_patches_{model_name}_ode.csv")
    out_df.to_csv(csv_path, index=False)
    logger.info(f"Saved simulation output to {csv_path}")

    plot_patch_subplots(t_range, out_ode, patches, plots_dir, model_name)
    logger.info(f"Saved all patch subplots to {plots_dir}/patch_timeseries_{model_name}_ode.png")
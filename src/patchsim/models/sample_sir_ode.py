"""
Sample SIR model implementation.
This serves as a reference implementation for creating new models.
"""
import os

import numpy as np
import pandas as pd
import yaml
from scipy.integrate import odeint

from patchsim.core.model import CompartmentalModel, NetworkModel
from patchsim.utils.logger import setup_logger
from patchsim.utils.viz import plot_patch_subplots


def create_sir_model(beta: float = 0.3, gamma: float = 0.1) -> CompartmentalModel:
    """Create a basic SIR model with given transmission and recovery rates."""
    return CompartmentalModel(
        compartments=["S", "I", "R"],
        parameters={"beta": beta, "gamma": gamma},
        transitions=[
            {"from": "S", "to": "I", "rate": "beta * S * lambda_i"},
            {"from": "I", "to": "R", "rate": "gamma * I"}
        ]
    )


def validate_seed_data(seed_data: dict, population: float) -> None:
    """Validate seed data for a single patch."""
    total = sum(seed_data[c] for c in ["S", "I", "R"])
    if not all(seed_data[c] >= 0 for c in ["S", "I", "R"]):
        raise ValueError("Seed values must be non-negative")
    if abs(total - population) >= 1e-6:
        raise ValueError(f"Seed values do not sum to population {population}")


def run_simulation(config_path: str) -> None:
    """Run the SIR model simulation with the given configuration."""
    # Load configuration
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Load patch data
    patch_df = pd.read_csv(config['PatchFile'])
    patches = patch_df['patch'].tolist()
    populations = patch_df.set_index('patch')['Population'].to_dict()

    for p, pop in populations.items():
        if pop <= 0:
            raise ValueError(f"Population for patch {p} must be positive")

    # Load and validate seed data
    seed_df = pd.read_csv(config['SeedFile'])
    for _, row in seed_df.iterrows():
        patch = row['patch']
        validate_seed_data(row, populations[patch])

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

    # Create model
    base_model = create_sir_model(
        beta=config.get('Beta', 0.3),
        gamma=config.get('Gamma', 0.1)
    )

    # Prepare initial conditions
    patch_idx = {p: i for i, p in enumerate(patches)}
    y0 = {}
    for _, row in seed_df.iterrows():
        for c in ["S", "I", "R"]:
            y0[f"{c}_{patch_idx[row['patch']]}"] = row[c]

    # Create network model
    net = NetworkModel(base_model=base_model, num_patches=num_patches, network_matrix=network_matrix)

    # Run simulation
    t_range = np.linspace(0, config['TMax']-1, int(config['TMax']))
    model_name = "sample-sir-ode"
    logger = setup_logger(model_name, config, num_patches, patches, base_model)

    # Create output directories
    for subdir in ['plots', 'runs']:
        dir_path = os.path.join(config['OutputDir'], subdir)
        os.makedirs(dir_path, exist_ok=True)
    plots_dir = os.path.join(config['OutputDir'], 'plots')
    runs_dir = os.path.join(config['OutputDir'], 'runs')

    # Run simulation
    _, out_ode = net.simulate_ode(y0, t_range, odeint)

    # Save results
    out_df = pd.DataFrame(out_ode)
    out_df['time'] = t_range
    cols = ['time'] + [c for c in out_df.columns if c != 'time']
    out_df = out_df[cols]
    csv_path = os.path.join(runs_dir, f"all_patches_{model_name}_ode.csv")
    out_df.to_csv(csv_path, index=False)
    logger.info(f"Saved simulation output to {csv_path}")

    # Plot results
    plot_patch_subplots(t_range, out_ode, patches, plots_dir, model_name)
    logger.info(f"Saved all patch subplots to {plots_dir}/patch_timeseries_{model_name}_ode.png")


if __name__ == "__main__":
    run_simulation("configs/sample-sir-ode.yaml")

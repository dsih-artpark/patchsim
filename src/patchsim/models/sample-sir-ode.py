import yaml
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from patchsim.core.model import CompartmentalModel, NetworkModel
from patchsim.utils.viz import plot_patch_subplots
from patchsim.utils.logger import setup_logger
import os

# Load config
with open('configs/sample-sir-ode.yaml') as f:
    config = yaml.safe_load(f)

# Read patch population
patch_df = pd.read_csv(config['PatchFile'])
patches = patch_df['patch'].tolist()
populations = patch_df.set_index('patch')['Population'].to_dict()

# Check populations are positive
for p, pop in populations.items():
    assert pop > 0, f"Population for patch {p} must be positive"

# Read seeds
seed_df = pd.read_csv(config['SeedFile'])

# Check seeds are non-negative and do not exceed population
for _, row in seed_df.iterrows():
    patch = row['patch']
    total = row['S'] + row['I'] + row['R']
    assert all(row[c] >= 0 for c in ['S', 'I', 'R']), f"Seed values must be non-negative for patch {patch}"
    assert abs(total - populations[patch]) < 1e-6, f"Seed values do not sum to population for patch {patch}"

# Read network matrix (static, day 0)
net_df = pd.read_csv(config['NetworkFile'])
net_df = net_df[net_df['day'] == 0]
patch_idx = {p: i for i, p in enumerate(patches)}
num_patches = len(patches)
network_matrix = np.zeros((num_patches, num_patches))
for _, row in net_df.iterrows():
    i = patch_idx[row['source'].strip('"')]
    j = patch_idx[row['target'].strip('"')]
    assert row['weight'] >= 0, f"Network weight must be non-negative between {row['source']} and {row['target']}"
    network_matrix[i, j] = row['weight']

# SIR parameters
beta = config['Beta']
gamma = config['Gamma']

# Initial conditions
y0 = {}
for _, row in seed_df.iterrows():
    for c in ['S', 'I', 'R']:
        y0[f"{c}_{patch_idx[row['patch']]}" ] = row[c]

# Model setup
base = CompartmentalModel(
    compartments=["S", "I", "R"],
    parameters={"beta": beta, "gamma": gamma},
    transitions=[
        {"from": "S", "to": "I", "rate": "beta * S * lambda_i"},
        {"from": "I", "to": "R", "rate": "gamma * I"},
    ]
)

net = NetworkModel(base_model=base, num_patches=num_patches, network_matrix=network_matrix)

t_range = np.linspace(0, config['TMax']-1, int(config['TMax']))

model_name = "sample-sir-ode"
logger = setup_logger(model_name, config, num_patches, patches, base)

# Output subdirectories
for subdir in ['plots', 'runs']:
    dir_path = os.path.join(config['OutputDir'], subdir)
    os.makedirs(dir_path, exist_ok=True)
plots_dir = os.path.join(config['OutputDir'], 'plots')
runs_dir = os.path.join(config['OutputDir'], 'runs')

_, out_ode = net.simulate_ode(y0, t_range, odeint)

out_df = pd.DataFrame(out_ode)
out_df['time'] = t_range
cols = ['time'] + [c for c in out_df.columns if c != 'time']
out_df = out_df[cols]
csv_path = os.path.join(runs_dir, f"all_patches_{model_name}_ode.csv")
out_df.to_csv(csv_path, index=False)
logger.info(f"Saved simulation output to {csv_path}")

plot_patch_subplots(t_range, out_ode, patches, plots_dir, model_name)
logger.info(f"Saved all patch subplots to {plots_dir}/patch_timeseries_{model_name}_ode.png")
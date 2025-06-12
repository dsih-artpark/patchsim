# PatchSim

> ⚠️ **WARNING: Active Development**  
> PatchSim is currently under **heavy development**. Features are experimental and subject to change. Not recommended for production use at this stage.

![PatchSim Banner](assets/patchsim-banner.png)

**PatchSim** is a modular metapopulation simulation framework for multi-disease epidemiological modelling.

---

## Vision

To develop a general-purpose, modular simulation framework for patch-based metapopulation epidemiology, enabling modellers and researchers to simulate disease transmission under diverse scenarios, diseases, and intervention strategies. The framework balances robust scientific modelling with flexibility for exploratory research and translational use cases.

---

## Core Features

PatchSim aims to support a range of modelling features commonly used in metapopulation disease simulations:

- 🗺️ **Spatial Networks**: Represent geographical units (e.g., subdistricts, regions) as interconnected patches with movement/contact matrices.
- 👥 **Stratification by Population Attributes**:
  - **Age groups**
  - **Species (e.g., cattle, buffalo)**
  - **Risk groups or occupations**
- 🧪 **Disease Agnostic Compartment Models**:
  - SIR, SEIR, SIRS and extensions
  - Supports both discrete timestep and ODE-based solvers
- 🛠️ **Scenario and Parameter Management**:
  - Batch simulations for scenario comparison
  - Sensitivity analysis and parameter sweeps
- 🧵 **Reproducibility**:
  - Random seed control and metadata logging
  - Version-tracked configurations
- 📦 **Modularity**:
  - Plug-in architecture for solvers, interventions, and input data pipelines

---

## Installation

To install PatchSim using [uv](https://github.com/astral-sh/uv):

```bash
# Clone the repository
git clone https://github.com/dsih-artpark/patchsim
cd patchsim

# Create a virtual environment and install dependencies
uv venv
source .venv/bin/activate
uv pip install -e .

# For development (with dev dependencies)
uv pip install -e .[dev]
```

---

## Usage

### Command Line Interface

PatchSim provides a command-line interface for running simulations:

```bash
# Show help and available options
patchsim --help

# Run simulation with sample SIR model
patchsim --config configs/sample-sir-ode.yaml

# Run simulation with custom model
patchsim --config path/to/your/config.yaml
```

The simulation outputs are saved in the following structure:
```
output/
├── logs/
│   └── cli_YYYYMMDD_HHMMSS.log  # Timestamped log files
└── sample-sir-ode/              # Model-specific output directory
    ├── plots/
    │   └── patch_timeseries_sample-sir-ode_ode.png
    └── runs/
        └── all_patches_sample-sir-ode_ode.csv
```

### Configuration

Simulations are configured using YAML files. The configuration file specifies:
- Model parameters (e.g., transmission rates)
- Input files (patch populations, seed data, network)
- Simulation settings (time horizon, output directory)

Example configuration:
```yaml
# Model parameters
Beta: 0.3
Gamma: 0.1

# Input files
PatchFile: data/patch/sample-sir-ode-patch-population.csv
SeedFile: data/seeds/sample-sir-ode-patchA-2.csv
NetworkFile: data/networks/sample-network-static.csv

# Simulation settings
OutputDir: output/sample-sir-ode
TMax: 50
```

---

## Contributing

We welcome contributions!

To contribute: fork the repo, create a branch, make your changes, and open a pull request.  
For major changes, please open an issue first to discuss.

Thanks for helping improve the framework!


## License

This project is licensed under the **GNU General Public License v3.0**.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

You may use, modify, and share this project under the same license terms. See the [LICENSE](./LICENSE) file for full details.

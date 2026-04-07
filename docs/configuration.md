# Configuration

PatchSim simulations are configured using YAML files.

## Model definition

```yaml
# Input files (required)
PatchFile: data/patch/patch-population.csv
NetworkFile: data/networks/network.csv
SeedFile: data/seeds/seed-initial.csv

# Model configuration
ModelName: sample-sir-ode

# Simulation parameters
TMax: 50
Tolerance: 1e-8
MaxIter: 10000
StartDate: 2020-01-01
EndDate: 2022-12-31
OutputDir: output/sample-sir-ode
compartments: ["S", "I", "R"]
Parameters:
  beta: 0.08
  gamma: 0.1
Transitions:
  "S -> I": "beta"
  "I -> R": "gamma * I"
```
# PatchSim

PatchSim is a modular, general-purpose metapopulation epidemiological simulation framework designed for research and translation use cases. It is a flexible and modular system for simulating epidemiological compartment models across single or multiple interacting patches.

The framework generalizes classical compartmental models (e.g. SIR) by allowing multiple subpopulations ("patches") to interact through a weighted directed network.

PatchSim emphasizes:
- Explicit separation between patch dynamics and network structure
- Continuous-time ODE-based simulation
- Configurable parameters via YAML
- Extensibility toward richer network and mobility models
 
## Features

- **Flexible Compartment Models**: Define custom SIR, SEIR, or other compartmental models that support multiple diseases and problem types
- **Multi-Patch Networks**: Simulate disease spread across connected geographical regions
- **CLI Integration**: Integrated CLI for running models
- **Multiple Simulation Methods**: Support for both ODE and discrete-time simulation modes
- **YAML Configuration**: Structured, readable configuration files
- **Per-Patch Parameters**: Built-in reproducibility and management of different parameters for each patch
- **Network Connectivity**: Define custom mixing matrices between patches
- **Extensibility**: Extensible and user-friendly

## Documentation
- [Core Module](reference/core.md)
- [Models](reference/models.md)
- [Utilities](reference/utils.md)


# Architecture

PatchSim follows a layered architecture:

## Core components

### Patch
Represents a single subpopulation with:
- Fixed population size
- Compartment state vector
- Parameter dictionary
- Validation logic

### Network
Defines interactions between patches via a weighted directed graph.

### Model
Orchestrates:
- ODE construction
- Numerical integration
- Visualization of results

This separation ensures that changes to network structure do not require
rewriting patch-level dynamics.

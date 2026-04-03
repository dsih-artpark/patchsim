# Model Implementations

This section documents concrete epidemiological models implemented using the PatchSim framework.

Unlike the core modules, models are **not reusable APIs** but
self-contained reference implementations that demonstrate how to
instantiate and run simulations using the framework.

---

## Design Philosophy

PatchSim separates:

- **Model definition** (compartments, transitions, parameters)
- **Simulation engine** (ODE solvers, network coupling)
- **Execution scripts** (models in this directory)

Files in `patchsim.models` serve as:
- Validation cases
- Reference implementations
- Templates for new models

---

## ka_fmd_sirsv_discrete

**File:**  
`patchsim/models/ka_fmd_sirsv_discrete.py`

**Description:**  
A discrete-time, multi-compartment SIRSV-style model designed for
foot-and-mouth disease (FMD) dynamics.

**Key features:**
- Discrete-time update rules
- Multiple infectious states
- Demonstrates non-ODE simulation capability
- Designed for Karnataka (KA) regional modeling

**Usage pattern:**
- Defines compartments and transitions explicitly
- Uses PatchSim core utilities for execution
- Can be adapted into YAML-driven workflows

---

## Extending Models

To create a new model:
1. Copy an existing model script
2. Modify compartments and transitions
3. Adjust parameters or network structure
4. Run via the `patchsim` CLI using `uv run`, for example `uv run patchsim run -c config.yaml`

For most use cases, users are encouraged to prefer
**YAML-based model definitions** over Python scripts.

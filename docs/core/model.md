# Model Module

The `model.py` module provides the foundation for defining epidemiological compartment models.

## CompartmentalModel

Base class for defining compartment-based disease models (e.g., SIR, SEIR, SIRS).

### Class Signature

```python
class CompartmentalModel:
    def __init__(
        self, 
        compartments: list[str], 
        parameters: dict[str, float], 
        transitions: list[dict[str, Any]]
    )
```

### Attributes

| Attribute      | Type                    | Description                                           |
|----------------|-------------------------|-------------------------------------------------------|
| `compartments` | `list[str]`             | Names of disease compartments (e.g., `['S', 'I', 'R']`) |
| `parameters`   | `dict[str, float]`      | Model parameters (e.g., `{'beta': 0.5, 'gamma': 0.1}`)  |
| `transitions`  | `list[dict[str, Any]]`  | Transition rules between compartments                 |

### Transition Format

Each transition is a dictionary with:
- `'from'`: Source compartment name
- `'to'`: Target compartment name  
- `'rate'`: Rate expression (string or float)

**Example:**
```python
transitions = [
    {'from': 'S', 'to': 'I', 'rate': 'beta * I / N'},
    {'from': 'I', 'to': 'R', 'rate': 'gamma'}
]
```

### Methods

#### `compute_rates(state: dict[str, float]) -> dict[str, float]`

Computes transition rates for the current state by evaluating rate expressions.

**Parameters:**

- `state` (dict): Current population in each compartment (e.g., `{'S': 990, 'I': 10, 'R': 0}`)

**Returns:**

- `dict`: Transition rates with keys formatted as `"{source}_to_{target}"`

**Example:**

```python
model = CompartmentalModel(
    compartments=['S', 'I', 'R'],
    parameters={'beta': 0.5, 'gamma': 0.1},
    transitions=[
        {'from': 'S', 'to': 'I', 'rate': 'beta * I / N'},
        {'from': 'I', 'to': 'R', 'rate': 'gamma'}
    ]
)

state = {'S': 990, 'I': 10, 'R': 0, 'N': 1000}
rates = model.compute_rates(state)
# Returns: {'S_to_I': 4.95, 'I_to_R': 1.0}
```

**Implementation Details:**

- String rate expressions are evaluated using parameter substitution
- Parameters and state variables are replaced with their numeric values
- The rate is multiplied by the source compartment population

---

## NetworkModel

Extends a base compartmental model to support multi-patch network simulations.

### Class Signature

```python
class NetworkModel:
    def __init__(
        self,
        base_model: CompartmentalModel,
        num_patches: int,
        network_matrix: list[list[float]]
    )
```

### Attributes

| Attribute          | Type                   | Description                                      |
|--------------------|------------------------|--------------------------------------------------|
| `base_model`       | `CompartmentalModel`   | Underlying compartment model                     |
| `num_patches`      | `int`                  | Number of patches in the network                 |
| `network`          | `list[list[float]]`    | NxN mixing matrix between patches                |
| `all_compartments` | `list[str]`            | All compartments across all patches (e.g., `['S_0', 'I_0', 'R_0', 'S_1', ...]`) |
| `patch_parameters` | `dict` (optional)      | Per-patch parameter overrides                    |

### Methods

#### `get_patch_state(full_state: Dict[str, float], patch_idx: int) -> Dict[str, float]`

Extracts the state variables for a specific patch from the full state.

**Parameters:**

- `full_state` (dict): Complete state with all compartments (e.g., `{'S_0': 990, 'I_0': 10, ...}`)
- `patch_idx` (int): Zero-based patch index

**Returns:**

- `dict`: State dictionary for the specified patch (e.g., `{'S': 990, 'I': 10, 'R': 0}`)

**Example:**

```python
full_state = {'S_0': 990, 'I_0': 10, 'R_0': 0, 'S_1': 480, 'I_1': 20, 'R_1': 0}
patch_0_state = network_model.get_patch_state(full_state, 0)
# Returns: {'S': 990, 'I': 10, 'R': 0}
```

---

#### `get_patch_population(state: Dict[str, float]) -> float`

Calculates total population for a single patch state.

**Parameters:**

- `state` (dict): Patch state dictionary (e.g., `{'S': 990, 'I': 10, 'R': 0}`)

**Returns:**

- `float`: Sum of all compartments

---

#### `compute_force_of_infection(full_state: dict[str, float], infected_compartment: str = "I") -> list[float]`

Calculates the infection pressure on each patch based on network connectivity.

**Parameters:**

- `full_state` (dict): Current state of all compartments across all patches
- `infected_compartment` (str): Name of infectious compartment (default: `"I"`)

**Returns:**

- `list[float]`: Force of infection for each patch

**Behavior:**

- **Single patch**: Returns local force of infection `beta * I / N`
- **Multi-patch**: Computes weighted sum using network matrix: `Σ(network[i][j] * I_j / N_j)`

**Example:**

```python
foi = network_model.compute_force_of_infection(full_state, infected_compartment='I')
# Returns: [0.005, 0.012]  # Force of infection for patch 0 and patch 1
```

---

#### `compute_derivatives(state: dict[str, float]) -> dict[str, float]`

Computes derivatives (rates of change) for all compartments based on transitions.

**Parameters:**

- `state` (dict): Current state of all compartments

**Returns:**

- `dict`: Derivatives for each compartment

**Implementation Details:**

- Processes each patch independently
- Supports per-patch parameter overrides via `patch_parameters` attribute
- Applies transition rules from the base model
- Returns rates for all compartments across all patches

---

#### `simulate_discrete(y0_dict: dict[str, float], t_range: list[float]) -> dict[str, list[float]]`

Runs a discrete-time simulation using Euler method.

**Parameters:**

- `y0_dict` (dict): Initial state for all compartments
- `t_range` (list): Time points for simulation

**Returns:**

- `dict`: History of each compartment over time

**Example:**

```python
initial_state = {
    'S_0': 990, 'I_0': 10, 'R_0': 0,
    'S_1': 480, 'I_1': 20, 'R_1': 0
}
t_range = list(range(0, 100))

history = network_model.simulate_discrete(initial_state, t_range)
# Returns: {'S_0': [...], 'I_0': [...], 'R_0': [...], ...}
```

---

#### `simulate_ode(y0_dict: dict[str, float], t_range: list[float], integrator: Callable = odeint) -> tuple[list[float], dict[str, list[float]]]`

Runs an ODE simulation using SciPy's `odeint` (or custom integrator).

**Parameters:**

- `y0_dict` (dict): Initial state for all compartments
- `t_range` (list): Time points for simulation
- `integrator` (Callable): Integration function (default: `scipy.integrate.odeint`)

**Returns:**

- `tuple`: 
  - Time points (list)
  - State history (dict) with arrays for each compartment

**Example:**

```python
from scipy.integrate import odeint

initial_state = {
    'S_0': 990, 'I_0': 10, 'R_0': 0,
    'S_1': 480, 'I_1': 20, 'R_1': 0
}
t_range = list(range(0, 100))

times, history = network_model.simulate_ode(initial_state, t_range)
# times: [0, 1, 2, ..., 99]
# history: {'S_0': array([...]), 'I_0': array([...]), ...}
```

---

## Complete Usage Example

```python
from patchsim.core.model import CompartmentalModel, NetworkModel

# Define SEIR model
seir = CompartmentalModel(
    compartments=['S', 'E', 'I', 'R'],
    parameters={'beta': 0.5, 'sigma': 0.2, 'gamma': 0.1},
    transitions=[
        {'from': 'S', 'to': 'E', 'rate': 'beta * I / N'},
        {'from': 'E', 'to': 'I', 'rate': 'sigma'},
        {'from': 'I', 'to': 'R', 'rate': 'gamma'}
    ]
)

# Create 2-patch network with mixing
network_matrix = [
    [0.95, 0.05],  # Patch 0: 95% local, 5% to patch 1
    [0.02, 0.98]   # Patch 1: 2% to patch 0, 98% local
]

network_model = NetworkModel(
    base_model=seir,
    num_patches=2,
    network_matrix=network_matrix
)

# Set per-patch parameters (optional)
network_model.patch_parameters = {
    'PatchA': {'beta': 0.5, 'sigma': 0.2, 'gamma': 0.1},
    'PatchB': {'beta': 0.4, 'sigma': 0.25, 'gamma': 0.12}
}

# Initial conditions
initial_state = {
    'S_0': 9900, 'E_0': 0, 'I_0': 100, 'R_0': 0,
    'S_1': 4980, 'E_1': 0, 'I_1': 20, 'R_1': 0
}

# Run ODE simulation
t_range = list(range(0, 365))
times, results = network_model.simulate_ode(initial_state, t_range)

# Access results
import matplotlib.pyplot as plt

plt.plot(times, results['I_0'], label='Patch 0 - Infected')
plt.plot(times, results['I_1'], label='Patch 1 - Infected')
plt.xlabel('Time (days)')
plt.ylabel('Number of Individuals')
plt.legend()
plt.show()
```

---

## Key Design Notes

1. **Compartment Naming**: In `NetworkModel`, compartments are suffixed with patch index (e.g., `S_0`, `S_1`)
2. **Parameter Overrides**: Use `patch_parameters` attribute to assign different parameters to each patch
3. **Rate Expressions**: Support string-based symbolic expressions evaluated at runtime
4. **Network Matrix**: Must be NxN square matrix where `network[i][j]` represents mixing from patch j to patch i

---

## See Also

- [Network Module](network.md) - Managing patch connectivity
- [Patch Module](patch.md) - Patch configuration
- [Simulation Module](simulation.md) - Running simulations
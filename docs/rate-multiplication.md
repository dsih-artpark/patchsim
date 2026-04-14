# Rate Multiplication Rule

## Overview

PatchSim uses a **rate multiplication** approach to scale transition rates based on compartment sizes or population dynamics. This section documents how rates are computed and applied across compartmental models.

## Basic Principle

In compartmental epidemiological models, transition rates represent the per-capita rate at which individuals move from one compartment to another. The actual number of individuals transitioning in a time step is computed by multiplying this per-capita rate by the number of individuals in the source compartment.

### Mathematical Formulation

For a transition from compartment A to compartment B:

$$\text{Flow}_{A \to B} = \text{rate}_{A \to B} \times n_A$$

where:
- $\text{Flow}_{A \to B}$ is the number of individuals transitioning per unit time
- $\text{rate}_{A \to B}$ is the per-capita transition rate (e.g., recovery rate, transmission rate)
- $n_A$ is the current population in compartment A

## Application in PatchSim

### Contact-Based Transmission

For disease transmission, the force of infection is computed as a weighted sum of infectious populations across patches:

$$\lambda_i(t) = \beta \sum_j W_{ij} \frac{I_j}{N_j}$$

where:
- $\lambda_i$ is the per-capita force of infection in patch i
- $\beta$ is the transmission rate
- $W_{ij}$ is the network weight from source patch j to target patch i
- $I_j$ is the number of infectious individuals in source patch j
- $N_j$ is the total population in source patch j

The actual number of new infections is then:

$$dS_i = -\lambda_i \times S_i \, dt$$

### Recovery and Other Transitions

For non-transmission transitions (recovery, loss of immunity, disease progression), the rate multiplication is simpler:

$$\text{Flow} = \text{rate} \times \text{compartment\_size}$$

For example:
$$dI = -\gamma \times I \, dt$$

where $\gamma$ is the per-capita recovery rate.

## Implementation Details

### Parameter-Based Rate Scaling

Parameters defined in the `Parameters` section of the configuration are used directly as per-capita rates:

```yaml
Parameters:
  beta: 0.5      # Transmission rate (contacts per individual per day)
  gamma: 0.1     # Recovery rate (1/disease duration)
  sigma: 0.2     # Rate of progression E → I
```

These rates multiply the compartment sizes during ODE integration.

### Validation in PatchSim

PatchSim validates that all transition expressions reference valid compartments and parameters. This ensures that rate multiplication is applied consistently across the model.

## Common Pitfalls

1. **Forgetting compartment size**: If you want a fixed number of individuals to transition (not rate-based), you must normalize by population size in the configuration or code.

2. **Mixing rates and absolute numbers**: All transition rates in PatchSim are assumed to be per-capita. Do not mix per-capita rates with fixed transition counts in the same model.

3. **Network weight interpretation**: Network weights are applied to the compartment proportion, not the absolute number. This ensures that transmission scales appropriately with patch population sizes.

## References

- [Mathematical Model](mathematical-model.md) for detailed network-coupled ODE formulations
- [Configuration Guide](configuration.md) for specifying parameters and transitions

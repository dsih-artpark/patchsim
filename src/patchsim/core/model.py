"""
Core model implementation for compartmental models.
"""

from typing import Any, Callable

import numpy as np
from scipy.integrate import odeint


class CompartmentalModel:
    """Base class for compartmental models."""

    def __init__(self, compartments: list[str], parameters: dict[str, float], transitions: list[dict[str, Any]]):
        """Initialize the model with compartments, parameters, and transitions."""
        self.compartments = compartments
        self.parameters = parameters
        self.transitions = transitions

    def compute_rates(self, state: dict[str, float]) -> dict[str, float]:
        """Compute transition rates for each compartment."""
        rates = {}
        for transition in self.transitions:
            source = transition['source']
            target = transition['target']
            rate = transition['rate']
            # Handle rate expressions
            if isinstance(rate, str):
                for param, value in self.parameters.items():
                    rate = rate.replace(param, str(value))
                for comp, value in state.items():
                    rate = rate.replace(comp, str(value))
                rate = eval(rate)
            rates[f"{source}_to_{target}"] = rate * state[source]
        return rates


class NetworkModel:
    """Network model for multi-patch simulations."""

    def __init__(self, base_model: CompartmentalModel, num_patches: int, network_matrix: list[list[float]]):
        """Initialize the network model."""
        self.base_model = base_model
        self.num_patches = num_patches
        self.network = network_matrix
        self.all_compartments = [f"{c}_{i}" for i in range(num_patches) for c in base_model.compartments]

    def compute_force_of_infection(self, full_state: dict[str, float]) -> list[float]:
        """Compute force of infection for each patch."""
        lambdas = []
        for i in range(self.num_patches):
            force = 0
            for j in range(self.num_patches):
                infected_j = full_state[f"I_{j}"]
                pop_j = sum(full_state[f"{c}_{j}"] for c in self.base_model.compartments)
                force += self.network[i][j] * infected_j / pop_j
            lambdas.append(force)
        return lambdas

    def simulate_discrete(self, y0_dict: dict[str, float], t_range: list[float]) -> dict[str, list[float]]:
        """Run discrete-time simulation."""
        state = y0_dict.copy()
        history = {c: [state[c]] for c in self.all_compartments}
        for _ in t_range[1:]:
            new_state = state.copy()
            lambdas = self.compute_force_of_infection(state)
            for i in range(self.num_patches):
                for c in self.base_model.compartments:
                    comp_key = f"{c}_{i}"
                    if c == 'S':
                        new_state[comp_key] -= lambdas[i] * state[comp_key]
                    elif c == 'I':
                        new_state[comp_key] += (
                            lambdas[i] * state[f"S_{i}"] - 
                            self.base_model.parameters['gamma'] * state[comp_key]
                        )
                    elif c == 'R':
                        new_state[comp_key] += self.base_model.parameters['gamma'] * state[f"I_{i}"]
            state = new_state
            for c in self.all_compartments:
                history[c].append(state[c])
        return history

    def simulate_ode(
        self, y0_dict: dict[str, float], t_range: list[float], integrator: Callable = odeint
    ) -> tuple[list[float], dict[str, list[float]]]:
        """Run ODE simulation."""
        if (
            hasattr(self, 'all_compartments') and 
            len(self.all_compartments) > 0 and 
            all('_' in c for c in self.all_compartments)
        ):
            y0 = [y0_dict[c] for c in self.all_compartments]

            def rhs(y, t):
                state = {c: y[i] for i, c in enumerate(self.all_compartments)}
                lambdas = self.compute_force_of_infection(state)
                dydt = np.zeros_like(y)
                for i in range(self.num_patches):
                    for c in self.base_model.compartments:
                        idx = self.all_compartments.index(f"{c}_{i}")
                        if c == 'S':
                            dydt[idx] = -lambdas[i] * state[f"S_{i}"]
                        elif c == 'I':
                            dydt[idx] = (
                                lambdas[i] * state[f"S_{i}"] - 
                                self.base_model.parameters['gamma'] * state[f"I_{i}"]
                            )
                        elif c == 'R':
                            dydt[idx] = self.base_model.parameters['gamma'] * state[f"I_{i}"]
                return dydt

            sol = integrator(rhs, y0, t_range)
            out = {c: sol[:, i] for i, c in enumerate(self.all_compartments)}
            return t_range, out

        y0 = [y0_dict[c] for c in self.base_model.compartments]

        def rhs(y, t):
            state = {c: y[i] for i, c in enumerate(self.base_model.compartments)}
            rates = self.base_model.compute_rates(state)
            dydt = np.zeros_like(y)
            for i, c in enumerate(self.base_model.compartments):
                for transition in self.base_model.transitions:
                    if transition['source'] == c:
                        dydt[i] -= rates[f"{c}_to_{transition['target']}"]
                    if transition['target'] == c:
                        dydt[i] += rates[f"{transition['source']}_to_{c}"]
            return dydt

        sol = integrator(rhs, y0, t_range)
        out = {c: sol[:, i] for i, c in enumerate(self.base_model.compartments)}
        return t_range, out

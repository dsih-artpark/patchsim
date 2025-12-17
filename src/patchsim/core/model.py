"""
Core model implementation for compartmental models.
"""
from typing import Any, Callable, Dict, List
from scipy.integrate import odeint
import numpy as np

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
            source = transition['from']
            target = transition['to']
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

    def get_patch_state(self, full_state: Dict[str, float], patch_idx: int) -> Dict[str, float]:
        """Get state for a specific patch."""
        return {c: full_state[f"{c}_{patch_idx}"] for c in self.base_model.compartments}

    def get_patch_population(self, state: Dict[str, float]) -> float:
        """Get total population for a patch."""
        return sum(state[c] for c in self.base_model.compartments)

    def compute_force_of_infection(self, full_state: dict[str, float], infected_compartment: str = "I") -> list[float]:
        """Compute force of infection for each patch.

        Args:
            full_state: Current state of all compartments
            infected_compartment: Name of the compartment representing infected individuals
        """
        lambdas = []
        for i in range(self.num_patches):
            if self.num_patches == 1:
                # Single patch case: local force of infection
                patch_state = self.get_patch_state(full_state, 0)
                infected = patch_state[infected_compartment]
                total_pop = self.get_patch_population(patch_state)
                force = self.base_model.parameters['beta'] * infected / total_pop
            else:
                # Multi-patch case: network-based force of infection
                force = 0
                for j in range(self.num_patches):
                    patch_state_j = self.get_patch_state(full_state, j)
                    infected_j = patch_state_j[infected_compartment]
                    pop_j = self.get_patch_population(patch_state_j)
                    force += self.network[i][j] * infected_j / pop_j
            lambdas.append(force)
        return lambdas

    def compute_derivatives(self, state: dict[str, float]) -> dict[str, float]:
        """Compute derivatives for all compartments based on transitions."""
        derivatives = {c: 0.0 for c in self.all_compartments}

        # Process each patch
        for i in range(self.num_patches):
            # Get state for this patch
            patch_state = self.get_patch_state(state, i)

            # Get rates for this patch
            if hasattr(self, "patch_parameters"):
                # get patch name using index (like PatchA, PatchB, etc.)
                patch_names = list(self.patch_parameters.keys())
                patch_name = patch_names[i] if i < len(patch_names) else None
                patch_params = {
                    **self.base_model.parameters,
                    **self.patch_parameters.get(patch_name, {})
                }

                #patch_params = self.patch_parameters.get(patch_name, self.base_model.parameters)
            else:
                patch_params = self.base_model.parameters

            # Temporarily override parameters for this patch
            old_params = self.base_model.parameters
            self.base_model.parameters = patch_params
            rates = self.base_model.compute_rates(patch_state)
            self.base_model.parameters = old_params  # restore global parameters

            # Update derivatives based on rates
            for transition in self.base_model.transitions:
                source = transition['from']
                target = transition['to']
                rate = rates[f"{source}_to_{target}"]

                # Decrease source compartment
                derivatives[f"{source}_{i}"] -= rate
                # Increase target compartment
                derivatives[f"{target}_{i}"] += rate

        return derivatives

    def simulate_discrete(self, y0_dict: dict[str, float], t_range: list[float]) -> dict[str, list[float]]:
        """Run discrete-time simulation."""
        state = y0_dict.copy()
        history = {c: [state[c]] for c in self.all_compartments}

        for _ in t_range[1:]:
            derivatives = self.compute_derivatives(state)
            new_state = {
                c: state[c] + derivatives[c]
                for c in self.all_compartments
            }
            state = new_state
            for c in self.all_compartments:
                history[c].append(state[c])

        return history

    def simulate_ode(
        self, y0_dict: dict[str, float], t_range: list[float], integrator: Callable = odeint
    ) -> tuple[list[float], dict[str, list[float]]]:
        """Run ODE simulation."""
        y0 = [y0_dict[c] for c in self.all_compartments]

        def rhs(y, t):
            state = {c: y[i] for i, c in enumerate(self.all_compartments)}
            derivatives = self.compute_derivatives(state)
            return [derivatives[c] for c in self.all_compartments]

        sol = integrator(rhs, y0, t_range)
        out = {c: sol[:, i] for i, c in enumerate(self.all_compartments)}
        return t_range, out
    
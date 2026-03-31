"""
Core model implementation for compartmental models.
"""

import re
from typing import Any, Callable, Dict

from scipy.integrate import odeint


class CompartmentalModel:
    """Base class for compartmental models."""

    def __init__(self, compartments: list[str], parameters: dict[str, float], transitions: list[dict[str, Any]]):
        """Initialize the model with compartments, parameters, and transitions."""
        self.compartments = compartments
        self.parameters = parameters
        self.transitions = transitions

    def compute_rates(self, state: dict[str, float], parameters: dict[str, float] | None = None) -> dict[str, float]:
        """Compute transition rates for each compartment.

        Args:
            state: Current compartment state
            parameters: Optional parameter override (defaults to self.parameters)
        """
        params = parameters if parameters is not None else self.parameters
        rates = {}
        for transition in self.transitions:
            transition_label = transition["transition"]
            source, target = [p.strip() for p in transition_label.split("->")]
            rate = transition["rate"]
            rate_expr = rate
            # Handle rate expressions
            if isinstance(rate_expr, str):
                # Safe evaluation: build scope from parameters and state, disable builtins
                scope = {**params, **state}
                try:
                    rate_val = eval(rate_expr, {"__builtins__": {}}, scope)
                except (KeyError, NameError, ValueError, SyntaxError, TypeError, ZeroDivisionError) as e:
                    msg = f"Invalid rate expression '{rate_expr}' in transition '{transition_label}': {e}"
                    raise ValueError(msg) from e
            else:
                rate_val = rate_expr

            # If expression already includes the source compartment, don't multiply again.
            if isinstance(rate, str) and re.search(rf"\b{re.escape(source)}\b", rate):
                flow = rate_val
            else:
                flow = rate_val * state[source]

            rates[transition_label] = flow
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
        """Compute force of infection for each patch (per-capita rate, before beta scaling).

        Args:
            full_state: Current state of all compartments
            infected_compartment: Name of the compartment representing infected individuals

        Returns:
            List of per-capita forces of infection (model_runner applies beta * FOI * S)
        """
        lambdas = []
        for i in range(self.num_patches):
            if self.num_patches == 1:
                # Single patch case: infected proportion
                patch_state = self.get_patch_state(full_state, 0)
                infected = patch_state[infected_compartment]
                total_pop = self.get_patch_population(patch_state)
                force = infected / total_pop if total_pop > 0 else 0
            else:
                # Multi-patch case: network-weighted infected proportion
                force = 0
                for j in range(self.num_patches):
                    patch_state_j = self.get_patch_state(full_state, j)
                    infected_j = patch_state_j[infected_compartment]
                    pop_j = self.get_patch_population(patch_state_j)
                    force += self.network[i][j] * (infected_j / pop_j if pop_j > 0 else 0)
            lambdas.append(force)
        return lambdas

    def compute_derivatives(self, state: dict[str, float]) -> dict[str, float]:
        """Compute derivatives for all compartments based on transitions, incorporating network-mediated FOI."""
        derivatives = {c: 0.0 for c in self.all_compartments}

        # Compute network-mediated force of infection for each patch
        lambdas = self.compute_force_of_infection(state)

        # Process each patch
        for i in range(self.num_patches):
            # Get state for this patch
            patch_state = self.get_patch_state(state, i)

            # Resolve patch parameters using canonical patch ordering when available.
            if hasattr(self, "patch_parameters"):
                patch_name = None
                if hasattr(self, "patch_names") and i < len(self.patch_names):
                    patch_name = self.patch_names[i]
                patch_params = {**self.base_model.parameters, **self.patch_parameters.get(patch_name, {})}
            else:
                patch_params = self.base_model.parameters

            # Compute rates with patch-specific parameters without mutating shared state
            rates = self.base_model.compute_rates(patch_state, parameters=patch_params)

            # Update derivatives based on rates, applying network-mediated FOI to infection transitions
            for transition in self.base_model.transitions:
                transition_label = transition["transition"]
                source, target = [p.strip() for p in transition_label.split("->")]
                rate = rates[transition_label]
                original_rate_expr = transition.get("rate", "")

                # Apply network-mediated FOI to susceptible-to-infection transitions.
                # Allow model-level override via `infection_compartments` attribute.
                infection_compartments = set(getattr(self, "infection_compartments", {"I", "E"}))
                is_infection_transition = source == "S" and target in infection_compartments

                has_network = self.num_patches > 1 and hasattr(self, "network") and self.network is not None

                if is_infection_transition and has_network:
                    # Network case: Apply network FOI (lambdas already computed above)
                    # Check if original expression includes beta term
                    beta = patch_params.get("beta", 1.0)
                    if isinstance(original_rate_expr, str) and "beta" in original_rate_expr:
                        # Rate expression includes beta; apply network FOI correction
                        adjusted_rate = beta * patch_state["S"] * lambdas[i]
                    else:
                        # Rate is already computed; apply FOI scaling
                        adjusted_rate = rate * lambdas[i] if patch_state["S"] > 0 else 0
                else:
                    # Single patch or non-infection transition: use rate as-is
                    adjusted_rate = rate

                # Decrease source compartment
                derivatives[f"{source}_{i}"] -= adjusted_rate
                # Increase target compartment
                derivatives[f"{target}_{i}"] += adjusted_rate

        return derivatives

    def simulate_discrete(self, y0_dict: dict[str, float], t_range: list[float]) -> dict[str, list[float]]:
        """Run discrete-time simulation."""
        state = y0_dict.copy()
        history = {c: [state[c]] for c in self.all_compartments}

        for _ in t_range[1:]:
            derivatives = self.compute_derivatives(state)
            new_state = {c: state[c] + derivatives[c] for c in self.all_compartments}
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

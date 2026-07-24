"""
Core model implementation for compartmental models.
"""

import re
from typing import Any, Callable, Dict

import numpy as np
from scipy.integrate import odeint

from patchsim.core.expressions import evaluate as evaluate_expression

# Negative values within this fraction of the initial population are treated as rounding error.
_NEGATIVITY_TOLERANCE = 1e-9
_MAX_EXPRESSION_IN_ERROR = 80  # max expression length shown in error messages


def _validated_time_grid(t_range: list[float]) -> "np.ndarray":
    """Return ``t_range`` as a finite, strictly increasing 1-D array of time points."""
    times = np.asarray(t_range, dtype=float)
    if times.ndim != 1:
        raise ValueError(f"Time grid must be one-dimensional; received {times.ndim} dimensions")
    if not np.all(np.isfinite(times)):
        raise ValueError("Time grid must contain only finite values")
    steps = np.diff(times)
    if steps.size and np.any(steps <= 0):
        raise ValueError(
            "Discrete simulation requires strictly increasing time points; "
            f"received a step of {steps.min()}. A zero step never advances and a "
            "negative step integrates backwards."
        )
    return times


def _validated_initial_total(y0_dict: dict[str, float]) -> float:
    """Return the total initial population, rejecting non-finite or negative values.

    Also rejects a non-finite total, which can occur when finite compartments sum past
    the float range and would otherwise make the negativity tolerance infinite.
    """
    for compartment, value in y0_dict.items():
        if not np.isfinite(value):
            raise ValueError(f"Initial value for '{compartment}' must be finite; received {value}")
        if value < 0:
            raise ValueError(f"Initial value for '{compartment}' must be non-negative; received {value}")
    total = float(sum(y0_dict.values()))
    if not np.isfinite(total):
        raise ValueError("Total initial population is not finite")
    return total


def _validated_population(
    value: float, compartment: str, time: float, step_index: int, dt: float, tolerance: float
) -> float:
    """Check one compartment value after a step, tolerating rounding residue near zero."""
    if not np.isfinite(value):
        raise ValueError(
            f"Discrete simulation diverged: '{compartment}' became {value} at "
            f"t={time} (step {step_index}). Reduce the step size."
        )
    if value < -tolerance:
        raise ValueError(
            f"Discrete simulation produced a negative population: '{compartment}' "
            f"became {value} at t={time} (step {step_index}). "
            f"The step size {dt} is too large for these rates."
        )
    # Returned unchanged rather than clipped to zero, so total population is conserved.
    return value


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
                scope = {**params, **state}
                try:
                    rate_val = evaluate_expression(rate_expr, scope)
                except ValueError as e:
                    shown = (
                        rate_expr
                        if len(rate_expr) <= _MAX_EXPRESSION_IN_ERROR
                        else f"{rate_expr[: _MAX_EXPRESSION_IN_ERROR - 3]}..."
                    )
                    msg = f"Invalid rate expression '{shown}' in transition '{transition_label}': {e}"
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

    def _adjust_infection_rate(
        self,
        patch_params: dict[str, float],
        original_rate_expr: Any,
        rate: float,
        patch_state: dict[str, float],
        lambdas: list[float],
        patch_idx: int,
        is_infection_transition: bool,
        has_network: bool,
    ) -> float:
        """Adjust infection rate for network-mediated FOI.

        Args:
            patch_params: Parameters for the current patch
            original_rate_expr: Original rate expression from transition definition
            rate: Computed rate from base model
            patch_state: Current state for the patch
            lambdas: Force of infection for each patch
            patch_idx: Current patch index
            is_infection_transition: Whether this is an infection transition
            has_network: Whether multi-patch network exists

        Returns:
            Adjusted rate incorporating network FOI if applicable
        """
        if is_infection_transition and has_network:
            # Network case: Apply network FOI (lambdas already computed)
            # Check if original expression includes beta term
            beta = patch_params.get("beta", 1.0)
            if isinstance(original_rate_expr, str) and re.search(r"\bbeta\b", original_rate_expr):
                # Rate expression includes beta; apply network FOI correction
                adjusted_rate = beta * patch_state["S"] * lambdas[patch_idx]
            else:
                # Rate is already computed; apply FOI scaling
                adjusted_rate = rate * lambdas[patch_idx] if patch_state["S"] > 0 else 0
        else:
            # Single patch or non-infection transition: use rate as-is
            adjusted_rate = rate
        return adjusted_rate

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
                elif self.patch_parameters:
                    import logging

                    logger = logging.getLogger(__name__)
                    logger.warning(
                        "patch_parameters defined but patch_names not set; "
                        "patch-specific parameters will be ignored for patch %d",
                        i,
                    )
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

                has_network = self.num_patches > 1 and self.network is not None

                # Use shared helper to adjust infection rate for network FOI
                adjusted_rate = self._adjust_infection_rate(
                    patch_params,
                    original_rate_expr,
                    rate,
                    patch_state,
                    lambdas,
                    i,
                    is_infection_transition,
                    has_network,
                )

                # Decrease source compartment
                derivatives[f"{source}_{i}"] -= adjusted_rate
                # Increase target compartment
                derivatives[f"{target}_{i}"] += adjusted_rate

        return derivatives

    def simulate_discrete(self, y0_dict: dict[str, float], t_range: list[float]) -> dict[str, list[float]]:
        """Run a discrete-time forward simulation.

        Takes one explicit Euler step per interval in ``t_range``, using that interval's
        own width. The grid names the points to simulate and need not be evenly spaced.
        The method has no stability control: a step that drives a compartment to a
        non-finite value, or significantly below zero, raises. Small negative residue from
        floating-point error is returned unchanged so that total population is conserved.

        Args:
            y0_dict: Initial state mapping for all compartment variables.
            t_range: Increasing sequence of time points to simulate.

        Returns:
            History of each compartment variable over time.

        Raises:
            ValueError: If the initial state or time grid is invalid, or if a step produces
                a non-finite or significantly negative compartment value.
        """
        state = y0_dict.copy()
        history = {c: [state[c]] for c in self.all_compartments}

        times = _validated_time_grid(t_range)
        # Scale the negativity check to the initial population. Validating here also rejects
        # a bad initial state before the early return below, not only when steps are taken.
        tolerance = _NEGATIVITY_TOLERANCE * max(abs(_validated_initial_total(y0_dict)), 1.0)

        if times.size < 2:
            return history

        for step_index, (dt, time) in enumerate(zip(np.diff(times), times[1:], strict=True), start=1):
            derivatives = self.compute_derivatives(state)
            state = {c: state[c] + derivatives[c] * float(dt) for c in self.all_compartments}
            for c in self.all_compartments:
                state[c] = _validated_population(state[c], c, time, step_index, float(dt), tolerance)
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

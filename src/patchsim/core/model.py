"""
Core model implementation for compartmental models.
"""

import re
from typing import Any, Callable, Dict

import numpy as np
from scipy.integrate import odeint

from patchsim.core.expressions import evaluate as evaluate_expression

_NEGATIVITY_ROUNDOFF_FACTOR = 64
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
    the float range and would make population-based rates invalid.
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
            f"t={time} (step {step_index}). Reduce TimeStep or the supplied interval width."
        )
    if value < -tolerance:
        raise ValueError(
            f"Discrete simulation produced a negative population: '{compartment}' "
            f"became {value} at t={time} (step {step_index}). "
            f"TimeStep or the supplied interval width ({dt}) is too large for these rates."
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

    def __init__(
        self,
        base_model: CompartmentalModel,
        num_patches: int,
        network_matrix: list[list[float]],
        groups: list[str] | None = None,
        interaction_matrix: list[list[float]] | None = None,
    ):
        """Initialize the network model."""
        self.base_model = base_model
        self.num_patches = num_patches
        self.network = network_matrix
        self.groups = list(groups or [])
        self.num_groups = len(self.groups) if self.groups else 1
        if len(set(self.groups)) != len(self.groups):
            raise ValueError("Group labels must be unique.")
        if interaction_matrix is not None and not self.groups:
            raise ValueError("An interaction matrix requires group labels.")
        self.interaction = np.asarray(
            interaction_matrix if interaction_matrix is not None else [[1.0]],
            dtype=float,
        )
        if self.interaction.shape != (self.num_groups, self.num_groups):
            raise ValueError(
                f"Interaction matrix must have shape {(self.num_groups, self.num_groups)}; "
                f"received {self.interaction.shape}."
            )
        self.all_compartments = [
            self.state_key(c, patch_idx, group_idx)
            for patch_idx in range(num_patches)
            for group_idx in range(self.num_groups)
            for c in base_model.compartments
        ]

    def state_key(self, compartment: str, patch_idx: int, group_idx: int = 0) -> str:
        """Return the internal state key for one compartment stratum."""
        if self.groups:
            return f"{compartment}_{patch_idx}_{group_idx}"
        return f"{compartment}_{patch_idx}"

    def get_patch_state(self, full_state: Dict[str, float], patch_idx: int, group_idx: int = 0) -> Dict[str, float]:
        """Get compartment state for a patch and optional group."""
        return {c: full_state[self.state_key(c, patch_idx, group_idx)] for c in self.base_model.compartments}

    def get_patch_population(self, state: Dict[str, float]) -> float:
        """Get total population for a patch."""
        return sum(state[c] for c in self.base_model.compartments)

    def compute_force_of_infection(
        self, full_state: dict[str, float], infected_compartment: str = "I"
    ) -> list[float] | list[list[float]]:
        """Compute force of infection for each patch (per-capita rate, before beta scaling).

        Args:
            full_state: Current state of all compartments
            infected_compartment: Name of the compartment representing infected individuals

        Returns:
            Per-capita forces by patch, or by patch and group for grouped models.
        """
        forces = np.zeros((self.num_patches, self.num_groups), dtype=float)
        for i in range(self.num_patches):
            for group_i in range(self.num_groups):
                for j in range(self.num_patches):
                    spatial_weight = 1.0 if self.num_patches == 1 else self.network[i][j]
                    for group_j in range(self.num_groups):
                        interaction_weight = self.interaction[group_i][group_j]
                        contributor_state = self.get_patch_state(full_state, j, group_j)
                        infected = contributor_state[infected_compartment]
                        population = self.get_patch_population(contributor_state)
                        prevalence = infected / population if population > 0 else 0.0
                        forces[i, group_i] += spatial_weight * interaction_weight * prevalence
        if self.groups:
            return forces.tolist()
        return forces[:, 0].tolist()

    def _adjust_infection_rate(
        self,
        patch_params: dict[str, float],
        original_rate_expr: Any,
        rate: float,
        patch_state: dict[str, float],
        force_of_infection: float,
        is_infection_transition: bool,
        has_mixing: bool,
    ) -> float:
        """Adjust infection rate for network-mediated FOI.

        Args:
            patch_params: Parameters for the current patch
            original_rate_expr: Original rate expression from transition definition
            rate: Computed rate from base model
            patch_state: Current state for the patch
            force_of_infection: Force of infection for the current stratum
            is_infection_transition: Whether this is an infection transition
            has_mixing: Whether spatial or group mixing is active

        Returns:
            Adjusted rate incorporating network FOI if applicable
        """
        if is_infection_transition and has_mixing:
            # Network case: Apply network FOI (lambdas already computed)
            # Check if original expression includes beta term
            beta = patch_params.get("beta", 1.0)
            if isinstance(original_rate_expr, str) and re.search(r"\bbeta\b", original_rate_expr):
                # Rate expression includes beta; apply network FOI correction
                adjusted_rate = beta * patch_state["S"] * force_of_infection
            else:
                # Rate is already computed; apply FOI scaling
                adjusted_rate = rate * force_of_infection if patch_state["S"] > 0 else 0
        else:
            # Single patch or non-infection transition: use rate as-is
            adjusted_rate = rate
        return adjusted_rate

    def compute_derivatives(self, state: dict[str, float]) -> dict[str, float]:
        """Compute derivatives for all compartments based on transitions, incorporating network-mediated FOI."""
        derivatives = {c: 0.0 for c in self.all_compartments}

        # Compute network-mediated force of infection for each patch
        lambdas = self.compute_force_of_infection(state)

        # Process each patch and optional group.
        for i in range(self.num_patches):
            for group_idx in range(self.num_groups):
                patch_state = self.get_patch_state(state, i, group_idx)

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

                rates = self.base_model.compute_rates(patch_state, parameters=patch_params)
                force = lambdas[i][group_idx] if self.groups else lambdas[i]

                for transition in self.base_model.transitions:
                    transition_label = transition["transition"]
                    source, target = [p.strip() for p in transition_label.split("->")]
                    rate = rates[transition_label]
                    original_rate_expr = transition.get("rate", "")

                    infection_compartments = set(getattr(self, "infection_compartments", {"I", "E"}))
                    is_infection_transition = source == "S" and target in infection_compartments
                    has_mixing = self.num_patches > 1 or bool(self.groups)
                    adjusted_rate = self._adjust_infection_rate(
                        patch_params,
                        original_rate_expr,
                        rate,
                        patch_state,
                        force,
                        is_infection_transition,
                        has_mixing,
                    )

                    derivatives[self.state_key(source, i, group_idx)] -= adjusted_rate
                    derivatives[self.state_key(target, i, group_idx)] += adjusted_rate

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
        expected_compartments = set(self.all_compartments)
        missing = expected_compartments - state.keys()
        extra = state.keys() - expected_compartments
        if missing or extra:
            raise ValueError(
                f"Initial state keys do not match model compartments (missing={sorted(missing)}, extra={sorted(extra)})"
            )
        history = {c: [state[c]] for c in self.all_compartments}

        times = _validated_time_grid(t_range)
        # Validate before the early return, not only when steps are taken.
        _validated_initial_total(y0_dict)

        if times.size < 2:
            return history

        for step_index, (dt, time) in enumerate(zip(np.diff(times), times[1:], strict=True), start=1):
            derivatives = self.compute_derivatives(state)
            next_state = {}
            for c in self.all_compartments:
                delta = derivatives[c] * float(dt)
                value = state[c] + delta
                tolerance = _NEGATIVITY_ROUNDOFF_FACTOR * np.finfo(float).eps * max(abs(state[c]), abs(delta), 1.0)
                next_state[c] = _validated_population(value, c, time, step_index, float(dt), tolerance)
                history[c].append(next_state[c])
            state = next_state

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

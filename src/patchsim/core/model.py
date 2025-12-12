"""
Core model implementation for compartmental models.
"""

from typing import Any, Callable, Dict, List, Optional
from scipy.integrate import odeint

import numpy as np
import matplotlib.pyplot as plt
from patchsim.core.network import Network


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




class Model:
    """
    Unified model class that manages network-based compartmental simulations.
    This class owns a Network object and provides three main capabilities: construct_ode, solve, visualize.
    Attributes:
        network: Network object containing patches and mixing matrix
        compartmental_model: Base compartmental model (e.g., SIR, SEIR)
        results: Dictionary storing simulation results after solve()
        times: Array of time points after solve()
    """
    
    def __init__(self, network: Network, compartmental_model: CompartmentalModel):
        """
        Initialize the Model.
        Args:
            network: Network object with patches and mixing matrix
            compartmental_model: Base compartmental model defining transitions
        """
        self.network = network
        self.compartmental_model = compartmental_model
        self.results = None
        self.times = None
        
        # Validate compatibility
        self._validate()
        
        # Build compartment names for all patches
        self.all_compartments = [
            f"{c}_{i}" 
            for i in range(self.network.num_patches) 
            for c in self.compartmental_model.compartments
        ]
    
    def _validate(self):
        """Validate that network patches are compatible with compartmental model."""
        # Check that all patches have the same compartments as the base model
        for patch in self.network.patches:
            if set(patch.compartments) != set(self.compartmental_model.compartments):
                raise ValueError(
                    f"Patch {patch.name} has compartments {patch.compartments} "
                    f"but model expects {self.compartmental_model.compartments}"
                )
    
    def get_patch_state(self, full_state: Dict[str, float], patch_idx: int) -> Dict[str, float]:
        """
        Get state for a specific patch.
        Args:
            full_state: Full state dictionary with keys like 'S_0', 'I_0', etc.
            patch_idx: Patch index
        Returns:
            State dictionary for this patch with keys like 'S', 'I', 'R'
        """
        return {
            c: full_state[f"{c}_{patch_idx}"] 
            for c in self.compartmental_model.compartments
        }
    
    def get_patch_population(self, state: Dict[str, float]) -> float:
        """
        Get total population for a patch state.
        Args:
            state: Patch state (not full state) 
        Returns:
            Sum of all compartments
        """
        return sum(state[c] for c in self.compartmental_model.compartments)
    
    def compute_force_of_infection(
        self, 
        full_state: dict[str, float], 
        infected_compartment: str = "I"
    ) -> list[float]:
        """
        Compute force of infection for each patch using network mixing.
        Args:
            full_state: Current state of all compartments
            infected_compartment: Name of infectious compartment   
        Returns:
            List of FOI values for each patch
        """
        lambdas = []
        for i in range(self.network.num_patches):
            if self.network.num_patches == 1:
                # Single patch case: local force of infection
                patch_state = self.get_patch_state(full_state, 0)
                infected = patch_state[infected_compartment]
                total_pop = self.get_patch_population(patch_state)
                # Use beta from patch parameters
                patch = self.network.patches[0]
                beta = patch.params.get('beta', self.compartmental_model.parameters.get('beta', 0))
                force = beta * infected / total_pop
            else:
                # Multi-patch case: network-based force of infection
                force = 0
                for j in range(self.network.num_patches):
                    patch_state_j = self.get_patch_state(full_state, j)
                    infected_j = patch_state_j[infected_compartment]
                    pop_j = self.get_patch_population(patch_state_j)
                    
                    # Network matrix[i][j] = influence from j to i
                    force += self.network.matrix[i, j] * infected_j / pop_j
            
            lambdas.append(force)
        return lambdas
    
    def compute_derivatives(self, state: dict[str, float]) -> dict[str, float]:
        """
        Compute derivatives for all compartments based on transitions.
        This is the core ODE computation that:
        1. Gets patch-specific state
        2. Uses patch-specific parameters
        3. Computes transition rates
        4. Updates compartment derivatives
        Args:
            state: Full state dictionary   
        Returns:
            Dictionary of derivatives for all compartments
        """
        derivatives = {c: 0.0 for c in self.all_compartments}
        
        # Process each patch
        for i in range(self.network.num_patches):
            patch = self.network.patches[i]
            
            # Get state for this patch
            patch_state = self.get_patch_state(state, i)
            
            # Add population to state for rate expressions
            patch_state['N'] = self.get_patch_population(patch_state)
            
            # Use patch-specific parameters
            patch_params = patch.params
            
            # Temporarily override parameters for rate computation
            old_params = self.compartmental_model.parameters
            self.compartmental_model.parameters = patch_params
            rates = self.compartmental_model.compute_rates(patch_state)
            self.compartmental_model.parameters = old_params  # restore
            
            # Update derivatives based on transitions
            for transition in self.compartmental_model.transitions:
                source = transition['from']
                target = transition['to']
                rate = rates[f"{source}_to_{target}"]
                
                # Decrease source compartment
                derivatives[f"{source}_{i}"] -= rate
                # Increase target compartment
                derivatives[f"{target}_{i}"] += rate
        
        return derivatives
    
    def construct_ode(self) -> Callable:
        """
        Construct the ODE system (right-hand side function).
        Returns:
            Callable: Function with signature f(y, t) -> dy/dt
        """
        def ode_rhs(y: np.ndarray, t: float) -> np.ndarray:
            """Right-hand side of ODE system."""
            # Convert array to dictionary
            state = {comp: y[i] for i, comp in enumerate(self.all_compartments)}
            # Compute derivatives
            derivatives = self.compute_derivatives(state)
            # Convert back to array
            return np.array([derivatives[c] for c in self.all_compartments])
        return ode_rhs
    
    def solve(
        self, 
        t_span: tuple[float, float],
        dt: float = 1.0,
        method: str = 'ode',
        integrator: Optional[Callable] = None
    ) -> tuple[np.ndarray, Dict[str, np.ndarray]]:
        """
        Solve the model and store results.
        Args:
            t_span: Time span (t_start, t_end)
            dt: Time step size
            method: 'ode' or 'discrete'
            integrator: Custom ODE integrator (default: scipy.integrate.odeint)    
        Returns:
            Tuple of (times, results_dict)   
        Side effects:
            Sets self.times and self.results
        """
        # Build initial conditions from patches
        y0_dict = {}
        for patch_idx, patch in enumerate(self.network.patches):
            for comp in self.compartmental_model.compartments:
                y0_dict[f"{comp}_{patch_idx}"] = patch.initial_conditions[comp]
        # Create time array
        t_start, t_end = t_span
        times = np.arange(t_start, t_end, dt)
        if method == 'ode':
            # Use ODE solver (matches old NetworkModel.simulate_ode)
            if integrator is None:
                integrator = odeint
            y0 = [y0_dict[c] for c in self.all_compartments]
            def rhs(y, t):
                state = {c: y[i] for i, c in enumerate(self.all_compartments)}
                derivatives = self.compute_derivatives(state)
                return [derivatives[c] for c in self.all_compartments]
            sol = integrator(rhs, y0, times)
            results = {c: sol[:, i] for i, c in enumerate(self.all_compartments)} 
        elif method == 'discrete':
            # Discrete time stepping (matches old NetworkModel.simulate_discrete)
            state = y0_dict.copy()
            history = {c: [state[c]] for c in self.all_compartments}
            for _ in times[1:]:
                derivatives = self.compute_derivatives(state)
                new_state = {
                    c: state[c] + derivatives[c]
                    for c in self.all_compartments
                }
                state = new_state
                for c in self.all_compartments:
                    history[c].append(state[c])
            # Convert lists to arrays
            results = {c: np.array(v) for c, v in history.items()}
        else:
            raise ValueError(f"Unknown method: {method}. Use 'ode' or 'discrete'")
        # Store results
        self.times = times
        self.results = results
        return times, results
    
    def visualize(
        self,
        compartments: Optional[List[str]] = None,
        patch_names: Optional[List[str]] = None,
        figsize: tuple = (12, 8),
        save_path: Optional[str] = None,
        **kwargs
    ) -> plt.Figure:
        """
        Visualize simulation results.
        Args:
            compartments: List of compartments to plot (default: all)
            patch_names: List of patch names to plot (default: all)
            figsize: Figure size
            save_path: Path to save figure (optional)
            **kwargs: Additional arguments passed to plt.plot()    
        Returns:
            matplotlib Figure object
        """
        if self.results is None:
            raise ValueError("No results to plot. Run solve() first.")
        # Default to all compartments except Susceptible (often dominates scale)
        if compartments is None:
            compartments = [c for c in self.compartmental_model.compartments if c != 'S']
        # Default to all patches
        if patch_names is None:
            patch_names = [p.name for p in self.network.patches]
        # Create subplots (one per patch)
        num_patches = len(patch_names)
        fig, axes = plt.subplots(
            num_patches, 1, 
            figsize=figsize, 
            squeeze=False
        )
        # Get patch name to index mapping
        name_to_idx = {p.name: p.id for p in self.network.patches}
        for plot_idx, patch_name in enumerate(patch_names):
            patch_idx = name_to_idx[patch_name]
            ax = axes[plot_idx, 0]
            # Plot each compartment for this patch
            for comp in compartments:
                comp_key = f"{comp}_{patch_idx}"
                if comp_key in self.results:
                    ax.plot(
                        self.times, 
                        self.results[comp_key],
                        label=comp,
                        **kwargs
                    ) 
            ax.set_xlabel('Time')
            ax.set_ylabel('Population')
            ax.set_title(f'Patch: {patch_name}')
            ax.legend()
            ax.grid(True, alpha=0.3)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        return fig



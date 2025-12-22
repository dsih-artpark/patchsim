from typing import List, Dict
import numpy as np

from patchsim.core.patch import Patch

class Network:
    def __init__(self, patches: List[Patch], adjacency: np.ndarray):
        """
        patches: list of Patch objects
        adjacency: num_patches x num_patches matrix with infection weights
                   adjacency[i,j] = influence of patch j on patch i
        """
        self.patches = patches
        self.adj = adjacency
        self.num_patches = len(patches)

    def get_full_state(self) -> Dict[str, float]:
        """Flatten all patch states into a dict with keys 'Compartment_i'."""
        state = {}
        for p in self.patches:
            for c in p.compartments:
                state[f"{c}_{p.id}"] = p.state[c]
        return state

    def compute_force_of_infection(self, infected_comp: str = "I") -> List[float]:
        """Compute lambda_i for each patch i considering other patches."""
        lambdas = []
        for i, patch_i in enumerate(self.patches):
            force = 0.0
            for j, patch_j in enumerate(self.patches):
                if patch_j.population > 0:
                    force += self.adj[i, j] * patch_j.state[infected_comp] / patch_j.population
            lambdas.append(force)
        return lambdas

    def compute_derivatives(self) -> Dict[str, float]:
        """Compute dX/dt for all patches using transitions + network force of infection."""
        derivs = {}
        full_state = self.get_full_state()
        lambdas = self.compute_force_of_infection()
        for idx, patch in enumerate(self.patches):
            # Pass network force of infection to patch for S->I transitions
            rates = patch.compute_transition_rates(force_of_infection=lambdas[idx])
            for t in patch.transitions:
                s = t['from']
                r = t['to']
                key = f"{s}_to_{r}"
                rate = rates[key]
                derivs[f"{s}_{patch.id}"] = derivs.get(f"{s}_{patch.id}", 0) - rate
                derivs[f"{r}_{patch.id}"] = derivs.get(f"{r}_{patch.id}", 0) + rate
        return derivs

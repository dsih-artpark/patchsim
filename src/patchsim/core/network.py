import numpy as np
from typing import Dict
from patchsim.core.patch import Patch


class Network:
    def __init__(self, patches: list[Patch], adjacency: np.ndarray):
        self.patches = patches
        self.adjacency = adjacency
        self.num_patches = len(patches)

        self.validate()

    def validate(self):
        if self.adjacency.shape != (self.num_patches, self.num_patches):
            raise ValueError("Adjacency matrix dimension mismatch")

    def get_state_vector(self) -> Dict[str, float]:
        state = {}
        for p in self.patches:
            for c, v in p.compartments.items():
                state[f"{c}_{p.id}"] = v
        return state

    def compute_derivatives(self) -> Dict[str, float]:
        derivatives = {}

        for p in self.patches:
            rates = p.compute_transition_rates()

            for t in p.transitions:
                src = t["from"]
                tgt = t["to"]
                rate = rates[f"{src}_to_{tgt}"]

                derivatives[f"{src}_{p.id}"] = (
                    derivatives.get(f"{src}_{p.id}", 0) - rate
                )
                derivatives[f"{tgt}_{p.id}"] = (
                    derivatives.get(f"{tgt}_{p.id}", 0) + rate
                )

        return derivatives

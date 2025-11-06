from typing import List
import numpy as np
from patchsim.core.patch import Patch

class Network:
    """
    Represents a multi-patch interaction network, with weighted connections
    between patches and helper lookup methods.
    """

    def __init__(self, patches: List[Patch], matrix):
        self.patches = patches
        self.matrix = np.array(matrix, dtype=float)
        self.num_patches = len(patches)

        self.validate()
        # name → id map
        self.name_to_id = {p.name: p.id for p in patches}

    def validate(self):
        """Validate network structure and matrix consistency."""
        if self.matrix.shape != (self.num_patches, self.num_patches):
            raise ValueError(
                f"Network matrix must be {self.num_patches}×{self.num_patches}, "
                f"got {self.matrix.shape}."
            )
        if np.any(self.matrix < 0):
            raise ValueError("Network weights must be non-negative.")
        # Optional rule: diagonal = 0 (no self-contact mixing)
        # if np.any(np.diag(self.matrix) != 0):
        #     raise ValueError("Diagonal of network matrix should be zero.")
        return True

    def get_patch(self, name: str) -> Patch:
        """Get Patch object by its name."""
        return self.patches[self.name_to_id[name]]

    def get_by_id(self, pid: int) -> Patch:
        return self.patches[pid]

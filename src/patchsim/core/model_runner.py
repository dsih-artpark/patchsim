import numpy as np
from scipy.integrate import odeint
from patchsim.core.network import Network
from patchsim.utils.viz import plot_patch_subplots

class Model:
    def __init__(self, network: Network):
        self.network = network
        self.results = None
        self.t_range = None
    
    def construct_ode(self):
        all_vars = []
        for p in self.network.patches:
            for c in p.compartments:
                all_vars.append(f"{c}_{p.id}")
        def rhs(y, t):
            idx = {name: i for i, name in enumerate(all_vars)}

            # update patch compartments
            for p in self.network.patches:
                for c in p.compartments:
                    p.compartments[c] = y[idx[f"{c}_{p.id}"]]

            derivs = self.network.compute_derivatives()
            return [derivs.get(v, 0.0) for v in all_vars]
        return rhs, all_vars

    def solve(self, t_range):
        rhs, variables = self.construct_ode()
        self.t_range = t_range

        y0 = []
        for p in self.network.patches:
            for c in p.compartments:
                y0.append(p.compartments[c])

        sol = odeint(rhs, y0, t_range)
        self.results = {variables[i]: sol[:, i] for i in range(len(variables))}
        return self.results

    def visualize(self, output_dir: str, model_name: str):
        if self.results is None:
            raise RuntimeError("Run solve() before visualize()")

        patch_names = [p.name for p in self.network.patches]
        plot_patch_subplots(
            self.t_range, self.results, patch_names, output_dir, model_name
        )

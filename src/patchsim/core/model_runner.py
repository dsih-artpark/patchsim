import numpy as np
from scipy.integrate import odeint

class Model:
    """
    High-level simulation model.
    Owns the Network and builds/solves the ODE.
    """

    def __init__(self, network_model, compartments):
        self.network = network_model
        self.compartments = compartments
        self.all_vars = self.network.all_compartments

    def construct_ode(self):
        def rhs(y, t):
            state = {v: y[i] for i, v in enumerate(self.all_vars)}
            dydt = {v: 0.0 for v in self.all_vars}

            # network infection term
            lambdas = self.network.compute_force_of_infection(state)
            beta = self.network.base_model.parameters["beta"]

            for i in range(self.network.num_patches):
                S = state[f"S_{i}"]
                infection = beta * S * lambdas[i]
                dydt[f"S_{i}"] -= infection
                dydt[f"I_{i}"] += infection

            # local transitions
            for i in range(self.network.num_patches):
                patch_state = {
                    c: state[f"{c}_{i}"] for c in self.compartments
                }

                rates = self.network.base_model.compute_rates(patch_state)

                for key, rate in rates.items():
                    src, tgt = key.split("_to_")
                    dydt[f"{src}_{i}"] -= rate
                    dydt[f"{tgt}_{i}"] += rate

            return [dydt[v] for v in self.all_vars]

        return rhs

    def solve(self, y0, t_range):
        rhs = self.construct_ode()
        y0_vec = [y0[v] for v in self.all_vars]
        sol = odeint(rhs, y0_vec, t_range)
        return {v: sol[:, i] for i, v in enumerate(self.all_vars)}

    def visualize(self, t, results, patches, outdir, model_name):
        from patchsim.utils.viz import plot_patch_subplots
        plot_patch_subplots(t, results, patches, outdir, model_name)

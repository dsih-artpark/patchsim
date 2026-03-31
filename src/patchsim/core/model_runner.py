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

            # Compute network-based force of infection (per-capita, without beta)
            lambdas = self.network.compute_force_of_infection(state)
            beta = self.network.base_model.parameters.get("beta", 0.0)

            # Apply transitions for all patches
            for i in range(self.network.num_patches):
                patch_state = {c: state[f"{c}_{i}"] for c in self.compartments}
                rates = self.network.base_model.compute_rates(patch_state)

                for key, rate in rates.items():
                    src, tgt = [p.strip() for p in key.split("->")]

                    # For S→I transitions, scale by network force of infection
                    if src == "S" and tgt == "I":
                        # rate from compute_rates is beta*S (before network scaling)
                        # Apply network FOI: lambda_i = network-weighted infected proportion
                        adjusted_rate = beta * state[f"S_{i}"] * lambdas[i]
                    else:
                        adjusted_rate = rate

                    dydt[f"{src}_{i}"] -= adjusted_rate
                    dydt[f"{tgt}_{i}"] += adjusted_rate

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

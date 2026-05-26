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

            # Apply transitions for all patches
            for i in range(self.network.num_patches):
                patch_state = {c: state[f"{c}_{i}"] for c in self.compartments}
                rates = self.network.base_model.compute_rates(patch_state)
                patch_params = self.network.base_model.parameters

                for key, rate in rates.items():
                    src, tgt = [p.strip() for p in key.split("->")]

                    # Determine if this is an infection transition
                    infection_compartments = set(getattr(self.network, "infection_compartments", {"I", "E"}))
                    is_infection_transition = src == "S" and tgt in infection_compartments

                    # Use network helper to apply FOI adjustment consistently
                    has_network = self.network.num_patches > 1 and self.network.network is not None
                    adjusted_rate = self.network._adjust_infection_rate(
                        patch_params, key, rate, patch_state, lambdas, i, is_infection_transition, has_network
                    )

                    dydt[f"{src}_{i}"] -= adjusted_rate
                    dydt[f"{tgt}_{i}"] += adjusted_rate

            return [dydt[v] for v in self.all_vars]

        return rhs

    def solve(self, y0, t_range):
        rhs = self.construct_ode()
        # Validate all required variables are present in y0
        missing = [v for v in self.all_vars if v not in y0]
        if missing:
            raise ValueError(f"Missing initial values for: {missing}")
        y0_vec = [y0[v] for v in self.all_vars]
        sol = odeint(rhs, y0_vec, t_range)
        return {v: sol[:, i] for i, v in enumerate(self.all_vars)}

    def visualize(self, t, results, patches, outdir, model_name):
        from patchsim.utils.viz import plot_patch_subplots

        plot_patch_subplots(t, results, patches, outdir, model_name, compartments=self.compartments)

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
            dydt = self.network.compute_derivatives(state)
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

        plot_patch_subplots(
            t,
            results,
            patches,
            outdir,
            model_name,
            compartments=self.compartments,
            groups=self.network.groups,
        )

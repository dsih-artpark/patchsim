from scipy.integrate import odeint


class Model:
    """High-level ODE simulation wrapper.

    This class converts a :class:`~patchsim.core.model.NetworkModel` into a
    solvable system of ordinary differential equations and provides a small
    convenience layer around integration and visualization.
    """

    def __init__(self, network_model, compartments):
        """Create a solver wrapper for a network model.

        Args:
            network_model: Network model containing compartments, parameters,
                and patch coupling.
            compartments: Ordered list of compartment names for a single patch.
        """
        self.network = network_model
        self.compartments = compartments
        self.all_vars = self.network.all_compartments

    def construct_ode(self):
        """Build the right-hand side function used by the ODE integrator.

        Returns:
            Callable ``rhs(y, t)`` that computes derivatives for all state
            variables in the flattened multi-patch system.
        """
        def rhs(y, t):
            state = {v: y[i] for i, v in enumerate(self.all_vars)}
            dydt = {v: 0.0 for v in self.all_vars}

            # Compute network-based force of infection (per-capita, without beta)
            lambdas = self.network.compute_force_of_infection(state)

            # Apply transitions for all patches
            for i in range(self.network.num_patches):
                patch_state = {c: state[f"{c}_{i}"] for c in self.compartments}
                if hasattr(self.network, "patch_parameters"):
                    patch_name = None
                    if hasattr(self.network, "patch_names") and i < len(self.network.patch_names):
                        patch_name = self.network.patch_names[i]
                    patch_params = {
                        **self.network.base_model.parameters,
                        **self.network.patch_parameters.get(patch_name, {}),
                    }
                else:
                    patch_params = self.network.base_model.parameters

                rates = self.network.base_model.compute_rates(patch_state, parameters=patch_params)

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
        """Integrate the ODE system for the provided initial state.

        Args:
            y0: Initial state mapping for every flattened compartment variable.
            t_range: Sequence of time points to integrate over.

        Returns:
            Mapping from variable name to the simulated trajectory over time.

        Raises:
            ValueError: If the initial state is missing any required variables.
        """
        rhs = self.construct_ode()
        # Validate all required variables are present in y0
        missing = [v for v in self.all_vars if v not in y0]
        if missing:
            raise ValueError(f"Missing initial values for: {missing}")
        y0_vec = [y0[v] for v in self.all_vars]
        sol = odeint(rhs, y0_vec, t_range)
        return {v: sol[:, i] for i, v in enumerate(self.all_vars)}

    def visualize(self, t, results, patches, outdir, model_name):
        """Render and save the patch-level time series plot.

        Args:
            t: Time values corresponding to ``results``.
            results: Mapping of flattened compartment variables to trajectories.
            patches: Ordered list of patch names.
            outdir: Output directory for the generated figure.
            model_name: Model name used to name the output image.
        """
        from patchsim.utils.viz import plot_patch_subplots

        plot_patch_subplots(t, results, patches, outdir, model_name, compartments=self.compartments)

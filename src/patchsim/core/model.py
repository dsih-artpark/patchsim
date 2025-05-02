import numpy as np

class Population:
    def __init__(self, patch_population, compartments, patch_id=None):
        """
        Initializes the population with compartments and a patch population.
        
        :param patch_population: Total population for this patch at the start of the simulation.
        :param compartments: A dictionary of compartments and their counts (e.g., susceptible, infected, etc.).
        :param patch_id: An optional identifier for this specific patch (e.g., district name or ID).
        """
        self.patch_population = patch_population
        self.compartments = compartments
        self.patch_id = patch_id
        self._validate_population()

    def _validate_population(self):
        """
        Validates that the total population equals the sum of the compartments at the start.
        """
        total_compartment_sum = sum(self.compartments.values())
        if total_compartment_sum != self.patch_population:
            raise ValueError(f"Compartment sum ({total_compartment_sum}) must equal the patch population ({self.patch_population}).")

    def update_population(self, new_population):
        """
        Updates the patch population over time (could vary based on external factors).

        :param new_population: New total population for the patch.
        """
        self.patch_population = new_population

    def update_compartment(self, compartment_name, new_value):
        """
        Updates the value of a specific compartment.
        
        :param compartment_name: The name of the compartment to update.
        :param new_value: The new value for the compartment.
        """
        if compartment_name not in self.compartments:
            raise ValueError(f"Compartment {compartment_name} does not exist.")
        self.compartments[compartment_name] = new_value

    def get_compartment(self, compartment_name):
        """
        Retrieves the current value of a specific compartment.
        
        :param compartment_name: The name of the compartment to retrieve.
        :return: The current value of the compartment.
        """
        return self.compartments.get(compartment_name, None)

    def __str__(self):
        return f"Patch {self.patch_id}: Population: {self.patch_population}, Compartments: {self.compartments}"


class CompartmentalModel:
    def __init__(self, compartments, parameters, transitions):
        self.compartments = compartments
        self.parameters = parameters
        self.transitions = transitions

    def compute_transition_rates(self, state, extras=None):
        local_env = {**state, **self.parameters}
        if extras:
            local_env.update(extras)

        deltas = {c: 0.0 for c in self.compartments}
        for tr in self.transitions:
            src = tr['from']
            dst = tr['to']
            rate_expr = tr['rate']
            rate = eval(rate_expr, {}, local_env)

            if src:
                deltas[src] -= rate
            if dst:
                deltas[dst] += rate

        return deltas

    def ode_rhs(self, y, t, extras_fn=None):
        """
        Returns the ODE right-hand side for use with ODE solvers.
        y: list/array of compartment values in the order of self.compartments
        t: current time (passed by ODE solver)
        extras_fn: optional function to provide extra variables (e.g., force of infection)
        """
        state = {c: y[i] for i, c in enumerate(self.compartments)}
        extras = extras_fn(t, y) if extras_fn else None
        deltas = self.compute_transition_rates(state, extras)
        return [deltas[c] for c in self.compartments]


class NetworkModel:
    def __init__(self, base_model, num_patches, network_matrix):
        self.base_model = base_model
        self.num_patches = num_patches
        self.network = np.array(network_matrix)

        self.all_compartments = [
            f"{c}_{i}" for i in range(num_patches) for c in base_model.compartments
        ]

    def compute_force_of_infection(self, full_state):
        lambdas = []
        for i in range(self.num_patches):
            force = 0
            for j in range(self.num_patches):
                I_j = full_state[f"I_{j}"]
                N_j = sum(full_state[f"{c}_{j}"] for c in self.base_model.compartments)
                force += self.network[i][j] * I_j / N_j
            lambdas.append(force)
        return lambdas

    def simulate_discrete(self, y0_dict, t_range):
        state = y0_dict.copy()
        history = {c: [state[c]] for c in self.all_compartments}

        for t in t_range[1:]:
            new_state = state.copy()
            lambdas = self.compute_force_of_infection(state)

            for i in range(self.num_patches):
                patch_state = {c: state[f"{c}_{i}"] for c in self.base_model.compartments}
                extras = {f"lambda_i": lambdas[i]}

                deltas = self.base_model.compute_transition_rates(patch_state, extras)

                for c, delta in deltas.items():
                    key = f"{c}_{i}"
                    new_state[key] += delta

            # Clamp negatives
            state = {k: max(v, 0) for k, v in new_state.items()}
            for c in self.all_compartments:
                history[c].append(state[c])

        return t_range, history

    def simulate_ode(self, y0_dict, t_range, odesolver, extras_fn=None):
        """
        Simulate the network model using ODE integration.
        y0_dict: dict of initial compartment values (e.g., {'S_0': ..., 'I_0': ..., ...})
        t_range: time points
        odesolver: ODE solver function (e.g., scipy.integrate.odeint)
        extras_fn: optional function to provide extra variables (e.g., force of infection)
        """
        y0 = [y0_dict[c] for c in self.all_compartments]
        def rhs(y, t):
            state = {c: y[i] for i, c in enumerate(self.all_compartments)}
            lambdas = self.compute_force_of_infection(state)
            dydt = []
            for i in range(self.num_patches):
                patch_state = {c: state[f"{c}_{i}"] for c in self.base_model.compartments}
                extras = {"lambda_i": lambdas[i]}
                deltas = self.base_model.compute_transition_rates(patch_state, extras)
                for c in self.base_model.compartments:
                    dydt.append(deltas[c])
            return dydt
        sol = odesolver(rhs, y0, t_range)
        # Convert solution to history dict
        sol = np.array(sol)
        history = {c: sol[:, i] for i, c in enumerate(self.all_compartments)}
        return t_range, history
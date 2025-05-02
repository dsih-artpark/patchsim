class Simulation:
    """
    General simulation class for both discrete (difference equation) and ODE-based compartmental models.
    """
    def __init__(
        self,
        model: object = None,
        network: object = None,
        y0: object = None,
        t_range: list = None,
        mode: str = 'discrete',
        odesolver: object = None,
        interventions: list = None
    ):
        """
        Args:
            model: CompartmentalModel instance (for discrete) or callable (for ODE)
            network: NetworkModel instance (for discrete metapopulation)
            y0: Initial state (dict for discrete, list/array for ODE)
            t_range: List or array of time points
            mode: 'discrete' or 'ode'
            odesolver: Function to solve ODEs (e.g., scipy.integrate.odeint)
            interventions: List of intervention objects (optional)
        """
        self.model = model
        self.network = network
        self.y0 = y0
        self.t_range = t_range
        self.mode = mode
        self.odesolver = odesolver
        self.interventions = interventions or []

    def set_state(self, state):
        self.y0 = state

    def apply_interventions(self, state, t=None):
        for intervention in self.interventions:
            intervention.apply(state, t)

    def run(self, **kwargs):
        if self.mode == 'discrete':
            if self.network is not None:
                t, out = self.network.simulate_discrete(self.y0, self.t_range)
                return t, out
            else:
                raise NotImplementedError("Discrete single-patch not implemented.")
        elif self.mode == 'ode':
            if self.odesolver is not None:
                y0 = self.y0
                t = self.t_range
                args = kwargs.get('args', tuple())
                sol = self.odesolver(self.model, y0, t, *args)
                return t, sol
            else:
                raise ValueError("ODESolver must be provided for ODE mode.")
        else:
            raise ValueError(f"Unknown simulation mode: {self.mode}")

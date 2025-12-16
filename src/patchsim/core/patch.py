from typing import Dict, List

class Patch:
    def __init__(
        self,
        patch_id: int,
        name: str,
        population: int,
        compartments: List[str],
        initial_conditions: Dict[str, float],
        parameters: Dict[str, float],
        transitions: List[Dict]
    ):
        self.id = patch_id
        self.name = name
        self.population = population
        self.compartments = compartments
        self.state = initial_conditions.copy()
        self.parameters = parameters.copy()
        self.transitions = transitions.copy()
        self.validate_patch()

    def validate_patch(self):
        # Check population consistency
        total = sum(self.state[c] for c in self.compartments)
        if abs(total - self.population) > 1e-6:
            raise ValueError(f"Patch {self.name}: compartment counts {total} do not sum to population {self.population}")
        # Check non-negative compartments
        if any(self.state[c] < 0 for c in self.compartments):
            raise ValueError(f"Patch {self.name}: negative compartment values detected")

    def compute_transition_rates(self):
        rates = {}
        for t in self.transitions:
            source = t['from']
            target = t['to']
            rate_expr = t['rate']
            # Replace parameters and compartment values
            if isinstance(rate_expr, str):
                for p, v in self.parameters.items():
                    rate_expr = rate_expr.replace(p, str(v))
                for c, v in self.state.items():
                    rate_expr = rate_expr.replace(c, str(v))
                rate_val = eval(rate_expr)
            else:
                rate_val = rate_expr
            rates[f"{source}_to_{target}"] = rate_val * self.state[source]
        return rates

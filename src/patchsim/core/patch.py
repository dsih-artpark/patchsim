
from typing import Dict, Any

class Patch:
    def __init__(
        self,
        patch_id: int,
        name: str,
        population: float,
        compartments: Dict[str, float],
        parameters: Dict[str, float],
        transitions: list[Dict[str, Any]],
    ):
        self.id = patch_id
        self.name = name
        self.population = population
        self.compartments = compartments
        self.parameters = parameters
        self.transitions = transitions

        self.validate()

    def validate(self) -> None:
        if self.population <= 0:
            raise ValueError(f"Patch {self.name}: population must be positive")

        if any(v < 0 for v in self.compartments.values()):
            raise ValueError(f"Patch {self.name}: negative compartment value")

        total = sum(self.compartments.values())
        if abs(total - self.population) > 1e-6:
            raise ValueError(
                f"Patch {self.name}: compartments do not sum to population"
            )

        for t in self.transitions:
            if t["from"] not in self.compartments or t["to"] not in self.compartments:
                raise ValueError(
                    f"Patch {self.name}: invalid transition {t}"
                )

    def compute_transition_rates(self) -> Dict[str, float]:
        rates = {}
        for t in self.transitions:
            expr = t["rate"]
            if isinstance(expr, str):
                local = {**self.parameters, **self.compartments}
                rate_value = eval(expr, {}, local)
            else:
                rate_value = float(expr)

            rates[f"{t['from']}_to_{t['to']}"] = rate_value

        return rates

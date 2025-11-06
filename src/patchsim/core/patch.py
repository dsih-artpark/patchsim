from typing import List, Dict

class Patch:
    """
    Represents a single patch within a multi-patch epidemiological model.
    Stores population, compartments, parameters, transitions, and initial conditions.
    Provides built-in validation to guarantee model integrity.
    """
    def __init__(
        self,
        id: int,
        name: str,
        population: float,
        compartments: List[str],
        params: Dict,
        initial_conditions: Dict[str, float],
        transitions: List[Dict]
    ):
        self.id = id
        self.name = name
        self.population = population
        self.compartments = compartments
        self.params = params
        self.initial_conditions = initial_conditions
        self.transitions = transitions

        self.validate()

    def validate(self):
        """Validate all patch attributes."""
        # Population must be positive
        if self.population <= 0:
            raise ValueError(f"Patch {self.name}: population must be positive.")
        # Compartments must be a non-empty list
        if not isinstance(self.compartments, list) or len(self.compartments) == 0:
            raise ValueError(f"Patch {self.name}: compartments must be a non-empty list.")
        # Initial conditions provided for all compartments
        for c in self.compartments:
            if c not in self.initial_conditions:
                raise ValueError(
                    f"Patch {self.name}: missing initial condition for compartment '{c}'."
                )
            if self.initial_conditions[c] < 0:
                raise ValueError(
                    f"Patch {self.name}: negative initial value for compartment '{c}'."
                )
        # Sum of initial conditions must match population
        total = sum(self.initial_conditions[c] for c in self.compartments)
        if abs(total - self.population) > 1e-6:
            raise ValueError(
                f"Patch {self.name}: initial conditions sum to {total}, "
                f"but population is {self.population}."
            )
        # Transition checks
        for tr in self.transitions:
            if "from" not in tr or "to" not in tr or "rate" not in tr:
                raise ValueError(f"Patch {self.name}: invalid transition {tr}")

            if tr["from"] not in self.compartments or tr["to"] not in self.compartments:
                raise ValueError(
                    f"Patch {self.name}: transition {tr['from']}→{tr['to']} "
                    "uses compartments that are not defined."
                )
        return True
    
    def get_initial_compartment_vector(self):
        """Return a vector [S, I, R, ...] preserving compartment order."""
        return [self.initial_conditions[c] for c in self.compartments]

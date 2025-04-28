import numpy as np
from scipy.integrate import odeint


def initialize_compartments(patch, seeds):
    """
    Initialize compartments for the SIR model.

    Args:
        patch (dict): Patch data containing population (N0).
        seeds (list): Seeding infection data.

    Returns:
        dict: Initialized compartments (S, I, R).
    """
    N0 = patch["N0"]
    I0 = 0  # Default infected population
    for seed in seeds:
        if seed["region"] == patch["region"]:
            I0 = seed["seed_count"]
            break
    S0 = N0 - I0
    R0 = 0
    return {"S": S0, "I": I0, "R": R0}


def sir_ode(y, t, beta, gamma):
    """
    Simple SIR model ODEs.

    Args:
        y (list): Current state [S, I, R].
        t (float): Time step.
        beta (float): Infection rate.
        gamma (float): Recovery rate.

    Returns:
        list: Derivatives [dS/dt, dI/dt, dR/dt].
    """
    S, I, R = y
    N = S + I + R
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return [dSdt, dIdt, dRdt]


def run_model(initial_conditions, days, config):
    """
    Run the SIR model simulation.

    Args:
        initial_conditions (dict): Initial conditions for the compartments (S, I, R).
        days (int): Number of days to simulate.
        config (dict): Configuration parameters.

    Returns:
        tuple: Time points and simulation results for [S, I, R] over time.
    """
    beta = config["InfectionRate"]
    gamma = config["RecoveryRate"]

    S0 = initial_conditions["S"]
    I0 = initial_conditions["I"]
    R0 = initial_conditions["R"]

    t = np.linspace(0, days, days + 1)
    y0 = [S0, I0, R0]
    solution = odeint(sir_ode, y0, t, args=(beta, gamma))
    return t, [{"S": S, "I": I, "R": R} for S, I, R in solution]
import warnings

import numpy as np
import logging
from tqdm import tqdm
import random

warnings.filterwarnings('ignore')


def sirsv_model_with_weibull_random_vaccination(params, scenario, random_seed=42, diagnosis=None):
    np.random.seed(random_seed)
    random.seed(random_seed)

    # Extract parameters
    beta = params['beta']
    gamma = params['gamma']
    vax_rate = params['vax_rate']
    weibull_shape_vax = params['weibull_shape_vax']
    weibull_scale_vax = params['weibull_scale_vax']
    weibull_shape_rec = params['weibull_shape_rec']
    weibull_scale_rec = params['weibull_scale_rec']
    days = params['days']
    seed_rate = params['seed_rate']
    vax_period = params['vax_period']
    vax_duration = params['vax_duration']
    start_vax_day = params['start_vax_day']

    # Initial conditions
    S0, I0, R0, V0 = params['S0'], params['I0'], params['R0'], params['V0']
    N = S0 + I0 + R0 + V0

    S, I, R, V = [np.zeros(days) for _ in range(4)]
    S[0], I[0], R[0], V[0] = S0, I0, R0, V0

    decay_times_vax = []
    decay_times_rec = []

    # Seed initial vaccinated and recovered individuals' waning times
    if V0 > 0:
        decay_times_vax.extend((weibull_scale_vax * np.random.weibull(weibull_shape_vax, int(V0))).astype(int).tolist())
    if R0 > 0:
        decay_times_rec.extend((weibull_scale_rec * np.random.weibull(weibull_shape_rec, int(R0))).astype(int).tolist())

    round_counter = 0

    logging.info(f"Starting simulation for scenario: {scenario}")

    for t in tqdm(range(1, days), desc=f"Running {scenario} simulation", unit="day"):
        new_seeds = min(seed_rate, S[t-1])

        # VACCINATION ROUND
        if t == start_vax_day or (t > start_vax_day and (t - start_vax_day) % vax_period == 0):
            round_counter += 1
            to_vaccinate = min(vax_rate * S[t-1], S[t-1])
            logging.info(f"Round {round_counter} start day: {t}")

            # Calculate the number of vaccinations to reset, considering the vaccination period
            num_vax_to_reset = int(min(vax_rate * vax_period * V[t-1], V[t-1]))
            if num_vax_to_reset > 0 and len(decay_times_vax) > 0:
                num_vax_to_reset = min(num_vax_to_reset, len(decay_times_vax))

                # Method 1: Randomly select decay times to reset
                indices_to_reset = random.sample(range(len(decay_times_vax)), num_vax_to_reset)

                # Method 2: Sort decay times and select indices with the lowest decay times
                # sorted_indices = np.argsort(decay_times_vax)
                # indices_to_reset = sorted_indices[:num_vax_to_reset]

                for i in indices_to_reset:
                    decay_times_vax[i] = weibull_scale_vax * np.random.weibull(weibull_shape_vax)
                logging.info(f"Day {t}: Re-vaccination reset: {num_vax_to_reset} decay times reset")

            # if diagnosis:
                # plot_histogram(decay_times_vax, decay_times_rec, scenario, round_counter, start=True)

        # Check if it's within a vaccination period
        is_vax_period = (t >= start_vax_day) and ((t - start_vax_day) % vax_period < vax_duration)
        if is_vax_period:
            new_vaccinations = int(to_vaccinate)
            logging.info(f"Day {t}: Daily vaccinations: {new_vaccinations}")

            # Update compartments for new vaccinations
            if new_vaccinations > 0:
                new_susceptible_decay_times = (weibull_scale_vax * np.random.weibull(weibull_shape_vax, new_vaccinations)).astype(int).tolist()
                initial_len = len(decay_times_vax)
                decay_times_vax.extend(new_susceptible_decay_times)
                updated_len = len(decay_times_vax)
                logging.info(f"Day {t}: New vaccinations: {new_vaccinations}, Length of decay_times_vax: {initial_len}, Length of decay_times_vax: {updated_len}")
            else:
                S[t] = S[t-1]
                R[t] = R[t-1]
                V[t] = V[t-1]

        else:
            new_vaccinations = 0

        # Calculate transitions
        new_infections = beta * S[t-1] * I[t-1] / N + new_seeds
        new_recoveries = gamma * I[t-1]

        # Update compartments
        S[t] = S[t-1] - new_infections - new_vaccinations
        I[t] = I[t-1] + new_infections - new_recoveries
        R[t] = R[t-1] + new_recoveries
        V[t] = V[t-1] + new_vaccinations

        if new_recoveries > 0:
            new_recovered_decay_times = (weibull_scale_rec * np.random.weibull(weibull_shape_rec, int(new_recoveries))).astype(int).tolist()
            decay_times_rec.extend(new_recovered_decay_times)

        # IMMUNITY WANING
        logging.info(f"Day {t}: Before waning: Length of decay_times_vax={len(decay_times_vax)}, Length of decay_times_rec={len(decay_times_rec)}")

        decay_times_vax_before = len(decay_times_vax)
        decay_times_vax = [x - 1 for x in decay_times_vax if x > 0]
        num_waned_vax = decay_times_vax_before - len(decay_times_vax)

        decay_times_rec_before = len(decay_times_rec)
        decay_times_rec = [x - 1 for x in decay_times_rec if x > 0]
        num_waned_rec = decay_times_rec_before - len(decay_times_rec)

        logging.info(f"Day {t}: Waned vaccinated: {num_waned_vax}, Waned recovered: {num_waned_rec}")
        logging.info(f"Day {t}: After waning: Length of decay_times_vax={len(decay_times_vax)}, Length of decay_times_rec={len(decay_times_rec)}")

        # Move waned individuals back to susceptible compartment
        S[t] += num_waned_vax + num_waned_rec
        V[t] -= num_waned_vax
        R[t] -= num_waned_rec

        # DIAGNOSIS AND LOG
        logging.info(f"Day {t}: Length of decay_times_vax={len(decay_times_vax)}, V[{t}]={V[t]}, Difference={V[t] - len(decay_times_vax)}")
        logging.info(f"Day {t}: S[t]={S[t]}, I[t]={I[t]}, R[t]={R[t]}, V[t]={V[t]}, Waned_vax={num_waned_vax}")

        if len(decay_times_vax) != V[t]:
            logging.warning(f"Day {t}: Length discrepancy: Length of decay_times_vax={len(decay_times_vax)}, V[t]={V[t]}")

        total_population = S[t] + I[t] + R[t] + V[t]
        if not np.isclose(total_population, N):
            logging.error(f"Population not conserved on day {t}: Total={total_population}, Expected={N}")

        if S[t] < 0 or I[t] < 0 or R[t] < 0 or V[t] < 0:
            logging.error(f"Negative compartment values on day {t}: S={S[t]}, I={I[t]}, R={R[t]}, V={V[t]}")

        # if is_vax_period and ((t - start_vax_day) % vax_period == vax_duration - 1) and diagnosis:
            # plot_histogram(decay_times_vax, decay_times_rec, scenario, round_counter, start=False)

        if t % 30 == 0 or is_vax_period:
            logging.info(f"Day {t}: S={S[t]:.2f}, I={I[t]:.2f}, R={R[t]:.2f}, V={V[t]:.2f}, New Vaccinations={new_vaccinations if is_vax_period else 0}")

    logging.info(f"Simulation of the {scenario.capitalize()} model completed.")

    return S, I, R, V


def sirsv_model_with_weibull_targetted_vaccination(params, scenario, random_seed=42, diagnosis=None):
    np.random.seed(random_seed)
    random.seed(random_seed)

    # Extract parameters
    beta = params['beta']
    gamma = params['gamma']
    vax_rate = params['vax_rate']
    weibull_shape_vax = params['weibull_shape_vax']
    weibull_scale_vax = params['weibull_scale_vax']
    weibull_shape_rec = params['weibull_shape_rec']
    weibull_scale_rec = params['weibull_scale_rec']
    days = params['days']
    seed_rate = params['seed_rate']
    vax_period = params['vax_period']
    vax_duration = params['vax_duration']
    start_vax_day = params['start_vax_day']

    # Initial conditions
    S0, I0, R0, V0 = params['S0'], params['I0'], params['R0'], params['V0']
    N = S0 + I0 + R0 + V0

    S, I, R, V = [np.zeros(days) for _ in range(4)]
    S[0], I[0], R[0], V[0] = S0, I0, R0, V0

    decay_times_vax = []
    decay_times_rec = []

    # Seed initial vaccinated and recovered individuals' waning times
    if V0 > 0:
        decay_times_vax.extend((weibull_scale_vax * np.random.weibull(weibull_shape_vax, int(V0))).astype(int).tolist())
    if R0 > 0:
        decay_times_rec.extend((weibull_scale_rec * np.random.weibull(weibull_shape_rec, int(R0))).astype(int).tolist())

    round_counter = 0

    logging.info(f"Starting simulation for scenario: {scenario}")

    for t in tqdm(range(1, days), desc=f"Running {scenario} simulation", unit="day"):
        new_seeds = min(seed_rate, S[t-1])

        # VACCINATION ROUND
        if t == start_vax_day or (t > start_vax_day and (t - start_vax_day) % vax_period == 0):
            round_counter += 1
            to_vaccinate = min(vax_rate * S[t-1], S[t-1])
            logging.info(f"Round {round_counter} start day: {t}")

            # Calculate the number of vaccinations to reset, considering the vaccination period
            num_vax_to_reset = int(min(vax_rate * vax_period * V[t-1], V[t-1]))
            if num_vax_to_reset > 0 and len(decay_times_vax) > 0:
                num_vax_to_reset = min(num_vax_to_reset, len(decay_times_vax))

                # Method 1: Randomly select decay times to reset
                # indices_to_reset = random.sample(range(len(decay_times_vax)), num_vax_to_reset)

                # Method 2: Sort decay times and select indices with the lowest decay times
                sorted_indices = np.argsort(decay_times_vax)
                indices_to_reset = sorted_indices[:num_vax_to_reset]

                for i in indices_to_reset:
                    decay_times_vax[i] = weibull_scale_vax * np.random.weibull(weibull_shape_vax)
                logging.info(f"Day {t}: Re-vaccination reset: {num_vax_to_reset} decay times reset")

            # if diagnosis:
                # plot_histogram(decay_times_vax, decay_times_rec, scenario, round_counter, start=True)

        # Check if it's within a vaccination period
        is_vax_period = (t >= start_vax_day) and ((t - start_vax_day) % vax_period < vax_duration)
        if is_vax_period:
            new_vaccinations = int(to_vaccinate)
            logging.info(f"Day {t}: Daily vaccinations: {new_vaccinations}")

            # Update compartments for new vaccinations
            if new_vaccinations > 0:
                new_susceptible_decay_times = (weibull_scale_vax * np.random.weibull(weibull_shape_vax, new_vaccinations)).astype(int).tolist()
                initial_len = len(decay_times_vax)
                decay_times_vax.extend(new_susceptible_decay_times)
                updated_len = len(decay_times_vax)
                logging.info(f"Day {t}: New vaccinations: {new_vaccinations}, Length of decay_times_vax: {initial_len}, Length of decay_times_vax: {updated_len}")
            else:
                S[t] = S[t-1]
                R[t] = R[t-1]
                V[t] = V[t-1]

        else:
            new_vaccinations = 0

        # Calculate transitions
        new_infections = beta * S[t-1] * I[t-1] / N + new_seeds
        new_recoveries = gamma * I[t-1]

        # Update compartments
        S[t] = S[t-1] - new_infections - new_vaccinations
        I[t] = I[t-1] + new_infections - new_recoveries
        R[t] = R[t-1] + new_recoveries
        V[t] = V[t-1] + new_vaccinations

        if new_recoveries > 0:
            new_recovered_decay_times = (weibull_scale_rec * np.random.weibull(weibull_shape_rec, int(new_recoveries))).astype(int).tolist()
            decay_times_rec.extend(new_recovered_decay_times)

        # IMMUNITY WANING
        logging.info(f"Day {t}: Before waning: Length of decay_times_vax={len(decay_times_vax)}, Length of decay_times_rec={len(decay_times_rec)}")

        decay_times_vax_before = len(decay_times_vax)
        decay_times_vax = [x - 1 for x in decay_times_vax if x > 0]
        num_waned_vax = decay_times_vax_before - len(decay_times_vax)

        decay_times_rec_before = len(decay_times_rec)
        decay_times_rec = [x - 1 for x in decay_times_rec if x > 0]
        num_waned_rec = decay_times_rec_before - len(decay_times_rec)

        logging.info(f"Day {t}: Waned vaccinated: {num_waned_vax}, Waned recovered: {num_waned_rec}")
        logging.info(f"Day {t}: After waning: Length of decay_times_vax={len(decay_times_vax)}, Length of decay_times_rec={len(decay_times_rec)}")

        # Move waned individuals back to susceptible compartment
        S[t] += num_waned_vax + num_waned_rec
        V[t] -= num_waned_vax
        R[t] -= num_waned_rec

        # DIAGNOSIS AND LOG
        logging.info(f"Day {t}: Length of decay_times_vax={len(decay_times_vax)}, V[{t}]={V[t]}, Difference={V[t] - len(decay_times_vax)}")
        logging.info(f"Day {t}: S[t]={S[t]}, I[t]={I[t]}, R[t]={R[t]}, V[t]={V[t]}, Waned_vax={num_waned_vax}")

        if len(decay_times_vax) != V[t]:
            logging.warning(f"Day {t}: Length discrepancy: Length of decay_times_vax={len(decay_times_vax)}, V[t]={V[t]}")

        total_population = S[t] + I[t] + R[t] + V[t]
        if not np.isclose(total_population, N):
            logging.error(f"Population not conserved on day {t}: Total={total_population}, Expected={N}")

        if S[t] < 0 or I[t] < 0 or R[t] < 0 or V[t] < 0:
            logging.error(f"Negative compartment values on day {t}: S={S[t]}, I={I[t]}, R={R[t]}, V={V[t]}")

        # if is_vax_period and ((t - start_vax_day) % vax_period == vax_duration - 1) and diagnosis:
            # plot_histogram(decay_times_vax, decay_times_rec, scenario, round_counter, start=False)

        if t % 30 == 0 or is_vax_period:
            logging.info(f"Day {t}: S={S[t]:.2f}, I={I[t]:.2f}, R={R[t]:.2f}, V={V[t]:.2f}, New Vaccinations={new_vaccinations if is_vax_period else 0}")

    logging.info(f"Simulation of the {scenario.capitalize()} model completed.")

    return S, I, R, V
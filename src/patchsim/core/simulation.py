import logging
from datetime import datetime
from patchsim.utils.io import read_config, read_patch_population, read_network, read_seeding_infection, apply_seeding_infections
import os
import csv
import importlib
import numpy as np
from pathlib import Path
import sys


current_dir = Path(__file__).resolve().parent
src_dir = current_dir.parents[2] / "src"
sys.path.insert(0, str(src_dir))

def run_simulation(config_path, model_type=None):
    """
    Run a model-agnostic metapopulation network simulation.

    Args:
        config_path (str): Path to the configuration file.
        model_type (str, optional): The type of model to run (e.g., "sample_sir_ode").
            If not provided, it will be read from the configuration file.
    """
    logging.info("Starting simulation...")

    # Read configuration
    config = read_config(config_path)

    # Extract parameters
    start_date = config["StartDate"]
    end_date = config["EndDate"]

    # Convert to datetime.datetime for consistency
    start_date = datetime.combine(start_date, datetime.min.time())
    end_date = datetime.combine(end_date, datetime.min.time())

    # Determine the model type
    if model_type is None:
        model_type = config.get("ModelType", "sample-sir-ode")

    # Read input files
    patches = read_patch_population(config["PatchFile"])
    network = read_network(config["NetworkFile"])
    seeds = read_seeding_infection(config["SeedFile"], start_date=start_date, end_date=end_date)

    # Apply seeding infections
    current_date = start_date
    patches = apply_seeding_infections(patches, seeds, current_date=current_date)

    # Update compartments in patches based on seeding infections
    for patch in patches:
        patch["I0"] = patch.get("I0", 0)  # Ensure I0 is set
        patch["S0"] = patch["N0"] - patch["I0"]  # Update S0 based on N0 and I0
        patch["R0"] = patch.get("R0", 0)  # Ensure R0 is set

    # Dynamically import the model module
    try:
        model_module = importlib.import_module(f"patchsim.models.{model_type}")
    except ModuleNotFoundError as e:
        logging.error(f"Model module '{model_type}' not found: {e}")
        raise

    # Initialize compartments for each patch using the model's initialization logic
    for patch in patches:
        patch["compartments"] = model_module.initialize_compartments(patch, seeds)

    # Run the model for each patch
    results = {}
    for patch in patches:
        initial_conditions = patch["compartments"]
        days = (end_date - start_date).days

        # Call the model runner for disease-specific logic
        t, solution = model_module.run_model(initial_conditions, days, config)
        results[patch["region"]] = solution
    
    # Ensure results[region] is structured as a 2D NumPy array
    for region in results:
        results[region] = np.array([list(compartments.values()) for compartments in results[region]])

    # Use the network to simulate interactions between patches
    # The contact matrix defines the connectivity between patches
    contact_matrix = {}
    for edge in network:
        source = edge["source"]
        target = edge["target"]
        weight = edge["weight"]

        if source not in contact_matrix:
            contact_matrix[source] = {}
        contact_matrix[source][target] = weight

    # Apply movement and force of infection
    for source, targets in contact_matrix.items():
        if source in results:
            S_source = results[source][:, 0]  # Susceptible individuals in the source patch
            for target, weight in targets.items():
                if target in results:
                    # Movement of susceptible individuals
                    moved_population = weight * S_source
                    results[source][:, 0] -= moved_population  # Decrease S in source
                    results[target][:, 0] += moved_population  # Increase S in target

    # Calculate force of infection for each patch
    for patch in patches:
        region = patch["region"]
        if region in results:
            S = results[region][:, 0]  # Susceptible individuals
            I = results[region][:, 1]  # Infected individuals
            N = S + I  # Total population
            beta = config["InfectionRate"]

            # Force of infection λj = βj * (Ij / Nj)
            force_of_infection = beta * (I / N)
            results[region][:, 1] += force_of_infection * S  # Update infected individuals

    # Save results to output/runs
    os.makedirs("output/runs", exist_ok=True)
    output_file = "output/runs/simulation_results.csv"
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Region", "Day", "Susceptible", "Infected", "Recovered"])
        for region, data in results.items():
            for day, (S, I, R) in enumerate(data):
                writer.writerow([region, day, S, I, R])

    logging.info(f"Simulation results saved to {output_file}")

    # Generate plots for the simulation results
    from patchsim.utils.viz import plot_simulation_results

    # Create a versioned/modular file name for the plots
    plot_output_dir = "output/plots"
    plot_version = datetime.now().strftime("%Y%m%d-%H%M%S")
    plot_file_prefix = f"simulation_results_{plot_version}"

    plot_simulation_results(results_file=output_file, output_dir=plot_output_dir, file_prefix=plot_file_prefix)

    logging.info(f"Plots saved to {plot_output_dir} with prefix {plot_file_prefix}")


if __name__ == "__main__":
    run_simulation("configs/global.yaml")
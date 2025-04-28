import matplotlib.pyplot as plt
import os

def plot_simulation_results(results_file: str, output_dir: str = "output/plots", file_prefix: str = "simulation_results"):
    """
    Generate a dummy visualization of the simulation results.

    Args:
        results_file (str): Path to the simulation results CSV file.
        output_dir (str): Directory to save the plots.
        file_prefix (str): Prefix for the plot file names.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Read the results file
    data = {}
    with open(results_file, "r") as f:
        next(f)  # Skip header
        for line in f:
            region, day, S, I, R = line.strip().split(",")
            day = int(day)
            S, I, R = float(S), float(I), float(R)
            if region not in data:
                data[region] = {"days": [], "S": [], "I": [], "R": []}
            data[region]["days"].append(day)
            data[region]["S"].append(S)
            data[region]["I"].append(I)
            data[region]["R"].append(R)

    # Plot results for each region
    for region, values in data.items():
        plt.figure()
        plt.plot(values["days"], values["S"], label="Susceptible")
        plt.plot(values["days"], values["I"], label="Infected")
        plt.plot(values["days"], values["R"], label="Recovered")
        plt.xlabel("Days")
        plt.ylabel("Population")
        plt.title(f"SIR Model Simulation for {region}")
        plt.legend()
        plt.grid()

        # Save the plot
        plot_file = os.path.join(output_dir, f"{file_prefix}_{region}_sir_plot.png")
        plt.savefig(plot_file)
        plt.close()

    print(f"Plots saved to {output_dir}")
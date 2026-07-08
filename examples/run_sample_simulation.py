"""Run the bundled sample SIR ODE simulation end to end.

This script is intended as a concrete, reproducible example for reviewers:
it loads the repository sample configuration, runs the simulation, and
prints the generated CSV and plot paths.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import patchsim


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    config_path = repo_root / "configs" / "sample-sir-ode.yaml"

    config = patchsim.load_config(str(config_path))
    net, y0, patches, num_patches = patchsim.setup_simulation(config)
    summary = patchsim.run_simulation(config, config["ModelName"], net, y0, patches, num_patches)

    print("PatchSim sample simulation completed successfully.")
    print(f"CSV output: {summary['csv_path']}")
    print(f"Plot output: {summary['plot_path']}")

    preview = pd.read_csv(summary["csv_path"]).head()
    print("\nFirst rows of the generated time series:")
    print(preview.to_string(index=False))


if __name__ == "__main__":
    main()
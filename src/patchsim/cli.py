import argparse
import logging
import os
from datetime import datetime

from patchsim.models.sample_sir_ode import run_simulation


def main():
    """Command-line interface for running the PatchSim simulation."""
    parser = argparse.ArgumentParser(
        description=(
            "PatchSim: A modular metapopulation simulation framework for multi-disease "
            "epidemiological modelling."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run simulation with sample SIR model
  patchsim --config configs/sample-sir-ode.yaml

  # Run simulation with custom model
  patchsim --config path/to/your/config.yaml
        """
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to the configuration file (YAML) containing model parameters, input files, and simulation settings",
    )
    args = parser.parse_args()

    # Create output directories
    os.makedirs("output/logs", exist_ok=True)

    # Set up logging with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"output/logs/cli_{timestamp}.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )

    try:
        logging.info("Starting PatchSim simulation...")
        run_simulation(args.config)
        logging.info("Simulation completed successfully.")
    except Exception as e:
        logging.error(f"Simulation failed: {e}")
        raise


if __name__ == "__main__":
    main()

import argparse
import logging
from patchsim.core.simulation import run_simulation

def main():
    """
    Command-line interface for running the PatchSim simulation.
    """
    parser = argparse.ArgumentParser(description="Run PatchSim simulation.")
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to the configuration file (e.g., configs/global.yaml)",
    )
    parser.add_argument(
        "--model",
        type=str,
        required=False,
        default=None,
        help="Name of the model to use (e.g., sample_sir_ode). Overrides the model specified in the config file.",
    )

    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler("output/logs/cli.log"),
            logging.StreamHandler(),
        ],
    )

    # Run the simulation
    try:
        logging.info("Starting PatchSim simulation...")
        run_simulation(config_path=args.config, model_type=args.model)
        logging.info("Simulation completed successfully.")
    except Exception as e:
        logging.error(f"Simulation failed: {e}")
        raise

if __name__ == "__main__":
    main()
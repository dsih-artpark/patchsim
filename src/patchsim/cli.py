import argparse
import logging

from patchsim.core.simulation import load_config, run_simulation, setup_simulation


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
        """
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to the configuration file (YAML) containing model parameters, input files, and simulation settings",
    )
    args = parser.parse_args()

    try:
        logging.info("Starting PatchSim simulation...")
        config = load_config(args.config)

        # Set up simulation
        net, y0, patches, num_patches = setup_simulation(config)
        run_simulation(config, config['ModelName'], net, y0, patches, num_patches)
        logging.info("Simulation completed successfully.")
    except Exception as e:
        logging.error(f"Simulation failed: {e}")
        raise


if __name__ == "__main__":
    main()

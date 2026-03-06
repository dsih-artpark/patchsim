import argparse
import logging
from pathlib import Path
import sys
import textwrap

import yaml

from patchsim.core.simulation import load_config, run_simulation, setup_simulation


def _configure_logging() -> None:
    """Configure CLI logging with rich handler when available."""
    try:
        from rich.logging import RichHandler

        logging.basicConfig(
            level=logging.INFO,
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler(rich_tracebacks=True)],
        )
    except Exception:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")


def _cmd_run(config_path: str) -> None:
    logging.info("Starting PatchSim simulation...")
    config = load_config(config_path)
    net, y0, patches, num_patches = setup_simulation(config)
    run_simulation(config, config['ModelName'], net, y0, patches, num_patches)
    logging.info("Simulation completed successfully.")


def _cmd_validate(config_path: str) -> None:
    config = load_config(config_path)
    setup_simulation(config)
    logging.info("Configuration is valid: %s", config_path)


def _cmd_init(output: str, force: bool = False) -> None:
    out_path = Path(output)
    if out_path.exists() and not force:
        raise FileExistsError(f"Refusing to overwrite existing file: {out_path}. Use --force to overwrite.")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    template = {
        "PatchFile": "data/patch/sample-sir-ode-patch-population.csv",
        "NetworkFile": "data/networks/sample-network-static.csv",
        "SeedFile": "data/seeds/sample-sir-ode-patchA-2.csv",
        "ModelName": "sample-sir-ode",
        "TMax": 50,
        "OutputDir": "output/sample-sir-ode",
        "compartments": ["S", "I", "R"],
        "Parameters": {"beta": 0.08, "gamma": 0.1},
        "Transitions": {
            "S -> I": "beta",
            "I -> R": "gamma * I",
        },
    }
    with out_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(template, f, sort_keys=False)

    logging.info("Created starter config: %s", out_path)


def main():
    """Command-line interface for running the PatchSim simulation."""
    _configure_logging()

    parser = argparse.ArgumentParser(
        description=(
            "PatchSim: A modular metapopulation simulation framework for multi-disease "
            "epidemiological modelling."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
Examples:
  patchsim run --config configs/sample-sir-ode.yaml
  patchsim validate --config configs/sample-sir-ode.yaml
  patchsim init --output configs/sample-sir-ode.yaml

Legacy:
  patchsim --config configs/sample-sir-ode.yaml
        """)
    )
    # Backward-compatible path: `patchsim --config ...` == `patchsim run --config ...`
    parser.add_argument(
        "--config",
        type=str,
        required=False,
        help="Path to the configuration file (legacy mode; equivalent to run --config)",
    )
    subparsers = parser.add_subparsers(dest="command")

    run_p = subparsers.add_parser("run", help="Run a simulation")
    run_p.add_argument("--config", required=True, help="Path to simulation config YAML")

    validate_p = subparsers.add_parser("validate", help="Validate config and inputs")
    validate_p.add_argument("--config", required=True, help="Path to simulation config YAML")

    init_p = subparsers.add_parser("init", help="Create a starter config file")
    init_p.add_argument("--output", default="configs/sample-sir-ode.yaml", help="Output path for starter config")
    init_p.add_argument("--force", action="store_true", help="Overwrite output file if it exists")

    args = parser.parse_args()

    try:
        if args.command == "run":
            _cmd_run(args.config)
        elif args.command == "validate":
            _cmd_validate(args.config)
        elif args.command == "init":
            _cmd_init(args.output, force=args.force)
        elif args.config:
            # Legacy mode
            _cmd_run(args.config)
        else:
            parser.print_help()
            raise SystemExit(2)
    except Exception as e:
        logging.error(f"Simulation failed: {e}")
        raise


if __name__ == "__main__":
    main()

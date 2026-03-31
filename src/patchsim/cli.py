import argparse
import logging
import shutil
import textwrap
from importlib import resources
from pathlib import Path

from patchsim import __version__
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
    run_simulation(config, config["ModelName"], net, y0, patches, num_patches)
    logging.info("Simulation completed successfully.")


def _cmd_validate(config_path: str) -> None:
    config = load_config(config_path)
    setup_simulation(config)
    logging.info("Configuration is valid: %s", config_path)


def _copy_template_tree(template_node, target_path: Path) -> None:
    if template_node.is_dir():
        target_path.mkdir(parents=True, exist_ok=True)
        for child in template_node.iterdir():
            _copy_template_tree(child, target_path / child.name)
        return

    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.write_bytes(template_node.read_bytes())


def _cmd_init(name: str, force: bool = False) -> None:
    project_dir = Path(name)

    # Safety checks: prevent deleting cwd or non-directories
    resolved = project_dir.resolve()
    if project_dir.exists() and not project_dir.is_dir():
        raise NotADirectoryError(f"Target exists and is not a directory: {project_dir}")
    if force and resolved == Path.cwd().resolve():
        raise ValueError(f"Refusing to overwrite current working directory: {resolved}")

    if project_dir.exists() and any(project_dir.iterdir()) and not force:
        raise FileExistsError(f"Refusing to overwrite existing directory: {project_dir}. Use --force to overwrite.")

    if project_dir.exists() and force:
        shutil.rmtree(project_dir)

    template_root = resources.files("patchsim").joinpath("templates", "project")
    _copy_template_tree(template_root, project_dir)

    config_path = project_dir / "config.yaml"
    rendered = config_path.read_text(encoding="utf-8").replace("{{PROJECT_NAME}}", project_dir.name)
    config_path.write_text(rendered, encoding="utf-8")

    logging.info("Created project scaffold at: %s", project_dir)


def _list_builtin_models() -> list[str]:
    models_dir = Path(__file__).resolve().parent / "models"
    models = []
    for f in sorted(models_dir.glob("*.py")):
        if f.name == "__init__.py":
            continue
        models.append(f.stem)
    return models


def _cmd_list_models() -> None:
    models = _list_builtin_models()
    if not models:
        logging.info("No built-in models found.")
        return
    logging.info("Built-in models:")
    for model in models:
        logging.info("- %s", model)


def main():
    """Command-line interface for running the PatchSim simulation."""
    _configure_logging()

    parser = argparse.ArgumentParser(
        description=(
            "PatchSim: A modular metapopulation simulation framework for multi-disease epidemiological modelling."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
Examples:
    patchsim init my-project
    patchsim run -c my-project/config.yaml
    patchsim validate -c my-project/config.yaml
    patchsim list-models
        """),
    )
    parser.add_argument("--version", action="version", version=f"patchsim {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    init_p = subparsers.add_parser("init", help="Scaffold a new PatchSim project")
    init_p.add_argument("name", help="Directory name for the new project")
    init_p.add_argument("--force", action="store_true", help="Overwrite target directory if it already exists")

    run_p = subparsers.add_parser("run", help="Run a simulation")
    run_p.add_argument("-c", "--config", required=True, help="Path to simulation config YAML")

    validate_p = subparsers.add_parser("validate", help="Validate config and inputs")
    validate_p.add_argument("-c", "--config", required=True, help="Path to simulation config YAML")

    subparsers.add_parser("list-models", help="List available built-in models")

    args = parser.parse_args()

    try:
        if args.command == "init":
            _cmd_init(args.name, force=args.force)
        elif args.command == "run":
            _cmd_run(args.config)
        elif args.command == "validate":
            _cmd_validate(args.config)
        elif args.command == "list-models":
            _cmd_list_models()
        else:
            parser.print_help()
            raise SystemExit(2)
    except Exception as e:
        logging.error(f"Simulation failed: {e}")
        raise


if __name__ == "__main__":
    main()

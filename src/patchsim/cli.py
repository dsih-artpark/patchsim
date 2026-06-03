import argparse
import json
import logging
import shutil
import textwrap
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as _np
import pandas as _pd
import yaml

from patchsim import __version__
from patchsim.core.simulation import (
    get_available_template_names,
    get_config_schema,
    get_init_template_config,
    get_model_catalog,
    load_config,
    run_simulation,
    setup_simulation,
)
from patchsim.utils import (
    distance_matrix_km,
    gravity_contact_matrix,
    load_geojson_centroids,
)


def _configure_logging(*, json_output: bool = False) -> None:
    """Configure CLI logging with rich handler when available."""
    if json_output:
        logging.basicConfig(level=logging.WARNING)
        return

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


def _emit_json(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, sort_keys=True))


def _cmd_run(config_path: str, *, json_output: bool = False) -> dict[str, Any]:
    if not json_output:
        print("Starting PatchSim simulation...")
    config = load_config(config_path)
    net, y0, patches, num_patches = setup_simulation(config)
    summary = run_simulation(config, config["ModelName"], net, y0, patches, num_patches)
    if not json_output:
        print("Simulation completed successfully.")
    return {"ok": True, "config": config_path, **summary} if json_output else summary


def _cmd_validate(config_path: str, *, json_output: bool = False, schema: bool = False) -> dict[str, Any] | None:
    if schema:
        return get_config_schema()

    config = load_config(config_path)
    _net, _y0, patches, num_patches = setup_simulation(config)
    if not json_output:
        print(f"Configuration is valid: {config_path}")
    if json_output:
        return {
            "ok": True,
            "config": config_path,
            "model_name": config.get("ModelName"),
            "num_patches": num_patches,
            "patches": patches,
        }
    return None


def _copy_template_tree(template_node, target_path: Path) -> None:
    if template_node.is_dir():
        target_path.mkdir(parents=True, exist_ok=True)
        for child in template_node.iterdir():
            _copy_template_tree(child, target_path / child.name)
        return

    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.write_bytes(template_node.read_bytes())


def _write_seed_for_template(project_dir: Path, config: dict[str, Any]) -> None:
    """Write a seed CSV whose columns match the template's compartments.

    The scaffold ships a single patch-population file; seed every patch fully
    susceptible and place one infectious individual in the first patch.
    """
    import pandas as pd

    compartments = list(config.get("compartments") or ["S", "I", "R"])
    patch_df = pd.read_csv(project_dir / config["PatchFile"])
    patch_col = next(c for c in patch_df.columns if c.lower() == "patch")
    pop_col = next(c for c in patch_df.columns if c.lower() == "population")

    susceptible = "S" if "S" in compartments else compartments[0]
    infectious = "I" if "I" in compartments else compartments[-1]

    rows = []
    for idx, record in patch_df.iterrows():
        seeded = 1 if idx == 0 else 0
        row = {"patch": record[patch_col], **{c: 0 for c in compartments}}
        row[infectious] = seeded
        row[susceptible] = int(record[pop_col]) - seeded
        rows.append(row)

    seed_path = project_dir / config["SeedFile"]
    seed_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows)[["patch", *compartments]].to_csv(seed_path, index=False)


def _cmd_init(name: str, force: bool = False, template: str = "sir") -> None:
    project_dir = Path(name)

    # Safety checks: prevent deleting cwd, parent dirs, root, or non-directories
    resolved = project_dir.resolve()
    cwd = Path.cwd().resolve()

    if project_dir.exists() and not project_dir.is_dir():
        raise NotADirectoryError(f"Target exists and is not a directory: {project_dir}")

    # Block deletion of root, cwd, or any ancestor of cwd
    if force and (
        resolved == Path(resolved.anchor)  # Filesystem root (/ or C:\)
        or resolved == cwd  # Current working directory
        or cwd in resolved.parents  # resolved is ancestor of cwd (e.g., ..)
    ):
        raise ValueError(f"Refusing to overwrite unsafe target: {resolved}")

    if project_dir.exists() and any(project_dir.iterdir()) and not force:
        raise FileExistsError(f"Refusing to overwrite existing directory: {project_dir}. Use --force to overwrite.")

    if project_dir.exists() and force:
        shutil.rmtree(project_dir)

    template_root = resources.files("patchsim").joinpath("templates", "project")
    _copy_template_tree(template_root, project_dir)

    template_config = get_init_template_config(template, project_dir.name)
    config_path = project_dir / "config.yaml"
    config_path.write_text(yaml.safe_dump(template_config, sort_keys=False), encoding="utf-8")
    _write_seed_for_template(project_dir, template_config)

    print(f"Created project scaffold at: {project_dir}")


def _list_builtin_models() -> list[dict[str, str]]:
    return get_model_catalog()


def _cmd_list_models(*, json_output: bool = False) -> list[dict[str, str]]:
    models = _list_builtin_models()
    if not models:
        if not json_output:
            print("No built-in models found.")
        return []

    if json_output:
        return models

    print("Built-in models and templates:")
    for model in models:
        print(f"- {model['name']} ({model['kind']})")
    return models


def main() -> None:
    """Command-line interface for running the PatchSim simulation."""
    parser = argparse.ArgumentParser(
        description=(
            "PatchSim: A modular metapopulation simulation framework for multi-disease epidemiological modelling."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """
Examples:
    uv run patchsim init my-project
    uv run patchsim init my-project --template seir
    uv run patchsim run -c my-project/config.yaml
    uv run patchsim validate -c my-project/config.yaml
    uv run patchsim list-models
        """
        ),
    )
    parser.add_argument("--version", action="version", version=f"patchsim {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    init_p = subparsers.add_parser("init", help="Scaffold a new PatchSim project")
    init_p.add_argument("name", help="Directory name for the new project")
    init_p.add_argument("--force", action="store_true", help="Overwrite target directory if it already exists")
    init_p.add_argument(
        "--template",
        choices=get_available_template_names(),
        default="sir",
        help="Starter template to use for config.yaml",
    )

    run_p = subparsers.add_parser("run", help="Run a simulation")
    run_p.add_argument("-c", "--config", required=True, help="Path to simulation config YAML")
    run_p.add_argument("--json", action="store_true", help="Emit machine-readable JSON summary")

    validate_p = subparsers.add_parser("validate", help="Validate config and inputs")
    validate_p.add_argument("-c", "--config", help="Path to simulation config YAML")
    validate_p.add_argument("--schema", action="store_true", help="Print the configuration JSON Schema")
    validate_p.add_argument("--json", action="store_true", help="Emit machine-readable JSON summary")

    list_p = subparsers.add_parser("list-models", help="List available built-in models")
    list_p.add_argument("--json", action="store_true", help="Emit machine-readable JSON list")

    gen_p = subparsers.add_parser("generate-contacts", help="Generate contact matrix from GeoJSON or centroids CSV")
    group = gen_p.add_mutually_exclusive_group(required=True)
    group.add_argument("--geojson", help="Path to GeoJSON file containing regions")
    group.add_argument("--centroids", help="Path to CSV with precomputed centroids (id,lat,lon,population)")
    gen_p.add_argument("--out", default="contacts.csv", help="Output CSV path")
    gen_p.add_argument("--id-prop", default="id", help="Property name to use as feature id in GeoJSON")
    gen_p.add_argument("--pop-prop", default=None, help="Property name for population in GeoJSON (optional)")
    gen_p.add_argument("--model", choices=["gravity", "distance"], default="gravity", help="Generation method")
    gen_p.add_argument("--decay", type=float, default=2.0, help="Gravity model decay exponent")
    gen_p.add_argument("--scale", type=float, default=1.0, help="Gravity model scale factor")
    gen_p.add_argument(
        "--min-distance",
        type=float,
        default=1e-3,
        help="Minimum distance (km) floor to avoid div/0",
    )
    gen_p.add_argument(
        "--normalize",
        choices=["none", "row", "col", "symmetric"],
        default="none",
        help="Normalize contact matrix",
    )
    gen_p.add_argument("--format", choices=["edge", "matrix"], default="edge", help="Output format for contacts CSV")
    gen_p.add_argument("--force", action="store_true", help="Overwrite output file if it exists")

    args = parser.parse_args()
    json_mode = bool(getattr(args, "json", False) or getattr(args, "schema", False))
    _configure_logging(json_output=json_mode)

    try:
        if args.command == "init":
            _cmd_init(args.name, force=args.force, template=args.template)
        elif args.command == "run":
            result = _cmd_run(args.config, json_output=args.json)
            if args.json:
                _emit_json(result)
        elif args.command == "validate":
            if not args.schema and not args.config:
                parser.error("the following arguments are required: -c/--config")
            result = _cmd_validate(args.config, json_output=args.json, schema=args.schema)
            if result is not None:
                _emit_json(result)
        elif args.command == "list-models":
            result = _cmd_list_models(json_output=args.json)
            if args.json:
                _emit_json({"models": result})
        elif args.command == "generate-contacts":
            # Load centroids
            if getattr(args, "geojson", None):
                df = load_geojson_centroids(
                    args.geojson,
                    id_prop=args.id_prop,
                    pop_prop_candidates=[args.pop_prop] if args.pop_prop else None,
                )
            else:
                df = _pd.read_csv(args.centroids)

            if df.empty:
                raise ValueError("No features/centroids loaded from input")

            if args.model == "gravity":
                if "population" not in df.columns:
                    raise ValueError("Population column required for gravity model")
                mat = gravity_contact_matrix(
                    df,
                    pop_col="population",
                    decay=args.decay,
                    scale=args.scale,
                    min_distance_km=args.min_distance,
                )
            else:
                # distance-based: weights = 1/d^decay
                D = distance_matrix_km(df)
                D_safe = _np.maximum(D, args.min_distance)
                mat = (1.0 / (D_safe ** args.decay))
                _np.fill_diagonal(mat, 0.0)

            # Normalization
            if args.normalize == "row":
                row_sums = mat.sum(axis=1, keepdims=True)
                row_sums[row_sums == 0] = 1.0
                mat = mat / row_sums
            elif args.normalize == "col":
                col_sums = mat.sum(axis=0, keepdims=True)
                col_sums[col_sums == 0] = 1.0
                mat = mat / col_sums
            elif args.normalize == "symmetric":
                s = mat.sum()
                if s > 0:
                    mat = mat / s

            out_path = Path(args.out)
            if out_path.exists() and not args.force:
                raise FileExistsError(f"Output exists: {out_path} (use --force to overwrite)")

            # Determine ID column to use for labeling patches in output
            cols_map = {c.lower(): c for c in df.columns}
            id_candidates = [args.id_prop or "", "patch", "id", "name"]
            id_col = None
            for cand in id_candidates:
                if cand and cand.lower() in cols_map:
                    id_col = cols_map[cand.lower()]
                    break
            if id_col is None:
                raise ValueError(
                    "Could not determine identifier column in centroids; provide --id-prop "
                    "or include a 'patch'/'id' column",
                )

            # Normalize population column name for gravity model if present
            pop_candidates = ["population", "pop", "pop_total", "population_total"]
            pop_col = None
            for pc in pop_candidates:
                if pc in cols_map:
                    pop_col = cols_map[pc]
                    break
            if pop_col and pop_col != "population":
                # rename to expected lowercase 'population'
                df = df.rename(columns={pop_col: "population"})

            ids = list(df[id_col].astype(str))

            # Write as network CSV matching existing PatchSim conventions: day,source,target,weight
            if args.format == "edge":
                rows = []
                n = mat.shape[0]
                for i in range(n):
                    for j in range(n):
                        w = float(mat[i, j])
                        if w == 0.0:
                            continue
                        rows.append({"day": 0, "source": ids[i], "target": ids[j], "weight": w})
                out_df = _pd.DataFrame(rows)
                out_df = out_df[["day", "source", "target", "weight"]]
                out_df.to_csv(out_path, index=False)
            else:
                mdf = _pd.DataFrame(mat, index=ids, columns=ids)
                # matrix CSV is not used directly by PatchSim but keep for completeness
                mdf.to_csv(out_path)

            print(f"Wrote contacts to: {out_path}")
        else:
            parser.print_help()
            raise SystemExit(2)
    except Exception as e:
        logging.error(f"Simulation failed: {e}")
        raise


if __name__ == "__main__":
    main()

from __future__ import annotations

import hashlib
import json
import platform
import re
import tempfile
import time
from copy import deepcopy
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from numbers import Real
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import least_squares

from patchsim.core.model import NetworkModel
from patchsim.core.simulation import _simulate_prepared, get_run_settings, load_config, setup_simulation

_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_ARTIFACT_NAMES = ("estimates.csv", "fitted-seeds.csv", "attempts.csv", "residuals.csv")


@dataclass(frozen=True)
class Observable:
    name: str
    columns: tuple[str, ...]
    scale: float


@dataclass(frozen=True)
class FitVariable:
    kind: str
    name: str
    lower: float
    upper: float
    baseline: float
    patch: str | None = None
    group: str | None = None
    compartment: str | None = None
    state_key: str | None = None
    remainder_key: str | None = None
    patch_index: int | None = None
    group_index: int | None = None

    @property
    def identity(self) -> tuple[str, str | None, str | None, str]:
        return (self.kind, self.patch, self.group, self.name)


@dataclass(frozen=True)
class CalibrationPlan:
    name: str
    observations_path: Path
    observations: pd.DataFrame
    observables: tuple[Observable, ...]
    variables: tuple[FitVariable, ...]
    starts: tuple[tuple[float, ...], ...]
    max_evaluations: int
    warnings: tuple[str, ...]

    @property
    def n(self) -> int:
        return len(self.observations)

    @property
    def p(self) -> int:
        return len(self.variables)

    @property
    def start_count(self) -> int:
        return len(self.starts)

    @property
    def max_forward_simulations(self) -> int:
        return self.start_count * self.max_evaluations


class _BudgetExceeded(RuntimeError):
    pass


def _finite_float(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{label} must be a finite real number")
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{label} must be a finite real number")
    return result


def _bounds(value: Any, label: str, *, non_negative: bool = False) -> tuple[float, float]:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{label} bounds must be a two-element list")
    lower = _finite_float(value[0], f"{label} lower bound")
    upper = _finite_float(value[1], f"{label} upper bound")
    if lower >= upper:
        raise ValueError(f"{label} bounds must satisfy lower < upper")
    if non_negative and lower < 0:
        raise ValueError(f"{label} lower bound must be non-negative")
    return lower, upper


def _safe_name(value: Any) -> str:
    if not isinstance(value, str) or not _NAME_PATTERN.fullmatch(value) or value in {".", ".."}:
        raise ValueError("Calibration Name must be a safe path component")
    return value


def _patch_parameter_names(config: dict[str, Any]) -> set[str]:
    names: set[str] = set()
    for entry in config.get("PatchParameters", []):
        if isinstance(entry, dict) and isinstance(entry.get("parameters"), dict):
            names.update(entry["parameters"])
    return names


def _parse_observables(block: Any, output_columns: set[str]) -> tuple[Observable, ...]:
    if not isinstance(block, dict) or not block:
        raise ValueError("Calibration Observables must be a non-empty mapping")
    observables = []
    for name, definition in block.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Calibration observable names must be non-empty strings")
        if not isinstance(definition, dict):
            raise ValueError(f"Calibration observable {name!r} must be a mapping")
        columns = definition.get("Columns")
        if not isinstance(columns, list) or not columns or not all(isinstance(column, str) for column in columns):
            raise ValueError(f"Calibration observable {name!r} Columns must be a non-empty string list")
        if len(columns) != len(set(columns)):
            raise ValueError(f"Calibration observable {name!r} contains duplicate output columns")
        unknown = sorted(set(columns) - output_columns)
        if unknown:
            raise ValueError(f"Calibration observable {name!r} uses unknown output columns: {unknown}")
        scale = _finite_float(definition.get("Scale"), f"Calibration observable {name!r} Scale")
        if scale <= 0:
            raise ValueError(f"Calibration observable {name!r} Scale must be positive")
        observables.append(Observable(name=name, columns=tuple(columns), scale=scale))
    return tuple(observables)


def _load_observations(
    path: Path,
    observables: tuple[Observable, ...],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, list[str]]:
    frame = pd.read_csv(path, keep_default_na=False)
    required = ["time", "observable", "value"]
    if frame.columns.tolist() != required:
        raise ValueError(f"Calibration observations must have exactly these columns in order: {required}")
    if frame.empty:
        raise ValueError("Calibration observations must contain at least one row")

    for column in ("time", "value"):
        try:
            frame[column] = pd.to_numeric(frame[column], errors="raise")
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Calibration observations column {column!r} must contain finite numbers") from exc
        values = frame[column].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            rows = (np.flatnonzero(~np.isfinite(values)) + 2).tolist()
            raise ValueError(f"Calibration observations column {column!r} contains non-finite values at rows {rows}")

    frame["observable"] = frame["observable"].astype(str)
    empty_names = frame.index[frame["observable"].str.strip().eq("")].tolist()
    if empty_names:
        rows = [row + 2 for row in empty_names]
        raise ValueError(f"Calibration observations contain empty observable names at rows {rows}")
    known = {observable.name for observable in observables}
    unknown = sorted(set(frame["observable"]) - known)
    if unknown:
        raise ValueError(f"Calibration observations use unknown observable names: {unknown}")
    unused = sorted(known - set(frame["observable"]))
    if unused:
        raise ValueError(f"Calibration observables have no observation rows: {unused}")
    duplicate = frame.duplicated(["time", "observable"], keep=False)
    if duplicate.any():
        rows = (frame.index[duplicate] + 2).tolist()
        raise ValueError(f"Calibration observations contain duplicate (time, observable) rows: {rows}")

    _solver, t_max, time_step = get_run_settings(config)
    grid = np.arange(t_max, dtype=np.float64) * time_step
    if not np.isfinite(grid).all() or (len(grid) > 1 and not np.all(np.diff(grid) > 0)):
        raise ValueError("Calibration reporting grid must be finite, strictly increasing, and unique in float64")
    min_spacing = float(np.min(np.diff(grid))) if len(grid) > 1 else time_step
    epsilon = np.finfo(np.float64).eps
    indices: list[int] = []
    unmatched: list[tuple[int, float]] = []
    for row_index, observation_time in enumerate(frame["time"].to_numpy(dtype=float), start=2):
        tolerance = min(
            min_spacing / 4.0,
            max(8.0 * epsilon * max(1.0, abs(observation_time)), 1e-9 * time_step),
        )
        matches = np.flatnonzero(np.abs(grid - observation_time) <= tolerance)
        if len(matches) == 1:
            indices.append(int(matches[0]))
        else:
            indices.append(-1)
            unmatched.append((row_index, float(observation_time)))
    frame["grid_index"] = indices
    warnings = []
    if unmatched:
        examples = ", ".join(f"row {row}: {value:g}" for row, value in unmatched[:10])
        warnings.append(f"{len(unmatched)} unmatched observation times ({examples})")
    return frame, warnings


def _global_variables(config: dict[str, Any], block: Any) -> list[FitVariable]:
    if block is None:
        return []
    if not isinstance(block, dict):
        raise ValueError("Calibration Parameters must be a mapping")
    global_parameters = config.get("Parameters", {})
    if not isinstance(global_parameters, dict):
        raise ValueError("Top-level Parameters must be a mapping")
    patch_names = _patch_parameter_names(config)
    variables = []
    for name, bounds_value in block.items():
        if name not in global_parameters:
            raise ValueError(f"Unknown global calibration parameter: {name!r}")
        if name in patch_names:
            raise ValueError(f"Cannot fit global parameter {name!r}; it is also set in PatchParameters")
        lower, upper = _bounds(bounds_value, f"Calibration parameter {name!r}")
        baseline = _finite_float(global_parameters[name], f"Configured parameter {name!r}")
        if not lower <= baseline <= upper:
            raise ValueError(f"Configured parameter {name!r} must lie inside its calibration bounds")
        variables.append(FitVariable("parameter", name, lower, upper, baseline))
    return variables


def _initial_variables(
    net: NetworkModel,
    y0: dict[str, float],
    block: Any,
) -> list[FitVariable]:
    if block is None:
        return []
    if not isinstance(block, list):
        raise ValueError("Calibration InitialConditions must be a list")
    patch_indices = {name: index for index, name in enumerate(net.patch_names)}
    group_indices = {name: index for index, name in enumerate(net.groups)}
    compartments = list(net.base_model.compartments)
    seen_cells: set[tuple[int, int]] = set()
    variables = []

    for entry in block:
        if not isinstance(entry, dict):
            raise ValueError("Each Calibration InitialConditions entry must be a mapping")
        patch = entry.get("Patch")
        if patch not in patch_indices:
            raise ValueError(f"Calibration InitialConditions uses unknown patch: {patch!r}")
        if net.groups:
            group = entry.get("Group")
            if group not in group_indices:
                raise ValueError(f"Calibration InitialConditions uses unknown group: {group!r}")
            group_index = group_indices[group]
        else:
            if "Group" in entry:
                raise ValueError("Calibration InitialConditions Group is forbidden for an ungrouped model")
            group = None
            group_index = 0
        patch_index = patch_indices[patch]
        cell = (patch_index, group_index)
        if cell in seen_cells:
            raise ValueError(f"Calibration InitialConditions repeats patch/group cell: {patch!r}, {group!r}")
        seen_cells.add(cell)

        remainder = entry.get("Remainder")
        if remainder not in compartments:
            raise ValueError(f"Calibration InitialConditions uses unknown remainder compartment: {remainder!r}")
        fit = entry.get("Fit")
        if not isinstance(fit, dict) or not fit:
            raise ValueError("Calibration InitialConditions Fit must be a non-empty mapping")
        if remainder in fit:
            raise ValueError("Calibration initial-condition remainder cannot also be fitted")
        unknown = sorted(set(fit) - set(compartments))
        if unknown:
            raise ValueError(f"Calibration InitialConditions fits unknown compartments: {unknown}")

        remainder_key = net.state_key(remainder, patch_index, group_index)
        fitted_upper = 0.0
        fitted = []
        for compartment, bounds_value in fit.items():
            lower, upper = _bounds(
                bounds_value,
                f"Calibration initial condition {patch!r}/{group!r}/{compartment!r}",
                non_negative=True,
            )
            state_key = net.state_key(compartment, patch_index, group_index)
            baseline = _finite_float(y0[state_key], f"Seed value {patch!r}/{group!r}/{compartment!r}")
            if not lower <= baseline <= upper:
                raise ValueError(
                    f"Seed value {patch!r}/{group!r}/{compartment!r} must lie inside its calibration bounds"
                )
            fitted_upper += upper
            fitted.append(
                FitVariable(
                    "initial",
                    compartment,
                    lower,
                    upper,
                    baseline,
                    patch=str(patch),
                    group=str(group) if group is not None else None,
                    compartment=compartment,
                    state_key=state_key,
                    remainder_key=remainder_key,
                    patch_index=patch_index,
                    group_index=group_index,
                )
            )
        fixed_sum = sum(
            y0[net.state_key(compartment, patch_index, group_index)]
            for compartment in compartments
            if compartment not in fit and compartment != remainder
        )
        population = sum(y0[net.state_key(compartment, patch_index, group_index)] for compartment in compartments)
        if fitted_upper + fixed_sum > population:
            raise ValueError(
                f"Calibration initial-condition bounds for patch {patch!r}, group {group!r} "
                "cannot guarantee a non-negative remainder"
            )
        variables.extend(fitted)
    return variables


def _parse_starts(calibration: dict[str, Any], variables: tuple[FitVariable, ...]) -> tuple[tuple[float, ...], ...]:
    baseline = tuple(variable.baseline for variable in variables)
    starts = [baseline]
    additional = calibration.get("Starts", [])
    if not isinstance(additional, list):
        raise ValueError("Calibration Starts must be a list")

    parameter_variables = [variable for variable in variables if variable.kind == "parameter"]
    initial_variables = [variable for variable in variables if variable.kind == "initial"]
    initial_by_cell: dict[tuple[str, str | None], list[FitVariable]] = {}
    for variable in initial_variables:
        initial_by_cell.setdefault((variable.patch or "", variable.group), []).append(variable)

    for start_index, entry in enumerate(additional, start=1):
        if not isinstance(entry, dict):
            raise ValueError(f"Calibration start {start_index} must be a mapping")
        parameter_values = entry.get("Parameters", {})
        if not isinstance(parameter_values, dict) or set(parameter_values) != {v.name for v in parameter_variables}:
            raise ValueError(f"Calibration start {start_index} must provide every fitted global parameter")

        values: dict[tuple[str, str | None, str | None, str], float] = {}
        for variable in parameter_variables:
            values[variable.identity] = _finite_float(
                parameter_values[variable.name], f"Calibration start {start_index} parameter {variable.name!r}"
            )

        initial_entries = entry.get("InitialConditions", [])
        if not isinstance(initial_entries, list):
            raise ValueError(f"Calibration start {start_index} InitialConditions must be a list")
        seen_cells: set[tuple[str, str | None]] = set()
        for initial_entry in initial_entries:
            if not isinstance(initial_entry, dict):
                raise ValueError(f"Calibration start {start_index} InitialConditions entries must be mappings")
            cell = (str(initial_entry.get("Patch")), initial_entry.get("Group"))
            if cell not in initial_by_cell or cell in seen_cells:
                raise ValueError(f"Calibration start {start_index} has an unknown or duplicate initial-condition cell")
            seen_cells.add(cell)
            cell_values = initial_entry.get("Values")
            cell_variables = initial_by_cell[cell]
            if not isinstance(cell_values, dict) or set(cell_values) != {v.name for v in cell_variables}:
                raise ValueError(f"Calibration start {start_index} must provide every fitted initial condition")
            for variable in cell_variables:
                values[variable.identity] = _finite_float(
                    cell_values[variable.name],
                    f"Calibration start {start_index} initial condition {variable.name!r}",
                )
        if seen_cells != set(initial_by_cell):
            raise ValueError(f"Calibration start {start_index} must provide every fitted initial condition")

        vector = tuple(values[variable.identity] for variable in variables)
        for variable, value in zip(variables, vector, strict=True):
            if not variable.lower <= value <= variable.upper:
                raise ValueError(f"Calibration start {start_index} value for {variable.name!r} is outside bounds")
        if vector in starts:
            raise ValueError(f"Calibration start {start_index} duplicates an existing start")
        starts.append(vector)
    return tuple(starts)


def get_calibration_plan(
    config: dict[str, Any],
    net: NetworkModel,
    y0: dict[str, float],
    *,
    required: bool = False,
) -> CalibrationPlan | None:
    calibration = config.get("Calibration")
    if calibration is None:
        if required:
            raise ValueError("Calibration configuration is required")
        return None
    if not isinstance(calibration, dict):
        raise ValueError("Calibration must be a mapping")
    name = _safe_name(calibration.get("Name"))
    if calibration.get("Method") != "least_squares":
        raise ValueError("Calibration Method must be 'least_squares'")
    max_evaluations = calibration.get("MaxEvaluations")
    if isinstance(max_evaluations, bool) or not isinstance(max_evaluations, int) or max_evaluations <= 0:
        raise ValueError("Calibration MaxEvaluations must be a positive integer")

    observables = _parse_observables(calibration.get("Observables"), set(net.all_compartments))
    observations_value = calibration.get("Observations")
    if not isinstance(observations_value, str) or not observations_value:
        raise ValueError("Calibration Observations must be a file path")
    observations_path = Path(observations_value)
    observations, warnings = _load_observations(observations_path, observables, config)

    variables = tuple(
        [
            *_global_variables(config, calibration.get("Parameters")),
            *_initial_variables(net, y0, calibration.get("InitialConditions")),
        ]
    )
    if not variables:
        raise ValueError("Calibration must fit at least one global parameter or initial condition")
    if max_evaluations < len(variables) + 1:
        raise ValueError("Calibration MaxEvaluations must be at least p + 1")
    if len(observations) < len(variables):
        raise ValueError(
            f"Calibration is underdetermined: n={len(observations)} usable observations, p={len(variables)} variables"
        )
    if len(observations) == len(variables):
        warnings.append("Calibration has n == p; there is no residual redundancy")
    starts = _parse_starts(calibration, variables)
    if len(starts) == 1:
        warnings.append("Calibration uses one starting point for a bounded local method")
    return CalibrationPlan(
        name=name,
        observations_path=observations_path,
        observations=observations,
        observables=observables,
        variables=variables,
        starts=starts,
        max_evaluations=max_evaluations,
        warnings=tuple(warnings),
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False, default=str).encode("utf-8")


def _software_versions() -> dict[str, str]:
    try:
        patchsim_version = version("patchsim")
    except PackageNotFoundError:
        patchsim_version = "0.1.0"
    return {
        "patchsim": patchsim_version,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }


def _request_record(config: dict[str, Any], config_path: Path, plan: CalibrationPlan) -> dict[str, Any]:
    input_paths = {
        key: Path(config[key])
        for key in ("PatchFile", "SeedFile", "NetworkFile", "GroupFile", "InteractionFile")
        if isinstance(config.get(key), str)
    }
    return {
        "normalized_config": deepcopy(config),
        "source_config_sha256": _sha256(config_path),
        "observation_sha256": _sha256(plan.observations_path),
        "input_sha256": {key: _sha256(path) for key, path in input_paths.items()},
        "variable_order": [list(variable.identity) for variable in plan.variables],
        "starts": [list(start) for start in plan.starts],
        "method": {
            "name": "least_squares",
            "method": "trf",
            "jac": "2-point",
            "loss": "linear",
            "ftol": 1e-8,
            "xtol": 1e-8,
            "gtol": 1e-8,
            "max_forward_simulations_per_start": plan.max_evaluations,
        },
        "versions": _software_versions(),
    }


def _artifact_paths(target: Path) -> dict[str, str]:
    return {
        "output_dir": str(target),
        "estimates_path": str(target / "estimates.csv"),
        "fitted_seeds_path": str(target / "fitted-seeds.csv"),
        "attempts_path": str(target / "attempts.csv"),
        "residuals_path": str(target / "residuals.csv"),
        "manifest_path": str(target / "manifest.json"),
    }


def _reuse_existing(target: Path, fingerprint: str, plan: CalibrationPlan) -> dict[str, Any]:
    manifest_path = target / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FileExistsError(f"Calibration output {target} exists but is missing or modified") from exc
    request = manifest.get("request")
    if (
        manifest.get("schema_version") != 1
        or not isinstance(request, dict)
        or hashlib.sha256(_canonical_bytes(request)).hexdigest() != fingerprint
        or manifest.get("fingerprint") != fingerprint
        or manifest.get("n") != plan.n
        or manifest.get("p") != plan.p
        or manifest.get("start_count") != plan.start_count
    ):
        raise FileExistsError(f"Calibration output {target} belongs to a different study")
    artifacts = manifest.get("artifacts", {})
    for filename in _ARTIFACT_NAMES:
        path = target / filename
        if not path.is_file() or artifacts.get(filename, {}).get("sha256") != _sha256(path):
            raise FileExistsError(f"Calibration output {target} has a missing or modified artifact")
    return {
        **_artifact_paths(target),
        "reused": True,
        "n": manifest["n"],
        "p": manifest["p"],
        "start_count": manifest["start_count"],
        "selected_start": manifest["selected_start"],
        "forward_simulations": 0,
        "warnings": manifest.get("warnings", []),
    }


def _apply_vector(
    net: NetworkModel,
    baseline_y0: dict[str, float],
    variables: tuple[FitVariable, ...],
    vector: np.ndarray,
) -> dict[str, float]:
    y0 = baseline_y0.copy()
    remainder_variables: dict[str, FitVariable] = {}
    for variable, value in zip(variables, vector, strict=True):
        value = float(value)
        if variable.kind == "parameter":
            net.base_model.parameters[variable.name] = value
            for patch_parameters in net.patch_parameters.values():
                patch_parameters[variable.name] = value
        else:
            assert variable.state_key is not None and variable.remainder_key is not None
            y0[variable.state_key] = value
            remainder_variables.setdefault(variable.remainder_key, variable)
    for remainder_key, variable in remainder_variables.items():
        assert variable.patch_index is not None and variable.group_index is not None
        baseline_total = sum(
            baseline_y0[net.state_key(compartment, variable.patch_index, variable.group_index)]
            for compartment in net.base_model.compartments
        )
        non_remainder_sum = sum(
            y0[net.state_key(compartment, variable.patch_index, variable.group_index)]
            for compartment in net.base_model.compartments
            if net.state_key(compartment, variable.patch_index, variable.group_index) != remainder_key
        )
        remainder = baseline_total - non_remainder_sum
        if not np.isfinite(remainder) or remainder < 0:
            raise ValueError(
                f"Fitted initial state produced an invalid remainder for "
                f"patch {variable.patch!r}, group {variable.group!r}"
            )
        y0[remainder_key] = remainder
    return y0


def _prediction_arrays(frame: pd.DataFrame, plan: CalibrationPlan) -> tuple[np.ndarray, np.ndarray]:
    observable_map = {observable.name: observable for observable in plan.observables}
    predicted = []
    scales = []
    for row in plan.observations.itertuples(index=False):
        observable = observable_map[row.observable]
        predicted.append(float(frame.loc[row.grid_index, list(observable.columns)].sum(skipna=False)))
        scales.append(observable.scale)
    predictions = np.asarray(predicted, dtype=float)
    scale_values = np.asarray(scales, dtype=float)
    if not np.isfinite(predictions).all():
        raise ValueError("Calibration simulation produced non-finite predictions")
    return predictions, scale_values


def _seed_frame(net: NetworkModel, y0: dict[str, float]) -> pd.DataFrame:
    rows = []
    for patch_index, patch in enumerate(net.patch_names):
        for group_index, group in enumerate(net.groups or [None]):
            row: dict[str, Any] = {"patch": patch}
            if net.groups:
                row["group"] = group
            row.update(
                {
                    compartment: y0[net.state_key(compartment, patch_index, group_index)]
                    for compartment in net.base_model.compartments
                }
            )
            rows.append(row)
    return pd.DataFrame(rows)


def run_calibration(
    config_path: str | Path,
    *,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    started = time.monotonic()
    source_path = Path(config_path).expanduser().resolve()
    config = load_config(str(source_path))
    net, baseline_y0, _patches, _num_patches = setup_simulation(config)
    plan = get_calibration_plan(config, net, baseline_y0, required=True)
    assert plan is not None
    if progress:
        progress(f"Observations: {plan.n}; fitted variables: {plan.p}; starts: {plan.start_count}")
        progress(f"Maximum forward simulations: {plan.max_forward_simulations}")
        for warning in plan.warnings:
            progress(f"Warning: {warning}")
    if (plan.observations["grid_index"] < 0).any():
        raise ValueError("Calibration cannot start with unmatched observation times")

    request = _request_record(config, source_path, plan)
    fingerprint = hashlib.sha256(_canonical_bytes(request)).hexdigest()
    output_root = Path(config["OutputDir"]).resolve()
    calibration_root = (output_root / "calibration").resolve()
    if not calibration_root.is_relative_to(output_root):
        raise ValueError("Calibration output path escapes OutputDir")
    target = calibration_root / plan.name
    if target.is_symlink() or not target.resolve().is_relative_to(calibration_root):
        raise ValueError("Calibration output path escapes its study directory")
    if target.exists():
        result = _reuse_existing(target, fingerprint, plan)
        result["elapsed_seconds"] = time.monotonic() - started
        return result

    lower = np.asarray([variable.lower for variable in plan.variables], dtype=float)
    upper = np.asarray([variable.upper for variable in plan.variables], dtype=float)
    observed = plan.observations["value"].to_numpy(dtype=float)
    attempts: list[dict[str, Any]] = []
    successes: list[tuple[int, Any, int]] = []
    total_forward_simulations = 0

    for start_index, start_vector in enumerate(plan.starts):
        calls = 0

        def residual(vector):
            nonlocal calls, total_forward_simulations
            if calls >= plan.max_evaluations:
                raise _BudgetExceeded(f"forward-simulation budget {plan.max_evaluations} exhausted")
            calls += 1
            total_forward_simulations += 1
            y0 = _apply_vector(net, baseline_y0, plan.variables, np.asarray(vector, dtype=float))
            frame = _simulate_prepared(config, net, y0)
            predictions, scales = _prediction_arrays(frame, plan)
            values = (predictions - observed) / scales
            if not np.isfinite(values).all():
                raise ValueError("Calibration produced non-finite residuals")
            return values

        try:
            result = least_squares(
                residual,
                np.asarray(start_vector, dtype=float),
                jac="2-point",
                bounds=(lower, upper),
                method="trf",
                ftol=1e-8,
                xtol=1e-8,
                gtol=1e-8,
                x_scale=upper - lower,
                loss="linear",
                max_nfev=plan.max_evaluations,
            )
            finite = all(
                np.isfinite(value).all()
                for value in (np.asarray(result.x), np.asarray(result.fun), np.asarray(result.jac))
            ) and np.isfinite(result.cost)
            success = bool(result.success and finite)
            message = str(result.message).replace("\n", " ")[:500]
            attempts.append(
                {
                    "start": start_index,
                    "success": success,
                    "status": int(result.status),
                    "cost": float(result.cost) if np.isfinite(result.cost) else None,
                    "optimality": float(result.optimality) if np.isfinite(result.optimality) else None,
                    "nfev": int(result.nfev),
                    "njev": int(result.njev) if result.njev is not None else None,
                    "forward_simulations": calls,
                    "message": message,
                }
            )
            if success:
                successes.append((start_index, result, calls))
        except Exception as exc:
            attempts.append(
                {
                    "start": start_index,
                    "success": False,
                    "status": None,
                    "cost": None,
                    "optimality": None,
                    "nfev": None,
                    "njev": None,
                    "forward_simulations": calls,
                    "message": str(exc).replace("\n", " ")[:500],
                }
            )

    if not successes:
        details = "; ".join(f"start {attempt['start']}: {attempt['message']}" for attempt in attempts)
        raise RuntimeError(f"Calibration failed: no starting point terminated successfully ({details})")
    selected_start, selected, _selected_calls = min(successes, key=lambda item: (float(item[1].cost), item[0]))
    selected_vector = np.asarray(selected.x, dtype=float)
    selected_y0 = _apply_vector(net, baseline_y0, plan.variables, selected_vector)
    selected_residuals = np.asarray(selected.fun, dtype=float)
    observable_scales = {observable.name: observable.scale for observable in plan.observables}
    scales = plan.observations["observable"].map(observable_scales).to_numpy(dtype=float)
    predictions = observed + selected_residuals * scales

    singular_values = np.linalg.svd(np.asarray(selected.jac, dtype=float), compute_uv=False)
    rank_tolerance = (
        float(singular_values[0]) * max(plan.n, plan.p) * np.finfo(np.float64).eps if len(singular_values) else 0.0
    )
    rank = int(np.sum(singular_values > rank_tolerance))
    condition = float(singular_values[0] / singular_values[-1]) if rank == plan.p and singular_values[-1] > 0 else None

    estimate_rows = []
    for index, (variable, value) in enumerate(zip(plan.variables, selected_vector, strict=True)):
        estimate_rows.append(
            {
                "kind": variable.kind,
                "name": variable.name,
                "patch": variable.patch,
                "group": variable.group,
                "value": float(value),
                "lower": variable.lower,
                "upper": variable.upper,
                "active_bound": int(selected.active_mask[index]),
            }
        )
    estimates = pd.DataFrame(estimate_rows)
    fitted_seeds = _seed_frame(net, selected_y0)
    attempts_frame = pd.DataFrame(attempts)
    residuals = pd.DataFrame(
        {
            "time": plan.observations["time"].to_numpy(dtype=float),
            "observable": plan.observations["observable"].tolist(),
            "observed": observed,
            "prediction": predictions,
            "residual": predictions - observed,
            "scale": scales,
            "standardized_residual": selected_residuals,
        }
    )
    residual_summary = {
        "rmse": float(np.sqrt(np.mean(np.square(residuals["residual"].to_numpy(dtype=float))))),
        "standardized_rmse": float(np.sqrt(np.mean(np.square(selected_residuals)))),
        "by_observable": {},
    }
    for observable in plan.observables:
        subset = residuals.loc[residuals["observable"] == observable.name]
        raw = subset["residual"].to_numpy(dtype=float)
        standardized = subset["standardized_residual"].to_numpy(dtype=float)
        residual_summary["by_observable"][observable.name] = {
            "n": len(subset),
            "rmse": float(np.sqrt(np.mean(np.square(raw)))),
            "standardized_rmse": float(np.sqrt(np.mean(np.square(standardized)))),
        }

    calibration_root.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{plan.name}.", dir=calibration_root))
    try:
        frames = {
            "estimates.csv": estimates,
            "fitted-seeds.csv": fitted_seeds,
            "attempts.csv": attempts_frame,
            "residuals.csv": residuals,
        }
        for filename, frame in frames.items():
            frame.to_csv(temporary / filename, index=False, lineterminator="\n")
        artifacts = {filename: {"sha256": _sha256(temporary / filename)} for filename in _ARTIFACT_NAMES}
        manifest = {
            "schema_version": 1,
            "fingerprint": fingerprint,
            "request": request,
            "n": plan.n,
            "p": plan.p,
            "start_count": plan.start_count,
            "max_forward_simulations": plan.max_forward_simulations,
            "selected_start": selected_start,
            "selected_attempt": attempts[selected_start],
            "forward_simulations": total_forward_simulations,
            "warnings": list(plan.warnings),
            "residual_summary": residual_summary,
            "jacobian": {
                "singular_values": singular_values.tolist(),
                "rank_tolerance": rank_tolerance,
                "rank": rank,
                "rank_deficient": rank < plan.p,
                "condition_number": condition,
            },
            "artifacts": artifacts,
        }
        (temporary / "manifest.json").write_bytes(_canonical_bytes(manifest) + b"\n")
        temporary.replace(target)
    except BaseException:
        for child in temporary.iterdir():
            child.unlink()
        temporary.rmdir()
        raise

    return {
        **_artifact_paths(target),
        "reused": False,
        "n": plan.n,
        "p": plan.p,
        "start_count": plan.start_count,
        "selected_start": selected_start,
        "forward_simulations": total_forward_simulations,
        "warnings": list(plan.warnings),
        "elapsed_seconds": time.monotonic() - started,
    }

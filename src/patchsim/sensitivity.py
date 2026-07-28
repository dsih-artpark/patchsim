"""Sobol sensitivity studies for PatchSim configurations."""

from __future__ import annotations

import hashlib
import json
import re
import shutil
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

from patchsim.core.simulation import _simulate_prepared, load_config, setup_simulation

_NAME_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")
_INPUT_FIELDS = ("PatchFile", "SeedFile", "NetworkFile", "GroupFile", "InteractionFile")
_ARTIFACT_NAMES = ("samples.csv", "responses.csv", "indices.csv")
_MIN_SALIB_VERSION = (1, 5, 2)
_NUM_RESAMPLES = 100
_CONF_LEVEL = 0.95


@dataclass(frozen=True)
class Metric:
    name: str
    columns: tuple[str, ...]
    reducer: str


@dataclass(frozen=True)
class SensitivityPlan:
    name: str
    base_samples: int
    seed: int
    parameters: tuple[tuple[str, float, float], ...]
    metrics: tuple[Metric, ...]

    @property
    def evaluation_count(self) -> int:
        return self.base_samples * (len(self.parameters) + 2)


def _parse_version(value: str) -> tuple[int, int, int]:
    match = re.match(r"(\d+)\.(\d+)\.(\d+)", value)
    return tuple(map(int, match.groups())) if match else (0, 0, 0)


def _load_salib():
    try:
        salib_version = version("SALib")
        from SALib.analyze import sobol as sobol_analyze
        from SALib.sample import sobol as sobol_sample
    except (ImportError, PackageNotFoundError) as exc:
        raise RuntimeError(
            "Sensitivity analysis requires SALib 1.5.2 or newer. "
            'Install it with `python -m pip install "patchsim[analysis]"`.'
        ) from exc
    if _parse_version(salib_version) < _MIN_SALIB_VERSION:
        raise RuntimeError(
            f"Sensitivity analysis requires SALib 1.5.2 or newer; found {salib_version}. "
            'Upgrade with `python -m pip install --upgrade "patchsim[analysis]"`.'
        )
    return sobol_sample, sobol_analyze, salib_version


def _patch_parameter_names(config: dict[str, Any]) -> set[str]:
    names: set[str] = set()
    for entry in config.get("PatchParameters", []):
        if isinstance(entry, dict) and isinstance(entry.get("parameters", {}), dict):
            names.update(entry.get("parameters", {}))
    return names


def get_sensitivity_plan(
    config: dict[str, Any],
    output_columns: list[str],
    *,
    required: bool = True,
) -> SensitivityPlan | None:
    """Validate the optional Sensitivity block against a prepared model."""
    raw = config.get("Sensitivity")
    if raw is None:
        if required:
            raise ValueError("The configuration has no 'Sensitivity' block")
        return None
    if not isinstance(raw, dict):
        raise ValueError("'Sensitivity' must be a mapping")

    name = raw.get("Name")
    if not isinstance(name, str) or not _NAME_PATTERN.fullmatch(name) or name in {".", ".."}:
        raise ValueError("'Sensitivity.Name' must be one safe path component using letters, numbers, '.', '_', or '-'")
    if raw.get("Method") != "sobol":
        raise ValueError("'Sensitivity.Method' must be 'sobol'")

    base_samples = raw.get("BaseSamples")
    if (
        isinstance(base_samples, bool)
        or not isinstance(base_samples, int)
        or base_samples < 2
        or base_samples & (base_samples - 1)
    ):
        raise ValueError("'Sensitivity.BaseSamples' must be a power of two of at least 2")
    seed = raw.get("Seed")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("'Sensitivity.Seed' must be a non-negative integer")

    global_parameters = config.get("Parameters", {})
    if not isinstance(global_parameters, dict):
        raise ValueError("'Parameters' must be a mapping")
    patch_parameters = _patch_parameter_names(config)
    raw_parameters = raw.get("Parameters")
    if not isinstance(raw_parameters, dict) or not raw_parameters:
        raise ValueError("'Sensitivity.Parameters' must contain at least one parameter")

    parameters = []
    for parameter_name, bounds in raw_parameters.items():
        if not isinstance(parameter_name, str) or not parameter_name or parameter_name == "sample_id":
            raise ValueError("Sensitivity parameter names must be non-empty and may not be 'sample_id'")
        if parameter_name not in global_parameters:
            raise ValueError(f"Unknown global sensitivity parameter: {parameter_name!r}")
        if parameter_name in patch_parameters:
            raise ValueError(f"Cannot sample {parameter_name!r}; it is also set in PatchParameters")
        if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
            raise ValueError(f"Bounds for {parameter_name!r} must be [lower, upper]")
        lower, upper = bounds
        if any(isinstance(value, bool) or not isinstance(value, Real) for value in bounds):
            raise ValueError(f"Bounds for {parameter_name!r} must be finite real numbers")
        lower, upper = float(lower), float(upper)
        if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
            raise ValueError(f"Bounds for {parameter_name!r} must satisfy finite lower < upper")
        parameters.append((parameter_name, lower, upper))

    raw_metrics = raw.get("Metrics")
    if not isinstance(raw_metrics, dict) or not raw_metrics:
        raise ValueError("'Sensitivity.Metrics' must contain at least one metric")
    parameter_names = {name for name, _lower, _upper in parameters}
    available_columns = set(output_columns)
    metrics = []
    for metric_name, metric_config in raw_metrics.items():
        if (
            not isinstance(metric_name, str)
            or not metric_name
            or metric_name == "sample_id"
            or metric_name in parameter_names
        ):
            raise ValueError(
                "Sensitivity metric names must be non-empty and distinct from parameter columns and 'sample_id'"
            )
        if not isinstance(metric_config, dict):
            raise ValueError(f"Metric {metric_name!r} must be a mapping")
        columns = metric_config.get("Columns")
        if (
            not isinstance(columns, list)
            or not columns
            or not all(isinstance(column, str) for column in columns)
            or len(columns) != len(set(columns))
        ):
            raise ValueError(f"Metric {metric_name!r} Columns must be a non-empty list of unique names")
        unknown = sorted(set(columns) - available_columns)
        if unknown:
            raise ValueError(f"Metric {metric_name!r} references unknown output columns: {unknown}")
        reducer = metric_config.get("Reduce")
        if reducer not in {"max", "final"}:
            raise ValueError(f"Metric {metric_name!r} Reduce must be 'max' or 'final'")
        metrics.append(Metric(metric_name, tuple(columns), reducer))

    return SensitivityPlan(name, base_samples, seed, tuple(parameters), tuple(metrics))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False, default=str).encode()


def _versions(salib_version: str) -> dict[str, str]:
    result = {"SALib": salib_version}
    for package in ("patchsim", "numpy", "pandas", "scipy"):
        try:
            result[package] = version(package)
        except PackageNotFoundError:
            result[package] = "unknown"
    return result


def _request_record(
    config: dict[str, Any],
    salib_version: str,
    source_config_sha256: str,
) -> dict[str, Any]:
    inputs = {
        field: {"path": str(config[field]), "sha256": _sha256(Path(config[field]))}
        for field in _INPUT_FIELDS
        if config.get(field)
    }
    method = {
        "name": "sobol",
        "calc_second_order": False,
        "scramble": True,
        "skip_values": 0,
        "num_resamples": _NUM_RESAMPLES,
        "conf_level": _CONF_LEVEL,
    }
    return {
        "normalized_config": json.loads(json.dumps(config, default=str)),
        "source_config_sha256": source_config_sha256,
        "inputs": inputs,
        "method": method,
        "versions": _versions(salib_version),
    }


def _artifact_paths(target: Path) -> dict[str, str]:
    return {
        "output_dir": str(target),
        "samples_path": str(target / "samples.csv"),
        "responses_path": str(target / "responses.csv"),
        "indices_path": str(target / "indices.csv"),
        "manifest_path": str(target / "manifest.json"),
    }


def _reuse_existing(
    target: Path,
    fingerprint: str,
    evaluation_count: int,
) -> dict[str, Any]:
    manifest_path = target / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FileExistsError(f"Existing sensitivity target is incomplete: {target}") from exc
    request = manifest.get("request")
    if (
        manifest.get("schema_version") != 1
        or not isinstance(request, dict)
        or hashlib.sha256(_canonical_bytes(request)).hexdigest() != fingerprint
        or manifest.get("study_fingerprint") != fingerprint
        or manifest.get("evaluation_count") != evaluation_count
        or manifest.get("source_config", {}).get("sha256") != request.get("source_config_sha256")
    ):
        raise FileExistsError(f"Existing sensitivity target was produced by a different study: {target}")

    recorded = manifest.get("artifacts", {})
    for filename in _ARTIFACT_NAMES:
        path = target / filename
        expected = recorded.get(filename, {}).get("sha256")
        if not path.is_file() or expected != _sha256(path):
            raise FileExistsError(f"Existing sensitivity artifact is missing or modified: {path}")

    return {
        **_artifact_paths(target),
        "reused": True,
        "planned_evaluations": evaluation_count,
        "completed_evaluations": 0,
    }


def _metric_value(frame: pd.DataFrame, metric: Metric) -> float:
    series = frame.loc[:, metric.columns].sum(axis=1, skipna=False)
    return float(series.max() if metric.reducer == "max" else series.iloc[-1])


def run_sensitivity(
    config_path: str | Path,
    *,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Run or reuse one configured Sobol sensitivity study."""
    started = time.monotonic()
    source_path = Path(config_path).expanduser().resolve()
    config = load_config(str(source_path))
    request_config = deepcopy(config)
    net, y0, _patches, _num_patches = setup_simulation(config)
    plan = get_sensitivity_plan(config, list(net.all_compartments))
    assert plan is not None
    if progress:
        progress(f"Planned model evaluations: {plan.evaluation_count}")

    sobol_sample, sobol_analyze, salib_version = _load_salib()
    source_config_sha256 = _sha256(source_path)
    request = _request_record(request_config, salib_version, source_config_sha256)
    fingerprint = hashlib.sha256(_canonical_bytes(request)).hexdigest()
    output_root = Path(config["OutputDir"]).resolve()
    sensitivity_root = output_root / "sensitivity"
    resolved_sensitivity_root = sensitivity_root.resolve()
    if not resolved_sensitivity_root.is_relative_to(output_root):
        raise ValueError(f"Sensitivity output path escapes OutputDir: {sensitivity_root}")
    target = sensitivity_root / plan.name
    if target.is_symlink() or not target.resolve().is_relative_to(resolved_sensitivity_root):
        raise ValueError(f"Sensitivity output path escapes its study directory: {target}")
    if target.exists():
        summary = _reuse_existing(target, fingerprint, plan.evaluation_count)
        summary["elapsed_seconds"] = time.monotonic() - started
        return summary

    names = [name for name, _lower, _upper in plan.parameters]
    problem = {
        "num_vars": len(names),
        "names": names,
        "bounds": [[lower, upper] for _name, lower, upper in plan.parameters],
    }
    samples = sobol_sample.sample(
        problem,
        plan.base_samples,
        calc_second_order=False,
        scramble=True,
        skip_values=0,
        seed=plan.seed,
    )
    responses = {metric.name: [] for metric in plan.metrics}

    for sample_id, sample in enumerate(samples):
        values = dict(zip(names, map(float, sample), strict=True))
        net.base_model.parameters.update(values)
        for patch_parameters in net.patch_parameters.values():
            patch_parameters.update(values)
        try:
            frame = _simulate_prepared(config, net, y0)
            metric_values = {metric.name: _metric_value(frame, metric) for metric in plan.metrics}
        except Exception as exc:
            raise RuntimeError(f"Sensitivity evaluation {sample_id} failed for parameters {values}: {exc}") from exc
        non_finite = [name for name, value in metric_values.items() if not np.isfinite(value)]
        if non_finite:
            raise RuntimeError(
                f"Sensitivity evaluation {sample_id} produced non-finite metrics {non_finite} for parameters {values}"
            )
        for name, value in metric_values.items():
            responses[name].append(value)

    index_rows = []
    for metric in plan.metrics:
        values = np.asarray(responses[metric.name], dtype=float)
        if np.ptp(values) == 0:
            raise ValueError(f"Sensitivity metric {metric.name!r} is constant at {values[0]}")
        indices = sobol_analyze.analyze(
            problem,
            values,
            calc_second_order=False,
            num_resamples=_NUM_RESAMPLES,
            conf_level=_CONF_LEVEL,
            print_to_console=False,
            parallel=False,
            seed=plan.seed,
        )
        for index, parameter_name in enumerate(names):
            estimates = {
                "S1": float(indices["S1"][index]),
                "S1_conf": float(indices["S1_conf"][index]),
                "ST": float(indices["ST"][index]),
                "ST_conf": float(indices["ST_conf"][index]),
            }
            if not all(np.isfinite(value) for value in estimates.values()):
                raise ValueError(
                    f"Sobol analysis produced non-finite indices for metric "
                    f"{metric.name!r}, parameter {parameter_name!r}"
                )
            index_rows.append(
                {
                    "metric": metric.name,
                    "parameter": parameter_name,
                    **estimates,
                }
            )

    samples_frame = pd.DataFrame(samples, columns=names)
    samples_frame.insert(0, "sample_id", np.arange(len(samples), dtype=int))
    responses_frame = pd.DataFrame(responses)
    responses_frame.insert(0, "sample_id", np.arange(len(samples), dtype=int))
    indices_frame = pd.DataFrame(index_rows)

    parent = target.parent
    parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{plan.name}.", dir=parent))
    try:
        samples_frame.to_csv(temporary / "samples.csv", index=False, lineterminator="\n")
        responses_frame.to_csv(temporary / "responses.csv", index=False, lineterminator="\n")
        indices_frame.to_csv(temporary / "indices.csv", index=False, lineterminator="\n")
        artifacts = {filename: {"sha256": _sha256(temporary / filename)} for filename in _ARTIFACT_NAMES}
        manifest = {
            "schema_version": 1,
            "study_fingerprint": fingerprint,
            "evaluation_count": plan.evaluation_count,
            "request": request,
            "source_config": {"path": str(source_path), "sha256": source_config_sha256},
            "artifacts": artifacts,
        }
        (temporary / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        try:
            temporary.rename(target)
        except OSError:
            if not target.exists():
                raise
            summary = _reuse_existing(target, fingerprint, plan.evaluation_count)
            summary["elapsed_seconds"] = time.monotonic() - started
            return summary
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)

    return {
        **_artifact_paths(target),
        "reused": False,
        "planned_evaluations": plan.evaluation_count,
        "completed_evaluations": plan.evaluation_count,
        "elapsed_seconds": time.monotonic() - started,
    }

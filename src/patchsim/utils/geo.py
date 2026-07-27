"""Spatial contact-matrix generation with explicit units and validation."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd

EARTH_RADIUS_KM = 6371.0088
REPORT_SCHEMA_VERSION = 1
_SOURCE_SUFFIXES = {".csv", ".geojson", ".json", ".shp"}
_SHAPEFILE_SUFFIXES = {".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix", ".sbn", ".sbx"}

KernelName = Literal["distance", "gravity"]
Normalization = Literal["none", "row"]


def _as_finite_float(name: str, value: Any, *, positive: bool = False, non_negative: bool = False) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not np.isfinite(number):
        raise ValueError(f"{name} must be finite")
    if positive and number <= 0:
        raise ValueError(f"{name} must be greater than zero")
    if non_negative and number < 0:
        raise ValueError(f"{name} must be non-negative")
    return number


def _normalized_identifiers(values: pd.Series, column: str) -> pd.Series:
    null_rows = values.index[values.isna()].tolist()
    if null_rows:
        raise ValueError(f"Identifier column '{column}' contains null values at rows {null_rows}")

    identifiers = values.astype(str).str.strip()
    empty_rows = identifiers.index[identifiers.eq("")].tolist()
    if empty_rows:
        raise ValueError(f"Identifier column '{column}' contains empty values at rows {empty_rows}")

    duplicated = sorted(identifiers[identifiers.duplicated(keep=False)].unique())
    if duplicated:
        raise ValueError(f"Identifier column '{column}' contains duplicate values after normalization: {duplicated}")
    return identifiers


def _numeric_series(values: pd.Series, name: str, *, positive: bool = False) -> pd.Series:
    numbers = pd.to_numeric(values, errors="coerce").astype(float)
    bad_rows = numbers.index[~np.isfinite(numbers)].tolist()
    if bad_rows:
        raise ValueError(f"Column '{name}' must contain finite numeric values; invalid rows: {bad_rows}")
    if positive:
        bad_rows = numbers.index[numbers.le(0)].tolist()
        if bad_rows:
            raise ValueError(f"Column '{name}' must contain positive values; invalid rows: {bad_rows}")
    return numbers


def _validated_regions(
    identifiers: pd.Series,
    latitudes: pd.Series,
    longitudes: pd.Series,
    *,
    id_column: str,
    population: pd.Series | None = None,
    population_column: str | None = None,
) -> pd.DataFrame:
    ids = _normalized_identifiers(identifiers, id_column)
    lat = _numeric_series(latitudes, "lat")
    lon = _numeric_series(longitudes, "lon")

    bad_lat = lat.index[~lat.between(-90, 90)].tolist()
    bad_lon = lon.index[~lon.between(-180, 180)].tolist()
    if bad_lat:
        raise ValueError(f"Latitude must be within [-90, 90]; invalid rows: {bad_lat}")
    if bad_lon:
        raise ValueError(f"Longitude must be within [-180, 180]; invalid rows: {bad_lon}")
    if len(ids) < 2:
        raise ValueError("Contact generation requires at least two regions")

    result = pd.DataFrame({"id": ids.to_numpy(), "lat": lat.to_numpy(), "lon": lon.to_numpy()})
    if population is not None:
        if population_column is None:
            raise ValueError("A population column name is required when population values are provided")
        result["population"] = _numeric_series(population, population_column, positive=True).to_numpy()
    return result


def _required_column(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        raise ValueError(f"Required column '{column}' not found")
    return frame[column]


def _load_centroid_csv(
    source: Path,
    *,
    id_column: str,
    population_column: str | None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    frame = pd.read_csv(source, dtype={id_column: "string"}, keep_default_na=False)
    population = _required_column(frame, population_column) if population_column else None
    regions = _validated_regions(
        _required_column(frame, id_column),
        _required_column(frame, "lat"),
        _required_column(frame, "lon"),
        id_column=id_column,
        population=population,
        population_column=population_column,
    )
    return regions, {"source_crs": None, "centroid_crs": None}


def _load_vector(
    source: Path,
    *,
    id_column: str,
    population_column: str | None,
    centroid_crs: str | None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    try:
        import geopandas as gpd
        from pyproj import CRS
    except ImportError as exc:  # pragma: no cover - exercised in an isolated environment
        raise RuntimeError("Vector input requires the optional geo dependencies: install 'patchsim[geo]'") from exc

    frame = gpd.read_file(source)
    identifiers = _required_column(frame, id_column)
    normalized_ids = _normalized_identifiers(identifiers, id_column)
    population = _required_column(frame, population_column) if population_column else None

    if frame.crs is None:
        raise ValueError("Vector input must declare a coordinate reference system (CRS)")
    invalid_geometry = frame.geometry.isna() | frame.geometry.is_empty | ~frame.geometry.is_valid
    if invalid_geometry.any():
        bad_ids = normalized_ids[invalid_geometry].tolist()
        raise ValueError(f"Vector input contains null, empty, or invalid geometries for identifiers: {bad_ids}")

    geometry_types = set(frame.geometry.geom_type)
    if geometry_types == {"Point"}:
        if centroid_crs is not None:
            raise ValueError("--centroid-crs is only valid for Polygon or MultiPolygon input")
        points = frame.to_crs("EPSG:4326").geometry
        resolved_centroid_crs = None
    elif geometry_types.issubset({"Polygon", "MultiPolygon"}):
        if centroid_crs is None:
            raise ValueError("Polygon input requires --centroid-crs with a suitable projected CRS")
        target_crs = CRS.from_user_input(centroid_crs)
        if not target_crs.is_projected:
            raise ValueError("--centroid-crs must be a projected CRS")
        projected = frame.to_crs(target_crs)
        centroids = gpd.GeoSeries(projected.geometry.centroid, index=frame.index, crs=target_crs)
        points = centroids.to_crs("EPSG:4326")
        resolved_centroid_crs = target_crs.to_string()
    else:
        raise ValueError(
            "Vector input must contain only Point geometries or only Polygon/MultiPolygon geometries; "
            f"received {sorted(geometry_types)}"
        )

    regions = _validated_regions(
        normalized_ids,
        pd.Series(points.y, index=frame.index),
        pd.Series(points.x, index=frame.index),
        id_column=id_column,
        population=population,
        population_column=population_column,
    )
    return regions, {
        "source_crs": frame.crs.to_string(),
        "centroid_crs": resolved_centroid_crs,
    }


def load_contact_regions(
    source: str | Path,
    *,
    id_column: str,
    population_column: str | None = None,
    centroid_crs: str | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load validated contact regions from a centroid CSV or vector file."""
    source_path = Path(source)
    suffix = source_path.suffix.lower()
    if suffix not in _SOURCE_SUFFIXES:
        raise ValueError(f"Unsupported source format '{suffix}'; use CSV, GeoJSON, JSON, or Shapefile")
    if suffix == ".csv":
        if centroid_crs is not None:
            raise ValueError("--centroid-crs is only valid for polygon vector input")
        return _load_centroid_csv(
            source_path,
            id_column=id_column,
            population_column=population_column,
        )
    return _load_vector(
        source_path,
        id_column=id_column,
        population_column=population_column,
        centroid_crs=centroid_crs,
    )


def distance_matrix_km(regions: pd.DataFrame) -> np.ndarray:
    """Return pairwise haversine distances in kilometres."""
    lat = np.radians(regions["lat"].to_numpy(dtype=float))
    lon = np.radians(regions["lon"].to_numpy(dtype=float))
    dlat = lat[:, None] - lat[None, :]
    dlon = lon[:, None] - lon[None, :]
    haversine = np.sin(dlat / 2) ** 2 + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2) ** 2
    angular_distance = 2 * np.arcsin(np.sqrt(np.clip(haversine, 0.0, 1.0)))
    distances = EARTH_RADIUS_KM * angular_distance
    np.fill_diagonal(distances, 0.0)
    return distances


def _summary(values: np.ndarray) -> dict[str, float]:
    return {
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "median": float(np.median(values)),
    }


def generate_contact_matrix(
    regions: pd.DataFrame,
    *,
    kernel: KernelName,
    decay: float,
    min_distance_km: float,
    normalize: Normalization,
    scale: float | None = None,
    self_weight: float | None = None,
    self_share: float | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Generate and validate a spatial contact matrix."""
    if kernel not in {"distance", "gravity"}:
        raise ValueError(f"Unknown kernel '{kernel}'")
    if normalize not in {"none", "row"}:
        raise ValueError(f"Unknown normalization mode '{normalize}'")

    population = _required_column(regions, "population") if kernel == "gravity" else None
    regions = _validated_regions(
        _required_column(regions, "id"),
        _required_column(regions, "lat"),
        _required_column(regions, "lon"),
        id_column="id",
        population=population,
        population_column="population" if population is not None else None,
    )

    decay_value = _as_finite_float("decay", decay, positive=True)
    distance_floor = _as_finite_float("min-distance-km", min_distance_km, positive=True)

    if normalize == "none":
        if scale is None or self_weight is None:
            raise ValueError("--normalize none requires --scale and --self-weight")
        if self_share is not None:
            raise ValueError("--self-share is only valid with --normalize row")
        scale_value = _as_finite_float("scale", scale, positive=True)
        diagonal_value = _as_finite_float("self-weight", self_weight, non_negative=True)
    else:
        if self_share is None:
            raise ValueError("--normalize row requires --self-share")
        if scale is not None or self_weight is not None:
            raise ValueError("--scale and --self-weight are only valid with --normalize none")
        share = _as_finite_float("self-share", self_share, non_negative=True)
        if share >= 1:
            raise ValueError("self-share must be less than one so the selected kernel contributes")
        scale_value = 1.0
        diagonal_value = share

    if kernel == "gravity":
        populations = regions["population"].to_numpy(dtype=float)
        if not np.all(np.isfinite(populations)) or np.any(populations <= 0):
            raise ValueError("Gravity kernel populations must be finite and positive")
    else:
        populations = None

    distances = distance_matrix_km(regions)
    diagonal = np.eye(len(regions), dtype=bool)
    safe_distances = np.maximum(distances, distance_floor)
    safe_distances[diagonal] = 1.0

    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            denominator = np.power(safe_distances, decay_value)
            if populations is None:
                raw = scale_value / denominator
            else:
                raw = scale_value * np.outer(populations, populations) / denominator
    except FloatingPointError as exc:
        raise ValueError("Kernel arithmetic overflowed; revise decay, scale, population, or distance floor") from exc

    np.fill_diagonal(raw, 0.0)
    raw_off_diagonal = raw[~diagonal]
    if np.any(~np.isfinite(raw_off_diagonal)):
        raise ValueError("Kernel produced non-finite off-diagonal weights")
    if np.any(raw_off_diagonal <= 0):
        bad_pair = np.argwhere((raw <= 0) & ~diagonal)[0]
        source_id = regions.iloc[int(bad_pair[0])]["id"]
        target_id = regions.iloc[int(bad_pair[1])]["id"]
        raise ValueError(
            f"Kernel weight underflowed for '{source_id}' -> '{target_id}'; revise decay or distance floor"
        )

    if normalize == "row":
        row_sums = raw.sum(axis=1)
        if np.any(~np.isfinite(row_sums)) or np.any(row_sums <= 0):
            raise ValueError("Every raw kernel row must have a positive finite off-diagonal sum")
        matrix = raw / row_sums[:, None] * (1.0 - diagonal_value)
    else:
        matrix = raw.copy()
    np.fill_diagonal(matrix, diagonal_value)

    if matrix.shape != (len(regions), len(regions)):
        raise ValueError("Kernel output shape does not match the region count")
    if np.any(~np.isfinite(matrix)) or np.any(matrix < 0):
        raise ValueError("Final contact matrix must contain only finite, non-negative weights")

    pair_distances = distances[np.triu_indices(len(regions), k=1)]
    pair_weights = raw[np.triu_indices(len(regions), k=1)]
    dynamic_range = float(np.max(pair_weights) / np.min(pair_weights))
    if not np.isfinite(dynamic_range):
        raise ValueError("Kernel weight dynamic range exceeds floating-point capacity; revise parameters")
    diagnostics = {
        "distance": _summary(pair_distances),
        "distance_floor_pair_count": int(np.count_nonzero(pair_distances < distance_floor)),
        "off_diagonal_weight": {
            **_summary(pair_weights),
            "dynamic_range": dynamic_range,
        },
        "row_sum": {
            "min": float(np.min(matrix.sum(axis=1))),
            "max": float(np.max(matrix.sum(axis=1))),
        },
        "diagonal": {
            "min": float(np.min(np.diag(matrix))),
            "max": float(np.max(np.diag(matrix))),
        },
        "raw_symmetric": bool(np.allclose(raw, raw.T, rtol=1e-12, atol=0.0)),
        "final_symmetric": bool(np.allclose(matrix, matrix.T, rtol=1e-12, atol=0.0)),
    }
    return matrix, diagnostics


def _path_aliases(first: Path, second: Path) -> bool:
    if first == second:
        return True
    if first.exists() and second.exists():
        return os.path.samefile(first, second)
    return False


def _source_components(source: Path) -> list[Path]:
    if source.suffix.lower() != ".shp":
        return [source]
    return sorted(
        path.resolve()
        for path in source.parent.iterdir()
        if path.stem == source.stem and path.suffix.lower() in _SHAPEFILE_SUFFIXES
    )


def _write_temp(parent: Path, prefix: str, content: bytes) -> Path:
    descriptor, name = tempfile.mkstemp(prefix=prefix, dir=parent)
    temp_path = Path(name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise
    return temp_path


def _write_output_pair(
    output: Path,
    report_path: Path,
    csv_bytes: bytes,
    report_bytes: bytes,
) -> None:
    csv_temp = _write_temp(output.parent, f".{output.name}.", csv_bytes)
    report_temp: Path | None = None
    try:
        report_temp = _write_temp(report_path.parent, f".{report_path.name}.", report_bytes)
        os.replace(csv_temp, output)
        os.replace(report_temp, report_path)
    finally:
        csv_temp.unlink(missing_ok=True)
        if report_temp is not None:
            report_temp.unlink(missing_ok=True)


def generate_contacts(
    source: str | Path,
    output: str | Path,
    *,
    id_column: str,
    kernel: KernelName,
    decay: float,
    min_distance_km: float,
    normalize: Normalization,
    population_column: str | None = None,
    scale: float | None = None,
    self_weight: float | None = None,
    self_share: float | None = None,
    centroid_crs: str | None = None,
    force: bool = False,
) -> tuple[Path, Path, dict[str, Any]]:
    """Generate a runtime-compatible network CSV and its validation report."""
    source_path = Path(source).resolve(strict=True)
    output_path = Path(output).resolve()
    if output_path.suffix.lower() != ".csv":
        raise ValueError("Output path must end in .csv")
    report_path = Path(f"{output_path}.validation.json")

    for source_component in _source_components(source_path):
        if _path_aliases(source_component, output_path) or _path_aliases(source_component, report_path):
            raise ValueError(f"Output paths must not overwrite source data: {source_component}")
    if not force and (output_path.exists() or report_path.exists()):
        raise FileExistsError(
            f"Refusing to overwrite existing output pair: {output_path}, {report_path}. Use --force to overwrite."
        )

    if kernel == "gravity" and population_column is None:
        raise ValueError("Gravity kernel requires --population-column")

    regions, source_metadata = load_contact_regions(
        source_path,
        id_column=id_column,
        population_column=population_column,
        centroid_crs=centroid_crs,
    )
    matrix, diagnostics = generate_contact_matrix(
        regions,
        kernel=kernel,
        decay=decay,
        min_distance_km=min_distance_km,
        normalize=normalize,
        scale=scale,
        self_weight=self_weight,
        self_share=self_share,
    )

    identifiers = regions["id"].tolist()
    rows = [
        {
            "day": 0,
            "source": identifiers[source_index],
            "target": identifiers[target_index],
            "weight": float(matrix[source_index, target_index]),
        }
        for source_index in range(len(identifiers))
        for target_index in range(len(identifiers))
        if matrix[source_index, target_index] > 0
    ]
    csv_bytes = (
        pd.DataFrame(rows, columns=["day", "source", "target", "weight"])
        .to_csv(index=False, lineterminator="\n")
        .encode("utf-8")
    )

    report: dict[str, Any] = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "source": str(source_path),
        "id_column": id_column,
        "population_column": population_column,
        "region_count": len(identifiers),
        "identifiers": identifiers,
        **source_metadata,
        "coordinate_crs": "EPSG:4326",
        "distance": {
            "algorithm": "haversine",
            "earth_radius_km": EARTH_RADIUS_KM,
            "unit": "km",
            **diagnostics["distance"],
            "floor_pair_count": diagnostics["distance_floor_pair_count"],
        },
        "kernel": {
            "name": kernel,
            "decay": float(decay),
            "scale": float(scale) if scale is not None else 1.0,
            "min_distance_km": float(min_distance_km),
            "raw_weight_unit": (
                "scale * population_i * population_j / km**decay" if kernel == "gravity" else "scale / km**decay"
            ),
        },
        "normalization": {
            "mode": normalize,
            "self_weight": float(self_weight) if self_weight is not None else None,
            "self_share": float(self_share) if self_share is not None else None,
            "final_weight_unit": "dimensionless" if normalize == "row" else "raw kernel unit",
        },
        "matrix": {
            "off_diagonal_weight": diagnostics["off_diagonal_weight"],
            "row_sum": diagnostics["row_sum"],
            "diagonal": diagnostics["diagonal"],
            "raw_symmetric": diagnostics["raw_symmetric"],
            "final_symmetric": diagnostics["final_symmetric"],
        },
        "csv_sha256": hashlib.sha256(csv_bytes).hexdigest(),
    }
    report_bytes = (json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    _write_output_pair(output_path, report_path, csv_bytes, report_bytes)
    return output_path, report_path, report


__all__ = [
    "EARTH_RADIUS_KM",
    "distance_matrix_km",
    "generate_contact_matrix",
    "generate_contacts",
    "load_contact_regions",
]

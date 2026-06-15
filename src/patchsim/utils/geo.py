import json
from math import atan2, cos, radians, sin, sqrt
from typing import Any, Dict

import numpy as np
import pandas as pd


def load_geofile_centroids(path: str, id_prop: str = "id", pop_prop_candidates=None) -> pd.DataFrame:
    """Load GeoJSON or Shapefile and return centroids DataFrame with columns ['id','lat','lon','population'].

    For GeoJSON this delegates to load_geojson_centroids. For shapefiles requires geopandas to be
    installed and will compute geometry centroids.
    """
    if path.lower().endswith(".shp"):
        try:
            import geopandas as gpd
        except Exception as e:  # pragma: no cover - optional dependency
            raise RuntimeError("Reading shapefiles requires geopandas. Install it or provide GeoJSON.") from e

        gdf = gpd.read_file(path)
        # compute centroids
        gdf = gdf.to_crs(epsg=4326)
        cents = gdf.geometry.centroid
        records = []
        for idx, row in gdf.iterrows():
            identifier = row.get(id_prop) or row.get("name") or row.get("NAME") or idx
            population = None
            if pop_prop_candidates is None:
                pop_prop_candidates = ["population", "pop", "pop_total", "POPULATION"]
            for key in pop_prop_candidates:
                if key in row and isinstance(row[key], (int, float)):
                    population = row[key]
                    break
            lon = cents.iloc[idx].x if hasattr(cents.iloc[idx], "x") else None
            lat = cents.iloc[idx].y if hasattr(cents.iloc[idx], "y") else None
            records.append({"id": identifier, "lon": lon, "lat": lat, "population": population})
        return pd.DataFrame.from_records(records)
    else:
        return load_geojson_centroids(path, id_prop=id_prop, pop_prop_candidates=pop_prop_candidates)


def generate_kernels(
    df: pd.DataFrame,
    methods: list[str] | None = None,
    decay: float = 2.0,
    scale: float = 1.0,
    min_distance_km: float = 1e-3,
) -> dict:
    """Generate contact kernels for the requested methods.

    Returns a dict mapping method name -> numpy.ndarray contact matrix.
    Supported methods: 'gravity', 'distance'
    """
    if methods is None:
        methods = ["gravity"]

    results: dict[str, np.ndarray] = {}
    for m in methods:
        if m == "gravity":
            if "population" not in df.columns:
                raise ValueError("Population column 'population' required for gravity kernel")
            results[m] = gravity_contact_matrix(
                df,
                pop_col="population",
                decay=decay,
                scale=scale,
                min_distance_km=min_distance_km,
            )
        elif m == "distance":
            D = distance_matrix_km(df)
            D_safe = np.maximum(D, min_distance_km)
            W = 1.0 / (D_safe ** decay)
            np.fill_diagonal(W, 0.0)
            results[m] = W
        else:
            raise ValueError(f"Unknown kernel method: {m}")

    return results


def _haversine_km(lat1, lon1, lat2, lon2):
    # approximate radius of earth in km
    R = 6371.0
    dlat = radians(lat2 - lat1)
    dlon = radians(lon2 - lon1)
    a = sin(dlat / 2) ** 2 + cos(radians(lat1)) * cos(radians(lat2)) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    return R * c


def load_geojson_centroids(path: str, id_prop: str = "id", pop_prop_candidates=None) -> pd.DataFrame:
    """Load GeoJSON and return a DataFrame with columns: ['id','lat','lon','population'].

    The function tries to infer a population property from common names if not provided.
    """
    if pop_prop_candidates is None:
        pop_prop_candidates = ["population", "pop", "pop_total", "POPULATION"]

    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    records = []
    for feat in data.get("features", []):
        props: Dict[str, Any] = feat.get("properties", {})
        geom = feat.get("geometry") or {}
        coords = None
        # handle Point or Polygon/MultiPolygon (use centroid)
        if geom.get("type") == "Point":
            coords = geom.get("coordinates")
        else:
            # compute centroid from coordinates arrays
            try:
                # flatten all coordinate tuples
                all_coords = []
                coords_list = geom.get("coordinates", [])

                # iterative flatten to avoid nested function closure issues
                if coords_list:
                    stack = [coords_list]
                    while stack:
                        c = stack.pop()
                        if isinstance(c[0], (float, int)):
                            all_coords.append(tuple(c))
                        else:
                            for x in c:
                                stack.append(x)
                    lon_mean = sum(c[0] for c in all_coords) / len(all_coords)
                    lat_mean = sum(c[1] for c in all_coords) / len(all_coords)
                    coords = [lon_mean, lat_mean]
            except Exception:
                coords = None

        if not coords:
            continue

        identifier = props.get(id_prop) or props.get("name") or props.get("NAME")

        population = None
        for key in ([id_prop] if id_prop in props else []) + pop_prop_candidates:
            if key in props and isinstance(props[key], (int, float)):
                population = props[key]
                break

        records.append({"id": identifier, "lon": coords[0], "lat": coords[1], "population": population})

    df = pd.DataFrame.from_records(records)
    return df


def distance_matrix_km(df: pd.DataFrame) -> np.ndarray:
    """Return pairwise distance matrix in kilometers from DataFrame with `lat` and `lon` columns."""
    lats = df["lat"].to_numpy()
    lons = df["lon"].to_numpy()
    n = len(lats)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = _haversine_km(lats[i], lons[i], lats[j], lons[j])
            D[i, j] = d
            D[j, i] = d
    return D


def gravity_contact_matrix(
    df: pd.DataFrame,
    pop_col: str = "population",
    decay: float = 2.0,
    scale: float = 1.0,
    min_distance_km: float = 1e-3,
) -> np.ndarray:
    """Generate a contact matrix using a simple gravity model.

    C_ij = scale * (pop_i * pop_j) / (distance_ij ** decay)
    Diagonal is set to zero.
    """
    if pop_col not in df.columns:
        raise ValueError(f"Population column '{pop_col}' not found in DataFrame")

    pops = df[pop_col].to_numpy(dtype=float)
    D = distance_matrix_km(df)
    D_safe = np.maximum(D, min_distance_km)
    outer = np.outer(pops, pops)
    C = scale * outer / (D_safe ** decay)
    np.fill_diagonal(C, 0.0)
    return C


__all__ = [
    "distance_matrix_km",
    "generate_kernels",
    "gravity_contact_matrix",
    "load_geofile_centroids",
    "load_geojson_centroids",
]

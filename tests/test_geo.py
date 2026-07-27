import builtins
import hashlib
import json

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point, Polygon

from patchsim.core.simulation import setup_simulation
from patchsim.utils.geo import (
    EARTH_RADIUS_KM,
    distance_matrix_km,
    generate_contact_matrix,
    generate_contacts,
    load_contact_regions,
)


@pytest.fixture
def regions():
    return pd.DataFrame(
        {
            "id": ["A", "B"],
            "lat": [0.0, 0.0],
            "lon": [0.0, 1.0],
            "population": [100.0, 200.0],
        }
    )


def test_distance_matrix_uses_haversine_kilometres(regions):
    distances = distance_matrix_km(regions)
    expected = EARTH_RADIUS_KM * np.pi / 180
    np.testing.assert_allclose(distances, [[0.0, expected], [expected, 0.0]])


def test_distance_kernel_preserves_raw_units(regions):
    matrix, diagnostics = generate_contact_matrix(
        regions,
        kernel="distance",
        decay=1.0,
        scale=2.0,
        min_distance_km=0.001,
        normalize="none",
        self_weight=3.0,
    )
    expected_cross_weight = 2.0 / (EARTH_RADIUS_KM * np.pi / 180)
    np.testing.assert_allclose(matrix, [[3.0, expected_cross_weight], [expected_cross_weight, 3.0]])
    assert diagnostics["raw_symmetric"] is True
    assert diagnostics["final_symmetric"] is True


def test_gravity_kernel_uses_both_patch_populations(regions):
    matrix, _diagnostics = generate_contact_matrix(
        regions,
        kernel="gravity",
        decay=2.0,
        scale=0.5,
        min_distance_km=0.001,
        normalize="none",
        self_weight=1.0,
    )
    distance = EARTH_RADIUS_KM * np.pi / 180
    expected_cross_weight = 0.5 * 100.0 * 200.0 / distance**2
    np.testing.assert_allclose(matrix, [[1.0, expected_cross_weight], [expected_cross_weight, 1.0]])


def test_row_normalization_has_explicit_self_share(regions):
    matrix, diagnostics = generate_contact_matrix(
        regions,
        kernel="distance",
        decay=2.0,
        min_distance_km=0.001,
        normalize="row",
        self_share=0.8,
    )
    np.testing.assert_allclose(matrix, [[0.8, 0.2], [0.2, 0.8]])
    assert diagnostics["row_sum"] == pytest.approx({"min": 1.0, "max": 1.0})


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        (
            {"normalize": "none", "scale": None, "self_weight": 1.0},
            "requires --scale and --self-weight",
        ),
        (
            {"normalize": "row", "scale": 1.0, "self_share": 0.5},
            "only valid with --normalize none",
        ),
        (
            {"normalize": "row", "self_share": 1.0},
            "less than one",
        ),
        (
            {"normalize": "none", "scale": 1.0, "self_weight": -1.0},
            "must be non-negative",
        ),
    ],
)
def test_kernel_rejects_invalid_normalization_options(regions, kwargs, message):
    with pytest.raises(ValueError, match=message):
        generate_contact_matrix(
            regions,
            kernel="distance",
            decay=2.0,
            min_distance_km=0.001,
            **kwargs,
        )


def test_kernel_rejects_floating_point_underflow(regions):
    with pytest.raises(ValueError, match="underflowed"):
        generate_contact_matrix(
            regions,
            kernel="distance",
            decay=20.0,
            scale=1e-300,
            min_distance_km=0.001,
            normalize="none",
            self_weight=1.0,
        )


def test_public_matrix_generator_rejects_invalid_coordinates(regions):
    regions.loc[1, "lat"] = 100
    with pytest.raises(ValueError, match="Latitude"):
        generate_contact_matrix(
            regions,
            kernel="distance",
            decay=2.0,
            min_distance_km=0.001,
            normalize="row",
            self_share=0.8,
        )


def test_centroid_csv_validation_normalizes_identifiers(tmp_path):
    source = tmp_path / "regions.csv"
    pd.DataFrame({"patch": ["A", " A "], "lat": [0, 1], "lon": [0, 1]}).to_csv(source, index=False)
    with pytest.raises(ValueError, match="duplicate values after normalization"):
        load_contact_regions(source, id_column="patch")


def test_centroid_csv_preserves_text_identifiers(tmp_path):
    source = tmp_path / "regions.csv"
    source.write_text("patch,lat,lon\n001,0,0\nNA,0,1\n", encoding="utf-8")

    loaded, _metadata = load_contact_regions(source, id_column="patch")
    assert loaded["id"].tolist() == ["001", "NA"]


@pytest.mark.parametrize(
    ("column", "values", "message"),
    [
        ("lat", [0, 91], "Latitude"),
        ("lon", [0, -181], "Longitude"),
        ("population", [100, 0], "positive"),
    ],
)
def test_centroid_csv_rejects_invalid_physical_values(tmp_path, column, values, message):
    data = {
        "patch": ["A", "B"],
        "lat": [0, 1],
        "lon": [0, 1],
        "population": [100, 200],
    }
    data[column] = values
    source = tmp_path / "regions.csv"
    pd.DataFrame(data).to_csv(source, index=False)
    with pytest.raises(ValueError, match=message):
        load_contact_regions(source, id_column="patch", population_column="population")


def test_point_geojson_is_transformed_to_wgs84(tmp_path):
    source = tmp_path / "points.geojson"
    gpd.GeoDataFrame(
        {"patch": ["A", "B"]},
        geometry=[Point(0, 0), Point(111319.490793, 0)],
        crs="EPSG:3857",
    ).to_file(source, driver="GeoJSON")

    loaded, metadata = load_contact_regions(source, id_column="patch")
    assert loaded["lon"].tolist() == pytest.approx([0.0, 1.0])
    assert metadata == {"source_crs": "EPSG:3857", "centroid_crs": None}


def test_point_shapefile_is_supported(tmp_path):
    source = tmp_path / "points.shp"
    gpd.GeoDataFrame(
        {"patch": ["A", "B"]},
        geometry=[Point(0, 0), Point(1, 0)],
        crs="EPSG:4326",
    ).to_file(source)

    loaded, metadata = load_contact_regions(source, id_column="patch")
    assert loaded["id"].tolist() == ["A", "B"]
    assert metadata["source_crs"] == "EPSG:4326"


def test_polygon_geojson_requires_explicit_projected_centroid_crs(tmp_path):
    source = tmp_path / "polygons.geojson"
    gpd.GeoDataFrame(
        {"patch": ["A", "B"]},
        geometry=[
            Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
            Polygon([(2, 0), (3, 0), (3, 1), (2, 1)]),
        ],
        crs="EPSG:4326",
    ).to_file(source, driver="GeoJSON")

    with pytest.raises(ValueError, match="requires --centroid-crs"):
        load_contact_regions(source, id_column="patch")

    loaded, metadata = load_contact_regions(source, id_column="patch", centroid_crs="EPSG:3857")
    assert loaded["lon"].tolist() == pytest.approx([0.5, 2.5])
    assert metadata["centroid_crs"] == "EPSG:3857"


def test_vector_input_requires_declared_crs(tmp_path, monkeypatch):
    source = tmp_path / "points.geojson"
    source.write_text("{}", encoding="utf-8")
    frame = gpd.GeoDataFrame({"patch": ["A", "B"]}, geometry=[Point(0, 0), Point(1, 0)])
    monkeypatch.setattr(gpd, "read_file", lambda _source: frame)

    with pytest.raises(ValueError, match="must declare"):
        load_contact_regions(source, id_column="patch")


def test_vector_input_rejects_invalid_geometry(tmp_path):
    source = tmp_path / "invalid.geojson"
    gpd.GeoDataFrame(
        {"patch": ["A", "B"]},
        geometry=[
            Polygon([(0, 0), (1, 1), (1, 0), (0, 1), (0, 0)]),
            Polygon([(2, 0), (3, 0), (3, 1), (2, 1)]),
        ],
        crs="EPSG:4326",
    ).to_file(source, driver="GeoJSON")
    with pytest.raises(ValueError, match="invalid geometries"):
        load_contact_regions(source, id_column="patch", centroid_crs="EPSG:3857")


def test_vector_input_without_geo_extra_has_actionable_error(tmp_path, monkeypatch):
    source = tmp_path / "regions.geojson"
    source.write_text("{}", encoding="utf-8")
    real_import = builtins.__import__

    def import_without_geo(name, *args, **kwargs):
        if name in {"geopandas", "pyproj"}:
            raise ImportError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_without_geo)
    with pytest.raises(RuntimeError, match=r"patchsim\[geo\]"):
        load_contact_regions(source, id_column="patch")


def test_generate_contacts_writes_stable_csv_and_validation_report(tmp_path, regions):
    source = tmp_path / "regions.csv"
    output = tmp_path / "networks" / "contacts.csv"
    regions.to_csv(source, index=False)

    output_path, report_path, report = generate_contacts(
        source,
        output,
        id_column="id",
        population_column="population",
        kernel="gravity",
        decay=2.0,
        min_distance_km=0.001,
        normalize="row",
        self_share=0.75,
    )

    edge_frame = pd.read_csv(output_path)
    assert edge_frame[["source", "target"]].values.tolist() == [
        ["A", "A"],
        ["A", "B"],
        ["B", "A"],
        ["B", "B"],
    ]
    assert edge_frame["day"].tolist() == [0, 0, 0, 0]
    assert report["normalization"]["final_weight_unit"] == "dimensionless"
    assert report["distance"]["unit"] == "km"
    assert report["kernel"]["raw_weight_unit"] == "scale * population_i * population_j / km**decay"
    assert report["csv_sha256"] == hashlib.sha256(output_path.read_bytes()).hexdigest()
    assert json.loads(report_path.read_text(encoding="utf-8")) == report


def test_generate_contacts_refuses_overwrite_and_source_alias(tmp_path, regions):
    source = tmp_path / "regions.csv"
    output = tmp_path / "contacts.csv"
    regions.to_csv(source, index=False)
    kwargs = {
        "id_column": "id",
        "kernel": "distance",
        "decay": 2.0,
        "scale": 1.0,
        "min_distance_km": 0.001,
        "normalize": "none",
        "self_weight": 1.0,
    }

    generate_contacts(source, output, **kwargs)
    with pytest.raises(FileExistsError):
        generate_contacts(source, output, **kwargs)
    generate_contacts(source, output, force=True, **kwargs)
    with pytest.raises(ValueError, match="must not overwrite source"):
        generate_contacts(source, source, force=True, **kwargs)


def test_generated_csv_loads_through_existing_simulation_setup(tmp_path, regions):
    source = tmp_path / "regions.csv"
    network = tmp_path / "network.csv"
    patch_file = tmp_path / "patches.csv"
    seed_file = tmp_path / "seeds.csv"
    regions.assign(id=["001", "NA"]).to_csv(source, index=False)
    patch_file.write_text("patch,population\n001,100\nNA,200\n", encoding="utf-8")
    seed_file.write_text("patch,S,I,R\n001,99,1,0\nNA,200,0,0\n", encoding="utf-8")
    generate_contacts(
        source,
        network,
        id_column="id",
        kernel="distance",
        decay=2.0,
        min_distance_km=0.001,
        normalize="row",
        self_share=0.8,
    )

    config = {
        "PatchFile": str(patch_file),
        "SeedFile": str(seed_file),
        "NetworkFile": str(network),
        "OutputDir": str(tmp_path / "output"),
        "TMax": 2,
        "compartments": ["S", "I", "R"],
        "Parameters": {"beta": 0.08, "gamma": 0.1},
        "Transitions": {"S -> I": "beta", "I -> R": "gamma * I"},
    }
    net, _state, patches, _count = setup_simulation(config)
    assert patches == ["001", "NA"]
    np.testing.assert_allclose(net.network, [[0.8, 0.2], [0.2, 0.8]])

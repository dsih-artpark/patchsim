import json
import subprocess
from copy import deepcopy
from importlib.metadata import PackageNotFoundError
from pathlib import Path

import pandas as pd
import pytest
import yaml

import patchsim
import patchsim.sensitivity as sensitivity_module
from patchsim.core.simulation import load_config
from patchsim.sensitivity import get_sensitivity_plan, run_sensitivity


def _write_project(tmp_path: Path, *, base_samples: int = 32) -> Path:
    pd.DataFrame({"patch": ["A"], "population": [100.0]}).to_csv(tmp_path / "patches.csv", index=False)
    pd.DataFrame({"patch": ["A"], "S": [0.0], "I": [100.0], "R": [0.0]}).to_csv(tmp_path / "seeds.csv", index=False)
    config = {
        "ModelName": "sensitivity-test",
        "PatchFile": "patches.csv",
        "SeedFile": "seeds.csv",
        "NetworkFile": None,
        "OutputDir": "output",
        "TMax": 5,
        "Solver": "discrete",
        "TimeStep": 1.0,
        "compartments": ["S", "I", "R"],
        "Parameters": {"theta": 0.1, "unused": 1.0},
        "Transitions": {"I -> R": "theta * I"},
        "Sensitivity": {
            "Name": "theta-study",
            "Method": "sobol",
            "BaseSamples": base_samples,
            "Seed": 42,
            "Parameters": {"theta": [0.01, 0.2], "unused": [1.0, 2.0]},
            "Metrics": {
                "final_r": {"Columns": ["R_0"], "Reduce": "final"},
                "peak_r": {"Columns": ["R_0"], "Reduce": "max"},
            },
        },
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return path


def test_simulate_is_side_effect_free_and_accepts_global_overrides(tmp_path):
    config_path = _write_project(tmp_path)
    config = load_config(str(config_path))
    original = deepcopy(config)

    baseline = patchsim.simulate(config)
    changed = patchsim.simulate(config, parameter_overrides={"theta": 0.2})

    assert config == original
    assert not Path(config["OutputDir"]).exists()
    assert baseline.columns.tolist() == ["time", "S_0", "I_0", "R_0"]
    assert changed["R_0"].iloc[-1] > baseline["R_0"].iloc[-1]


def test_simulate_rejects_unknown_non_finite_and_patch_specific_overrides(tmp_path):
    config = load_config(str(_write_project(tmp_path)))

    with pytest.raises(ValueError, match="Unknown global"):
        patchsim.simulate(config, parameter_overrides={"missing": 1.0})
    with pytest.raises(ValueError, match="finite real"):
        patchsim.simulate(config, parameter_overrides={"theta": float("nan")})

    config["PatchParameters"] = [{"patch": "A", "parameters": {"theta": 0.05}}]
    with pytest.raises(ValueError, match="PatchParameters"):
        patchsim.simulate(config, parameter_overrides={"theta": 0.1})


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("Name", "../escape", "safe path"),
        ("BaseSamples", 1, "at least 2"),
        ("BaseSamples", 12, "power of two"),
        ("Seed", -1, "non-negative"),
        ("Parameters", {"theta": [0.2, 0.1]}, "lower < upper"),
    ],
)
def test_sensitivity_plan_rejects_invalid_contract(tmp_path, field, value, message):
    config = load_config(str(_write_project(tmp_path)))
    config["Sensitivity"][field] = value

    with pytest.raises(ValueError, match=message):
        get_sensitivity_plan(config, ["S_0", "I_0", "R_0"])


def test_sensitivity_plan_rejects_ambiguous_or_unknown_columns(tmp_path):
    config = load_config(str(_write_project(tmp_path)))
    config["PatchParameters"] = [{"patch": "A", "parameters": {"theta": 0.05}}]
    with pytest.raises(ValueError, match="PatchParameters"):
        get_sensitivity_plan(config, ["S_0", "I_0", "R_0"])

    config["PatchParameters"] = []
    config["Sensitivity"]["Metrics"]["final_r"]["Columns"] = ["R_9"]
    with pytest.raises(ValueError, match="unknown output"):
        get_sensitivity_plan(config, ["S_0", "I_0", "R_0"])


def test_sobol_study_writes_complete_deterministic_artifacts(tmp_path):
    config_path = _write_project(tmp_path)

    summary = run_sensitivity(config_path)

    assert summary["reused"] is False
    assert summary["planned_evaluations"] == summary["completed_evaluations"] == 128
    samples = pd.read_csv(summary["samples_path"])
    responses = pd.read_csv(summary["responses_path"])
    indices = pd.read_csv(summary["indices_path"])
    manifest = json.loads(Path(summary["manifest_path"]).read_text(encoding="utf-8"))

    assert samples.columns.tolist() == ["sample_id", "theta", "unused"]
    assert responses.columns.tolist() == ["sample_id", "final_r", "peak_r"]
    assert len(samples) == len(responses) == 128
    assert set(indices["metric"]) == {"final_r", "peak_r"}
    final = indices.set_index(["metric", "parameter"]).loc["final_r"]
    assert final.loc["theta", "ST"] > 0.95
    assert abs(final.loc["unused", "ST"]) < 1e-12
    assert manifest["evaluation_count"] == 128
    assert manifest["request"]["method"]["calc_second_order"] is False
    assert manifest["request"]["method"]["num_resamples"] == 100
    assert manifest["request"]["normalized_config"]["Parameters"] == {
        "theta": 0.1,
        "unused": 1.0,
    }
    assert set(manifest["artifacts"]) == {"samples.csv", "responses.csv", "indices.csv"}


def test_identical_existing_study_is_reused_without_solves(tmp_path, monkeypatch):
    config_path = _write_project(tmp_path, base_samples=8)
    first = run_sensitivity(config_path)

    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("reused study must not run the model")

    monkeypatch.setattr(sensitivity_module, "_simulate_prepared", fail_if_called)
    second = run_sensitivity(config_path)

    assert second["reused"] is True
    assert second["completed_evaluations"] == 0
    assert second["samples_path"] == first["samples_path"]


@pytest.mark.parametrize("change", ["config", "artifact", "manifest"])
def test_existing_study_rejects_changed_inputs_or_artifacts(tmp_path, change):
    config_path = _write_project(tmp_path, base_samples=8)
    summary = run_sensitivity(config_path)

    if change == "config":
        config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        config["Sensitivity"]["Seed"] = 43
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        message = "different study"
    else:
        if change == "artifact":
            Path(summary["samples_path"]).write_text("modified\n", encoding="utf-8")
            message = "missing or modified"
        else:
            manifest_path = Path(summary["manifest_path"])
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["evaluation_count"] = 1
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            message = "different study"

    with pytest.raises(FileExistsError, match=message):
        run_sensitivity(config_path)


def test_sensitivity_rejects_output_symlink_escape(tmp_path):
    config_path = _write_project(tmp_path, base_samples=8)
    outside = tmp_path / "outside"
    outside.mkdir()
    output = tmp_path / "output"
    output.mkdir()
    (output / "sensitivity").symlink_to(outside, target_is_directory=True)

    with pytest.raises(ValueError, match="escapes OutputDir"):
        run_sensitivity(config_path)
    assert not list(outside.iterdir())


def test_failed_or_constant_study_publishes_nothing(tmp_path, monkeypatch):
    config_path = _write_project(tmp_path, base_samples=8)

    def fail(*_args, **_kwargs):
        raise ValueError("solver failed")

    monkeypatch.setattr(sensitivity_module, "_simulate_prepared", fail)
    with pytest.raises(RuntimeError, match=r"evaluation 0.*theta.*solver failed"):
        run_sensitivity(config_path)
    assert not (tmp_path / "output").exists()

    monkeypatch.undo()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["Sensitivity"]["Metrics"] = {"susceptible": {"Columns": ["S_0"], "Reduce": "final"}}
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    with pytest.raises(ValueError, match="constant at 0"):
        run_sensitivity(config_path)
    assert not (tmp_path / "output").exists()


def test_multicolumn_metric_does_not_skip_non_finite_values(tmp_path, monkeypatch):
    config_path = _write_project(tmp_path, base_samples=8)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["Sensitivity"]["Metrics"] = {"combined": {"Columns": ["I_0", "R_0"], "Reduce": "final"}}
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    def non_finite(*_args, **_kwargs):
        return pd.DataFrame({"time": [0.0], "S_0": [0.0], "I_0": [1.0], "R_0": [float("nan")]})

    monkeypatch.setattr(sensitivity_module, "_simulate_prepared", non_finite)
    with pytest.raises(RuntimeError, match="non-finite metrics.*combined"):
        run_sensitivity(config_path)
    assert not (tmp_path / "output").exists()


def test_recomputation_is_byte_identical(tmp_path):
    config_path = _write_project(tmp_path, base_samples=8)
    first = run_sensitivity(config_path)
    first_bytes = {
        name: (Path(first["output_dir"]) / name).read_bytes()
        for name in (*sensitivity_module._ARTIFACT_NAMES, "manifest.json")
    }
    target = Path(first["output_dir"])
    for path in target.iterdir():
        path.unlink()
    target.rmdir()

    second = run_sensitivity(config_path)

    assert second["reused"] is False
    for name, expected in first_bytes.items():
        assert (Path(second["output_dir"]) / name).read_bytes() == expected


def test_validate_and_sensitivity_cli_json(tmp_path):
    config_path = _write_project(tmp_path, base_samples=8)
    validation = subprocess.run(
        ["patchsim", "validate", "-c", str(config_path), "--json"],
        capture_output=True,
        text=True,
        timeout=300,
    )
    result = subprocess.run(
        ["patchsim", "sensitivity", "-c", str(config_path), "--json"],
        capture_output=True,
        text=True,
        timeout=300,
    )

    assert validation.returncode == 0, validation.stderr
    assert json.loads(validation.stdout)["sensitivity"]["evaluation_count"] == 32
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["completed_evaluations"] == 32
    assert Path(payload["indices_path"]).is_file()
    assert "Planned model evaluations: 32" in result.stderr


def test_salib_dependency_errors_are_actionable(monkeypatch):
    def missing(_package):
        raise PackageNotFoundError

    monkeypatch.setattr(sensitivity_module, "version", missing)
    with pytest.raises(RuntimeError, match=r"patchsim\[analysis\]"):
        sensitivity_module._load_salib()

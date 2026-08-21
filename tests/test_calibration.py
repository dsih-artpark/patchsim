import json
import math
import shutil
import subprocess
from datetime import date
from pathlib import Path

import pandas as pd
import pytest
import yaml

import patchsim
import patchsim.calibration as calibration_module
from patchsim.calibration import get_calibration_plan, run_calibration
from patchsim.core.simulation import load_config, setup_simulation


def _write_project(tmp_path: Path, *, observations: pd.DataFrame | None = None) -> Path:
    pd.DataFrame({"patch": ["north"], "population": [100.0]}).to_csv(tmp_path / "patches.csv", index=False)
    pd.DataFrame({"patch": ["north"], "S": [92.0], "I": [8.0], "R": [0.0]}).to_csv(tmp_path / "seeds.csv", index=False)
    if observations is None:
        times = [1.0, 2.0, 3.0, 4.0, 5.0]
        observations = pd.DataFrame(
            {
                "time": times,
                "observable": ["removed"] * len(times),
                "value": [20.0 * (1.0 - 0.8**time) for time in times],
            }
        )
    observations.to_csv(tmp_path / "observations.csv", index=False)

    config = {
        "ModelName": "calibration-test",
        "PatchFile": "patches.csv",
        "SeedFile": "seeds.csv",
        "NetworkFile": None,
        "OutputDir": "output",
        "TMax": 6,
        "Solver": "discrete",
        "TimeStep": 1.0,
        "compartments": ["S", "I", "R"],
        "Parameters": {"theta": 0.08},
        "Transitions": {"I -> R": "theta * I"},
        "Calibration": {
            "Name": "theta-i0",
            "Method": "least_squares",
            "Observations": "observations.csv",
            "MaxEvaluations": 200,
            "Observables": {"removed": {"Columns": ["R_0"], "Scale": 1.0}},
            "Parameters": {"theta": [0.01, 0.5]},
            "InitialConditions": [
                {
                    "Patch": "north",
                    "Remainder": "S",
                    "Fit": {"I": [1.0, 50.0]},
                }
            ],
            "Starts": [
                {
                    "Parameters": {"theta": 0.4},
                    "InitialConditions": [
                        {"Patch": "north", "Values": {"I": 40.0}},
                    ],
                }
            ],
        },
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return path


def _plan(config_path: Path):
    config = load_config(str(config_path))
    net, y0, _patches, _num_patches = setup_simulation(config)
    return get_calibration_plan(config, net, y0, required=True)


def test_calibration_plan_aligns_observations_and_declares_fit_variables(tmp_path):
    plan = _plan(_write_project(tmp_path))

    assert plan.name == "theta-i0"
    assert plan.n == 5
    assert plan.p == 2
    assert plan.start_count == 2
    assert plan.max_forward_simulations == 400
    assert plan.warnings == ()
    assert [variable.name for variable in plan.variables] == ["theta", "I"]
    assert plan.observations["grid_index"].tolist() == [1, 2, 3, 4, 5]


def test_time_mismatch_warns_during_planning_and_blocks_fitting(tmp_path):
    observations = pd.DataFrame({"time": [1.0, 1.5, 2.0], "observable": ["removed"] * 3, "value": [4.0, 5.0, 7.2]})
    config_path = _write_project(tmp_path, observations=observations)

    plan = _plan(config_path)

    assert len(plan.warnings) == 1
    assert "row 3" in plan.warnings[0]
    with pytest.raises(ValueError, match="unmatched observation times"):
        run_calibration(config_path)
    assert not (tmp_path / "output").exists()


def test_observation_duplicates_and_unsafe_initial_bounds_are_rejected(tmp_path):
    duplicate = pd.DataFrame({"time": [1.0, 1.0], "observable": ["removed", "removed"], "value": [4.0, 4.0]})
    with pytest.raises(ValueError, match="duplicate"):
        _plan(_write_project(tmp_path, observations=duplicate))

    config_path = _write_project(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["Calibration"]["InitialConditions"][0]["Fit"]["I"] = [1.0, 101.0]
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    with pytest.raises(ValueError, match="non-negative remainder"):
        _plan(config_path)


def test_calibration_recovers_parameter_and_initial_state_and_reuses_artifacts(tmp_path, monkeypatch):
    config_path = _write_project(tmp_path)

    summary = run_calibration(config_path)

    assert summary["reused"] is False
    assert summary["n"] == 5
    assert summary["p"] == 2
    assert summary["start_count"] == 2
    estimates = pd.read_csv(summary["estimates_path"]).set_index("name")
    assert estimates.loc["theta", "value"] == pytest.approx(0.2, abs=1e-6)
    assert estimates.loc["I", "value"] == pytest.approx(20.0, abs=1e-5)
    fitted_seeds = pd.read_csv(summary["fitted_seeds_path"])
    assert fitted_seeds.loc[0, "S"] == pytest.approx(80.0, abs=1e-5)
    assert fitted_seeds.loc[0, "I"] == pytest.approx(20.0, abs=1e-5)
    attempts = pd.read_csv(summary["attempts_path"])
    assert len(attempts) == 2
    assert attempts["success"].all()
    residuals = pd.read_csv(summary["residuals_path"])
    assert residuals["residual"].abs().max() < 1e-5
    manifest = json.loads(Path(summary["manifest_path"]).read_text(encoding="utf-8"))
    assert manifest["n"] == 5
    assert manifest["p"] == 2
    assert manifest["max_forward_simulations"] == 400
    assert manifest["request"]["method"]["max_forward_simulations_per_start"] == 200
    assert manifest["selected_start"] in {0, 1}
    assert manifest["selected_attempt"]["success"] is True
    assert manifest["selected_attempt"]["status"] > 0
    assert manifest["residual_summary"]["rmse"] < 1e-5
    assert manifest["residual_summary"]["standardized_rmse"] < 1e-5
    assert manifest["residual_summary"]["by_observable"]["removed"]["n"] == 5

    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("reused calibration must not run the model")

    monkeypatch.setattr(calibration_module, "_simulate_prepared", fail_if_called)
    reused = patchsim.run_calibration(config_path)
    assert reused["reused"] is True
    assert reused["forward_simulations"] == 0


def test_calibration_continues_after_one_start_fails(tmp_path, monkeypatch):
    config_path = _write_project(tmp_path)
    simulate = calibration_module._simulate_prepared

    def fail_high_start(config, net, y0):
        if net.base_model.parameters["theta"] > 0.3:
            raise ValueError("deliberate high-start failure")
        return simulate(config, net, y0)

    monkeypatch.setattr(calibration_module, "_simulate_prepared", fail_high_start)
    summary = run_calibration(config_path)
    attempts = pd.read_csv(summary["attempts_path"])

    assert attempts["success"].tolist() == [True, False]
    assert "deliberate high-start failure" in attempts.loc[1, "message"]


def test_all_failed_starts_publish_nothing(tmp_path, monkeypatch):
    config_path = _write_project(tmp_path)

    def fail(*_args, **_kwargs):
        raise ValueError("solver failed")

    monkeypatch.setattr(calibration_module, "_simulate_prepared", fail)
    with pytest.raises(RuntimeError, match="no starting point.*solver failed"):
        run_calibration(config_path)
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("change", ["artifact", "manifest", "config"])
def test_existing_calibration_rejects_changed_request_or_artifacts(tmp_path, change):
    config_path = _write_project(tmp_path)
    summary = run_calibration(config_path)
    if change == "artifact":
        Path(summary["estimates_path"]).write_text("modified\n", encoding="utf-8")
        message = "missing or modified"
    elif change == "manifest":
        manifest_path = Path(summary["manifest_path"])
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["n"] = 1
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        message = "different study"
    else:
        config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        config["Calibration"]["MaxEvaluations"] = 201
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        message = "different study"

    with pytest.raises(FileExistsError, match=message):
        run_calibration(config_path)


def test_grouped_initial_condition_uses_semantic_identifiers_and_preserves_population(tmp_path):
    pd.DataFrame({"patch": ["north"], "population": [100.0]}).to_csv(tmp_path / "patches.csv", index=False)
    pd.DataFrame({"patch": ["north", "north"], "group": ["adult", "child"], "population": [60.0, 40.0]}).to_csv(
        tmp_path / "groups.csv", index=False
    )
    pd.DataFrame(
        {
            "patch": ["north", "north"],
            "group": ["adult", "child"],
            "S": [55.0, 40.0],
            "I": [5.0, 0.0],
            "R": [0.0, 0.0],
        }
    ).to_csv(tmp_path / "seeds.csv", index=False)
    pd.DataFrame(
        {
            "focal_group": ["adult", "adult", "child", "child"],
            "contributor_group": ["adult", "child", "adult", "child"],
            "weight": [1.0, 0.0, 0.0, 1.0],
        }
    ).to_csv(tmp_path / "interactions.csv", index=False)
    pd.DataFrame({"time": [0.0, 1.0, 2.0], "observable": ["adult_i"] * 3, "value": [5.0, 4.5, 4.05]}).to_csv(
        tmp_path / "observations.csv", index=False
    )
    config = {
        "ModelName": "grouped-calibration",
        "PatchFile": "patches.csv",
        "GroupFile": "groups.csv",
        "InteractionFile": "interactions.csv",
        "InteractionUnits": "contacts/person/day",
        "SeedFile": "seeds.csv",
        "NetworkFile": None,
        "OutputDir": "output",
        "TMax": 3,
        "Solver": "discrete",
        "TimeStep": 1.0,
        "compartments": ["S", "I", "R"],
        "Parameters": {"gamma": 0.1},
        "Transitions": {"I -> R": "gamma * I"},
        "Calibration": {
            "Name": "adult-i0",
            "Method": "least_squares",
            "Observations": "observations.csv",
            "MaxEvaluations": 20,
            "Observables": {"adult_i": {"Columns": ["I_0_0"], "Scale": 1.0}},
            "InitialConditions": [{"Patch": "north", "Group": "adult", "Remainder": "S", "Fit": {"I": [1.0, 10.0]}}],
            "Starts": [
                {
                    "Parameters": {},
                    "InitialConditions": [
                        {"Patch": "north", "Group": "adult", "Values": {"I": 8.0}},
                    ],
                }
            ],
        },
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    loaded = load_config(str(config_path))
    net, y0, _patches, _num_patches = setup_simulation(loaded)
    plan = get_calibration_plan(loaded, net, y0, required=True)
    assert plan is not None

    trial = calibration_module._apply_vector(net, y0, plan.variables, pd.Series([10.0]).to_numpy())

    assert plan.variables[0].state_key == "I_0_0"
    assert trial["S_0_0"] == 50.0
    assert trial["I_0_0"] == 10.0
    assert trial["S_0_0"] + trial["I_0_0"] + trial["R_0_0"] == 60.0
    assert trial["S_0_1"] == 40.0


def test_validate_and_calibrate_cli_json(tmp_path):
    config_path = _write_project(tmp_path)
    validation = subprocess.run(
        ["patchsim", "validate", "-c", str(config_path), "--json"],
        capture_output=True,
        text=True,
        timeout=300,
    )
    result = subprocess.run(
        ["patchsim", "calibrate", "-c", str(config_path), "--json"],
        capture_output=True,
        text=True,
        timeout=300,
    )

    assert validation.returncode == 0, validation.stderr
    validation_payload = json.loads(validation.stdout)
    assert validation_payload["calibration"] == {
        "name": "theta-i0",
        "method": "least_squares",
        "n": 5,
        "p": 2,
        "start_count": 2,
        "max_forward_simulations": 400,
        "warnings": [],
    }
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["n"] == 5
    assert Path(payload["manifest_path"]).is_file()
    assert "Maximum forward simulations: 400" in result.stderr


def test_observation_and_parameter_contract_errors_are_explicit(tmp_path):
    config_path = _write_project(tmp_path)
    observations = pd.read_csv(tmp_path / "observations.csv")
    observations["notes"] = "ignored?"
    observations.to_csv(tmp_path / "observations.csv", index=False)
    with pytest.raises(ValueError, match="exactly these columns"):
        _plan(config_path)

    observations.drop(columns="notes").to_csv(tmp_path / "observations.csv", index=False)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["Calibration"]["Observables"]["removed"]["Columns"] = ["R_9"]
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown output columns.*R_9"):
        _plan(config_path)

    config["Calibration"]["Observables"]["removed"]["Columns"] = ["R_0"]
    config["PatchParameters"] = [{"patch": "north", "parameters": {"theta": 0.1}}]
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    with pytest.raises(ValueError, match="also set in PatchParameters"):
        _plan(config_path)


def test_every_configured_observable_must_have_data(tmp_path):
    config_path = _write_project(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["Calibration"]["Observables"]["infectious"] = {"Columns": ["I_0"], "Scale": 1.0}
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    with pytest.raises(ValueError, match="no observation rows.*infectious"):
        _plan(config_path)


def test_start_sample_size_and_budget_contracts(tmp_path):
    config_path = _write_project(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    del config["Calibration"]["Starts"][0]["InitialConditions"]
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    with pytest.raises(ValueError, match="every fitted initial condition"):
        _plan(config_path)

    config_path = _write_project(tmp_path)
    observations = pd.read_csv(tmp_path / "observations.csv").iloc[:1]
    observations.to_csv(tmp_path / "observations.csv", index=False)
    with pytest.raises(ValueError, match="underdetermined: n=1.*p=2"):
        _plan(config_path)

    observations = pd.DataFrame({"time": [1.0, 2.0], "observable": ["removed"] * 2, "value": [4.0, 7.2]})
    observations.to_csv(tmp_path / "observations.csv", index=False)
    plan = _plan(config_path)
    assert "n == p" in plan.warnings[0]

    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["Calibration"]["MaxEvaluations"] = 3
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    with pytest.raises(RuntimeError, match="budget 3 exhausted"):
        run_calibration(config_path)
    assert not (tmp_path / "output").exists()


def test_rank_deficient_fit_is_reported_without_identifiability_claim(tmp_path):
    config_path = _write_project(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["Parameters"]["unused"] = 1.0
    config["Calibration"]["Parameters"]["unused"] = [0.5, 1.5]
    config["Calibration"]["Starts"][0]["Parameters"]["unused"] = 1.2
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    summary = run_calibration(config_path)
    manifest = json.loads(Path(summary["manifest_path"]).read_text(encoding="utf-8"))

    assert manifest["jacobian"]["rank"] < 3
    assert manifest["jacobian"]["condition_number"] is None


def test_ode_calibration_recovers_parameter_with_fixed_initial_state(tmp_path):
    config_path = _write_project(tmp_path)
    pd.DataFrame({"patch": ["north"], "S": [80.0], "I": [20.0], "R": [0.0]}).to_csv(tmp_path / "seeds.csv", index=False)
    times = [1.0, 2.0, 3.0, 4.0, 5.0]
    pd.DataFrame(
        {
            "time": times,
            "observable": ["removed"] * len(times),
            "value": [20.0 * (1.0 - math.exp(-0.2 * time)) for time in times],
        }
    ).to_csv(tmp_path / "observations.csv", index=False)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["Solver"] = "ode"
    config["Parameters"]["theta"] = 0.05
    del config["Calibration"]["InitialConditions"]
    config["Calibration"]["Starts"] = [{"Parameters": {"theta": 0.4}}]
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    summary = run_calibration(config_path)
    estimates = pd.read_csv(summary["estimates_path"]).set_index("name")

    assert summary["p"] == 1
    assert estimates.loc["theta", "value"] == pytest.approx(0.2, abs=1e-6)


def test_clean_recalibration_produces_identical_artifacts(tmp_path):
    config_path = _write_project(tmp_path)
    first = run_calibration(config_path)
    artifact_names = ("estimates.csv", "fitted-seeds.csv", "attempts.csv", "residuals.csv", "manifest.json")
    first_bytes = {name: (Path(first["output_dir"]) / name).read_bytes() for name in artifact_names}

    shutil.rmtree(first["output_dir"])
    second = run_calibration(config_path)

    assert {name: (Path(second["output_dir"]) / name).read_bytes() for name in artifact_names} == first_bytes


def test_yaml_date_is_serialized_in_calibration_provenance(tmp_path):
    config_path = _write_project(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["StartDate"] = date(2020, 1, 1)
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    summary = run_calibration(config_path)
    manifest = json.loads(Path(summary["manifest_path"]).read_text(encoding="utf-8"))

    assert manifest["request"]["normalized_config"]["StartDate"] == "2020-01-01"


def test_interrupted_publication_removes_temporary_study(tmp_path, monkeypatch):
    config_path = _write_project(tmp_path)

    def interrupt(*_args, **_kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(pd.DataFrame, "to_csv", interrupt)
    with pytest.raises(KeyboardInterrupt):
        run_calibration(config_path)

    calibration_root = tmp_path / "output" / "calibration"
    assert not list(calibration_root.iterdir())


def test_calibration_rejects_output_symlink_escape(tmp_path):
    config_path = _write_project(tmp_path)
    outside = tmp_path / "outside"
    outside.mkdir()
    output = tmp_path / "output"
    output.mkdir()
    (output / "calibration").symlink_to(outside, target_is_directory=True)

    with pytest.raises(ValueError, match="escapes OutputDir"):
        run_calibration(config_path)
    assert not list(outside.iterdir())

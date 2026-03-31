import subprocess


def test_patchsim_import():
    import patchsim

    assert hasattr(patchsim, "__version__")


def test_patchsim_sdk_exports():
    import patchsim

    assert hasattr(patchsim, "CompartmentalModel")
    assert hasattr(patchsim, "NetworkModel")
    assert hasattr(patchsim, "load_config")
    assert hasattr(patchsim, "setup_simulation")
    assert hasattr(patchsim, "run_simulation")
    assert hasattr(patchsim, "plot_patch_subplots")


def test_core_imports():
    pass


def test_cli_help():
    result = subprocess.run(
        ["patchsim", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "patchsim" in result.stdout.lower()


def test_cli_version():
    result = subprocess.run(
        ["patchsim", "--version"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "patchsim" in result.stdout.lower()


def test_cli_list_models():
    result = subprocess.run(
        ["patchsim", "list-models"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    output = (result.stdout + result.stderr).lower()
    assert "ka_fmd_sirsv_discrete" in output


def test_cli_init_scaffold(tmp_path):
    project_dir = tmp_path / "demo-project"
    result = subprocess.run(
        ["patchsim", "init", str(project_dir)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0

    assert (project_dir / "config.yaml").exists()
    assert (project_dir / "data" / "patch" / "patch-population.csv").exists()
    assert (project_dir / "data" / "networks" / "network-static.csv").exists()
    assert (project_dir / "data" / "seeds" / "seed-initial.csv").exists()
    assert (project_dir / "output").exists()


def test_yaml_model_loading():
    from patchsim.core.simulation import load_config

    config = load_config("configs/sample-sir-ode.yaml")
    assert "compartments" in config
    assert "Parameters" in config
    assert "Transitions" in config


def test_minimal_sir_simulation():
    import numpy as np

    from patchsim.core.model import CompartmentalModel, NetworkModel
    from patchsim.core.model_runner import Model

    base_model = CompartmentalModel(
        compartments=["S", "I", "R"],
        parameters={"beta": 0.3, "gamma": 0.1},
        transitions=[{"transition": "S->I", "rate": "beta"}, {"transition": "I->R", "rate": "gamma * I"}],
    )
    # Create single-patch network model
    network_model = NetworkModel(base_model=base_model, num_patches=1, network_matrix=[[1.0]])
    # Initial conditions
    y0 = {"S_0": 999.0, "I_0": 1.0, "R_0": 0.0}
    t_range = np.linspace(0, 10, 11)

    # Use Model runner
    model = Model(network_model, compartments=["S", "I", "R"])
    results = model.solve(y0, t_range)

    # Smoke-level assertions
    assert results is not None
    assert "S_0" in results
    assert "I_0" in results
    assert "R_0" in results

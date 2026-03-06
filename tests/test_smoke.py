import subprocess
import sys

def test_patchsim_import():
    import patchsim
    assert hasattr(patchsim, "__version__")


def test_core_imports():
    from patchsim.core.model import CompartmentalModel, NetworkModel
    from patchsim.core.simulation import load_config, setup_simulation, run_simulation
    from patchsim.core.network import Network
    
def test_cli_help():
    result = subprocess.run(
        ["patchsim", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "patchsim" in result.stdout.lower()
    
    
def test_yaml_model_loading():
    from patchsim.core.simulation import load_config
    config = load_config("configs/sample-sir-ode.yaml")
    assert "compartments" in config
    assert "Parameters" in config
    assert "Transitions" in config


def test_minimal_sir_simulation():
    from patchsim.core.model import CompartmentalModel, NetworkModel
    from patchsim.core.model_runner import Model
    import numpy as np
    
    base_model = CompartmentalModel(
        compartments=["S", "I", "R"],
        parameters={"beta": 0.3, "gamma": 0.1},
        transitions=[
            {"from": "S", "to": "I", "rate": "beta"},
            {"from": "I", "to": "R", "rate": "gamma * I"}
        ]
    )
    # Create single-patch network model
    network_model = NetworkModel(
        base_model=base_model,
        num_patches=1,
        network_matrix=[[1.0]]
    )
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

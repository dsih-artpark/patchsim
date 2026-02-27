import subprocess
import sys
#from patchsim.core.model import 
#from patchsim.utils.loader import 

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
    
#def test_minimal_sir_simulation():
    # model = Model(
    #     compartments=["S", "I", "R"],
    #     parameters={"beta": 0.3, "gamma": 0.1},
    #     transitions={
    #         "S -> I": "beta * S * I / N",
    #         "I -> R": "gamma * I",
    #     },
    #     population=1000,
    #     initial_conditions={"S": 999, "I": 1, "R": 0},
    # )

    #results = model.simulate(t_max=10)

    # Smoke-level assertions
    #assert results is not None
    #assert len(results) == 11  # includes t=0
    #assert "S" in results.columns
    #assert "I" in results.columns
    #assert "R" in results.columns

#def test_yaml_model_loading():
    #config = load_model_config("examples/sir.yaml")

    #assert "compartments" in config
    #assert "parameters" in config
    #assert "transitions" in config
import subprocess
import sys
#from patchsim.core.model import Model
#from patchsim.utils.loader import load_model_config

def test_patchsim_import():
    import patchsim
    assert hasattr(patchsim, "__version__")


def test_core_imports():
    from patchsim.core.model import CompartmentalModelodel, NetworkModel
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
    

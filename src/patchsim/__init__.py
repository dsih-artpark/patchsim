"""PatchSim public SDK interface."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from patchsim.core.model import CompartmentalModel, NetworkModel
from patchsim.core.simulation import load_config, run_simulation, setup_simulation
from patchsim.utils.viz import plot_patch_subplots

try:
    __version__ = version("patchsim")
except PackageNotFoundError:  # pragma: no cover - local editable fallback
    __version__ = "0.0.0"

__all__ = [
    "__version__",
    "CompartmentalModel",
    "NetworkModel",
    "load_config",
    "setup_simulation",
    "run_simulation",
    "plot_patch_subplots",
]

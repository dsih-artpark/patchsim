"""PatchSim public SDK interface."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from patchsim.calibration import run_calibration
from patchsim.core.model import CompartmentalModel, NetworkModel
from patchsim.core.simulation import load_config, run_simulation, setup_simulation, simulate
from patchsim.utils.viz import plot_patch_subplots

try:
    __version__ = version("patchsim")
except PackageNotFoundError:  # pragma: no cover - local editable fallback
    __version__ = "0.1.0"

__all__ = [
    "CompartmentalModel",
    "NetworkModel",
    "__version__",
    "load_config",
    "plot_patch_subplots",
    "run_calibration",
    "run_simulation",
    "setup_simulation",
    "simulate",
]

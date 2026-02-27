
from __future__ import annotations

try:
	from importlib.metadata import version as _pkg_version
except Exception:  # pragma: no cover - fallback for old environments
	_pkg_version = None  # type: ignore[assignment]

try:
	__version__ = _pkg_version("patchsim") if _pkg_version else "0.0.0"
except Exception:  # pragma: no cover - package metadata missing
	__version__ = "0.0.0"

__all__ = ["__version__"]

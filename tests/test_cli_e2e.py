import json
import subprocess
from pathlib import Path

import pytest

TEMPLATES = ["sir", "sis", "sirs", "seir"]


def _patchsim(*args, cwd=None):
    return subprocess.run(["patchsim", *args], capture_output=True, text=True, cwd=cwd)


@pytest.mark.parametrize("template", TEMPLATES)
def test_init_then_run_produces_existing_artifacts(tmp_path, template):
    """Every built-in template must scaffold AND run end-to-end, and the artifact
    paths reported by `run --json` must actually exist on disk."""
    proj = tmp_path / f"{template}-proj"
    init = _patchsim("init", str(proj), "--template", template)
    assert init.returncode == 0, init.stderr

    run = _patchsim("run", "-c", "config.yaml", "--json", cwd=str(proj))
    assert run.returncode == 0, f"`run` failed for template {template}:\n{run.stderr}"

    summary = json.loads(run.stdout)
    assert Path(summary["csv_path"]).is_file(), f"reported csv missing: {summary['csv_path']}"
    assert Path(summary["plot_path"]).is_file(), f"reported plot missing: {summary['plot_path']}"


def test_list_models_output_is_clean():
    """Human-facing `list-models` output matches the documented format and does not
    leak logging decoration (timestamps / source-file references)."""
    res = _patchsim("list-models")
    assert res.returncode == 0
    assert "- sir (yaml-template)" in res.stdout
    assert "- sis (yaml-template)" in res.stdout
    combined = res.stdout + res.stderr
    assert "cli.py" not in combined
    assert "INFO" not in combined


def test_list_models_excludes_unimplemented_models():
    """The catalog only advertises models that actually exist (no dangling entries)."""
    res = _patchsim("list-models", "--json")
    assert res.returncode == 0
    names = {m["name"] for m in json.loads(res.stdout)["models"]}
    assert names == {"seir", "sir", "sirs", "sis"}

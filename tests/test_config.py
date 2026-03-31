import pytest

from patchsim.core.simulation import load_config


def test_load_config_success(tmp_path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text(
        "PatchFile: foo\n"
        "SeedFile: bar\n"
        "OutputDir: out\n"
        "TMax: 10\n"
        "Transitions:\n"
        '  "S -> I": "beta"\n'
    )
    cfg = load_config(str(cfg_file))
    assert "PatchFile" in cfg
    assert cfg["OutputDir"] == "out"


def test_load_config_missing_field(tmp_path):
    cfg_file = tmp_path / "bad.yaml"
    cfg_file.write_text("PatchFile: foo\n")  # missing SeedFile and OutputDir
    with pytest.raises(ValueError):
        load_config(str(cfg_file))

"""Empirical sanity checks on simulation dynamics.

These lock in epidemiologically-correct behaviour for a multi-patch SIRS model
(waning immunity) so future changes can't silently break the maths.
"""

import numpy as np
import pandas as pd
import yaml

from patchsim.core.model_runner import Model
from patchsim.core.simulation import load_config, setup_simulation


def _solve_sirs(tmp_path, *, waning=0.02, tmax=160, npatches=5, seed_patch=0):
    """Build a multi-patch SIRS scenario with a ring network and solve it."""
    data = tmp_path / "data"
    for sub in ("patch", "seeds", "networks"):
        (data / sub).mkdir(parents=True, exist_ok=True)

    patches = [f"P{i}" for i in range(npatches)]
    pops = [1000 * (i + 1) for i in range(npatches)]
    pd.DataFrame({"patch": patches, "Population": pops}).to_csv(data / "patch" / "p.csv", index=False)

    infected = [5 if i == seed_patch else 0 for i in range(npatches)]
    pd.DataFrame(
        {
            "patch": patches,
            "S": [p - inf for p, inf in zip(pops, infected, strict=True)],
            "I": infected,
            "R": [0] * npatches,
        }
    ).to_csv(data / "seeds" / "s.csv", index=False)

    rows = []
    for i, src in enumerate(patches):
        rows.append((0, src, src, 0.8))
        rows.append((0, src, patches[(i + 1) % npatches], 0.1))
        rows.append((0, src, patches[(i - 1) % npatches], 0.1))
    pd.DataFrame(rows, columns=["day", "source", "target", "weight"]).to_csv(data / "networks" / "n.csv", index=False)

    cfg = {
        "ModelName": "sanity",
        "PatchFile": str(data / "patch" / "p.csv"),
        "SeedFile": str(data / "seeds" / "s.csv"),
        "NetworkFile": str(data / "networks" / "n.csv"),
        "OutputDir": str(tmp_path / "out"),
        "TMax": tmax,
        "compartments": ["S", "I", "R"],
        "Parameters": {"beta": 0.35, "gamma": 0.1, "waning": waning},
        "Transitions": {"S -> I": "beta", "I -> R": "gamma * I", "R -> S": "waning * R"},
    }
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    config = load_config(str(cfg_path))
    net, y0, patches_out, num_patches = setup_simulation(config)
    model = Model(net, compartments=list(net.base_model.compartments))
    t = np.arange(tmax, dtype=float)
    return model.solve(y0, t), num_patches, t


def _total(out, num_patches, comp):
    return np.sum([out[f"{comp}_{i}"] for i in range(num_patches)], axis=0)


def test_population_is_conserved(tmp_path):
    out, n, t = _solve_sirs(tmp_path)
    total = sum(_total(out, n, c) for c in ("S", "I", "R"))
    assert np.max(np.abs(total - total[0])) < 1e-3 * total[0]


def test_no_compartment_goes_negative(tmp_path):
    out, n, t = _solve_sirs(tmp_path)
    for key, series in out.items():
        assert np.min(series) > -1e-6, f"{key} went negative"


def test_epidemic_grows_from_seed(tmp_path):
    out, n, t = _solve_sirs(tmp_path)
    total_I = _total(out, n, "I")
    assert total_I.max() > 5 * total_I[0]


def test_infection_spreads_through_network_to_unseeded_patch(tmp_path):
    out, n, t = _solve_sirs(tmp_path, seed_patch=0)
    assert out["I_2"][0] == 0  # P2 starts with no infections
    assert out["I_2"].max() > 1.0  # but gets infected via the network


def test_waning_immunity_sustains_endemic_infection(tmp_path):
    with_waning, n, t = _solve_sirs(tmp_path, waning=0.05, tmax=200)
    without_waning, n2, t2 = _solve_sirs(tmp_path, waning=0.0, tmax=200)
    final_with = sum(with_waning[f"I_{i}"][-1] for i in range(n))
    final_without = sum(without_waning[f"I_{i}"][-1] for i in range(n2))
    assert final_with > final_without  # waning replenishes S -> infection persists
    assert final_with > 1.0

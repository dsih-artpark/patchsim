from patchsim.core.simulation import setup_simulation, load_config
import pandas as pd
import numpy as np
import pytest

def test_setup_simulation_returns_objects(tmp_data_dir):
    cfg = load_config(tmp_data_dir["config"])
    net, y0, patches, num_patches = setup_simulation(cfg)
    # basic assertions:
    assert num_patches == 2
    assert isinstance(y0, dict)
    # check that y0 keys include S_0 and I_1 etc.
    assert "S_0" in y0 and "I_1" in y0
    assert patches == ["A", "B"]

def test_setup_simulation_population_check(tmp_data_dir):
    # mutate seed so it doesn't sum to population to trigger error
    import pandas as pd
    seed_df = pd.read_csv(tmp_data_dir["seed_csv"])
    seed_df.loc[0, "S"] = 500  # break conservation for PatchA
    seed_df.to_csv(tmp_data_dir["seed_csv"], index=False)
    cfg = load_config(tmp_data_dir["config"])
    with pytest.raises(ValueError):
        setup_simulation(cfg)

import pytest

from patchsim.core.model import CompartmentalModel, NetworkModel


def test_parameter_and_compartment_names_must_be_distinct():
    with pytest.raises(ValueError, match=r"shared names: \['S'\]"):
        CompartmentalModel(compartments=["S", "I"], parameters={"S": 0.5}, transitions=[])


def test_compute_rates_with_parameters():
    compartments = ["S", "I", "R"]
    transitions = [
        {"transition": "S->I", "rate": "beta * S * I / (S + I + R)"},
        {"transition": "I->R", "rate": "gamma * I"},
    ]
    params = {"beta": 0.5, "gamma": 0.1}
    model = CompartmentalModel(compartments=compartments, parameters=params, transitions=transitions)

    state = {"S": 1000.0, "I": 10.0, "R": 0.0}
    rates = model.compute_rates(state)
    # expected infection rate per susceptible: beta * S * I / N  => numeric value
    # check keys present
    assert "S->I" in rates
    assert "I->R" in rates
    # numeric sanity: rates positive
    assert rates["S->I"] > 0
    assert rates["I->R"] == pytest.approx(0.1 * state["I"])


def test_networkmodel_derivatives_conserve_population():
    compartments = ["S", "I", "R"]
    transitions = [
        {"transition": "S->I", "rate": "beta * S * I / (S + I + R)"},
        {"transition": "I->R", "rate": "gamma * I"},
    ]
    params = {"beta": 0.5, "gamma": 0.1}
    base = CompartmentalModel(compartments=compartments, parameters=params, transitions=transitions)
    net_matrix = [[0, 0], [0, 0]]
    nm = NetworkModel(base_model=base, num_patches=2, network_matrix=net_matrix)
    # attach patch params to keep code path consistent (optional)
    nm.patch_parameters = {"PatchA": params, "PatchB": params}
    state = {"S_0": 990.0, "I_0": 10.0, "R_0": 0.0, "S_1": 995.0, "I_1": 5.0, "R_1": 0.0}
    deriv = nm.compute_derivatives(state)
    # population conservation: sum deriv across all compartments ~ 0
    total_change = sum(deriv.values())
    assert total_change == pytest.approx(0.0)

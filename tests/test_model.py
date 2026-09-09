import logging

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


def _two_patch_model(rate_expr: str, params: dict[str, float]) -> NetworkModel:
    base = CompartmentalModel(
        compartments=["S", "I", "R"],
        parameters=params,
        transitions=[{"transition": "S->I", "rate": rate_expr}, {"transition": "I->R", "rate": "gamma * I"}],
    )
    return NetworkModel(base_model=base, num_patches=2, network_matrix=[[0.9, 0.1], [0.1, 0.9]])


_TWO_PATCH_STATE = {"S_0": 990.0, "I_0": 10.0, "R_0": 0.0, "S_1": 995.0, "I_1": 5.0, "R_1": 0.0}


def test_network_infection_flow_is_local_flow_times_infectious_pressure():
    deriv = _two_patch_model("beta", {"beta": 0.4, "gamma": 0.1}).compute_derivatives(_TWO_PATCH_STATE)
    lambda_0 = 0.9 * (10.0 / 1000.0) + 0.1 * (5.0 / 1000.0)
    assert deriv["S_0"] == pytest.approx(-0.4 * 990.0 * lambda_0)


def test_network_infection_flow_keeps_every_factor_of_the_rate_expression():
    params = {"beta": 0.4, "gamma": 0.1}
    full = _two_patch_model("beta", params).compute_derivatives(_TWO_PATCH_STATE)
    halved = _two_patch_model("beta * 0.5", params).compute_derivatives(_TWO_PATCH_STATE)
    assert halved["S_0"] == pytest.approx(0.5 * full["S_0"])
    assert halved["S_1"] == pytest.approx(0.5 * full["S_1"])


def test_network_infection_flow_does_not_depend_on_the_parameter_name():
    named_beta = _two_patch_model("beta * 0.5", {"beta": 0.4, "gamma": 0.1}).compute_derivatives(_TWO_PATCH_STATE)
    named_b = _two_patch_model("b * 0.5", {"b": 0.4, "gamma": 0.1}).compute_derivatives(_TWO_PATCH_STATE)
    assert named_b["S_0"] == pytest.approx(named_beta["S_0"])


def test_recovery_flow_is_not_scaled_by_infectious_pressure_under_mixing():
    deriv = _two_patch_model("beta", {"beta": 0.4, "gamma": 0.1}).compute_derivatives(_TWO_PATCH_STATE)
    assert deriv["R_0"] == pytest.approx(0.1 * 10.0)


def test_warns_when_coupled_infection_expression_names_an_infectious_compartment(caplog):
    with caplog.at_level(logging.WARNING):
        _two_patch_model("beta * I / (S + I + R)", {"beta": 0.4, "gamma": 0.1})
    assert "infectious pressure" in caplog.text


def test_no_warning_for_per_capita_infection_expression(caplog):
    with caplog.at_level(logging.WARNING):
        _two_patch_model("beta", {"beta": 0.4, "gamma": 0.1})
    assert "infectious pressure" not in caplog.text

"""Tests for explicit-Euler stepping and numerical safeguards."""

import pytest

from patchsim.core.model import CompartmentalModel, NetworkModel


def _decay_network(gamma: float) -> NetworkModel:
    """Single patch, pure decay I -> R, so each Euler step has a closed form."""
    base = CompartmentalModel(
        compartments=["I", "R"],
        parameters={"gamma": gamma},
        transitions=[{"transition": "I->R", "rate": "gamma * I"}],
    )
    return NetworkModel(base_model=base, num_patches=1, network_matrix=[[0.0]])


def test_step_size_scales_the_update():
    net = _decay_network(gamma=0.2)

    history = net.simulate_discrete({"I_0": 100.0, "R_0": 0.0}, [0.0, 0.5])

    # one explicit Euler step of size dt: I = I0 - gamma * I0 * dt
    assert history["I_0"][-1] == pytest.approx(100.0 - 0.2 * 100.0 * 0.5)


def test_halving_the_step_size_halves_the_first_move():
    net = _decay_network(gamma=0.2)

    coarse = net.simulate_discrete({"I_0": 100.0, "R_0": 0.0}, [0.0, 1.0])
    fine = net.simulate_discrete({"I_0": 100.0, "R_0": 0.0}, [0.0, 0.5])

    coarse_move = 100.0 - coarse["I_0"][-1]
    fine_move = 100.0 - fine["I_0"][-1]
    assert fine_move == pytest.approx(coarse_move / 2.0)


def test_each_interval_uses_its_own_step_size():
    """A time grid names the points to simulate; it need not be evenly spaced."""
    net = _decay_network(gamma=0.2)

    history = net.simulate_discrete({"I_0": 100.0, "R_0": 0.0}, [0.0, 1.0, 1.5])

    # first interval dt=1.0:  100 - 0.2*100*1.0 = 80
    # second interval dt=0.5:  80 - 0.2*80*0.5  = 72
    assert history["I_0"][1] == pytest.approx(80.0)
    assert history["I_0"][2] == pytest.approx(72.0)


def test_negative_initial_population_is_rejected():
    net = _decay_network(gamma=0.2)

    with pytest.raises(ValueError, match="non-negative"):
        net.simulate_discrete({"I_0": -1.0, "R_0": 0.0}, [0.0, 1.0])


def test_initial_state_is_validated_even_when_no_steps_are_taken():
    """A single-point grid takes no steps but must still reject an invalid initial state."""
    net = _decay_network(gamma=0.2)

    with pytest.raises(ValueError, match="non-negative"):
        net.simulate_discrete({"I_0": -1.0, "R_0": 0.0}, [0.0])


def test_missing_initial_state_key_is_reported_clearly():
    net = _decay_network(gamma=0.2)

    with pytest.raises(ValueError, match=r"missing=\['R_0'\]"):
        net.simulate_discrete({"I_0": 100.0}, [0.0, 1.0])


def test_extra_initial_state_key_is_reported_clearly():
    net = _decay_network(gamma=0.2)

    with pytest.raises(ValueError, match=r"extra=\['X_0'\]"):
        net.simulate_discrete({"I_0": 100.0, "R_0": 0.0, "X_0": 1.0}, [0.0, 1.0])


def test_initial_total_overflowing_to_infinity_is_rejected():
    """Each compartment is finite but their sum is not."""
    import sys

    half_max = sys.float_info.max
    net = _decay_network(gamma=0.2)

    with pytest.raises(ValueError, match="finite"):
        net.simulate_discrete({"I_0": half_max, "R_0": half_max}, [0.0, 1.0])


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_non_finite_initial_state_is_rejected(bad):
    """A compartment unused by transitions must still be validated."""
    base = CompartmentalModel(
        compartments=["I", "R", "X"],
        parameters={"gamma": 0.2},
        transitions=[{"transition": "I->R", "rate": "gamma * I"}],
    )
    net = NetworkModel(base_model=base, num_patches=1, network_matrix=[[0.0]])

    with pytest.raises(ValueError, match="finite"):
        net.simulate_discrete({"I_0": 100.0, "R_0": 0.0, "X_0": bad}, [0.0, 1.0])


@pytest.mark.parametrize("times", [[0.0, float("nan"), 2.0], [0.0, float("inf")]])
def test_non_finite_time_points_are_rejected(times):
    net = _decay_network(gamma=0.2)

    with pytest.raises(ValueError, match="finite"):
        net.simulate_discrete({"I_0": 100.0, "R_0": 0.0}, times)


@pytest.mark.parametrize("times", [[0.0, 0.0, 0.0], [0.0, -1.0, -2.0]])
def test_non_positive_time_steps_are_rejected(times):
    """A zero step never advances; a negative step integrates backwards in time."""
    net = _decay_network(gamma=0.2)

    with pytest.raises(ValueError, match="strictly increasing"):
        net.simulate_discrete({"I_0": 100.0, "R_0": 0.0}, times)


def test_negative_population_is_raised_not_returned():
    # gamma * dt > 1 makes explicit Euler overshoot below zero
    net = _decay_network(gamma=5.0)

    with pytest.raises(ValueError, match="TimeStep"):
        net.simulate_discrete({"I_0": 100.0, "R_0": 0.0}, [0.0, 1.0, 2.0])


def test_floating_point_noise_near_zero_is_not_flagged():
    """A compartment emptying to zero lands just below it; that is arithmetic, not divergence."""
    net = _decay_network(gamma=1.0 + 1e-15)

    history = net.simulate_discrete({"I_0": 100.0, "R_0": 0.0}, [0.0, 1.0])

    assert history["I_0"][-1] == pytest.approx(0.0, abs=1e-9)


def test_roundoff_allowance_scales_with_large_compartment():
    net = _decay_network(gamma=1.0 + 1e-15)

    history = net.simulate_discrete({"I_0": 1e12, "R_0": 0.0}, [0.0, 1.0])

    assert history["I_0"][-1] == pytest.approx(0.0, abs=0.01)


def test_large_compartment_still_rejects_material_overshoot():
    net = _decay_network(gamma=1.0 + 1e-12)

    with pytest.raises(ValueError, match="negative"):
        net.simulate_discrete({"I_0": 1e12, "R_0": 0.0}, [0.0, 1.0])


def test_roundoff_allowance_tracks_compartment_growth():
    base = CompartmentalModel(
        compartments=["I", "B", "C"],
        parameters={"gamma": 1.0 + 1e-15},
        transitions=[
            {"transition": "I->B", "rate": "I"},
            {"transition": "B->C", "rate": "gamma * B"},
        ],
    )
    net = NetworkModel(base_model=base, num_patches=1, network_matrix=[[0.0]])

    history = net.simulate_discrete(
        {"I_0": 1e12, "B_0": 0.0, "C_0": 0.0},
        [0.0, 1.0, 2.0],
    )

    assert history["B_0"][-1] == pytest.approx(0.0, abs=0.01)


def test_total_population_is_conserved_exactly():
    """Clipping residue to zero would invent mass, breaking a core model invariant.

    Tolerating a tiny negative without raising is the point of the tolerance; rewriting
    the value is a different thing, and it silently violates conservation.
    """
    net = _decay_network(gamma=1.0 + 1e-15)

    history = net.simulate_discrete({"I_0": 100.0, "R_0": 0.0}, [0.0, 1.0])

    # Clipping would fabricate ~1.1e-13 here, about 8x the float64 precision at this
    # magnitude (~1.4e-14), so the bound below distinguishes the two behaviours.
    assert history["I_0"][-1] + history["R_0"][-1] == pytest.approx(100.0, abs=5e-14)


def test_large_other_compartment_does_not_hide_material_negative_value():
    base = CompartmentalModel(
        compartments=["I", "R", "X"],
        parameters={"gamma": 1.0 + 1e-8},
        transitions=[{"transition": "I->R", "rate": "gamma * I"}],
    )
    net = NetworkModel(base_model=base, num_patches=1, network_matrix=[[0.0]])

    with pytest.raises(ValueError, match="negative"):
        net.simulate_discrete(
            {"I_0": 1.0, "R_0": 0.0, "X_0": 1e18},
            [0.0, 1.0],
        )


def test_slightly_irregular_grid_is_integrated_not_rejected():
    """Grids from real data are rarely exactly even; each interval is honoured as given."""
    net = _decay_network(gamma=0.2)

    history = net.simulate_discrete({"I_0": 100.0, "R_0": 0.0}, [0.0, 1.0, 2.000005])

    assert history["I_0"][2] == pytest.approx(80.0 - 0.2 * 80.0 * 1.000005)


def test_unit_step_behaviour_is_unchanged():
    net = _decay_network(gamma=0.2)

    history = net.simulate_discrete({"I_0": 100.0, "R_0": 0.0}, [0.0, 1.0, 2.0])

    assert history["I_0"][1] == pytest.approx(80.0)
    assert history["I_0"][2] == pytest.approx(64.0)

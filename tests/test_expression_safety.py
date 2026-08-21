"""Rate expressions come from user-supplied YAML, so they are untrusted input.

Configs are the unit PatchSim asks users to share, so evaluating one must never be able
to run arbitrary code.
"""

import warnings  # noqa: F401  - loads catch_warnings, which the escape payload reaches

import numpy as np
import pytest

from patchsim.core import expressions
from patchsim.core.expressions import evaluate
from patchsim.core.model import CompartmentalModel


def _model(rate: str) -> CompartmentalModel:
    return CompartmentalModel(
        compartments=["S", "I"],
        parameters={"beta": 0.5},
        transitions=[{"transition": "S->I", "rate": rate}],
    )


def test_rate_expression_cannot_write_to_the_filesystem(tmp_path):
    marker = tmp_path / "escaped.txt"
    payload = (
        "[c for c in ().__class__.__bases__[0].__subclasses__() "
        "if c.__name__ == 'catch_warnings'][0]()._module.__builtins__"
        f"['open']({str(marker)!r}, 'w').write('escaped')"
    )

    try:
        _model(payload).compute_rates({"S": 1.0, "I": 1.0})
    except ValueError:
        pass  # rejecting the expression is the desired outcome

    assert not marker.exists(), "a rate expression executed arbitrary code"


@pytest.mark.parametrize(
    "expr",
    [
        "().__class__",
        "S.__class__",
        "().__class__.__bases__[0].__subclasses__()",
        "__import__('os')",
        "open('x', 'w')",
    ],
)
def test_disallowed_constructs_are_rejected(expr):
    with pytest.raises(ValueError, match="not allowed"):
        _model(expr).compute_rates({"S": 1.0, "I": 1.0})


def test_unknown_name_is_reported_clearly():
    with pytest.raises(ValueError, match="undefined_param"):
        _model("undefined_param * S").compute_rates({"S": 1.0, "I": 1.0})


@pytest.mark.parametrize(
    ("expr", "reason"),
    [
        ("1e308 * 1e308", "overflow to infinity"),
        ("0 - 1e308 * 1e308", "overflow to negative infinity"),
        ("9**9**9**9", "exponent overflow"),
    ],
)
def test_non_finite_results_are_rejected(expr, reason):
    """A rate of inf silently poisons every downstream compartment."""
    with pytest.raises(ValueError, match="finite"):
        _model(expr).compute_rates({"S": 1.0, "I": 1.0})


def test_complex_results_are_rejected():
    """A negative base with a fractional exponent yields a complex number in Python."""
    with pytest.raises(ValueError, match="real"):
        _model("(0 - 1) ** 0.5").compute_rates({"S": 1.0, "I": 1.0})


def test_division_by_zero_is_reported_as_an_expression_error():
    with pytest.raises(ValueError, match="division"):
        _model("beta / 0").compute_rates({"S": 1.0, "I": 1.0})


@pytest.mark.parametrize("expr", ["S // 2", "S % 2"])
def test_discontinuous_binary_operators_are_rejected(expr):
    with pytest.raises(ValueError, match="not allowed"):
        _model(expr).compute_rates({"S": 4.0, "I": 1.0})


def test_deeply_nested_expression_is_rejected_not_crashed():
    """Nesting deep enough to exhaust the interpreter stack must surface as a config error.

    Parentheses collapse during parsing, but repeated unary operators build real AST depth.
    This stays under the length cap so it exercises the depth guard rather than that one.
    """
    expr = "-" * 100 + "beta"
    with pytest.raises(ValueError, match="too deeply nested"):
        _model(expr).compute_rates({"S": 1.0, "I": 1.0})


@pytest.mark.parametrize("expr", ["1e400", "0 - 1e400"])
def test_non_finite_literals_are_rejected(expr):
    """A literal large enough to overflow reaches float() as inf without any operator."""
    with pytest.raises(ValueError, match="finite"):
        _model(expr).compute_rates({"S": 1.0, "I": 1.0})


def test_non_finite_parameter_is_rejected():
    model = CompartmentalModel(
        compartments=["S", "I"],
        parameters={"beta": float("inf")},
        transitions=[{"transition": "S->I", "rate": "beta"}],
    )
    with pytest.raises(ValueError, match="finite"):
        model.compute_rates({"S": 1.0, "I": 1.0})


def test_oversized_integer_literal_is_rejected():
    """An integer too large for float must not escape as a raw OverflowError."""
    with pytest.raises(ValueError, match="too large"):
        _model("9" * 400).compute_rates({"S": 1.0, "I": 1.0})


def test_numpy_scalars_are_accepted():
    """Solver state arrives as numpy scalars; only float64 subclasses Python float."""
    assert evaluate("beta * S", {"beta": np.float32(0.5), "S": np.int64(4)}) == pytest.approx(2.0)


def test_large_exponent_is_rejected_quickly():
    """Operands are coerced to float before the operator, so a huge power overflows and is
    rejected instead of triggering an expensive big-integer computation."""
    with pytest.raises(ValueError, match="finite"):
        _model("2 ** 100000").compute_rates({"S": 1.0, "I": 1.0})


def test_oversized_expression_is_rejected_before_parsing():
    """Bound the input before handing it to the parser, not only the AST depth after."""
    with pytest.raises(ValueError, match="too long"):
        _model("1+" * 10000 + "1").compute_rates({"S": 1.0, "I": 1.0})


def test_parsed_expression_is_reused(monkeypatch):
    expressions._parse.cache_clear()
    parse = expressions.ast.parse
    parse_calls = 0

    def counting_parse(*args, **kwargs):
        nonlocal parse_calls
        parse_calls += 1
        return parse(*args, **kwargs)

    monkeypatch.setattr(expressions.ast, "parse", counting_parse)

    assert evaluate("beta * S", {"beta": 0.5, "S": 4.0}) == pytest.approx(2.0)
    assert evaluate("beta * S", {"beta": 0.25, "S": 8.0}) == pytest.approx(2.0)
    assert parse_calls == 1


def test_arithmetic_expressions_still_evaluate():
    rates = _model("beta * S * I").compute_rates({"S": 10.0, "I": 2.0})
    assert rates["S->I"] == pytest.approx(0.5 * 10.0 * 2.0)


def test_numeric_literals_and_precedence_are_preserved():
    rates = _model("beta * S * I / (S + I) + 1").compute_rates({"S": 3.0, "I": 1.0})
    assert rates["S->I"] == pytest.approx(0.5 * 3.0 * 1.0 / 4.0 + 1)

"""Evaluate transition rate expressions from user config files.

Evaluation is restricted to arithmetic over named values so that a config file cannot
execute arbitrary code. Any other construct is rejected.
"""

import ast
import math
import numbers
import operator
from functools import lru_cache
from typing import Any, Mapping

_MAX_DEPTH = 64  # maximum expression nesting depth
_MAX_LENGTH = 1000  # maximum expression length in characters

_BINARY_OPS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.FloorDiv: operator.floordiv,
    ast.Mod: operator.mod,
    ast.Pow: operator.pow,
}

_UNARY_OPS = {
    ast.UAdd: operator.pos,
    ast.USub: operator.neg,
}


@lru_cache(maxsize=256)
def _parse(expression: str) -> ast.AST:
    if len(expression) > _MAX_LENGTH:
        raise ValueError(f"is too long to evaluate ({len(expression)} characters, limit {_MAX_LENGTH})")

    try:
        tree = ast.parse(expression, mode="eval")
    except SyntaxError as exc:
        raise ValueError(f"could not be parsed: {exc.msg}") from exc
    except (RecursionError, MemoryError) as exc:
        raise ValueError("is too deeply nested to parse") from exc

    return tree.body


def evaluate(expression: str, scope: Mapping[str, Any]) -> float:
    """Evaluate an arithmetic rate expression against a scope of named values.

    Args:
        expression: Arithmetic expression over parameter and compartment names.
        scope: Mapping of names available to the expression.

    Returns:
        The numeric value of the expression.

    Raises:
        ValueError: If the expression uses an unsupported construct, references an
            unknown name, or cannot be evaluated numerically.
    """
    return _evaluate_node(_parse(expression), scope, depth=0)


def _evaluate_node(node: ast.AST, scope: Mapping[str, Any], depth: int) -> float:
    if depth > _MAX_DEPTH:
        raise ValueError(f"is too deeply nested to evaluate (limit {_MAX_DEPTH} levels)")

    if isinstance(node, ast.Constant):
        if isinstance(node.value, bool) or not isinstance(node.value, numbers.Real):
            raise ValueError(f"{type(node.value).__name__} literals are not allowed in rate expressions")
        return _coerce(node.value, "a numeric literal")

    if isinstance(node, ast.Name):
        if node.id not in scope:
            raise ValueError(f"unknown name '{node.id}'; expected a parameter or compartment")
        return _coerce(scope[node.id], f"name '{node.id}'")

    if isinstance(node, ast.BinOp) and type(node.op) in _BINARY_OPS:
        left = _evaluate_node(node.left, scope, depth + 1)
        right = _evaluate_node(node.right, scope, depth + 1)
        return _apply(_BINARY_OPS[type(node.op)], left, right)

    if isinstance(node, ast.UnaryOp) and type(node.op) in _UNARY_OPS:
        return _apply(_UNARY_OPS[type(node.op)], _evaluate_node(node.operand, scope, depth + 1))

    raise ValueError(
        f"{_describe(node)} is not allowed in rate expressions, "
        "which may only combine parameters and compartments arithmetically"
    )


def _coerce(value: Any, description: str) -> float:
    """Convert a value to a finite real float, or raise ValueError.

    Accepts ``numbers.Real`` rather than ``float`` because solver state can be a numpy
    scalar and only ``numpy.float64`` subclasses ``float``. Complex results (from a
    negative base to a fractional power) and non-finite results (from overflow) are
    rejected.
    """
    if isinstance(value, complex):
        raise ValueError(f"{description} is a complex number; rates must be real")
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ValueError(f"{description} is not a real number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{description} is too large to represent as a floating point number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{description} is {result}; rates must be finite")
    return result


def _apply(op, *operands: float) -> float:
    """Apply an arithmetic operator and validate its result."""
    try:
        value = op(*operands)
    except ZeroDivisionError as exc:
        raise ValueError("attempted division by zero") from exc
    except OverflowError as exc:
        raise ValueError("overflowed to a value that is not finite") from exc
    except ValueError as exc:
        raise ValueError(f"could not be evaluated: {exc}") from exc

    return _coerce(value, "the expression")


def _describe(node: ast.AST) -> str:
    """Return a readable name for a disallowed construct, for the error message."""
    descriptions = {
        ast.Attribute: "attribute access",
        ast.Call: "function calls",
        ast.Subscript: "indexing",
        ast.ListComp: "comprehensions",
        ast.GeneratorExp: "comprehensions",
        ast.DictComp: "comprehensions",
        ast.SetComp: "comprehensions",
        ast.Lambda: "lambdas",
        ast.List: "list literals",
        ast.Tuple: "tuple literals",
        ast.Dict: "dict literals",
    }
    for node_type, description in descriptions.items():
        if isinstance(node, node_type):
            return description
    return f"{type(node).__name__} expressions"

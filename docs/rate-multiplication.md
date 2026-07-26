# Rate multiplication

Transition expressions can be written either as per-capita rates or as complete
flows. PatchSim distinguishes the two forms by checking whether the source
compartment appears as an identifier in the expression.

## Source rule

For `A -> B`:

- if the expression does not name `A`, PatchSim multiplies it by `A`;
- if the expression names `A`, PatchSim uses the evaluated expression as the
  flow.

| Transition | Expression | Flow |
| --- | --- | --- |
| `I -> R` | `"gamma"` | $\gamma I$ |
| `I -> R` | `"gamma * I"` | $\gamma I$ |
| `E -> I` | `"sigma"` | $\sigma E$ |
| `S -> I` | `"beta * I / N"` | $\beta SI/N$ if `N` is a parameter |
| `S -> I` | `"beta * S * I / N"` | $\beta SI/N$ |

This is an identifier check, not dimensional analysis. An expression that
mentions the source for another reason is still treated as a complete flow.

## Multi-patch infection rule

The CLI ODE runner applies one additional rule to transitions from `S` to `I`
or `E` when more than one patch is present: it multiplies the local flow by the
network infectious pressure

$$
\lambda_i = \sum_j W_{ij}\frac{I_j}{N_j}.
$$

The built-in multi-patch templates use:

```yaml
"S -> I": "beta"
```

The source rule produces $\beta S_i$, then network coupling produces
$\beta S_i\lambda_i$.

Do not add a second infectious proportion to this template expression unless
that extra factor is part of the intended model.

## Single-patch infection rule

A one-patch model does not apply network coupling. To represent
frequency-dependent infection, write the infectious proportion explicitly:

```yaml
"S -> I": "beta * I / (S + I + R)"
```

PatchSim multiplies this by `S`, producing $\beta SI/N$.

## Expression safety

Expressions are parsed as arithmetic, not executed as Python. They may use
finite real numbers, known names, and the documented arithmetic operators.
Calls, attributes, indexing, comparisons, boolean logic, and comprehensions are
rejected. Results must be finite and real.

See [Configuration](configuration.md) for the complete operator list and
[Mathematical model](mathematical-model.md) for the equations.

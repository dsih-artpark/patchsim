# Mathematical model

PatchSim represents each patch with continuous-valued compartments and
transitions between them. The current CLI solves the resulting ordinary
differential equations (ODEs).

## Transition flow

For a transition from compartment $A$ to compartment $B$, let $F_{A\to B}$ be
the flow evaluated from the YAML expression. PatchSim contributes

$$
\frac{dA}{dt} = -F_{A\to B},
\qquad
\frac{dB}{dt} = F_{A\to B}.
$$

The expression-to-flow rule is:

$$
F_{A\to B} =
\begin{cases}
A\,r(\mathbf{x}, \mathbf{p}), & \text{if the expression does not name } A,\\
r(\mathbf{x}, \mathbf{p}), & \text{if the expression names } A.
\end{cases}
$$

For example, both recovery expressions below produce $\gamma I$:

```yaml
"I -> R": "gamma"
"I -> R": "gamma * I"
```

See [Rate multiplication](rate-multiplication.md) for the exact expression
rules.

## Single-patch SIR

Network coupling is not applied when the model has one patch. A
frequency-dependent SIR model can be written as:

```yaml
compartments: [S, I, R]
Parameters:
  beta: 0.2
  gamma: 0.1
Transitions:
  "S -> I": "beta * I / (S + I + R)"
  "I -> R": "gamma"
```

With $N=S+I+R$, this gives:

$$
\frac{dS}{dt} = -\beta S\frac{I}{N},
$$

$$
\frac{dI}{dt} = \beta S\frac{I}{N} - \gamma I,
$$

$$
\frac{dR}{dt} = \gamma I.
$$

An explicit density-dependent expression such as `"beta * S * I"` is also
accepted, but it represents a different model and requires a correspondingly
scaled $\beta$.

## Multi-patch infectious pressure

For patch $i$, let $I_i$ be infectious population and $N_i$ total population.
The current runtime computes

$$
\lambda_i(t) = \sum_{j=1}^{n} W_{ij}\frac{I_j(t)}{N_j(t)}.
$$

`W[i, j]` is loaded from the day-zero CSV row whose `source` is patch $i$ and
whose `target` is patch $j$. Thus row $i$ selects the infectious patches that
contribute to focal patch $i$. See [Network design](network-design.md) for the
file convention.

For an `S -> I` or `S -> E` transition in a multi-patch ODE run, PatchSim
multiplies the local flow by $\lambda_i$. The built-in SIR template uses:

```yaml
"S -> I": "beta"
"I -> R": "gamma * I"
```

The local source rule first gives $\beta S_i$, and network coupling gives:

$$
\frac{dS_i}{dt} = -\beta S_i\lambda_i,
$$

$$
\frac{dI_i}{dt} = \beta S_i\lambda_i - \gamma I_i,
$$

$$
\frac{dR_i}{dt} = \gamma I_i.
$$

Weights are not probabilities unless the input author makes them so. PatchSim
does not normalize $W$.

## Other built-in templates

The CLI can scaffold these transition structures:

| Template | Compartments | Additional transition |
| --- | --- | --- |
| `sir` | `S, I, R` | `I -> R` |
| `seir` | `S, E, I, R` | `E -> I` |
| `sirs` | `S, I, R` | `R -> S` |
| `sis` | `S, I` | `I -> S` |

They are YAML templates using the same transition engine, not separate
hard-coded model classes.

## Numerical integration

The two built-in solvers use the same derivative function and differ only in
numerical integration.

`Solver: ode` uses `scipy.integrate.odeint`, SciPy's interface to LSODA. LSODA
chooses internal step sizes adaptively. `Solver: discrete` uses deterministic
explicit Euler with one step per reporting interval:

$$
t_k = k\,\mathrm{TimeStep}, \qquad
x_{k+1} = x_k + \mathrm{TimeStep}\,f(x_k).
$$

The discrete method can fail or be inaccurate when `TimeStep` is too large.
Compare results at successively smaller intervals over the same horizon.
`Tolerance` and `MaxIter` are not passed to `odeint`.

## Assumptions and limits

- Compartments are continuous real-valued states.
- Each patch is well mixed within its compartment definitions.
- Transition expressions have no explicit time variable.
- Only compartment transfers are represented; births, deaths, or imports require
  explicit model structure.
- Patch-specific parameter records are applied by both built-in solvers.
- Network rows after `day == 0` are currently ignored.

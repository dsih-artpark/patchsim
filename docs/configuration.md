# Configuration

PatchSim reads one YAML mapping plus CSV files for patches, seeds, and optionally
a network. Run `patchsim validate -c CONFIG` before `patchsim run`.

## Conventions

- Field names are case-sensitive. Use the spelling shown on this page.
- Relative file paths are resolved from the directory containing the YAML file.
- Patch identifiers must match exactly across all CSV files.
- Patch identifiers should be unique. Output indices follow `PatchFile` row order.
- Transition keys use `SOURCE -> TARGET`.
- Quote transition expressions. This keeps YAML parsing predictable.
- Use a separate `OutputDir` for each retained scenario because the CSV and PNG
  are overwritten on rerun.

## Active fields

| Field | Required | Meaning |
| --- | --- | --- |
| `PatchFile` | Yes | Patch identifiers and positive populations |
| `SeedFile` | Yes | Initial compartment values for each seeded patch |
| `OutputDir` | Yes | Root directory for `runs`, `plots`, and `logs` |
| `Transitions` | Yes | Non-empty arrow-to-expression mapping |
| `TMax` | Yes | Positive integer count of reporting points |
| `Solver` | No | `ode` (default) or deterministic explicit Euler (`discrete`) |
| `TimeStep` | No | Positive reporting-grid interval; defaults to `1.0` |
| `ModelName` | For `run` | Identifier used in output filenames |
| `NetworkFile` | No | Day-zero network CSV; `null` creates a zero matrix |
| `GroupFile` | Grouped runs | Patch-by-group populations |
| `InteractionFile` | Grouped runs | Shared focal-by-contributor group matrix |
| `InteractionUnits` | Grouped runs | Non-empty description of interaction-weight units |
| `compartments` | No | Ordered compartment names; inferred from `SeedFile` if absent |
| `Parameters` | No | Global names available to transition expressions |
| `PatchParameters` | No | Per-patch parameter overrides |
| `Sensitivity` | For `sensitivity` | Sobol study name, seed, bounds, and scalar metrics |

`Compartments` with an uppercase `C` is also accepted for compatibility.
`compartments` is the documented spelling.

The three grouped-run fields are optional as a set and invalid when only partly
configured. See [Group stratification](group-stratification.md).

## File paths and output

```yaml
ModelName: baseline
PatchFile: data/patch/patch-population.csv
SeedFile: data/seeds/seed-initial.csv
NetworkFile: data/networks/network-static.csv
OutputDir: output/baseline
TMax: 60
Solver: ode
TimeStep: 1.0
```

If this YAML file is `/work/study/config.yaml`, `OutputDir: output/baseline`
resolves to `/work/study/output/baseline`, regardless of the shell's current
directory.

## Solver and time grid

The reporting grid contains `TMax` points:

$$
t_k = k\,\mathrm{TimeStep}, \qquad k = 0,\ldots,\mathrm{TMax}-1.
$$

`TMax` is a point count, not the final simulated time. PatchSim does not infer
calendar units. If rates are per day, `TimeStep` must be in days.

```yaml
Solver: discrete
TMax: 241
TimeStep: 0.25
```

`Solver: ode` uses LSODA with adaptive internal steps; `TimeStep` controls only
the reported times. `Solver: discrete` takes one deterministic explicit-Euler
step per reporting interval. It uses floating-point compartment states and is
not a stochastic or integer-state method.

Euler accuracy is not automatic. Compare a run using `dt` with one using
`dt / 2` over the same horizon. When halving `TimeStep`, use
`TMax_fine = 2 * (TMax_coarse - 1) + 1`. A run that remains finite and
non-negative has not thereby demonstrated convergence.

## Patch file

`PatchFile` must contain patch and population columns. Their capitalization is
matched case-insensitively.

```csv
patch,population
A,1000
B,500
```

Populations must be positive. Row order is significant: patch `A` above becomes
output suffix `_0`, and patch `B` becomes `_1`.

## Seed file

`SeedFile` must use the exact lowercase column name `patch` followed by one
column per compartment:

```csv
patch,S,I,R
A,999,1,0
B,500,0,0
```

For every seed row:

- the patch must exist in `PatchFile`;
- all compartment values must be non-negative; and
- the compartment sum must equal that patch's population within `1e-6`.

If `compartments` is provided, the seed compartment columns must match it
exactly. If it is omitted, PatchSim infers compartments from all seed columns
other than `patch`.

Grouped runs instead use one row per patch/group pair and exclude both `patch`
and `group` when inferring compartments. Group populations and the full file
contract are documented in [Group stratification](group-stratification.md).

## Network file

`NetworkFile` uses:

```csv
day,source,target,weight
0,A,A,0.9
0,A,B,0.1
0,B,A,0.1
0,B,B,0.9
```

Current behavior is precise:

- only rows with `day == 0` are loaded;
- weights must be non-negative;
- weights are not normalized;
- omitted matrix entries remain zero; and
- a repeated `(source, target)` pair is overwritten by its last day-zero row.

The matrix is stored as `W[source, target]`. The runtime computes the infectious
pressure for row patch $i$ as

$$
\lambda_i = \sum_j W_{ij}\frac{I_j}{N_j}.
$$

Therefore, in the current implementation, `source` selects the focal matrix row
and `target` selects the infectious patch contributing to that row. See
[Network design](network-design.md) before constructing a custom matrix.

If `NetworkFile` is missing or `null`, PatchSim creates a zero matrix. In a
multi-patch run this means no network infectious pressure. A single-patch run
does not apply network coupling.

## Model definition

```yaml
compartments: [S, I, R]

Parameters:
  beta: 0.08
  gamma: 0.1

Transitions:
  "S -> I": "beta"
  "I -> R": "gamma * I"
```

Parameter and compartment names become identifiers in transition expressions.
Names are case-sensitive, and a parameter cannot have the same name as a
compartment.

### Transition expression rules

The current evaluator accepts:

- finite real numeric literals;
- parameter and compartment names;
- binary `+`, `-`, `*`, `/`, `//`, `%`, and `**`; and
- unary `+` and `-`.

Function calls, attribute access, indexing, comparisons, booleans, conditionals,
and comprehensions are rejected. `//` and `%` are accepted by the current
evaluator, but they introduce discontinuities and should be avoided in rate
expressions.

An expression that does not mention its source compartment is treated as a
per-capita rate and multiplied by the source:

```yaml
"I -> R": "gamma"       # flow = gamma * I
```

If the source compartment appears as an identifier, the expression is already a
flow and is used as written:

```yaml
"I -> R": "gamma * I"   # same flow
```

See [Rate multiplication](rate-multiplication.md) for edge cases and examples.

### Infection transitions

For a multi-patch run with either built-in solver, transitions from `S` to `I`
or `E` are additionally scaled by the network infectious pressure. The built-in
templates therefore use:

```yaml
"S -> I": "beta"
```

This produces $\beta S_i\lambda_i$ in a multi-patch run.

For a one-patch frequency-dependent SIR model, make the infectious proportion
explicit because network coupling is not applied:

```yaml
"S -> I": "beta * I / (S + I + R)"
```

Because this expression does not contain `S`, PatchSim multiplies it by `S`,
giving $\beta SI/N$.

## Patch parameters

The accepted shape is a list of patch records:

```yaml
PatchParameters:
  - patch: A
    parameters:
      beta: 0.12
  - patch: B
    parameters:
      beta: 0.05
```

Patch names are validated, and each override is merged with the global parameter
mapping during setup. Both built-in solvers use these merged parameters.

## Sensitivity block

`Sensitivity` is optional for normal runs and required by `patchsim
sensitivity`. `patchsim validate` checks it when present. Version 1 supports
first-order and total-order Sobol indices for independent uniform bounds on
global parameters. It does not sample patch-specific, group-specific, integer,
categorical, correlated, network, or interaction-matrix inputs.

`BaseSamples` must be a power of two of at least 2, and `Seed` is required. `Name`
must be a single safe path component. Metrics sum exact output columns and use
either `max` or `final`; names may not collide with sampled parameter columns or
`sample_id`.

See {ref}`Simulation workflow <sensitivity-study>` for the complete
configuration and interpretation guidance.

## Accepted compatibility fields

The generated scaffold also contains:

```yaml
Logging: false
Tolerance: 1.0e-8
MaxIter: 10000
StartDate: "2020-01-01"
EndDate: "2022-12-31"
```

These fields are accepted by the schema but do not currently change the CLI ODE
solver:

- `Tolerance` and `MaxIter` are not passed to `odeint`;
- `StartDate` and `EndDate` do not change the numeric time grid; and
- `Logging: false` does not disable the timestamped run log.

`TMax` and `TimeStep` define the reporting grid.

## Complete example

```yaml
ModelName: baseline
PatchFile: data/patch/patch-population.csv
SeedFile: data/seeds/seed-initial.csv
NetworkFile: data/networks/network-static.csv
OutputDir: output/baseline
TMax: 60
Solver: ode
TimeStep: 1.0

compartments: [S, I, R]

Parameters:
  beta: 0.08
  gamma: 0.1

Transitions:
  "S -> I": "beta"
  "I -> R": "gamma * I"
```

Validate it with:

```bash
patchsim validate -c config.yaml
```

Use `patchsim validate --schema` when integrating an editor or another config
tool.

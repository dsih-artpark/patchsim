# Architecture

PatchSim currently has one configuration-driven CLI path and a small Python API.
The runtime is organized around explicit model, network, solver, and output
modules.

## Runtime flow

```text
CLI
  -> load_config
  -> setup_simulation
  -> NetworkModel + initial state
  -> NetworkModel.simulate_ode | NetworkModel.simulate_discrete
  -> CSV + PNG + log
```

## Components

### `patchsim.cli`

Defines `init`, `validate`, `run`, and `list-models`. It owns argument parsing,
human or JSON command output, and process exit behavior.

### `patchsim.core.simulation`

Owns the configuration schema and built-in template data. It resolves
config-relative paths, validates CSV inputs, constructs the dense network matrix,
creates the compartment/network model, and writes run artifacts.

### `patchsim.core.model`

Contains:

- `CompartmentalModel`, which evaluates transition expressions and applies the
  source-compartment rule; and
- `NetworkModel`, which stores patch-indexed state, computes infectious pressure,
  and provides discrete and ODE helpers.

Rate expressions are parsed by `patchsim.core.expressions` through a restricted
arithmetic evaluator.

### `patchsim.core.model_runner`

Provides the compatibility `Model` Python API. It solves the coupled ODE with
`scipy.integrate.odeint`. The CLI dispatches through `NetworkModel` so its ODE
and discrete paths share the same derivative function.

### `patchsim.utils`

`logger` writes the timestamped run log. `viz` writes the multi-panel PNG.

## Data ownership

- YAML owns model and run configuration.
- `PatchFile` owns patch order and total populations.
- `SeedFile` owns initial compartment values.
- `NetworkFile` owns the day-zero dense matrix entries.
- `OutputDir` owns user-facing run artifacts.

There is no hidden project state directory and no persistent cache in the current
implementation.

## Current extension boundary

Built-in model names are YAML template data interpreted by the same runtime.
PatchSim exposes configuration-selectable built-in solvers but no plugin
registry. New solver behavior requires Python changes; new compartmental models
can use the existing transition language.

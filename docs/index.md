---
title: PatchSim documentation
---

# PatchSim documentation

```{toctree}
:maxdepth: 2

getting-started.md
configuration.md
mathematical-model.md
network-design.md
simulation-workflow.md
architecture.md
results.md
cli-reference.md
rate-multiplication.md
```
# PatchSim

PatchSim is a modular metapopulation simulation framework for compartmental epidemiology.
It supports single-patch and multi-patch models, YAML-based configuration, and ODE-based
simulation with network-coupled transmission.

## What PatchSim provides

- **Compartment models**: SIR, SEIR, SIRS, SIS, and custom variants
- **Multi-patch networks**: Weighted directed patch-to-patch transmission
- **YAML configuration**: Model structure, parameters, seed data, and network inputs
- **CLI workflows**: Scaffold, validate, run, inspect models, and export summaries
- **Results outputs**: CSV time series, plots, and logs

## Start here

If you are new to the project, follow the docs in this order:

1. [Getting Started](getting-started.md) — install, initialize, validate, and run
2. [Configuration](configuration.md) — every YAML field and a worked SEIR example
3. [CLI Reference](cli-reference.md) — command and flag reference
4. [Mathematical Model](mathematical-model.md) — equations and the rate rule
5. [Network Design](network-design.md) — how multi-patch mixing works
6. [Results](results.md) — how to interpret CSVs, plots, and logs

## Core concepts

### Compartments

Compartments represent model states such as `S`, `E`, `I`, and `R`.

### Parameters

Parameters are per-capita rates like `beta`, `sigma`, `gamma`, and `rho`.

### Transitions

Transitions define movement between compartments using YAML expressions.
See [Mathematical Model](mathematical-model.md) for the rate multiplication rule.

### Network structure

Patch-to-patch transmission is represented by a weighted directed network.
See [Network Design](network-design.md) for the expected CSV format.

## Common workflows

- Use [Getting Started](getting-started.md) for a fast first run
- Use [Configuration](configuration.md) to edit model behavior
- Use [CLI Reference](cli-reference.md) to see flags and outputs
- Use [Results](results.md) to interpret output files

## Output overview

PatchSim writes results to the configured `OutputDir`:

- `runs/` — time series CSV output
- `plots/` — visualization images
- `logs/` — run logs with configuration and solver details

## Documentation map

- [Getting Started](getting-started.md)
- [Configuration](configuration.md)
- [CLI Reference](cli-reference.md)
- [Mathematical Model](mathematical-model.md)
- [Rate Multiplication](rate-multiplication.md)
- [Network Design](network-design.md)
- [Simulation Workflow](simulation-workflow.md)
- [Architecture](architecture.md)
- [Results](results.md)
